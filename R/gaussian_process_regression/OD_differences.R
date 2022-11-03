library(here)
library(tidyverse)
library(withr)
library(performance)
library(modelbased)
library(see)
library(parameters)
library(lme4)

source(here::here("R", "gaussian_process_regression", "lgpr_funs.R"))

# Read data ---------------------------------------------------------------

opd <- read_rds(here::here("data", "gaussian_process_regression", "formatted_for_cluster.rds"))[['opd']] %>%
  mutate(opd=10^y) %>%
  filter(day < 40) %>%
  mutate(treatment = paste(pred, evo, sep="_")) %>%
  mutate(treatment = factor(treatment, levels = c("none_none", "none_evolution", "predation_none", "predation_evolution")))


# First look at ATP data --------------------------------------------------

ggplot(opd) +
  geom_point(aes(y=opd, x=day, group=id, color=evo)) +
  geom_line(aes(y=opd, x=day, group=id, color=evo)) + 
  facet_grid(~pred) +
  scale_y_log10()


# modeling ----------------------------------------------------------------

manymodels <- function(data){
  list(
  m01 = lmer(opd ~ day + (day|id), data), # Just time. does not model the effect of treatment
  m02 = lmer(opd ~ treatment * day + (day|id), data),
  m03 = lmer(opd ~ treatment * poly(day, 2) + (poly(day, 2)|id), data),
  m04 = lmer(opd ~ treatment * day + (1|id) + (0 + day|id), data),
  m05 = lmer(opd ~ treatment * poly(day, 2) + (1 + treatment|id) + (0 + poly(day, 2)|id), data)
  )
}

with_seed(1237853,
          result <- manymodels(opd)
)


# manymodels_glm <- function(data){
#   list(
#     m01 = glm(opd ~ treatment * day, family=gaussian, data=data),
#     m02 = glm(opd ~ treatment * poly(day, 2), family=gaussian, data=data),
#     m03 = glm(opd ~ treatment * day, family = Gamma(link="inverse"), data=data),
#     m04 = glm(opd ~ treatment * poly(day, 2), family = Gamma(link="inverse"), data=data)
#   )
# }

# with_seed(1237853,
#   result_glm <- manymodels_glm(opd)
# )

# m02 and m03 are linear mixed models in which the fixed effects are the
# representative intercept and slope for the population and the random effects
# are the deviations in intercept and slope associated with mesocosm id.

# We also consider a model with uncorrelated random effects. This is model m To express this we
# use two random-effects terms with the same grouping factor and different
# left-hand sides. In the formula for an lmer model, distinct random effects
# terms are modeled as being independent. Thus we specify the model with two
# distinct random effects terms, each of which has mesocosm id as the grouping
# factor. The model matrix for one term is intercept only (1) and for the other
# term is the column for Days only, which can be written 0+Days. (The expression
# Days generates a column for Days and an intercept. To suppress the intercept
# we add 0+ to the expression; -1 also works.)

# We also include a polynomial term (m03, m05) since there is some curvature in
# the data over time

# Comparing models --------------------------------------------------------

compare_performance(result$m01,
                    result$m02,
                    result$m03,
                    result$m04,
                    result$m05,
                    rank = FALSE)

# compare_performance(result_glm$m01,
#                     result_glm$m02,
#                     result_glm$m03,
#                     result_glm$m04,
#                     rank = TRUE)

# compare_performance(result_glm$m01,
#                     result_glm$m02,
#                     result_glm$m03,
#                     result_glm$m04,
#                     result$m01,
#                     result$m02,
#                     result$m03,
#                     result$m04,
#                     result$m05,
#                     rank = FALSE)

# these results suggest that m04 is best fairly closely followed by m02

# Model m02 contains model m04 in the sense that if the parameter values for
# model m02 were constrained so as to force the correlation, and hence the
# covariance, to be zero, and the model were re-fit, we would get model m04. The
# value 0, to which the correlation is constrained, is not on the boundary of
# the allowable parameter values. In these circumstances a likelihood ratio test
# and a reference distribution of a χ 2 on 1 degree of freedom is suitable.

anova(result$m04, result$m02)

# Because the large p-value indicates that we would not reject m04 in favor of
# m02, we prefer the more parsimonious m04. This is consistent with the AIC and
# BIC values from the compare_performance test (where smaller - more negative -
# is better)

# Plotting ----------------------------------------------------------------

plot(estimate_expectation(result$m05, data = "grid"))

check_model(result$m03)

parameters(result$m05, effects = "fixed")

estimate_contrasts(result$m05, contrast = "treatment")

p1 <- estimate_means(result$m05, at = c("treatment")) %>%
  plot(show_data = "boxplot") + 
  scale_y_log10(limits = c(0.03, 1.5), breaks = c(1.0, 0.3, 0.1, 0.03)) +
  theme_bw() +
  mytheme()

p1

ggsave(plot=p1, filename=here::here("figs", "fig_opd_mod.svg"), device="svg", height = 2.5, width = 2.5)

# Report ------------------------------------------------------------------

library(report)

report(result$m05)
