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

atp <- read_rds(here::here("data", "gaussian_process_regression", "formatted_for_cluster.rds"))[['atp']] %>%
  mutate(atp=10^y) %>%
  filter(day < 40) %>%
  mutate(treatment = paste(pred, evo, sep="_")) %>%
  mutate(treatment = factor(treatment, levels = c("none_none", "none_evolution", "predation_none", "predation_evolution")))


# First look at ATP data --------------------------------------------------

ggplot(atp) +
  geom_point(aes(y=atp, x=day, group=id, color=evo)) +
  geom_line(aes(y=atp, x=day, group=id, color=evo)) + 
  facet_grid(~pred) +
  scale_y_log10()


# modeling ----------------------------------------------------------------

manymodels <- function(data){
  list(
  m01 = lmer(y ~ day + (day|id), data), # Just time. does not model the effect of treatment
  m02 = lmer(y ~ treatment * day + (day|id), data),
  m03 = lmer(y ~ treatment * poly(day, 2) + (poly(day, 2)|id), data),
  m04 = lmer(y ~ treatment * day + (1|id) + (0 + day|id), data),
  m05 = lmer(y ~ treatment * poly(day, 2) + (1 + treatment|id) + (0 + poly(day, 2)|id), data)
  )
}

with_seed(1237853,
  result <- manymodels(atp)
)

# m02 and m03 are linear mixed models in which the fixed effects are the
# representative intercept and slope for the population and the random effects
# are the deviations in intercept and slope associated with mesocosm id.

# We also consider a model with uncorrelated random effects. This is model m To express this we
# use two random-effects terms with the same grouping factor and different
# left-hand sides. In the formula for an lmer model, distinct random effects
# terms are modeled as being independent. Thus we specify the model with two
# distinct random effects terms, each of which has mesocosm id as the grouping
# factor. The model matrix for one term is intercept only (1) and for the other
# term is the column for day only, which can be written 0+day. (The expression
# day generates a column for day and an intercept. To suppress the intercept
# we add 0+ to the expression; -1 also works.)

# We also include a polynomial term (m03, m05) since there is some curvature in
# the data over time

# Comparing models --------------------------------------------------------

compare_performance(result$m01,
                    result$m02,
                    result$m03,
                    result$m04,
                    result$m05,
                    rank = TRUE)

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

check_model(result$m05)

parameters(result$m05, effects = "fixed")

estimate_contrasts(result$m05, contrast = "treatment")

p1 <- estimate_means(result$m05, at = c("treatment")) %>%
  plot(show_data = "boxplot") + 
  scale_y_continuous(limits = c(4.5, 6.5)) + #breaks = c(4.5, 5, 5.5, 6, 6.5), ) +
  theme_bw() +
  mytheme()

p1

ggsave(plot=p1, filename=here::here("figs", "fig_atp_mod.svg"), device="svg", height = 2.5, width = 2.5)


# Report ------------------------------------------------------------------

library(report)

report(result$m05)
