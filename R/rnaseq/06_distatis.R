library(here)
library(tidyverse)
library(DistatisR)
library(withr)

source(here::here("R", "rnaseq", "rnaseq_funs.R"))

# Read data ---------------------------------------------------------------
load(file = here::here("data", "rnaseq", "deseqinputs.RData"))

rlog_04 <- read_rds(here::here("data", "rnaseq", "rlog_04_STC.rds"))
rlog_45 <- read_rds(here::here("data", "rnaseq", "rlog_45_STC.rds"))

# Format ------------------------------------------------------------------
# convert each transcriptome to distance matrix
rlog_04_dist <- map(rlog_04, dist) %>% map(., as.matrix)
rlog_45_dist <- map(rlog_45, dist) %>% map(., as.matrix)

# Now combine the arrays in each list into multi-dimensional arrays
rlog_04_dist_arr <- abind::abind(rlog_04_dist, along = 3)
rlog_45_dist_arr <- abind::abind(rlog_45_dist, along = 3)

# Run DiSTATIS ------------------------------------------------------------
with_seed(123542,
          distatis_res_04 <-
            DistatisR::distatis(rlog_04_dist_arr, nfact2keep = 10))

with_seed(123542,
          distatis_res_45 <-
            DistatisR::distatis(rlog_45_dist_arr, nfact2keep = 10))

partial_df <- bind_rows(get_partial(distatis_res_04, meta_04),
                        get_partial(distatis_res_45, meta_45)) %>%
  mutate(genome = factor(table, levels = genome_order_full)) %>%
  select(sample, genome, `Factor 1`, `Factor 2`, days:summary_cat) %>%
  mutate(score = "partial")

# Get dataframe of comprimise values aka factor scores
compromise_df_04 <-
  as_tibble(distatis_res_04$res4Splus$F, rownames = "sample") %>%
  left_join(., rownames_to_column(meta_04, var = "sample"))

compromise_df_45 <-
  as_tibble(distatis_res_45$res4Splus$F, rownames = "sample") %>%
  left_join(., rownames_to_column(meta_45, var = "sample"))

compromise_df <- bind_rows(compromise_df_04, compromise_df_45) %>%
  select(sample, `Factor 1`, `Factor 2`, days:last_col()) %>%
  mutate(score = "compromise") %>%
  select(sample:summary_cat, score)

# Bootstrap of compromise
with_seed(123542,
  distatis_comp_boot_04 <-
    DistatisR::BootFromCompromise(rlog_04_dist_arr, nfact2keep = 10, 
                                  niter = 100)
)

fboot_04 <- boot2df(distatis_comp_boot_04, meta_04)

with_seed(123542,
  distatis_comp_boot_45 <-
    DistatisR::BootFromCompromise(rlog_45_dist_arr, nfact2keep = 10, 
                                  niter = 100)
)

fboot_45 <- boot2df(distatis_comp_boot_45, meta_45)

fboot <- bind_rows(fboot_04, fboot_45) %>%
  select(sample:summary_cat)

# clustering of compromise
library(NbClust)

with_seed(123542,
  clust_04 <-
    NbClust(
      distatis_res_04$res4Splus$F,
      method = "kmeans",
      index = "kl",
      min.nc = 2,
      max.nc = 8
    )
)
with_seed(123542,
  clust_45 <-
    NbClust(
      distatis_res_45$res4Splus$F,
      method = "kmeans",
      index = "kl",
      min.nc = 2,
      max.nc = 8
    )
)

clusters <-
  bind_rows(
    clust_04$Best.partition %>% 
      enframe(name = "sample", value = "cluster") %>% mutate(days = 4),
    clust_45$Best.partition %>% 
      enframe(name = "sample", value = "cluster") %>% mutate(days = 45)
  ) %>%
  left_join(., compromise_df) %>%
  mutate(cluster = factor(cluster))

# Get eigenvalues and tau (% inertia explained)
l04 <-
  enframe(distatis_res_04$res4Cmat$eigValues, name = "dimension") %>%
  dplyr::rename(lambda = value) %>%
  dplyr::mutate(day = 4)

t04 <- enframe(distatis_res_04$res4Cmat$tau, name = "dimension") %>%
  dplyr::rename(tau = value) %>%
  dplyr::mutate(day = 4)

l45 <-
  enframe(distatis_res_45$res4Cmat$eigValues, name = "dimension") %>%
  dplyr::rename(lambda = value) %>%
  dplyr::mutate(day = 45)

t45 <- enframe(distatis_res_45$res4Cmat$tau, name = "dimension") %>%
  dplyr::rename(tau = value) %>%
  dplyr::mutate(day = 45)

bind_rows(left_join(l04, t04), left_join(l45, t45)) %>%
  filter(dimension %in% c("dim 1", "dim 2"))

# dimension   lambda   day   tau
# dim 1       12.7     4     91
# dim 2       0.591    4     4
# dim 1       11.5     45    77
# dim 2       1.10     45    7

# On day 04 91% of variance explained on dim 1 and 04% on dim 2
# On day 45 77% of variance explained on dim 1 and 07% on dim 2

# DiSTATIS plot -----------------------------------------------------------

# plot
library(ggforce)

fullplot <- bind_rows(partial_df, compromise_df)

p <- ggplot(fullplot) +
  stat_ellipse(
    data = fboot,
    type = "norm",
    geom = "polygon",
    aes(x = `Factor 1`, y = `Factor 2`, fill = summary_cat)
  ) +
  geom_point(aes(x = `Factor 1`, y = `Factor 2`, color = genome),
             alpha = 1,) +  #size = 1.5
  geom_jitter(
    data = compromise_df,
    aes(x = `Factor 1`, y = `Factor 2`, fill = summary_cat),
    size = 4,
    width = 0.025,
    height = 0.025,
    shape = 22
  ) +
  geom_mark_ellipse(data = clusters, aes(x = `Factor 1`,
                                         y = `Factor 2`, 
                                         linetype = cluster)) +
  coord_fixed() +
  scale_color_manual(values = mycols19, na.value = "black") +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  labs(
    x = "Factor 1",
    y = "Factor 2",
    shape = "",
    color = ""
  ) +
  facet_grid(days ~ score) +
  theme_bw() +
  mytheme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical"
  )

p

ggsave(here("figs", "fig4.svg"), p, device="svg", width=7, height=7, units="in")

# PERMANOVA of compromise factors -----------------------------------------
library(vegan)

# dimensions 1:5 explain 100% of variance
distatis_res_04$res4Cmat$tau

# run the permanova
with_seed(1234,
  adonis2(
    distatis_res_04$res4Splus$F[, 1:5] ~ predation * pseudomonas_hist,
    data = meta_04,
    method = "euclidean",
    by = "terms"
  )
)

# dimensions 1:9 explain 100% of variance
distatis_res_45$res4Cmat$tau

with_seed(1234,
  adonis2(
    distatis_res_45$res4Splus$F[, 1:9] ~ predation * pseudomonas_hist,
    data = meta_45,
    method = "euclidean",
    by = "terms"
  )
) 


# check for homogeneity of dispersion
with_seed(12354,
          betadisp_04 <-
            betadisper(
              vegdist(distatis_res_04$res4Splus$F[, 1:5], method = "euclidean"),
              meta_04$summary_cat
            ))

# To test if one or more groups is more variable than the others, 
# ANOVA of the distances to group centroids can be performed and parametric 
# theory used to interpret the significance of F. *Not Significant*
print(anova(betadisp_04)) %>% xtable::xtable()
# anova(mod0)%>% xtable::xtable()

# permutest permutes model residuals to generate a permutation distribution of 
# F under the null hypothesis of no difference in dispersion between groups.
# *Not Significant*
print(permutest(betadisp_04, permutations = 999))

with_seed(12354,
          betadisp_45 <-
            betadisper(
              vegdist(distatis_res_45$res4Splus$F[, 1:9] , method = "euclidean"),
              meta_45$summary_cat
            ))

print(anova(betadisp_45)) %>% xtable::xtable()

# permanova on one combined dataset of all species
# produces same result as with the factor scores

perm_df_04 <-
  map(rlog_04, t) %>% map(., as.data.frame) %>% bind_rows() %>% t()

with_seed(1234,
  adonis2(
    perm_df_04 ~ predation * pseudomonas_hist,
    data = meta_04,
    method = "euclidean",
    by = "terms"
  )
)

perm_df_45 <-
  map(rlog_45, t) %>% 
  map(., as.data.frame) %>% 
  bind_rows() %>% 
  t()

with_seed(1234,
  adonis2(
    perm_df_45 ~ predation * pseudomonas_hist,
    data = meta_45,
    method = "euclidean",
    by = "terms"
  )
)

