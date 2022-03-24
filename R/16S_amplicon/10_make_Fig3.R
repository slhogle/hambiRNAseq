library(here)
library(tidyverse)

source(here::here("R", "16S_amplicon", "amplicon_funs.R"))

tax <- read_tsv(here::here("data", "species_taxonomy.tsv"))

# Fig. 3A -----------------------------------------------------------------
# Shannon diversity

shannonDivNet <- read_rds(here::here("data", "16S_amplicon", "shannon_divnet.rds")) 
shannonPlugin <- read_rds(here::here("data", "16S_amplicon", "shannon_plugin.rds"))

p_shannon <- ggplot() +
  geom_jitter(
    data = shannonPlugin,
    aes(
      x = predation,
      y = estimate,
      color = expcombo,
      group = expcombo
    ),
    size = 0.5,
    position = position_dodge2(width = 0.5)
  ) +
  geom_linerange(
    data = shannonDivNet,
    aes(
      x = predation,
      ymin = lower,
      ymax = upper,
      color = expcombo,
      group = expcombo,
      shape = days
    ),
    position = position_dodge2(width = 0.5)
  ) +
  geom_point(
    data = shannonDivNet,
    aes(
      x = predation,
      y = estimate,
      fill = expcombo,
      group = expcombo,
      shape = days
    ),
    position = position_dodge2(width = 0.5)
  ) +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_color_brewer(type = "qual", palette = "Paired") +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  facet_wrap(~ days) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    #legend.position = "none"
  )

#p_shannon

# Fig. 3B -----------------------------------------------------------------
# NMDS spider plot
nmds_ord <- read_rds(here::here("data", "16S_amplicon", "ordination.rds"))
metadf <- read_tsv(here("data", "sample_metadata.tsv"), col_types = "cfffff") %>%
  mutate(dayexpcombo = interaction(days, pseudomonas_hist, predation)) %>%
  mutate(expcombo = interaction(pseudomonas_hist, predation)) %>%
  mutate(
    predation = factor(predation, levels = c("no", "yes")),
    expcombo = factor(
      expcombo,
      levels = c(
        "ancestral.no",
        "coevolved.no",
        "ancestral.yes",
        "coevolved.yes"
      )
    )
  )

nmds2plot <- data.frame(nmds_ord$points) %>%
  rownames_to_column(var = "sample") %>%
  left_join(., metadf) %>%
  mutate(method = "NMDS") %>%
  dplyr::rename(axis_1 = MDS1, axis_2 = MDS2)

centroids <- nmds2plot %>%
  group_by(expcombo, days) %>%
  dplyr::mutate(
    NMDS1 = mean(axis_1),
    NMDS2 = mean(axis_2),
    NMDS1sd = sd(axis_1),
    NMDS2sd = sd(axis_2)
  ) %>%
  ungroup() 

spiders <- nmds2plot %>%
  select(sample, fromx = axis_1, fromy = axis_2) %>%
  left_join(., centroids) %>%
  dplyr::rename(tox = NMDS1, toy = NMDS2)

p_nmds <- ggplot() +
  geom_point(data = nmds2plot, aes(
    x = axis_1,
    y = axis_2,
    shape = days,
    color = expcombo
  )) +
  geom_segment(data = spiders,
               aes(
                 x = fromx,
                 xend = tox,
                 y = fromy,
                 yend = toy,
                 group = dayexpcombo,
                 color = expcombo
               )) +
  geom_point(data = centroids,
             aes(
               x = NMDS1,
               y = NMDS2,
               fill = expcombo,
               shape = days
             ),
             size = 2) +
  coord_fixed() +
  scale_shape_manual(values = c(21, 22, 24)) +
  scale_color_brewer(type = "qual", palette = "Paired") +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none"
  )

#p_nmds

# Fig. 3C -----------------------------------------------------------------
# beta binomial model coefficients
library(corncob)
modInteraction <- read_rds(here::here("data", "16S_amplicon", "corncob_interaction.rds"))
modPredation   <- read_rds(here::here("data", "16S_amplicon", "corncob_predation.rds"))
modEvolution   <- read_rds(here::here("data", "16S_amplicon", "corncob_evolution.rds"))

predPval <- tibble(strainID=names(modPredation$p_fdr[modPredation$p_fdr <= 0.1]),
       p_fdr=modPredation$p_fdr[modPredation$p_fdr <= 0.1]) %>% drop_na

predCoef <- plot(modPredation, data_only = T) %>%
  mutate(effect = "predation") %>%
  filter(variable == "predationyes\nDifferential Abundance") %>%
  separate(taxa, into = c("genus", "species")) %>%
  left_join(., tax) %>%
  left_join(., predPval)

evolPval <- tibble(strainID = names(modEvolution$p_fdr[modEvolution$p_fdr <= 0.1]),
                  p_fdr = modEvolution$p_fdr[modEvolution$p_fdr <= 0.1]) %>% drop_na 

evolCoef <- plot(modEvolution, data_only = T) %>% 
  mutate(effect = "evolution") %>%
  filter(variable == "pseudomonas_histcoevolved\nDifferential Abundance") %>%
  separate(taxa, into = c("genus", "species")) %>%
  left_join(., tax) %>%
  left_join(., evolPval)

interactionPval <- tibble(strainID=names(modInteraction$p_fdr[modInteraction$p_fdr <= 0.1]),
                  p_fdr=modInteraction$p_fdr[modInteraction$p_fdr <= 0.1]) %>% drop_na

interactionCoef <- plot(modInteraction, data_only = T) %>% 
  mutate(effect = "interaction") %>%
  filter(variable == "pseudomonas_histcoevolved:predationyes\nDifferential Abundance") %>%
  separate(taxa, into = c("genus", "species")) %>%
  left_join(., tax) %>%
  left_join(., interactionPval)

cc.coeffs <- bind_rows(predCoef, evolCoef, interactionCoef) %>%
  mutate(effect=factor(effect, levels=c("predation", "evolution", "interaction")))

p_cccf <- 
  ggplot(cc.coeffs, aes(
    x = x,
    xmin = xmin,
    xmax = xmax,
    y = strainID
  )) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    size = 1,
    color = "grey80"
  ) +
  geom_point(size = 1) +
  geom_linerange() +
  facet_grid( ~ effect) +
  theme_bw() +
  theme(legend.position = "none")

#p_cccf

# Fig. 3D -----------------------------------------------------------------
# SBW25 abundance

pf <- read_rds(here::here("data", "16S_amplicon", "corncob_SBW25.rds")) 

pf_qm <- pf[[2]]
pf_dat <- pf[[1]]$dat

p_sbw25 <- ggplot(pf_qm) + 
  geom_jitter(data=pf_dat, aes(x=predation, color=expcombo, y=W/M), size=0.5,
              position=position_dodge2(width=0.5)) +
  geom_linerange(aes(x=predation, color=expcombo, group=expcombo, 
                     y=q50, ymin=q2.5, ymax=q97.5),
                 position=position_dodge2(width=0.5)) +
  geom_point(aes(x=predation, fill=expcombo, group=expcombo, y=q50, shape=days),
             position=position_dodge2(width=0.5)) +
  scale_color_brewer(type="qual", palette="Paired") +
  scale_fill_brewer(type="qual", palette="Paired") +
  scale_shape_manual(values = c(21, 22, 24)) +
  facet_wrap(~days) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

#p_sbw25

# Fig. 3 ------------------------------------------------------------------
# Make final figure
library(patchwork)

pfinal <- (p_shannon / p_sbw25 / p_cccf) | p_nmds + plot_layout(guides = "collect")
pfinal

ggsave(plot=pfinal, here("figs", "fig3.svg"), width=17.8, height=14, units="cm")