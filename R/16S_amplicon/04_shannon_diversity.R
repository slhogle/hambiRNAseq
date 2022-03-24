library(here)
library(tidyverse)
library(phyloseq)
library(DivNet)
library(withr)

physeq <- read_rds(here::here("data", "16S_amplicon", "physeq.rds"))
metadf <- rownames_to_column(data.frame(sample_data(physeq)), var = "sample")

## plugin estimate of shannon diversity 
plugin_sha <- estimate_richness(physeq, measures = "Shannon") %>%
  rownames_to_column(var = "sample") %>%
  dplyr::rename(estimate = Shannon) %>%
  mutate(div = "Shannon", inference = "plugin") %>%
  left_join(., metadf)

write_rds(plugin_sha, here::here("data", "16S_amplicon", "shannon_plugin.rds"))

## divnet estimate of shannon diversity. Takes a long time
with_seed(1234,
  dv <- divnet(
    physeq,
    ncores = 1,
    tuning = "fast",
    X = model.matrix(~ dayexpcombo, data = metadf)
))

## make a summary dataframe
dvnt_sha <- summary(dv$shannon) %>%
  dplyr::select(sample = sample_names, estimate, lower, upper) %>%
  mutate(div = "Shannon", inference = "divnet") %>%
  left_join(., metadf) %>%
  select(-sample, -replicate) %>%
  distinct()

write_rds(dvnt_sha, here::here("data", "16S_amplicon", "shannon_divnet.rds"))

# Plugin and divnet estimates visually agree
divcheck  <- ggplot() +
  geom_pointrange(
    data = dvnt_sha,
    aes(
      x = days,
      y = estimate,
      ymin = lower,
      ymax = upper,
      color = expcombo,
      group = expcombo
    ),
    fatten = 0.5,
    position = position_dodge2(width = 0.1)
  ) +
  geom_point(
    data = plugin_sha,
    aes(
      x = days,
      y = estimate,
      color = expcombo,
      group = expcombo
    ),
    position = position_dodge2(width = 0.1)
  ) +
  facet_grid(~ inference) +
  scale_fill_viridis_c(option = "B", name = NULL) +
  labs(y = "Shannon Diversity", x = "Day") +
  theme_bw() +
  theme(legend.position = "none")

divcheck

## differences in shannon diversity across day and treatment using betta
with_seed(123784,
          ba <- breakaway::betta(
            sapply(dv[["shannon"]], function(x)
              x$estimate),
            sapply(dv[["shannon"]], function(x)
              x$error),
            X = model.matrix(~ days + pseudomonas_hist * predation, data = metadf)
))

print(ba)

# ba$table %>%
#   xtable::xtable()