library(here)
library(tidyverse)
library(phyloseq)
library(DivNet)
library(withr)

physeq = readRDS(here::here("data", "physeq.rds"))
metadf = rownames_to_column(data.frame(sample_data(physeq)), var = "sample")

## plugin estimate of shannon diversity 
plugin.sha = estimate_richness(physeq, measures = "Shannon") %>%
  rownames_to_column(var = "sample") %>%
  dplyr::rename(estimate = Shannon) %>%
  mutate(div = "Shannon", inference = "plugin") %>%
  left_join(., metadf)

saveRDS(plugin.sha, here::here("data", "shannonPlugin.rds"))

## divnet estimate of shannon diversity. Takes a long time
with_preserve_seed(dv <- divnet(
  physeq,
  ncores = 1,
  tuning = "fast",
  X = model.matrix(~ dayexpcombo, data = metadf)
))

## make a summary dataframe
dvnt.sha <- summary(dv$shannon) %>%
  dplyr::select(sample = sample_names, estimate, lower, upper) %>%
  mutate(div = "Shannon", inference = "divnet") %>%
  left_join(., metadf) %>%
  select(-sample, -replicate) %>%
  distinct()

saveRDS(dvnt.sha, here::here("data", "shannonDivNet.rds"))

# Plugin and divnet estimates visually agree
divcheck = ggplot() +
  geom_pointrange(
    data = dvnt.sha,
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
    data = plugin.sha,
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
with_preserve_seed(ba <- breakaway::betta(
  sapply(dv[["shannon"]], function(x)
    x$estimate),
  sapply(dv[["shannon"]], function(x)
    x$error),
  X = model.matrix(~ days + pseudomonas_hist * predation, data = metadf)
))

print(ba)

# ba$table %>%
#   xtable::xtable()