library(here)
library(tidyverse)
library(phyloseq)
library(corncob)
library(withr)

physeq = readRDS(here::here("data", "physeq.rds"))

# Fit a model just for SBW25
with_preserve_seed(pfmod <- bbdml(data = physeq,
                                  formula = PsFluSBW25 ~ days + pseudomonas_hist * predation, 
                                  phi.formula =        ~ days + pseudomonas_hist * predation))

with_preserve_seed(pf.q <- grouped_bbdml_quantiles95(pfmod, B=1000))

pf.qm = pf.q %>%
  pivot_wider(values_from=x, names_from=quantile) %>%
  mutate(expcombo=interaction(pseudomonas_hist, predation))

saveRDS(list(pfmod, pf.qm), file=here::here("data", "corncobSBW25.rds"))
