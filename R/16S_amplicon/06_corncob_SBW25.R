library(here)
library(tidyverse)
library(phyloseq)
library(corncob)
library(withr)

source(here::here("R", "16S_amplicon", "amplicon_funs.R"))

# Read data ---------------------------------------------------------------

physeq <- read_rds(here::here("data", "16S_amplicon", "physeq.rds"))

# Fit a model just for SBW25
with_seed(1237853, pfmod <- bbdml(data = physeq,
                                  formula = PsFluSBW25 ~ days + pseudomonas_hist * predation, 
                                  phi.formula =        ~ days + pseudomonas_hist * predation))

with_seed(123894, pf_q <- grouped_bbdml_quantiles95(pfmod, B=1000))

pf_qm <- pf_q %>%
  pivot_wider(values_from=x, names_from=quantile) %>%
  mutate(expcombo=interaction(pseudomonas_hist, predation))

write_rds(list(pfmod, pf_qm), file=here::here("data", "16S_amplicon", "corncob_SBW25.rds"))
