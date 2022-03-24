library(here)
library(tidyverse)
library(phyloseq)
library(corncob)
library(withr)

physeq <- readRDS(here::here("data", "16S_amplicon", "physeq.rds"))

with_seed(1237853,
  modEvolution <- differentialTest(
    formula =      ~ days + pseudomonas_hist + predation + pseudomonas_hist:predation,
    formula_null = ~ days + predation,
    
    # include this if interested in test differential dispersion
    phi.formula =      ~ days + pseudomonas_hist * predation,
    phi.formula_null = ~ days + pseudomonas_hist * predation,
    
    data = physeq,
    test = "LRT",
    boot = TRUE,
    B = 1000,
    fdr_cutoff = 0.1,
    full_output = TRUE
  )
)

write_rds(modEvolution, here::here("data", "corncob_evolution.rds"))