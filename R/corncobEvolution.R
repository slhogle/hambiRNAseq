library(here)
library(tidyverse)
library(phyloseq)
library(corncob)
library(withr)

physeq = readRDS(here::here("data", "physeq.rds"))

with_preserve_seed(
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

saveRDS(modEvolution, here::here("data", "corncobEvolution.rds"))