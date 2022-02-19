library(here)
library(tidyverse)
library(phyloseq)
library(corncob)
library(withr)

physeq = readRDS(here::here("data", "physeq.rds"))

# test for effect of predation term (in absence of interaction)
# days * pseudomonas_hist * predation
# days + pseudomonas_hist + predation + days:pseudomonas_hist + days:predation

# days + pseudomonas_hist * predation
# days + pseudomonas_hist + predation + pseudomonas_hist:predation

with_preserve_seed(
  modPredation <- differentialTest(
    formula =      ~ days + pseudomonas_hist + predation + pseudomonas_hist:predation,
    formula_null = ~ days + pseudomonas_hist,
    
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

saveRDS(modPredation, here::here("data", "corncobPredation.rds"))