library(here)
library(tidyverse)
library(phyloseq)
library(corncob)
library(withr)

source(here::here("R", "16S_amplicon", "amplicon_funs.R"))

# Read data ---------------------------------------------------------------

physeq <- read_rds(here::here("data", "16S_amplicon", "physeq.rds"))


# Subset ------------------------------------------------------------------


p4y <- subset_samples(physeq,  predation == "yes")

with_seed(1237853, p4ymod <- bbdml(data = p4y,
                                  formula = PsFluSBW25 ~ days + pseudomonas_hist, 
                                  phi.formula =        ~ days + pseudomonas_hist))

p4ymod
