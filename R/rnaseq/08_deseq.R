library(here)
library(tidyverse)
library(withr)
library(DESeq2)

source(here::here("R", "rnaseq", "rnaseq_funs.R"))

# Load data ---------------------------------------------------------------
load(file=here::here("data", "rnaseq", "deseqinputs.RData"))

# use seed for reproducibility
with_seed(125234,
  dds_04 <-
    imap(
      counts_formatted_filt_04_split_w[genome_order_04],
      deseqifySTC,
      meta_04,
      rnatotals_04
    )
)

write_rds(dds_04, here::here("data", "rnaseq", "dds_04.rds"))

with_seed(125234,
  dds_45 <-
    imap(
      counts_formatted_filt_45_split_w[genome_order_45],
      deseqifySTC,
      meta_45,
      rnatotals_45
    )
)

write_rds(dds_45, here::here("data", "rnaseq", "dds_45.rds"))
