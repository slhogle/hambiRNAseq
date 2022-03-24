library(here)
library(tidyverse)
library(withr)
library(DESeq2)

source(here::here("R", "rnaseq", "rnaseq_funs.R"))

# Load data ---------------------------------------------------------------
load(file = here::here("data", "rnaseq", "deseqinputs.RData"))

# rlog normalized by amplicon relative abundance --------------------------
with_seed(125234,
  rlog_04_ASC <-
    imap(
      counts_formatted_filt_04_split_w[genome_order_04],
      rlogifyASC,
      meta_04,
      rnatotals_04
    )
)

write_rds(rlog_04_ASC, here::here("data", "rnaseq", "rlog_04_ASC.rds"))

with_seed(125234,
  rlog_45_ASC <-
    imap(
      counts_formatted_filt_45_split_w[genome_order_45],
      rlogifyASC,
      meta_45,
      rnatotals_45
    )
)

write_rds(rlog_45_ASC, here::here("data", "rnaseq", "rlog_45_ASC.rds"))


# rlog normalized by sum coding RNA per species ---------------------------

with_seed(125234,
  rlog_04_STC <-
    imap(
      counts_formatted_filt_04_split_w[genome_order_04],
      rlogifySTC,
      meta_04,
      rnatotals_04
    )
)

write_rds(rlog_04_STC, here::here("data", "rnaseq", "rlog_04_STC.rds"))

with_seed(125234,
  rlog_45_STC <-
    imap(
      counts_formatted_filt_45_split_w[genome_order_45],
      rlogifySTC,
      meta_45,
      rnatotals_45
    )
)

write_rds(rlog_45_STC, here::here("data", "rnaseq", "rlog_45_STC.rds"))
