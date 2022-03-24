library(here)
library(tidyverse)
library(withr)
library(DESeq2)

source(here::here("R", "rnaseq", "rnaseq_funs.R"))

# Load data ---------------------------------------------------------------

dds_04_overall  <- read_rds(here::here("data", "rnaseq", "dds_04_overall.rds"))
dds_45_overall  <- read_rds(here::here("data", "rnaseq", "dds_45_overall.rds"))

# Normalized counts and Z-scores ------------------------------------------

with_seed(1521,
  dds_04_df_counts <-
    imap_dfr(dds_04_overall, get_deseq_counts) %>%
    mutate(days = 4)
)

with_seed(1521,
  dds_45_df_counts <-
    imap_dfr(dds_45_overall, get_deseq_counts) %>%
    mutate(days = 45)
)

compiled_results <- bind_rows(dds_04_df_counts, dds_45_df_counts)

write_rds(compiled_results,
          here::here("data", "rnaseq", "compiled_deseq_norm_counts.rds"))
