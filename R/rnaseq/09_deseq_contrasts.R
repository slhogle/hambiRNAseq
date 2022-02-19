library(here)
library(tidyverse)
library(withr)
library(DESeq2)

source(here::here("R", "rnaseq", "rnaseq_funs.R"))

# Load data ---------------------------------------------------------------
load(file=here::here("data", "rnaseq", "deseqinputs.RData"))

# deseq results
dds_04 <- read_rds(here::here("data", "rnaseq", "dds_04.rds"))
dds_45 <- read_rds(here::here("data", "rnaseq", "dds_45.rds"))


# Deseq results with contrasts --------------------------------------------

# anc vs evo Pseuomonas in presence of predation
mycontrast_pred <- c("summary_cat", "COEVO_COEVO", "ANC_COEVO")

with_seed(1521,
  dds_04_df_pred <-
    imap_dfr(dds_04, get_deseq_results_shrunk, mycontrast_pred) %>%
    mutate(days = 4, predation = "yes")
)

with_seed(5672,
  dds_45_df_pred <-
    imap_dfr(dds_45, get_deseq_results_shrunk, mycontrast_pred) %>%
    mutate(days = 45, predation = "yes")
)

# anc vs evo Pseduomonas in absence of predation
mycontrast_nopred <- c("summary_cat", "COEVO_NONE", "ANC_NONE")

with_seed(452,
  dds_04_df_nopred <-
    imap_dfr(
      dds_04,
      get_deseq_results_shrunk,
      mycontrast_nopred,
      contrast2coef = TRUE
    ) %>%
    mutate(days = 4, predation = "no")
)

with_seed(56334,
  dds_45_df_nopred <-
    imap_dfr(
      dds_45,
      get_deseq_results_shrunk,
      mycontrast_nopred,
      contrast2coef = TRUE
    ) %>%
    mutate(days = 45, predation = "no")
)

# Combine results and save ------------------------------------------------

compiled_results <- bind_rows(dds_04_df_pred,
                              dds_45_df_pred,
                              dds_04_df_nopred,
                              dds_45_df_nopred)

write_rds(compiled_results,
          here::here("data", "rnaseq", "compiled_deseq_results.rds"))