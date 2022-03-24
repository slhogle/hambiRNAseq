library(here)
library(tidyverse)
library(withr)
library(DESeq2)

source(here::here("R", "rnaseq", "rnaseq_funs.R"))

# Load data ---------------------------------------------------------------
load(file=here::here("data", "rnaseq", "deseqinputs.RData"))

dds_04_overall  <- read_rds(here::here("data", "rnaseq", "dds_04_overall.rds"))
dds_45_overall  <- read_rds(here::here("data", "rnaseq", "dds_45_overall.rds"))

dds_04_specific  <- read_rds(here::here("data", "rnaseq", "dds_04_specific.rds"))
dds_45_specific  <- read_rds(here::here("data", "rnaseq", "dds_45_specific.rds"))


# DESeq2 contrasts overall effects ----------------------------------------

# Important to shrink the log fold change (effect size) so that we don't get
# lots of significant genes with very small counts and likely insignificant
# biological consequence


# overall effect of predation
mycontrast_overall_pred <- c("predation", "yes", "no")


with_seed(1521,
          dds_04_overall_pred <-
            imap_dfr(dds_04_overall, get_deseq_results_shrunk, mycontrast_overall_pred) %>%
            mutate(days = 4)
)

with_seed(5672,
          dds_45_overall_pred <-
            imap_dfr(dds_45_overall, get_deseq_results_shrunk, mycontrast_overall_pred) %>%
            mutate(days = 45)
)

# overall effect of SBW25 coevolution
mycontrast_overall_coevo <- c("pseudomonas_hist", "coevolved", "ancestral")

with_seed(1521,
          dds_04_overall_coevo <-
            imap_dfr(dds_04_overall, get_deseq_results_shrunk, mycontrast_overall_coevo) %>%
            mutate(days = 4)
)

with_seed(5672,
          dds_45_overall_coevo <-
            imap_dfr(dds_45_overall, get_deseq_results_shrunk, mycontrast_overall_coevo) %>%
            mutate(days = 45)
)


# DESeq2 specific effects -------------------------------------------------


# anc vs evo Pseudomonas in presence of predation
mycontrast_specific_pred <- c("summary_cat", "COEVO_COEVO", "ANC_COEVO")

with_seed(1521,
  dds_04_specific_pred <-
    imap_dfr(dds_04_specific, get_deseq_results_shrunk, mycontrast_specific_pred) %>%
    mutate(days = 4)
)

with_seed(5672,
  dds_45_specific_pred <-
    imap_dfr(dds_45_specific, get_deseq_results_shrunk, mycontrast_specific_pred) %>%
    mutate(days = 45)
)

# anc vs evo Pseduomonas in absence of predation
mycontrast_specific_coevo <- c("summary_cat", "COEVO_NONE", "ANC_NONE")

with_seed(452,
  dds_04_specific_coevo <-
    imap_dfr(
      dds_04_specific,
      get_deseq_results_shrunk,
      mycontrast_specific_coevo,
      contrast2coef = TRUE
    ) %>%
    mutate(days = 4)
)

with_seed(56334,
  dds_45_specific_coevo <-
    imap_dfr(
      dds_45_specific,
      get_deseq_results_shrunk,
      mycontrast_specific_coevo,
      contrast2coef = TRUE
    ) %>%
    mutate(days = 45)
)

# Combine results and save ------------------------------------------------

compiled_results <- bind_rows(
  dds_04_overall_pred,
  dds_45_overall_pred,
  dds_04_overall_coevo,
  dds_45_overall_coevo,
  dds_04_specific_pred,
  dds_45_specific_pred,
  dds_04_specific_coevo,
  dds_45_specific_coevo
)

write_rds(compiled_results,
          here::here("data", "rnaseq", "compiled_deseq_results.rds"))
