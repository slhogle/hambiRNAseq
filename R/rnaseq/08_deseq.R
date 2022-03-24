library(here)
library(tidyverse)
library(withr)
library(DESeq2)

source(here::here("R", "rnaseq", "rnaseq_funs.R"))

# Load data ---------------------------------------------------------------
load(file=here::here("data", "rnaseq", "deseqinputs.RData"))

# DESeq2 overall effects --------------------------------------------------

# Here we use the design formula ~species_amplicon_proportion + pseudomonas_hist
# + predation. There is no interaction term here so we only look for the overall
# effects of predation controlling for differences due to SBW25 coevolution or
# the overall effects of SBW25 coevolution controlling for predation. We are not
# looking specifically for genes whose predation effect is not consistent across
# SBW coevolution treatments or vice versa.

with_seed(125234,
          dds_04_overall <-
            imap(
              counts_formatted_filt_04_split_w[genome_order_04],
              deseqifySTC,
              meta_04,
              rnatotals_04,
              overall = TRUE
            )
)

write_rds(dds_04_overall, here::here("data", "rnaseq", "dds_04_overall.rds"))

with_seed(125234,
          dds_45_overall <-
            imap(
              counts_formatted_filt_45_split_w[genome_order_45],
              deseqifySTC,
              meta_45,
              rnatotals_45,
              overall = TRUE
            )
)

write_rds(dds_45_overall, here::here("data", "rnaseq", "dds_45_overall.rds"))


# DESeq2 coevolution-specific predation effects ---------------------------

# Instead of using interaction terms (e.g. ~species_amplicon_proportion +
# pseudomonas_hist + predation + pseudomonas_hist:predation) this approach
# simply combines the predation and coevolution factors into a single summary
# factor. Then differences between specific levels in the summary factor can be
# examined later using contrasts. Generally, this is approach is easiest. For example,
# recommended here:
# https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions

# use seed for reproducibility
with_seed(125234,
  dds_04_specific <-
    imap(
      counts_formatted_filt_04_split_w[genome_order_04],
      deseqifySTC,
      meta_04,
      rnatotals_04,
      overall = FALSE
    )
)

write_rds(dds_04_specific, here::here("data", "rnaseq", "dds_04_specific.rds"))

with_seed(125234,
  dds_45_specific <-
    imap(
      counts_formatted_filt_45_split_w[genome_order_45],
      deseqifySTC,
      meta_45,
      rnatotals_45,
      overall = FALSE
    )
)

write_rds(dds_45_specific, here::here("data", "rnaseq", "dds_45_specific.rds"))


# tmp <- results(dds_04_frm$HAMBI1287,
#         name = "predation_yes_vs_no") %>%
#   as.data.frame() %>%
#   rownames_to_column(var="Geneid") %>%
#   tibble() %>%
#   filter(padj <= 0.1 & abs(log2FoldChange) >= 2) %>%
#   arrange(Geneid)
# 
# resultsNames(dds_04_frm$PsFluSBW25)
# 
# ?counts
# 
# plotCounts(dds_04_frm$HAMBI1287, gene="HAMBI-1287_04355", intgroup="predation", 
#            normalized=TRUE, returnData = TRUE) 
