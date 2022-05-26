library(here)
library(tidyverse)
library(withr)

source(here::here("R", "rnaseq", "rnaseq_funs.R"))

# Load data ---------------------------------------------------------------
# counts
counts_formatted <-
  read_rds(here::here("data", "rnaseq", "counts_formatted.rds"))

# metadata
metadata <- read_tsv(here::here("data", "sample_metadata.tsv"))

# tetrahymena annotations
tet_annotation <-
  read_tsv(here::here("data_raw", "rnaseq", "tetrahymena_annotations.tsv.xz"))

# Tetrahymena: filter and format ------------------------------------------
counts_formatted_filt_tet <-
  left_join(counts_formatted, metadata) %>%
  filter(rnatype == "coding") %>%
  filter(genome == "Tthermophila")

# day 4
counts_formatted_filt_04_tet <- counts_formatted_filt_tet %>%
  filter(str_detect(sample, "...Y.")) %>%
  filter(days == 4) %>%
  select(Geneid, sample, mappedreads) %>%
  arrange(Geneid, sample) %>%
  pivot_wider(names_from = sample, values_from = mappedreads) %>%
  column_to_rownames(var = "Geneid") %>%
  as.matrix()

# day 45
counts_formatted_filt_45_tet <- counts_formatted_filt_tet %>%
  filter(str_detect(sample, "...Y.")) %>%
  filter(days == 45) %>%
  select(Geneid, sample, mappedreads) %>%
  arrange(Geneid, sample) %>%
  pivot_wider(names_from = sample, values_from = mappedreads) %>%
  column_to_rownames(var = "Geneid") %>%
  as.matrix()

meta_04 <- metadata %>%
  filter(days == 4, predation == "yes",!(sample %in% c("04AYC", "04CNC"))) %>%
  arrange(sample) %>%
  column_to_rownames(var = "sample") %>%
  as.data.frame()

meta_45 <- metadata %>%
  filter(days == 45, predation == "yes",!(sample %in% c("04AYC", "04CNC"))) %>%
  arrange(sample) %>%
  column_to_rownames(var = "sample") %>%
  as.data.frame()

# DeSeq2 differential expression analysis ---------------------------------
library(DESeq2)

# with this contrast genes upregulated with progenitor Pseudomonas have a positive
# log2 fold change

mycontrast <- c("summary_cat", "COEVO_COEVO", "ANC_COEVO")

with_seed(123784,
  diff_genes_04_tet <-
    deseqify(counts_formatted_filt_04_tet, meta_04, mycontrast)
)

diff_genes_04_tet_ann <-
  left_join(diff_genes_04_tet, tet_annotation)

write_rds(
  diff_genes_04_tet_ann,
  here::here("data", "rnaseq", "diff_genes_04_tet_annotated.rds")
)

with_seed(123784,
  diff_genes_45_tet <-
    deseqify(counts_formatted_filt_45_tet, meta_45, mycontrast)
)

diff_genes_45_tet_ann <-
  left_join(diff_genes_45_tet, tet_annotation)

# nothing differentially regulated at day 45
#write_rds(diff_genes_04_tet_ann, here::here("data", "rnaseq", "diff_genes_45_tet_annotated.rds"))

diff_genes_04_tet_ann %>%
  filter(abs(log2FoldChange) > 2) %>% nrow()

# Check for enriched genes ------------------------------------------------
library(clusterProfiler)

diff_genes_04_tet_ann <-
  read_rds(here::here("data", "rnaseq", "diff_genes_04_tet_annotated.rds"))

# genes that are upregulated with co-evolved Pseudomonas with fold change > 2
tetpos <- diff_genes_04_tet_ann %>%
  filter(log2FoldChange > 2) %>%
  pull(gene_symbol)

# genes that are downregulated with co-evolved Pseudomonas fold change < -2
tetneg <- diff_genes_04_tet_ann %>%
  filter(log2FoldChange < -2) %>%
  pull(gene_symbol)

with_seed(12367,
          kk_pos <- enrichKEGG(
            gene = tetpos,
            organism = 'tet',
            pvalueCutoff = 0.05
          )
        )

head(kk_pos)

with_seed(12367,
          kk_neg <- enrichKEGG(
            gene = tetneg,
            organism = 'tet',
            pvalueCutoff = 0.05
          )
        )

head(kk_neg)

# ID               Description GeneRatio BgRatio     pvalue   p.adjust     qvalue                          geneID Count
# tet03015 tet03015 mRNA surveillance pathway       2/6 53/1871 0.01098134 0.04602062 0.01614759 TTHERM_00295340/TTHERM_00051730     2
# tet03018 tet03018           RNA degradation       2/6 63/1871 0.01534021 0.04602062 0.01614759 TTHERM_00295340/TTHERM_00051730
