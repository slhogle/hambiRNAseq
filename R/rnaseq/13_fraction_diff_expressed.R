library(here)
library(tidyverse)

source(here::here("R", "rnaseq", "rnaseq_funs.R"))

# Load data ---------------------------------------------------------------

# deseq input data
load(file = here::here("data", "rnaseq", "deseqinputs.RData"))

# compiled deseq results
compiled_results <-
  read_rds(here::here("data", "rnaseq", "compiled_deseq_results.rds"))

# community gene annotations
comm_genes <-
  read_tsv(here::here("dataRaw", "rnaseq", "hambi_annotations_prokka.tsv.xz")) %>%
  mutate(locus_tag = str_replace(locus_tag, "-", ""))

# Percent of transcriptomes differentially regulated ----------------------

gene_counts <- comm_genes %>%
  separate(locus_tag,
           c("genome", "number"),
           sep = "_",
           remove = FALSE) %>%
  filter(ftype == "CDS") %>%
  group_by(genome) %>%
  count(genome, name = "gene_num") %>%
  ungroup()

exp_counts <- compiled_results %>%
  group_by(genome, days, predation) %>%
  filter(abs(log2FoldChange) > 2) %>%
  filter(padj < 0.1) %>%
  count(genome, name = "diff_exp") %>%
  ungroup() %>%
  complete(genome, days, predation, fill = list(diff_exp = 0))

left_join(exp_counts, gene_counts) %>%
group_by(predation, days) %>%
  mutate(f = diff_exp / gene_num * 100) %>%
  summarize(
    mn = mean(f, na.rm = T),
    min = min(f, na.rm = T),
    max = max(f, na.rm = T)
  )

exp_counts %>%
  group_by(predation, days) %>%
  filter(diff_exp > 0) %>%
  count()
