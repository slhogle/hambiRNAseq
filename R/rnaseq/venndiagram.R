library(here)
library(tidyverse)
library(patchwork)
library(ggVennDiagram)

# Load data ---------------------------------------------------------------

# compiled deseq results
diff_exp <-
  read_rds(here::here("data", "rnaseq", "compiled_deseq_results.rds")) %>%
  mutate(genome = str_replace(genome, "HAMBI", "HAMBI-"))


# Format data -------------------------------------------------------------

gene_list_04 <- diff_exp %>%
  filter(days==4) %>%
  mutate(padj = if_else(is.na(padj), 0.99, padj)) %>%
  filter(padj < 0.1 & abs(log2FoldChange) > 2) %>%
  group_by(days, coefficient) %>%
  group_split() %>%
  map(., pull, Geneid)

gene_list_45 <- diff_exp %>%
  filter(days==45) %>%
  mutate(padj = if_else(is.na(padj), 0.99, padj)) %>%
  filter(padj < 0.1 & abs(log2FoldChange) > 2) %>%
  group_by(days, coefficient) %>%
  group_split() %>%
  map(., pull, Geneid)

# predation_yes_vs_no
# pseudomonas_hist_coevolved_vs_ancestral
# summary_cat_COEVO_COEVO_vs_ANC_COEVO
# summary_cat_COEVO_NONE_vs_ANC_NONE


# Plot --------------------------------------------------------------------

pd4 <- ggVennDiagram(gene_list_04, 
              category.names = c("pred","coevo","coevo:pred","coevo:nopred")) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_distiller(palette = "YlGn", direction = 1, limits = c(0, 1400))

pd45 <- ggVennDiagram(gene_list_45, 
                     category.names = c("pred","coevo","coevo:pred","coevo:nopred")) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_distiller(palette = "YlGn", direction = 1, limits = c(0, 1400))

pf <- pd4 + pd45 +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels="A")

pf
  
ggsave(
  here("figs", "venn.svg"),
  plot = pf,
  device = "svg",
  width = 10,
  height = 10,
  units = "in"
)
