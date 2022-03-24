library(here)
library(tidyverse)
library(patchwork)
library(scales)
library(rcartocolor)
library(ggrepel)

source(here::here("R", "rnaseq", "rnaseq_funs.R"))

# Read data ---------------------------------------------------------------
load(file = here::here("data", "rnaseq", "deseqinputs.RData"))

compiled_results <-
  read_rds(here::here("data", "rnaseq", "compiled_deseq_results.rds")) %>%
  #mutate(genome = str_replace(genome, "HAMBI", "HAMBI-")) %>%
  filter(coefficient %in% c("summary_cat_COEVO_COEVO_vs_ANC_COEVO", "summary_cat_COEVO_NONE_vs_ANC_NONE")) %>%
  mutate(predation=if_else(coefficient == "summary_cat_COEVO_COEVO_vs_ANC_COEVO", "yes", "no")) %>%
  mutate(
    days = factor(days, levels = c(4, 45)),
    predation = factor(predation, levels = c("no", "yes")),
    genome = factor(genome, levels = genome_order_full)
  ) %>%
  mutate(padj = if_else(is.na(padj), 0.99, padj)) %>%
  mutate(lfchi = if_else(padj < 0.1, "yes", "no")) %>%
  mutate(lfchi = if_else(abs(log2FoldChange) < 2, "no", "yes"))

volcano_annotations <- read_rds(here::here("data", "rnaseq", "volcano_annotations.rds"))

compiled_results <-   left_join(compiled_results, volcano_annotations)

# Volcano plot ------------------------------------------------------------

myvolcano <- function(.data) {
  ggplot(.data,) +
    geom_point(
      data = filter(.data, lfchi == "yes"),
      aes(x = baseMean, y = log2FoldChange),
      shape = 16,
      size = 0.25,
      color="red",
      alpha = 0.5
    ) +
    geom_point(
      data = filter(.data, lfchi == "no"),
      aes(x = baseMean, y = log2FoldChange),
      size = 0.1,
      alpha = 0.01
    ) +
    geom_point(
      data = filter(filter(.data, lfchi == "yes"), !is.na(group)),
      aes(x = baseMean, y = log2FoldChange, color = group), 
      size = 1
    ) +
    # geom_text_repel(
    #   data = filter(filter(.data, lfchi == "yes"), !is.na(group)),
    #   aes(x = baseMean, y = log2FoldChange, label=genename, color=group),
    #   max.overlaps = Inf,
    #   min.segment.length = 0,
    #   size=1
    # ) +
    labs(x = "Mean of normalized counts", y = "log2 fold change") +
    facet_wrap( ~ genome, nrow = 2) +
    scale_x_continuous(
      trans = "log10",
      breaks = c(0.1, 10, 1000),
      labels = trans_format('log10', math_format(10 ^ .x))
    ) +
    #scale_color_manual(values = mycols19, na.value = "black") +
    #scale_color_carto_d(palette="Vivid") +
    scale_color_manual(values=c("#6a3d9a", "#cab2d6", "#33a02c", "#a6cee3", "#1f78b4")) +
    #scale_fill_viridis_d(begin=0.4, end=1) +
    theme_bw() +
    mytheme()
}

volcano_wrap <- function(.data) {
  .data %>%
    group_by(days, predation) %>%
    group_split() %>%
    map(., myvolcano) %>%
    wrap_plots(ncol = 1) +
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = "A") &
    theme(plot.margin = unit(c(0, 5, 0, 5), "pt"))
}

ggsave(
  here("figs", "figVolcano4.png"),
  plot = volcano_wrap(filter(compiled_results, days==4)),
  device = "png",
  dpi = 300,
  width = 8,
  height = 5,
  units = "in"
)
