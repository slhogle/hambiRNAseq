library(here)
library(tidyverse)
library(patchwork)
library(scales)

source(here::here("R", "rnaseq", "rnaseq_funs.R"))

# Read data ---------------------------------------------------------------
load(file = here::here("data", "rnaseq", "deseqinputs.RData"))

compiled_results <-
  read_rds(here::here("data", "rnaseq", "compiled_deseq_results.rds")) %>%
  mutate(
    days = factor(days, levels = c(4, 45)),
    predation = factor(predation, levels = c("no", "yes")),
    genome = factor(genome, levels = genome_order_full)
  ) %>%
  mutate(padj = if_else(is.na(padj), 0.99, padj)) %>%
  mutate(lfchi = if_else(padj < 0.1, "yes", "no")) %>%
  mutate(lfchi = if_else(abs(log2FoldChange) < 2, "no", "yes"))

# Volcano plot ------------------------------------------------------------

myvolcano <- function(.data) {
  ggplot(.data,) +
    geom_point(
      data = filter(.data, lfchi == "yes"),
      aes(x = baseMean, y = log2FoldChange, color = genome),
      size = 0.25,
      alpha = 0.7
    ) +
    geom_point(
      data = filter(.data, lfchi == "no"),
      aes(x = baseMean, y = log2FoldChange),
      size = 0.1,
      alpha = 0.01
    ) +
    labs(x = "Mean of normalized counts", y = "log2 fold change") +
    facet_wrap( ~ genome, nrow = 2) +
    scale_x_continuous(
      trans = "log10",
      breaks = c(0.1, 10, 1000),
      labels = trans_format('log10', math_format(10 ^ .x))
    ) +
    scale_color_manual(values = mycols19, na.value = "black") +
    theme_bw() +
    mytheme(legend.position = "none")
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
  here("figs", "figS5.png"),
  plot = volcano_wrap(compiled_results),
  device = "png",
  dpi = 300,
  width = 7,
  height = 10,
  units = "in"
)
