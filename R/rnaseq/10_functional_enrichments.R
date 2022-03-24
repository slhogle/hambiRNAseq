library(here)
library(tidyverse)
library(clusterProfiler)
library(withr)
library(colorspace)
library(pheatmap)
library(Polychrome)
library(patchwork)

source(here::here("R", "rnaseq", "10_enrich_funs.R"))

# Load data ---------------------------------------------------------------

# deseq input data
load(file = here::here("data", "rnaseq", "deseqinputs.RData"))
names(mycols19) <- str_replace(names(mycols19), "HAMBI", "HAMBI-")
names(mycols19) <- str_replace(names(mycols19), "JE2571", "JE2571-")

# compiled deseq results
diff_exp <-
  read_rds(here::here("data", "rnaseq", "compiled_deseq_results.rds")) %>%
  mutate(genome = str_replace(genome, "HAMBI", "HAMBI-"))

#deseq normalized counts
gene_counts <-
  read_rds(here::here("data", "rnaseq", "compiled_deseq_norm_counts.rds"))

# Load kegg pathways and their descriptions
path2ko <- read_tsv(here("dataRaw", "rnaseq", "path2ko.tsv"), 
                    col_names = c("PATH", "KO"), 
                    col_types="cc")
pathdesc <- read_tsv(here("dataRaw", "rnaseq", "path_desc.tsv"), 
                     col_names = c("PATH", "DESC"), 
                     col_types="cc")

# EggNogMapper results for the HAMBI genomes
eggnog <- read_tsv(here("dataRaw", "rnaseq", "hambi_annotations_eggnog.tsv"),
                   col_names = c("Geneid",
                                 "seed.eggNOG.ortholog",
                                 "seed.ortholog.evalue",
                                 "seed.ortholog.score",
                                 "predicted.tax.group",
                                 "predicted.protein",
                                 "go.terms",
                                 "EC.number",
                                 "KEGG.ko",
                                 "KEGG.Pathway",
                                 "KEGG.Module",
                                 "KEGG.Reaction",
                                 "KEGG.rclass",
                                 "KEGG.BRITE",
                                 "KEGG.TC",
                                 "CAZy",
                                 "BiGG.rxn",
                                 "tax.scope", #eggNOG taxonomic level used for annotation
                                 "eggNOG.OGs",
                                 "bestOG", #deprecated, use smallest from eggnog OGs
                                 "COG.category",
                                 "description"),
                   col_types = "ccddcccccccccccccccccc")



# Formatting --------------------------------------------------------------

# convert gene ID to genome
keggpath2gene <- eggnog %>%
  select(PATH = KEGG.Pathway, Geneid) %>%
  separate_rows(PATH) %>%
  filter(!str_detect(PATH, "^ko")) %>%
  mutate(dum = "path:") %>%
  unite(PATH, dum, PATH, sep = '', remove = TRUE) %>%
  mutate(PATH = ifelse(PATH == "path:NA", NA, PATH))

keggpath2name <- pathdesc

# Get differentially expressed genes --------------------------------------

diff_exp_uniq <- bind_rows(get_condition_specific_genes(diff_exp, 4),
                           get_condition_specific_genes(diff_exp, 45))

with_seed(12378,
          hyper_tests <- diff_exp_uniq %>%
            group_by(genome, days, coefficient) %>%
            group_split() %>%
            map_dfr(., enrichfun, keggpath2gene, keggpath2name)
)

# with_seed(12378,
#           hyper_tests <- diff_exp %>%
#             mutate(padj = if_else(is.na(padj), 0.99, padj)) %>%
#             filter(padj < 0.1 & abs(log2FoldChange) > 2) %>%
#             group_by(genome, days, coefficient) %>%
#             group_split() %>%
#             map_dfr(., enrichfun, keggpath2gene, keggpath2name)
# )



d4order <- c("04ANA", "04ANB", "04ANC", "04CNA", "04CNB", "04AYA", "04AYB", "04CYA", "04CYB", "04CYC")
d45order <- c("45ANA", "45ANB", "45ANC", "45CNA", "45CNB", "45CNC", "45AYA", "45AYB", "45AYC", "45CYA", "45CYB", "45CYC")

# predation

pp4 <- customheatmap(4, "predation_yes_vs_no", diff_exp_uniq, hyper_tests,
              eggnog, gene_counts, meta_04, mycols19, d4order, -2, 2)

ggsave(here::here("figs", "pred_heat_04.svg"), plot = pp4, device = "svg",
       width = 7, height = 10, units = "in")

pp45 <- customheatmap(45, "predation_yes_vs_no", diff_exp_uniq, hyper_tests,
              eggnog, gene_counts, meta_45, mycols19, d45order, -2, 2)

ggsave(here::here("figs", "pred_heat_45.svg"), plot = pp45, device = "svg",
  width = 5.5, height = 10, units = "in")

# SBW25 coevolution
pc4 <- customheatmap(4, "pseudomonas_hist_coevolved_vs_ancestral", diff_exp_uniq, hyper_tests,
              eggnog, gene_counts, meta_04, mycols19, 
              c("04ANA", "04ANB", "04ANC",
                "04AYA", "04AYB", 
                "04CNA", "04CNB",
                "04CYA", "04CYB", "04CYC"), -2, 2)

ggsave(here::here("figs", "coevo_heat_04.svg"), plot = pc4, device = "svg",
       width = 6, height = 5.5, units = "in")

pc45 <- customheatmap(45, "pseudomonas_hist_coevolved_vs_ancestral", diff_exp_uniq, hyper_tests,
              eggnog, gene_counts, meta_45, mycols19, 
              c("45ANA", "45ANB", "45ANC", 
                "45AYA", "45AYB", "45AYC",
                "45CNA", "45CNB", "45CNC", 
                "45CYA", "45CYB", "45CYC"), -2, 2)

ggsave(here::here("figs", "coevo_heat_45.svg"), plot = pc45, device = "svg",
       width = 5, height = 10, units = "in")

# SBW25 coevolution no predation
psc4 <- customheatmap(4, "summary_cat_COEVO_NONE_vs_ANC_NONE", diff_exp_uniq, hyper_tests,
              eggnog, gene_counts, meta_04, mycols19, d4order, 2, -2)

psc45 <- customheatmap(45, "summary_cat_COEVO_NONE_vs_ANC_NONE", diff_exp_uniq, hyper_tests,
              eggnog, gene_counts, meta_45, mycols19, d45order, 2, -2)

# SBW25 coevolution with predation

psp4 <- customheatmap(4, "summary_cat_COEVO_COEVO_vs_ANC_COEVO", diff_exp_uniq, hyper_tests,
              eggnog, gene_counts, meta_04, mycols19, d4order,
              -2, 2)

ggsave(here::here("figs", "summary_cat_pred_heat_04.svg"), plot = psp4, device = "svg",
       width = 5.5, height = 10, units = "in")

psp45 <- customheatmap(45, "summary_cat_COEVO_COEVO_vs_ANC_COEVO", diff_exp_uniq, hyper_tests,
              eggnog, gene_counts, meta_45, mycols19, d45order, -2, 2)




# Predation Log Fold Change Plot ------------------------------------------

dexp     <- diff_subset(diff_exp_uniq, 4, "predation_yes_vs_no", genomecutoff=0)
genomes  <- pull(distinct(select(dexp, genome)), genome)
enriched <- hyper_tests_subset(hyper_tests, eggnog, 4, "predation_yes_vs_no", genomes)
top_fams <- enriched %>% count(Description) %>% slice_max(order_by=n, n=5, with_ties=FALSE) %>% pull(Description)
p4 <- plotlfc(dexp, enriched, top_fams)

dexp     <- diff_subset(diff_exp_uniq, 45, "predation_yes_vs_no", genomecutoff=0)
genomes  <- pull(distinct(select(dexp, genome)), genome)
enriched <- hyper_tests_subset(hyper_tests, eggnog, 45, "predation_yes_vs_no", genomes) %>%
  mutate(Description = if_else(Description == "Photosynthesis", "Oxidative phosphorylation", Description))
top_fams <- enriched %>% count(Description) %>% slice_max(order_by=n, n=5, with_ties=FALSE) %>% pull(Description)
p45 <- plotlfc(dexp, enriched, top_fams)

pf <- p4 / p45 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels="A")

pf
ggsave(here::here("figs", "pred_lfc.svg"), plot = pf, device = "svg",
       width = 4, height = 5.5, units = "in")

# Coevolution log fold change plot ----------------------------------------

dexp     <- diff_subset(diff_exp, 4, "pseudomonas_hist_coevolved_vs_ancestral", genomecutoff=0)
genomes  <- pull(distinct(select(dexp, genome)), genome)
enriched <- hyper_tests_subset(hyper_tests, eggnog, 4, "pseudomonas_hist_coevolved_vs_ancestral", genomes)
top_fams <- enriched %>% count(Description) %>% slice_max(order_by=n, n=5, with_ties=FALSE) %>% pull(Description)
p4 <- plotlfc(dexp, enriched, top_fams)

dexp     <- diff_subset(diff_exp, 45, "pseudomonas_hist_coevolved_vs_ancestral", genomecutoff=0)
genomes  <- pull(distinct(select(dexp, genome)), genome)
enriched <- hyper_tests_subset(hyper_tests, eggnog, 45, "pseudomonas_hist_coevolved_vs_ancestral", genomes) %>%
  mutate(Description = if_else(Description == "Photosynthesis", "Oxidative phosphorylation", Description))
top_fams <- enriched %>% count(Description) %>% slice_max(order_by=n, n=5, with_ties=FALSE) %>% pull(Description)
p45 <- plotlfc(dexp, enriched, top_fams)

pf <- p4 / p45 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels="A")

pf

ggsave(here::here("figs", "coevo_lfc.svg"), plot = pf, device = "svg",
       width = 4, height = 5.5, units = "in")


# summary_cat_COEVO_COEVO_vs_ANC_COEVO log fold change plot -----------------

dexp     <- diff_subset(diff_exp_uniq, 4, "summary_cat_COEVO_COEVO_vs_ANC_COEVO", genomecutoff=0)
genomes  <- pull(distinct(select(dexp, genome)), genome)
enriched <- hyper_tests_subset(hyper_tests, eggnog, 4, "summary_cat_COEVO_COEVO_vs_ANC_COEVO", genomes)
top_fams <- enriched %>% count(Description) %>% slice_max(order_by=n, n=5, with_ties=FALSE) %>% pull(Description)
p4 <- plotlfc(dexp, enriched, top_fams)

dexp     <- diff_subset(diff_exp_uniq, 45, "summary_cat_COEVO_COEVO_vs_ANC_COEVO", genomecutoff=0)
genomes  <- pull(distinct(select(dexp, genome)), genome)
enriched <- hyper_tests_subset(hyper_tests, eggnog, 45, "summary_cat_COEVO_COEVO_vs_ANC_COEVO", genomes) #%>%
  mutate(Description = if_else(Description == "Photosynthesis", "Oxidative phosphorylation", Description))
top_fams <- enriched %>% count(Description) %>% slice_max(order_by=n, n=5, with_ties=FALSE) %>% pull(Description)
p45 <- plotlfc(dexp, enriched, top_fams)

pf <- p4 / p45 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels="A")

pf

# summary_cat_COEVO_NONE_vs_ANC_NONE log fold change plot -----------------















# Coevo vs ancestral in presence of predation

diff_exp1 <- diff_exp %>%
  mutate(padj = if_else(is.na(padj), 0.99, padj)) %>%
  filter(padj < 0.1 & abs(log2FoldChange) > 2) %>%
  filter(days == 4 & coefficient == "summary_cat_COEVO_COEVO_vs_ANC_COEVO")


dexp     <- diff_subset(diff_exp, 4, "summary_cat_COEVO_COEVO_vs_ANC_COEVO", genomecutoff=0)
genomes  <- pull(distinct(select(dexp, genome)), genome)
enriched <- hyper_tests_subset(hyper_tests, eggnog, 4, "summary_cat_COEVO_COEVO_vs_ANC_COEVO", genomes)


fams <- c("ABC transporters", "Sulfur metabolism", "Quorum sensing",  "Photosynthesis", 
          "Taurine and hypotaurine metabolism", "Selenocompound metabolism")

plfc <- left_join(enriched, dexp) %>%
  filter(Description %in% fams) %>%
  mutate(Description=factor(Description, levels=fams)) %>%
  ggplot(aes(x=Description, y=log2FoldChange, color=genome, size=baseMean)) +
  geom_jitter(width=0.2, height=0.5) +
  geom_hline(yintercept=0) +
  scale_color_manual(values=mycols19) +
  labs(x="Log Fold Change", y="", size="Normalized mean") + 
  scale_size_continuous(breaks=c(10, 100, 1000, 10000), range=c(1,3)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) + 
  #scale_x_discrete(guide = guide_axis(angle = 90)) + #labels = abbreviate
  theme_bw() +
  mytheme()

plfc

ggsave(here::here("figs", "lfc_04.svg"), plot = plfc, device = "svg",
       width = 4, height = 5.5, units = "in")
