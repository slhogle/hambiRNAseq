library(here)
library(tidyverse)
library(clusterProfiler)
library(withr)
library(colorspace)

source(here::here("R", "rnaseq", "rnaseq_funs.R"))

# Load data ---------------------------------------------------------------

# deseq input data
load(file = here::here("data", "rnaseq", "deseqinputs.RData"))

# compiled deseq results
compiled_results <-
  read_rds(here::here("data", "rnaseq", "compiled_deseq_results.rds")) %>%
  mutate(genome = str_replace(genome, "HAMBI", "HAMBI-"))

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



# Calculate enrichments for individual genomes ----------------------------

# convert gene ID to genome
keggpath2gene <- eggnog %>%
  select(PATH = KEGG.Pathway, Geneid) %>%
  separate_rows(PATH) %>%
  filter(!str_detect(PATH, "^ko")) %>%
  mutate(dum = "path:") %>%
  unite(PATH, dum, PATH, sep = '', remove = TRUE) %>%
  mutate(PATH = ifelse(PATH == "path:NA", NA, PATH))

keggpath2name <- pathdesc

# run hypergeometric test
# Run for each day with and without predation

with_seed(12378,
  hyper_tests <- compiled_results %>%
    mutate(padj = if_else(is.na(padj), 0.99, padj)) %>%
    filter(padj <= 0.1 & abs(log2FoldChange) > 2) %>%
    group_by(genome, days, predation) %>%
    group_split() %>%
    map_dfr(., enrichfun, keggpath2gene, keggpath2name)
)

# this could potentially be a supplementary table?
# supptable <- hyper_tests %>%
#   select(-Geneid) %>%
#   distinct() %>%
#   arrange(genome, predation, days) %>%
#   group_by(Description) %>%
#   mutate(total=sum(Count))


# Plotting ----------------------------------------------------------------

hyper_tests_filt <- left_join(hyper_tests, eggnog, by=c("Geneid")) %>%
  group_by(Geneid, days, predation) %>%
  filter(enrichFactor == max(enrichFactor)) %>%
  group_by(genome) %>%
  add_count(genome) %>%
  ungroup() %>% 
  filter(n > 1) %>%
  separate(KEGG.ko, c("ko1", "ko2", "ko3", "ko4", "ko5", "ko6", "ko7", "ko8", "ko9"), sep = ",", remove=FALSE)

# get pathways to plot
pathways <- hyper_tests_filt %>%
  select(ko1, Description) %>%
  distinct(ko1, Description) %>%
  group_by(ko1) %>%
  add_count() %>%
  filter(n == 1)

# these KOs have multiple descriptions/pathways. Reduce to the most common ones
pathways_add <- tibble(
    ko1 = c(
      "ko:K00244",
      "ko:K00955",
      "ko:K00956",
      "ko:K01738",
      "ko:K02035",
      "ko:K03119"
    ),
    Description = c(
      "Two-component system",
      "Sulfur metabolism",
      "Sulfur metabolism",
      "Biosynthesis of antibiotics",
      "Quorum sensing",
      "Taurine and hypotaurine metabolism"
  )
)
  
# order KO ids for plotting
ko_levels <- bind_rows(pathways, pathways_add) %>%
  select(-n) %>%
  group_by(Description) %>%
  add_count(Description) %>%
  ungroup() %>%
  arrange(desc(n), Description, ko1) %>% 
  pull(ko1)  

# plot the pathway descriptions
ppath <- bind_rows(pathways, pathways_add) %>%
  mutate(ko1=factor(ko1, levels=ko_levels)) %>%
  ggplot() +
  geom_tile(aes(x=ko1, y="path", fill=Description)) +
  theme_bw() +
  mytheme(axis.text.x = element_text(angle = 90, hjust = 1, size=6),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.text = element_text(size=7),
          legend.key.size = unit(0.5, 'cm'))

ppath

ggsave(here("figs", "ppath2.svg"), ppath, device="svg", width=14, height=5, units="in")

genomeorderaxis <- rev(c("PsFluSBW25", "HAMBI-1287", "HAMBI-1992", "HAMBI-1972", "HAMBI-0403",
                         "HAMBI-2659", "HAMBI-2160", "HAMBI-0105", "HAMBI-2443", "HAMBI-3031", "HAMBI-1299"))

pheat <- left_join(hyper_tests_filt, compiled_results, by = c("Geneid", "genome", "predation", "days")) %>%
  complete(ko1, genome, days, predation) %>%
  mutate(ko1=factor(ko1, levels=ko_levels)) %>%
  mutate(genome = factor(genome, levels=genomeorderaxis)) %>%
  ggplot() +
  geom_tile(aes(x=ko1, y=genome, fill=log2FoldChange)) +
  facet_wrap(days ~ predation, ncol=1) +
  scale_fill_continuous_diverging(palette = "Purple-Green", na.value = "grey95") +
  theme_bw() +
  mytheme(axis.text.x = element_text(angle = 90, hjust = 1, size=4))

pheat

ggsave(here("figs", "pheat.svg"), pheat, device="svg", width=7, height=7, units="in")
