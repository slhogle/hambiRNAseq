library(here)
library(tidyverse)
library(withr)
library(Polychrome)
library(patchwork)

source(here::here("R", "rnaseq", "rnaseq_funs.R"))

# Read data ---------------------------------------------------------------

# counts
counts_formatted <- read_rds(here::here("data", "rnaseq", "counts_formatted.rds"))

# metadata
metadata <- read_tsv(here::here("data", "amplicon_metadata.tsv"))

# number of noncoding rna per genome
ncrnas <- read_tsv(here::here("dataRaw", "rnaseq", "hambi_noncoding_rna_ids.txt"), col_names = "gene") %>%
  separate(gene, c("genome", "number"), sep = "_", remove=FALSE) %>%
  count(genome, name = "n_ncrna" ) %>%
  mutate(rnatype="noncoding")

hambi_anno <- read_tsv(here::here("dataRaw", "rnaseq", "hambi_annotations_prokka.tsv")) %>%
  filter(locus_tag!="locus_tag") %>%
  separate(locus_tag, c("genome", "number"), sep = "_", remove=FALSE) %>%
  group_by(genome) %>%
  count(ftype) %>%
  pivot_wider(id_cols=genome, names_from=ftype, values_from=n) %>%
  mutate(ncrna=sum(rRNA, tmRNA, tRNA, na.rm=T)) %>%
  select(genome, coding=CDS, noncoding=ncrna) %>%
  pivot_longer(cols=-genome, values_to = "count", names_to="rnatype") 
  
# Intermediate results for future plotting --------------------------------

# get 30 species order by relative abundance in total RNA
sp_order <- 
  counts_formatted %>%
  filter(!str_detect(genome, "Tthermo")) %>%
  group_by(genome) %>%
  summarize(s=sum(mappedreads)) %>%
  arrange(s) %>%
  pull(genome)

# get species proportions
# while normalizing for genome size

sp_comp <- counts_formatted %>%
  filter(!str_detect(genome, "Tthermo")) %>%
  group_by(sample, genome, rnatype) %>%
  summarize(s=sum(mappedreads)) %>%
  ungroup() %>%
  left_join(., hambi_anno) %>%
  mutate(s=s/count) %>%
  group_by(sample, rnatype) %>%
  mutate(total=sum(s)) %>%
  ungroup() %>%
  mutate(f=(s/total)*100) %>%
  left_join(., metadata) %>%
  mutate(summary_cat=factor(summary_cat, levels=c("ANC_NONE", "COEVO_NONE", "ANC_COEVO", "COEVO_COEVO")),
         genome=factor(genome, levels=sp_order))

# get species, lumping together all species < 1%
counts_formatted_rnatype <- counts_formatted %>%
  filter(!str_detect(genome, "Tthermo")) %>%
  group_split(rnatype)

names(counts_formatted_rnatype) <- c("coding", "noncoding")

sp_comp_lumped <- map(counts_formatted_rnatype, comp_lump) %>%
  bind_rows(.id="rnatype") %>%
  ungroup() %>%
  left_join(., metadata) %>%
  rename(mappedreads=n)

# how many genomes left
distinct(sp_comp_lumped, genome2) %>%
  nrow()

# order species by relative abundance
sp_order_lumped <- 
  sp_comp_lumped %>%
  group_by(genome) %>%
  summarize(s=sum(mappedreads)) %>%
  arrange(s) %>%
  filter(genome !="Other") %>%
  pull(genome)

# add "Other" category last
sp_order_lumped <- c(sp_order_lumped, "Other")

# set up the tibble for plotting
sp_comp_lumped_fct <- sp_comp_lumped %>%
  mutate(summary_cat=factor(summary_cat, levels=c("ANC_NONE", "COEVO_NONE", "ANC_COEVO", "COEVO_COEVO")),
         genome=factor(genome, levels=sp_order_lumped))

# create color palette
with_seed(4000,
          mypal <- unname(createPalette(16, c("#F3874AFF", "#FCD125FF"), M=5000)))

mypal <- c("#D3D3D3", mypal)
names(mypal) <- rev(sp_order_lumped)
swatch(mypal)

# Plotting ----------------------------------------------------------------
plot_04 <- sp_comp_lumped_fct %>%
  filter(days==4) %>%
  ggplot() +
  geom_bar(aes(y = f, x=replicate, fill = genome), 
           color="black", size=0.25, stat="identity") +
  facet_grid(rnatype ~ summary_cat, scales = "free_x") +
  scale_fill_manual(values = mypal) + 
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  labs(x="Replicate", y="% abundance", fill="") + 
  theme_bw() + 
  mybartheme(legend.position="bottom")

plot_45 <- sp_comp_lumped_fct %>%
  filter(days==45) %>%
  ggplot() +
  geom_bar(aes(y = f, x=replicate, fill = genome), 
           color="black", size=0.25, stat="identity") +
  facet_grid(rnatype ~ summary_cat, scales = "free_x") +
  scale_fill_manual(values = mypal) + 
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  labs(x="Replicate", y="% abundance", fill="") + 
  theme_bw() + 
  mybartheme(legend.position="bottom")

pfinal <- plot_04 + plot_45 + 
  plot_layout(nrow=1, guides="collect", ) + 
  plot_annotation(tag_levels="A") & theme(legend.position = 'bottom')

# Save --------------------------------------------------------------------

ggsave(here::here("figs", "figS4.svg"), pfinal, 
       width=7, height=5, units="in", device="svg")
