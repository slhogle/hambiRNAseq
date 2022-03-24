library(here)
library(tidyverse)
source(here::here("R", "rnaseq", "rnaseq_funs.R"))

# Load data ---------------------------------------------------------------

# counts
counts_formatted <- read_rds(here::here("data", "rnaseq", "counts_formatted.rds")) %>%
  mutate(genome=str_replace(genome, "-", ""))

# metadata
metadata <- read_tsv(here::here("data", "amplicon_metadata.tsv"))

# species relative abundances from amplicon
sp_abund <-  left_join(read_tsv(here("data", "amplicon", "speciesCounts.tsv"), col_types = "ccccd"),
                       read_tsv(here("data", "amplicon_metadata.tsv"), col_types = "cfffff")) %>%
  mutate(strainID=str_replace(strainID, "-", ""))

# Rank species ------------------------------------------------------------

# including only coding RNAs and exclude tetrahymena
counts_formatted_filt <- left_join(counts_formatted, metadata) %>% 
  filter(rnatype=="coding") %>%
  filter(genome != "Tthermophila")

# selecting which genomes to include for each day
# Perform the transformations for each day separately
counts_formatted_filt_lump <- counts_formatted_filt %>%
  group_by(days, genome) %>%
  summarize(n=sum(mappedreads)) %>%
  mutate(f=n/sum(n)*100) %>%
  mutate(genome2=fct_lump(genome, prop = 0.005, w = n)) %>%
  mutate(genome2=as.character(genome2))

genome_order_04 <- counts_formatted_filt_lump %>%
  filter(days==4, genome2 != "Other") %>%
  arrange(desc(f)) %>%
  pull(genome2)

genome_order_45 <- counts_formatted_filt_lump %>%
  filter(days==45, genome2 != "Other") %>%
  arrange(desc(f)) %>%
  pull(genome2)

genome_order_full <- counts_formatted_filt_lump %>%
  arrange(days, desc(f)) %>%
  filter(genome2 != "Other") %>%
  group_by(genome2) %>%
  filter(f==max(f)) %>%
  pull(genome2)


# Get amplicon abundances -------------------------------------------------

sp_abund.f <-  sp_abund %>%
  group_by(sample) %>%
  mutate(count1 = sum(count)) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(freq = count / count1) %>% select(-count1) %>%
  ungroup() %>%
  select(sample, strainID, days, freq)

sp_abund_04 <- sp_abund.f %>%
  filter(days == 4, !(sample %in% c("04AYC", "04CNC"))) %>%
  select(-days) %>%
  arrange(sample) %>%
  pivot_wider(names_from=strainID, values_from = freq)

sp_abund_45 <- sp_abund.f %>%
  filter(days == 45, !(sample %in% c("04AYC", "04CNC"))) %>%
  select(-days) %>%
  arrange(sample) %>%
  pivot_wider(names_from=strainID, values_from = freq)

# Species abundances RNA --------------------------------------------------

# Counts of coding and noncoding RNA from bacterial community
# Save for later incorporation in DESeq2 normalization factors
# using taxa sum-scaling

rnatotals_04 <- counts_formatted %>%
  filter(str_detect(sample, "04")) %>%
  rnadfformat

names(rnatotals_04) <- c("coding", "noncoding")

rnatotals_45 <- counts_formatted %>%
  filter(str_detect(sample, "45")) %>%
  rnadfformat

names(rnatotals_45) <- c("coding", "noncoding")

# Filter and format -------------------------------------------------------

# splitting df into list by genome
# day 4
counts_formatted_filt_04_split <- counts_formatted_filt %>%
  filter(days==4) %>%
  arrange(genome) %>%
  group_by(genome) %>%
  group_split()

names(counts_formatted_filt_04_split) <- pull(arrange(distinct(filter(counts_formatted_filt, days==4), genome)))

# day 45
counts_formatted_filt_45_split <- counts_formatted_filt %>%
  filter(days==45) %>%
  arrange(genome) %>%
  group_by(genome) %>%
  group_split()

names(counts_formatted_filt_45_split) <- pull(arrange(distinct(filter(counts_formatted_filt, days==45), genome)))

# deseq requires matrix. it indexs by column row "index" not name so need to
# make sure that these are always in the same order for every sample

counts_formatted_filt_04_split_w <- split_wide_mat(counts_formatted_filt_04_split)
counts_formatted_filt_45_split_w <- split_wide_mat(counts_formatted_filt_45_split)

# format metadata for deseq. Again because indexing is numeric the column/row order matters
meta_04 <- metadata %>%
  filter(days == 4, !(sample %in% c("04AYC", "04CNC"))) %>%
  left_join(., sp_abund_04) %>%
  arrange(sample) %>%
  mutate(
    days = factor(days, levels = c(4, 45)),
    pseudomonas_hist = factor(pseudomonas_hist, levels = c("ancestral", "coevolved")),
    predation = factor(predation, levels = c("no", "yes")),
    summary_cat = factor(summary_cat)
  ) %>% 
  column_to_rownames(var = "sample") %>%
  as.data.frame()

meta_45 <- metadata %>%
  filter(days == 45, !(sample %in% c("04AYC", "04CNC"))) %>%
  left_join(., sp_abund_45) %>%
  arrange(sample) %>%
  mutate(
    days = factor(days, levels = c(4, 45)),
    pseudomonas_hist = factor(pseudomonas_hist, levels = c("ancestral", "coevolved")),
    predation = factor(predation, levels = c("no", "yes")),
    summary_cat = factor(summary_cat)
  ) %>% 
  column_to_rownames(var = "sample") %>%
  as.data.frame()

# setup colors for subsequent plotting
mycols19 <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(genome_order_full)-1)
mycols19 <- append(mycols19, "#000000", after=1)
names(mycols19) <- genome_order_full

# SAVE --------------------------------------------------------------------

save(counts_formatted_filt_04_split_w, counts_formatted_filt_45_split_w,
     genome_order_04, genome_order_45, genome_order_full,
     meta_04, meta_45,
     rnatotals_04, rnatotals_45,
     mycols19, # colors for distatis and other plots
     file=here::here("data", "rnaseq", "deseqinputs.RData"))
