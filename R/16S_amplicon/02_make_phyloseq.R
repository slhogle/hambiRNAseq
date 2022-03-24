library(here)
library(tidyverse)
library(phyloseq)

# Read data ---------------------------------------------------------------

counts_metadata <- read_tsv(here("data", "sample_metadata.tsv"), col_types = "cfffff")
counts <- read_tsv(here("data", "16S_amplicon", "species_counts.tsv"), col_types = "ccccd")

# make otu matrix
otumat <- left_join(counts, counts_metadata) %>%
  dplyr::select(sample, strainID, count) %>%
  arrange(sample) %>%
  pivot_wider(names_from=sample, values_from=count) %>%
  column_to_rownames(var="strainID") %>%
  as.matrix()

taxmat <- counts %>% 
  dplyr::select(strainID, Genus=genus, Species=species) %>%
  distinct(strainID, Genus, Species) %>%
  column_to_rownames(var="strainID") %>%
  as.matrix()

metadf <- counts_metadata %>%
  mutate(
    pseudomonas_hist = factor(pseudomonas_hist, levels = c("ancestral", "coevolved")),
    predation = factor(predation, levels = c("no", "yes")),
    days = factor(days, levels = c("4", "41", "45"))
  ) %>%
  mutate(dayexpcombo = interaction(days, pseudomonas_hist, predation)) %>%
  mutate(expcombo = interaction(pseudomonas_hist, predation)) %>%
  arrange(sample) %>%
  column_to_rownames(var = "sample")

physeq <- phyloseq(otu_table(otumat, taxa_are_rows = TRUE), 
                         tax_table(taxmat),
                         sample_data(metadf))

write_rds(physeq, here::here("data", "16S_amplicon", "physeq.rds"))
