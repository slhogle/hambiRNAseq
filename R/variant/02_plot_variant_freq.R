library(here)
library(tidyverse)
library(patchwork)
library(scales)
library(rcartocolor)

source(here::here("R", "variant", "parallelism_utils.R"))


# Load data ---------------------------------------------------------------

pfanc <- read_tsv(here::here("data", "variant", "Pf_ANC.ASM922.final.tsv"))

pfE1 <- left_join(
  read_tsv(here::here("data", "variant", "Pf_E1.ASM922.final.tsv")), 
  read_delim(here::here("data", "variant", "Pf_E1.ASM922.hapfreq.tsv"), 
             col_names = c("CHROM", "FILTER", "POS", "ALLELENUM", "REFALLELE",
                           "ALTALLELE", "PHASEGROUP","GT", "ALLELEDEPTH", 
                           "ALLELEFREQ", "MAP_HAPFREQ"))) %>% 
  mutate(replicate="pfE1")

pfE2 <- left_join(
  read_tsv(here::here("data", "variant", "Pf_E2.ASM922.final.tsv")), 
  read_delim(here::here("data", "variant", "Pf_E2.ASM922.hapfreq.tsv"),
             col_names = c("CHROM", "FILTER", "POS", "ALLELENUM", "REFALLELE",
                           "ALTALLELE", "PHASEGROUP","GT", "ALLELEDEPTH", 
                           "ALLELEFREQ", "MAP_HAPFREQ"))) %>% 
  mutate(replicate="pfE2")

pfE3 <- left_join(
  read_tsv(here::here("data", "variant", "Pf_E3.ASM922.final.tsv")), 
  read_delim(here::here("data", "variant", "Pf_E3.ASM922.hapfreq.tsv"),
             col_names = c("CHROM", "FILTER", "POS", "ALLELENUM", "REFALLELE",
                           "ALTALLELE", "PHASEGROUP","GT", "ALLELEDEPTH", 
                           "ALLELEFREQ", "MAP_HAPFREQ"))) %>% 
  mutate(replicate="pfE3")

# write this as supplementary dataset
df_final <- bind_rows(pfE1, pfE2, pfE3) %>%
  # exclude calls present in ancestor strain
  anti_join(., dplyr::select(pfanc, CHROM, POS)) %>%
  mutate(SPECIES="Pseudomonas fluorescens SBW25") %>%
  separate(ALLELEFREQ, into=c("REF.FREQ", "ALT.FREQ"), sep=",") %>%
  separate(MAP_HAPFREQ, into=c("HAP1.FREQ", "HAP2.FREQ", "HAP3.FREQ"), sep=",") %>%
  select(SPECIES, REPLICATE=replicate, POS, REF, ALT, FILTER, ALLELENUM, GENOTYPE=GT, 
         REF.FREQ, ALT.FREQ, PHASEGROUP, HAP1.FREQ, HAP2.FREQ, HAP3.FREQ, FTYPE:PRODUCT) %>%
  arrange(POS) %>%
  mutate(ALT.FREQ=as.numeric(ALT.FREQ))

write_tsv(df_final, here::here("tables", "final_variants.tsv"))

# read significant parallel evolved genes
sig_genes <- read_tsv(here::here("tables", "significant_genes.tsv")) %>%
  rename(LOCUS_TAG=Locus)

annotate_genes <- semi_join(df_final, sig_genes) %>%
  left_join(., sig_genes)

p <- ggplot(df_final) + 
  geom_vline(data=annotate_genes, aes(xintercept=POS), linetype="dashed") + 
  geom_point(aes(x=POS, y=ALT.FREQ)) +
  geom_point(data=annotate_genes, aes(x=POS, y=ALT.FREQ, fill=LOCUS_TAG), size=4, shape=21) + 
  facet_grid(REPLICATE~.) +
  labs(x="SBW25 genome position", y="Alt. allele frequency") +
  scale_x_continuous(labels = label_number_si(accuracy = 1, unit = "bp", sep = NULL)) + 
  scale_y_continuous(labels = label_percent()) + 
  scale_fill_carto_d(palette="Vivid") + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom"
  )

p

ggsave(here::here("figs", "fig2_variants.svg"), p, device="svg", width=7, height=4, units="in")
