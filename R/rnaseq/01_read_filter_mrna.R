library(here)
library(tidyverse)

# Read data ---------------------------------------------------------------

filelist <- list.files(here::here("dataRaw", "rnaseq", "counts"), full.names = TRUE)
files <- set_names(filelist, str_extract(filelist, regex("(?<=[/])([^/]+)(?=\\.[^.]+)")))

counts <- map_df(files, read_tsv, 
                 comment = "#", 
                 col_names = TRUE,
                 col_types = "ccddcdd",
                 .id = "sample")

# Read ids of all noncoding RNAs including 16S, 23S, 5S, 18S rRNAs and tRNAs
ncrnas <- read_tsv(here("dataRaw", "rnaseq", "hambi_noncoding_rna_ids.txt"), 
                   col_names="Geneid",
                   col_types="c") %>%
  mutate(rnatype="noncoding")


# rename and filter -------------------------------------------------------

# exclude all genes that are non-coding rnas

# the The A3_CO3_d4 library size is <5% of the other libraries 
# and there are other weird things about it in the quality control fastp files
# so we'll filter that sample out here as well.

# The Tetrahymena thermophila assembly is a transcriptome assembly - GCF_000189635.1
# according to here - https://www.ncbi.nlm.nih.gov/data-hub/gene/table/taxon/312017/
# all these genes are protein coding

counts_formatted <- left_join(counts, ncrnas, by="Geneid") %>%
  select(sample, Geneid, mappedreads, rnatype) %>%
  separate(Geneid, c("genome", "number"), sep = "_", remove=FALSE) %>%
  mutate(genome=ifelse(str_detect(genome, "gene"), "Tthermophila", genome)) %>%
  mutate(rnatype=if_else(genome=="Tthermophila", "coding", rnatype)) %>%
  mutate(rnatype=if_else(is.na(rnatype), "coding", rnatype)) %>%
  select(-number) %>%
  filter(sample != "04AYC") %>%
  arrange(genome, sample, Geneid) 

# the The A3_CO3_d4 library size is <5% of the other libraries 
# and there are other weird things about it in the quality control fastp files
# so we'll filter that sample out here as well.

write_rds(counts_formatted, here::here("data", "rnaseq", "counts_formatted.rds"))
