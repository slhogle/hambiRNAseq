library(here)
library(tidyverse)
source(here::here("R", "genericUtils.R"))

mapdir <- here::here("dataRaw", "16SAmplicon", "bbmapRPKM")

files <- set_names(list.files(mapdir, full.names = TRUE),
                   str_extract(
                     list.files(mapdir, full.names = TRUE),
                     regex("(?<=[/])([^/]+)(?=\\.[^.]+)")
                   ))

counts <- map_df(
  files,
  read_tsv,
  comment = "#",
  col_names = c(
    "strainID",
    "Length",
    "Bases",
    "Coverage",
    "count",
    "RPKM",
    "Frags",
    "FPKM"
  ),
  .id = "sample"
) %>%
  left_join(., tax) %>%
  select(sample, strainID, genus, species, count)

### accounting for multiple 16S copies in genomes.

#### accounting for multiple 16S copies in genomes.
## HAMBI-0216 - 3 copies
## HAMBI-1842 - 4 copies
## HAMBI-2659 - 2 copies
## JE2571-RP4 - 7 copies
## PsFluSBW25 - 5 copies

counts1 <- counts %>%
  mutate(
    count = case_when(
      # trunc ensures we get integers. important for some alpha div measures
      .$strainID == "HAMBI-0216" ~ trunc(.$count/3),
      .$strainID == "HAMBI-1842" ~ trunc(.$count/4),
      .$strainID == "HAMBI-2659" ~ trunc(.$count/2),
      .$strainID == "JE2571RP4"  ~ trunc(.$count/7),
      .$strainID == "PsFluSBW25" ~ trunc(.$count/5),
      TRUE ~ as.numeric(.$count)
    ))

counts1_wide <- counts1 %>% 
  group_by_at(vars(-count)) %>%
  mutate(row_id=1:n()) %>% ungroup() %>%  # build group index
  spread(key=sample, value=count) %>% # spread
  select(-row_id) %>% # drop the index 
  drop_na()

write_tsv(counts1, here::here("data", "speciesCounts.tsv"))
write_tsv(counts1_wide, here::here("data", "speciesCountsWide.tsv"))