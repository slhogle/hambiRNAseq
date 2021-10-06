library(here)
library(tidyverse)

full = left_join(read_tsv(here("data", "speciesCounts.tsv"), col_types = "ccccd"),
                 read_tsv(here("data", "ampliconMetadata.tsv"), col_types = "cfffff"))

# get species in terms of frequencies
full.f = full %>%
  group_by(sample) %>%
  mutate(count1 = sum(count)) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(freq = count / count1) %>% select(-count1) %>%
  ungroup()

# get all species above 1% frequency and pull them to a vector ordered by abundance
full.99 = full.f %>%
  group_by(strainID) %>%
  mutate(sum = sum(count)) %>%
  distinct(strainID, sum) %>%
  ungroup() %>%
  mutate(sum1 = sum(sum)) %>%
  mutate(rank = sum / sum1 * 100) %>%
  filter(rank >= 1) %>%
  arrange(desc(rank)) %>%
  pull(strainID)

# make tibble with < 1% species collapsed
full.r = full.f %>%
  mutate(genus = ifelse(strainID %in% full.99, genus, "Other")) %>%
  mutate(species = ifelse(strainID %in% full.99, species, "Other")) %>%
  mutate(strainID = ifelse(strainID %in% full.99, strainID, "Other")) %>%
  group_by(strainID, genus, species, replicate, days, summary_cat) %>%
  mutate(freq = sum(freq)) %>%
  ungroup() %>%
  distinct(freq,
           strainID,
           genus,
           species,
           replicate,
           days,
           summary_cat,
           .keep_all = TRUE) %>%
  mutate(strainID = factor(strainID, levels = c(full.99, "Other"))) %>%
  mutate(replicate = factor(
    replicate,
    levels = c("a", "b", "c"),
    labels = c("Rep A", "Rep B", "Rep C")
  )) %>%
  mutate(summary_cat = factor(
    summary_cat,
    levels = c(
      "ANC_ANC",
      "EVO_ANC",
      "COEVO_ANC",
      "ANC_COEVO",
      "EVO_COEVO",
      "COEVO_COEVO",
      "ANC_NONE",
      "EVO_NONE",
      "COEVO_NONE"
    ),
    labels = c(
      "Pseudo.anc + Tetra.anc",
      "Pseudo.evo + Tetra.anc",
      "Pseudo.coevo + Tetra.anc",
      "Pseudo.anc + Tetra.coevo",
      "Pseudo.evo + Tetra.coevo",
      "Pseudo.coevo + Tetra.coevo",
      "Pseudo.anc + No pred.",
      "Pseudo.evo + No pred.",
      "Pseudo.coevo + No pred."
    )
  ))

cols = colorRampPalette(RColorBrewer::brewer.pal(9, "Set3"))(11)
cols[8] = "#000000" # for Pseudomonas


f.freq = ggplot(data = full.r) +
  geom_bar(aes(
    x = interaction(replicate),
    y = freq,
    fill = strainID
  ),
  stat = "identity") + #color="black",
  facet_grid(as.factor(days) ~ summary_cat) +
  scale_fill_manual(values = cols) +
  labs(x = "Treatment", y = "Scaled relative abundance", fill = "Strain") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )

#f.freq

ggsave(
  plot = f.freq,
  here("figs", "figS2.svg"),
  width = 17.8,
  height = 12.8,
  units = "cm"
)
