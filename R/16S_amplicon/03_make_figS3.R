library(here)
library(tidyverse)
library(Polychrome)
library(withr)

source(here::here("R", "16S_amplicon", "amplicon_funs.R"))
# Read data ---------------------------------------------------------------

full <- left_join(read_tsv(here::here("data", "16S_amplicon", "species_counts.tsv"), col_types = "ccccd"),
                 read_tsv(here::here("data", "sample_metadata.tsv"), col_types = "cfffff"))

# get species in terms of frequencies
full_f <- full %>%
  group_by(sample) %>%
  mutate(count1 = sum(count)) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(freq = (count / count1)*100) %>% select(-count1) %>%
  ungroup()

# get all species above 1% frequency and pull them to a vector ordered by abundance
full_99 <- full_f %>%
  group_by(strainID) %>%
  mutate(sum = sum(count)) %>%
  distinct(strainID, sum) %>%
  ungroup() %>%
  mutate(sum1 = sum(sum)) %>%
  mutate(rank = sum / sum1 * 100) %>%
  filter(rank >= 0.1) %>%
  arrange(desc(rank)) %>%
  pull(strainID)

# make tibble with < 1% species collapsed
full_r <- full_f %>%
  mutate(genus = ifelse(strainID %in% full_99, genus, "Other")) %>%
  mutate(species = ifelse(strainID %in% full_99, species, "Other")) %>%
  mutate(strainID = ifelse(strainID %in% full_99, strainID, "Other")) %>%
  group_by(strainID, genus, species, replicate, days, summary_cat) %>%
  mutate(freq = sum(freq)) %>%
  ungroup() %>%
  distinct(freq,
           strainID,
           genus,
           replicate,
           days,
           summary_cat,
           .keep_all = TRUE) %>%
  mutate(strainID = factor(strainID, levels = c(full_99, "Other"))) %>%
  mutate(replicate = factor(
    replicate,
    levels = c("a", "b", "c"),
    labels = c("Rep A", "Rep B", "Rep C")
  )) %>%
  mutate(summary_cat = factor(
    summary_cat,
    levels = c("ANC_NONE", "COEVO_NONE", "ANC_COEVO", "COEVO_COEVO")
  ))

# create color palette
with_seed(4785,
          mypal <- unname(createPalette(18, c("#F3874AFF", "#FCD125FF"), M=500)))

mypal <- c("#D3D3D3", mypal)
names(mypal) <- c("Other", rev(full_99))
swatch(mypal)

f_freq <- ggplot(data = full_r) +
  geom_bar(aes(
    x = interaction(replicate),
    y = freq,
    fill = strainID
  ),
  stat = "identity") + #color="black",
  facet_grid(as.factor(days) ~ summary_cat) +
  scale_fill_manual(values = mypal) +
  labs(x = "Treatment", y = "Scaled relative abundance", fill = "Strain") +
  scale_y_continuous(limits = c(0,101), expand = c(0, 0)) +
  labs(x="Replicate", y="% abundance", fill="") + 
  theme_bw() + 
  mybartheme(legend.position="bottom")

f_freq

ggsave(
  plot = f_freq,
  here("figs", "figS3.svg"),
  width = 7,
  height = 7,
  units = "in"
)
