library(here)
library(tidyverse)
library(withr)

# Read data ---------------------------------------------------------------

metadata <- read_tsv(here("data_raw", "atp_cell_density", "metadata.tsv"))

atp <-
  read_tsv(here("data_raw", "atp_cell_density", "ATP_concentration.tsv"),
           col_types = "cddddddddddd") %>%
  pivot_longer(-sample, names_to = "day", values_to = "value") %>%
  mutate(day = as.numeric(day),
         var = "atp")

# scale response var there is a detection limit for ciliates at about 1000 per
# mL. Any lower than that and noise takes over so you get lots of zeros and
# weird low values. Here we are just randomly sampling from normal distribution
# for all counts below 1000

with_seed(345761,
          cil <-
            read_tsv(here("data_raw", "atp_cell_density", "cilliate_counts.tsv"),
                     col_types = "cddddddddddd") %>%
            mutate(`1` = rnorm(18, mean = 5000, sd = 1000)) %>%
            pivot_longer(-sample, names_to = "day", values_to = "value") %>%
            rowwise() %>%
            mutate(value = ifelse(
              value < 1000, rnorm(1, mean = 1000, sd = 500), value
            )) %>%
            ungroup() %>%
            mutate(value = ifelse(value < 0, 300, value)) %>%
            mutate(day = as.numeric(day),
                   var = "cil"))

with_seed(123785,
          cfu <-
            read_tsv(here("data_raw", "atp_cell_density", "ppy_colony_counts.tsv"),
                     col_types = "ccdddddd") %>%
            dplyr::select(-platetype) %>%
            rowwise() %>%
            mutate(`0` = rnorm(1, mean = 2.9e+07 + 1e+07, sd = 1e6)) %>%
            ungroup() %>%
            pivot_longer(-sample, names_to = "day", values_to = "value") %>%
            mutate(day = as.numeric(day),
                   var = "cfu"))

opd <-
  read_tsv(here("data_raw", "atp_cell_density", "OD600.tsv"),
           col_types = "cddddddddddd") %>%
  pivot_longer(-sample, names_to = "day", values_to = "value") %>%
  mutate(day = as.numeric(day),
         var = "opd")

# Format data -------------------------------------------------------------

fmt <- bind_rows(opd, cfu, cil, atp) %>%
  left_join(., metadata, by = "sample") %>%
  filter(pseudomonas_hist %in% c("ancestral", "coevolved")) %>%
  filter(tetra_hist %in% c("coevolved", "no_tetra")) %>%
  mutate(
    pred = ifelse(str_detect(tetra_hist, "coevolved"), "predation", "none"),
    evo = ifelse(str_detect(pseudomonas_hist, "coevolved"), "evolution", "none"),
    pred_evo = ifelse(str_detect(tetra_hist, "coevolved") & str_detect(pseudomonas_hist, "coevolved"), "interaction", "none"),
    y=log10(value), 
    recoveryonset=day-43) %>%
  group_by(var) %>%
  group_split() %>%
  map(., make_id) %>%
  map(., factorize_my_vars)

names(fmt) <- map(fmt, select, var) %>%
  bind_rows() %>%
  distinct() %>% 
  pull()

fmt <- map(fmt, select, -var) %>%
  map(., as.data.frame)

# Save data ---------------------------------------------------------------

write_rds(fmt, here::here("data", "gaussian_process_regression", "formatted_for_cluster.rds"))
