library(here)
library(tidyverse)

# Read data ---------------------------------------------------------------

# counts
counts_formatted <- read_rds(here::here("data", "rnaseq", "counts_formatted.rds"))

# Summarize ---------------------------------------------------------------

# how much of total is coding RNA (mRNA) excluding Tetrahymena
summary_stats <- counts_formatted %>%
  filter(!str_detect(genome, "Tthermo")) %>%
  group_by(sample, rnatype) %>%
  summarize(s=sum(mappedreads)) %>%
  pivot_wider(id_cols=c("sample"), names_from = rnatype, values_from=s) %>%
  mutate(percent_mrna=coding/(coding+noncoding)*100) %>%
  ungroup()

summary_stats %>%
  summarize(medf=median(percent_mrna), 
            minf=min(percent_mrna),
            maxf=max(percent_mrna))

# percentage coding RNA
# medf  minf  maxf
# 2.33  0.82  5.06

# how much coding RNA belongs to each species.
sp_comp <- counts_formatted %>%
  filter(rnatype=="coding") %>%
  group_by(sample, genome) %>%
  summarize(s=sum(mappedreads)) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(total=sum(s)) %>%
  mutate(f=round((s/total)*100, digits=2)) %>%
  ungroup()

# how much RNA is tetrahymena
sp_comp %>%
  filter(str_detect(genome, "Tthermo")) %>%
  filter(str_detect(sample, "04.Y.")) %>%
  summarize(medf=median(f), 
            minf=min(f),
            maxf=max(f))

# percentage tetrahymena coding RNA in predation samples
# medf  minf  maxf
# 28.0  23.0  84.2


tmp <- counts_formatted %>%
  filter(str_detect(genome, "Tthermo")) %>%
  filter(str_detect(sample, "04.Y.")) %>%
  group_by(Geneid) %>%
  summarize(mmmapped=mean(mappedreads, na.rm=T)) %>%
  mutate(f=mmmapped/sum(mmmapped)*100)
  