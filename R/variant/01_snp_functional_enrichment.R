library(here)
library(tidyverse)
library(latex2exp)
library(clusterProfiler)
library(future.apply)
library(qvalue)
library(patchwork)

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
  arrange(POS)


# Multiplicity ------------------------------------------------------------

# Crude first look at which genes are mutated in more than one replicate
# population

mult1 <-  df_final %>%
  group_by(REPLICATE, LOCUS_TAG) %>%
  slice(1) %>%
  ungroup() %>%
  add_count(LOCUS_TAG) %>%
  arrange(desc(n), LOCUS_TAG, REPLICATE) %>%
  select(REPLICATE, POS, REF, ALT, LOCUS_TAG, PRODUCT)


# Parallelism at nucleotide level -----------------------------------------

## Nucleotide multiplicity.
# For each site, we defined them multiplicity, $m_{i}$, as the number of
# populations with a point mutation detected at that site.  We calculated the
# multiplicity separately for both the mutator and nonmutator populations, so
# that the multiplicity could range from 1 to 3. Each mutation was then assigned
# a multiplicity score according to the site in which it occurred.

df_final %>%
  group_by(POS, REF, ALT) %>%
  summarize(mi=n()) %>%
  ungroup() %>%
  count(mi)

# So there are 5 sites where the same mutation occurred in 3 independent lines,
# 11 sites where the same mutation occurred in 2 independent lines, and 632
# sites where the mutation was unique in one population. To put these
# observations in context, we can compare them to a null model in which
# mutations are uniformly distributed across the sites in the genome.

Ltot <- 6722539

ntot  <- df_final %>%
  group_by(POS, REF, ALT) %>%
  summarize(mi=n()) %>%
  ungroup() %>%
  count(mi) %>% summarize(sum(n)) %>% pull

poisson_pmf <- function(mi, ntot, Ltot=6722539){
  dpois(mi, lambda = ntot/Ltot)
}

allele_mult <- df_final %>%
  group_by(POS, REF, ALT) %>%
  summarize(mi=n()) %>%
  ungroup() %>%
  count(mi, name="Survival_Observed") %>% 
  mutate(ntot = sum(Survival_Observed)) %>%
  mutate(poisson_pmf = poisson_pmf(mi-1, ntot)) %>%
  arrange(desc(mi)) %>%
  mutate(poisson_survive=cumsum(poisson_pmf)) %>%
  mutate(Survival_Expected=ntot*poisson_survive) %>%
  select(mi, Survival_Expected, Survival_Observed) %>%
  pivot_longer(-mi) %>%
  mutate(name=factor(name, levels=c("Survival_Observed", "Survival_Expected"), 
                     labels= c("Survival Observed", "Survival Expected")))

# plot survival function
mutations_survive_plot <- ggplot(allele_mult, aes(x=mi, y=value, color=name)) + 
  geom_step() +
  scale_y_log10( breaks = scales::trans_breaks("log10", function(x) 10^x),
                 labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x=TeX("Nucleotide multiplicitiy, \\textit{m}"), 
       y=TeX("Total mutations $\\geq \\textit{m}$")) +
  scale_color_manual(values = c("blue", "grey")) +
  mytheme() + 
  theme(legend.position = c(.25, .25)) +
  annotation_logticks(sides = "l", outside=T) +
  coord_cartesian(clip = "off")

mutations_survive_plot

# Parallelism at the gene level -------------------------------------------


# If selection pressures and mutation rates did not vary between genes, the
# number of mutations in each gene should be proportional to the target size.
# While it is difficult to estimate the local target size for beneficial,
# deleterious, and neutral mutations in any particular gene, we assume that gene
# length is a primary driver of the target size.  Similar to our
# nucleotide-level analysis above, we then define a multiplicity for each gene
# according to
#
# $ m_{i} = n_{i} \times \frac{\bar{L}}{L_{i}} $
#
# where $n_{i}$ is  the  number  of  mutations  in  $gene_{i}$ across  all
# replicate  populations  (including  indels and structural variants, but
# excluding synonymous mutations), $L_{i}$ is the total number of nonsynonymous
# and noncoding sites in $gene_{i}$, and $\bar{L}$ is the average value of
# $L_{i}$ across all genes in the genome. This definition ensures that under the
# null hypothesis, all genes have the same expected multiplicity $\bar{m} =
# \frac{n_{tot}}{n_{genes}}$.
#
# Using the observed and expected values, we can quantify the net increase of
# the log-likelihood of the alternative hypothesis relative to the null.
#
# $ \Delta l = \sum\limits_{i}n_{i} \log \left( \frac{m_{i}}{\bar{m}} \right) $
#
# where significance is assessed using permutation tests. Because this measure
# can be sensitive to $n_{tot}$  for comparisons across different strains and
# treatments we randomly sub-sampled mutations as a multinomial distribution,
# where the probability of sampling a mutation at gene i was given by $ p_{i} =
# {n_{i}} / {n_{tot}}$  Multinomial sampling was performed 10,000 times with a
# sub-sampled n tot set to 50.

lns <- read_tsv(here::here("data", "variant", "gene_lengths.tsv"), 
                col_names=c("LOCUS_TAG", "length"))

obs_gene_multiplicity <- df_final %>%
  dplyr::filter(!is.na(LOCUS_TAG)) %>% 
  #dplyr::filter(!str_detect(EFFECT, "synonymous_variant")) %>%
  group_by(LOCUS_TAG) %>%
  summarize(n_i = n(), n_populations=n_distinct(REPLICATE)) %>%
  full_join(., lns)  %>%
  rename(l_i=length) %>%
  arrange(LOCUS_TAG) %>%
  mutate(
    across(everything(), ~replace_na(.x, 0))
  )

# Observed total G_score
G_obs <- G_score(obs_gene_multiplicity)
G_obs


# Now we want to check whether this G score value is statistically significant
# so we use permutation tests

# required input stats about the different populations

df_pop_summary <- df_final %>%
  dplyr::filter(!is.na(LOCUS_TAG)) %>% 
  #dplyr::filter(!str_detect(EFFECT, "synonymous_variant")) %>%
  count(REPLICATE) %>%
  mutate(sum=sum(n), p=n/sum)

nmutations       <- df_pop_summary$sum[1]
population_num   <- nrow(df_pop_summary)
population_probs <- pull(df_pop_summary, p)
gene_names       <- lns$LOCUS_TAG

# equal probability independent of length
# gene_probs       = rep(1/nrow(lns), nrow(lns))

# probability depends on gene length
gene_probs       <- lns$length / sum(lns$length)

# Simulate 1 draw

G_score(G_prep(nmutations, population_num, population_probs, gene_probs, gene_names, lns))


# Simulate 10000 draws in parallel
parallelly::availableCores()

plan(multisession, workers = 15)

g_reps <- future_replicate(10000, G_score(G_prep(nmutations, population_num, population_probs, gene_probs, gene_names, lns)))

future:::ClusterRegistry("stop")

# Calculate P-value
length(g_reps[g_reps > G_obs])/length(g_reps)
# 0.0104

# plot

G_score_plot <- ggplot(tibble(g_reps), aes(x=g_reps)) + 
  geom_histogram() +
  #geom_vline(xintercept=G_obs, color="red") +
  geom_segment(x=G_obs, xend=G_obs, y=0, yend=150, color="red", size=2) +
  annotate("text", x = 2.38, y = 750, label = "P value = 0.01") +
  annotate("text", x = 2.38, y = 600, label = "Observed G\nscore = 2.37") +
  labs(x="G score", y="Count") +
  mytheme()

G_score_plot


# Calculate P value for each gene
gene_qvalues <- calculate_parallelism_qvalues(obs_gene_multiplicity)

# Calculate log-P value for each gene
gene_logpvalues <- calculate_parallelism_logpvalues(obs_gene_multiplicity)

# Collect P values for genes with at least nmin mutations across all populations
nmin <- 2
pooled_pvalues <- sort(gene_logpvalues[from_parallelism_statistics(obs_gene_multiplicity)$ns >= nmin])

# Unnormalized poisson survival vector
tsrv <- calculate_unnormalized_survival_from_vector(pooled_pvalues)

observed_ps <- tsrv[[1]]
observed_pvalue_survival <- tsrv[[2]]

# Estimate final survival function for genes
nglpsf_out <- NullGeneLogpSurvivalFunction(obs_gene_multiplicity, rep(0:399), 0.1)


# Prepare output for plotting

df2plot <- tibble(
  obs_p = observed_ps,
  Survival_Expected = nglpsf_out$survivals,
  Survival_Observed = observed_pvalue_survival
) %>%
  pivot_longer(-obs_p)  %>%
  mutate(name=factor(name, levels=c("Survival_Observed", "Survival_Expected"), 
                     labels= c("Survival Observed", "Survival Expected")))

segs <- tibble(
  x = c(min(observed_ps), nglpsf_out$pstar),
  y = c(nglpsf_out$num_significant, 0),
  xend = c(nglpsf_out$pstar, nglpsf_out$pstar),
  yend = c(nglpsf_out$num_significant, nglpsf_out$num_significant)
)

anns <- tibble(x = nglpsf_out$pstar, y = nglpsf_out$num_significant)

# Survival function plot

genes_survive_plot <- ggplot() +
  geom_segment(data = segs,
    aes(x = x, y = y, xend = xend, yend = yend), 
    linetype = "dashed", size = 0.25) +
  geom_step(data = df2plot, aes(x = obs_p, y = value, color = name)) +
  annotate("point", x = anns$x, y = anns$y) +
  annotate("text", x = anns$x + 3, y = anns$y + 3.5,
    label = TeX("Critical P*\nthreshold = $1.6 \\times 10^{-4}$")) +
  labs(x = TeX("$-log_{10}P$"), y = "Number of genes") +
  scale_color_manual(values = c("blue", "grey")) +
  scale_y_sqrt() +
  mytheme(legend.position = c(.8, .75))

genes_survive_plot

# make the final plot
pfinal <- mutations_survive_plot + G_score_plot + genes_survive_plot + plot_annotation(tag_levels = 'A')
pfinal

ggsave(here::here("figs", "figS1.svg"), pfinal, device="svg", width=7, height=2.5, units="in")


# Which genes have more parallel mutations than expected by chance?

gene_table_out <- data.frame(cbind(pooled_pvalues, nglpsf_out$survivals)) %>%
  rownames_to_column(var="LOCUS_TAG") %>%
  tibble() %>%
  left_join(., obs_gene_multiplicity) %>%
  filter(pooled_pvalues >= nglpsf_out$pstar) %>%
  left_join(., distinct(select(df_final, LOCUS_TAG, GENE, PRODUCT))) %>%
  rename(neg.log10P=pooled_pvalues, Expected=V2) %>%
  select(Locus=LOCUS_TAG, Gene=GENE, Product=PRODUCT, Length=l_i, Populations_detected=n_populations, Observed=n_i, Expected, neg.log10P)

write_tsv(gene_table_out, here::here("tables", "significant_genes.tsv"))

gene_table_out #%>% xtable::xtable()

obs_gene_multiplicity %>%
  mutate(mult=n_i*(mean(l_i)/l_i)) %>%
  filter(LOCUS_TAG %in% pull(gene_table_out, Locus)) %>%
  rename(Locus=LOCUS_TAG) %>%
  left_join(., gene_table_out ) %>%
  arrange(neg.log10P) %>%
  select(mult)


# Functional enrichment ---------------------------------------------------

# Load KEGG pathways and their descriptions

path2ko <- readr::read_tsv(here("data", "variant", "path2ko.tsv"), 
                          col_names = c("PATH", "KO"), 
                          col_types="cc")
pathdesc <- readr::read_tsv(here("data", "variant", "path_desc.tsv"), 
                           col_names = c("PATH", "DESC"), 
                           col_types="cc")

# Read in results of eggnogmapper for SBW25
eggnog <- read_tsv(here("data", "variant", "SBW25.annotations"),
                  col_names = c("query.name",
                                "seed.eggNOG.ortholog",
                                "seed.ortholog.evalue",
                                "seed.ortholog.score",
                                "predicted.tax.group",
                                "predicted.protein",
                                "go.terms",
                                "EC.number",
                                "KEGG.ko",
                                "KEGG.Pathway",
                                "KEGG.Module",
                                "KEGG.Reaction",
                                "KEGG.rclass",
                                "KEGG.BRITE",
                                "KEGG.TC",
                                "CAZy",
                                "BiGG.rxn",
                                "tax.scope", #eggNOG taxonomic level used for annotation
                                "eggNOG.OGs",
                                "bestOG", #deprecated, use smallest from eggnog OGs
                                "COG.category",
                                "description"),
                  col_types = "ccddcccccccccccccccccc")


# Formatting for use w/ clusterprofiler

keggpath2gene <- eggnog %>%
  select(PATH=KEGG.Pathway, Geneid=query.name) %>%
  separate_rows(PATH) %>%
  filter(!str_detect(PATH, "^ko")) %>%
  mutate(dum="path:") %>%
  unite(PATH, dum, PATH, sep='', remove=TRUE) %>% 
  mutate(PATH=ifelse(PATH=="path:NA", NA, PATH))

keggpath2name <- pathdesc

# GENES <- pfE3 %>% 
#   #dplyr::filter(FILTER=="PASS") %>% 
#   dplyr::filter(!is.na(LOCUS_TAG)) %>% 
#   dplyr::filter(!str_detect(EFFECT, "synonymous_variant")) %>%
#   #dplyr::filter(PRODUCT !="hypothetical protein") %>%
#   pull(LOCUS_TAG)

GENES <- names(pooled_pvalues[pooled_pvalues > nglpsf_out$pstar])

x <-  try(enricher(GENES, 
                 TERM2GENE=keggpath2gene, 
                 TERM2NAME=keggpath2name, 
                 pAdjustMethod = "fdr"), silent = TRUE)

x

#https://www.genome.jp/pathway/map02020
