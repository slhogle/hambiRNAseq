library(here)
library(tidyverse)
library(vegan)
library(withr)
library(ape)
library(compositions)
source(here::here("R", "16S_amplicon", "amplicon_funs.R"))


# Read data ---------------------------------------------------------------

counts <- read_tsv(here("data", "16S_amplicon", "species_counts.tsv"), col_types = "ccccd")

metadf <- read_tsv(here("data", "sample_metadata.tsv"), col_types = "cfffff") %>%
  mutate(
    pseudomonas_hist = factor(pseudomonas_hist, levels = c("ancestral", "coevolved")),
    predation = factor(predation, levels = c("no", "yes")),
    days = factor(days, levels = c("4", "41", "45"))
  ) %>%
  mutate(dayexpcombo = interaction(days, pseudomonas_hist, predation)) %>%
  mutate(expcombo = interaction(pseudomonas_hist, predation))

full <- left_join(counts, metadf)


# Format ------------------------------------------------------------------


## transform. replace 0s with smallest non-zero count
mymat <- full %>%
  select(sample, strainID, count) %>%
  group_by(sample) %>%
  mutate(count = ifelse(count == 0, minnz(count), count)) %>%
  ungroup() %>%
  pivot_wider(names_from = "strainID", values_from = "count") %>%
  column_to_rownames(var = "sample") %>%
  data.frame()

comp <- acomp(mymat)
# balclr <- clr(comp)

# distance matrices -------------------------------------------------------

# Aitchison dissimilarity matrix
# with_preserve_seed(dist <- compositions::dist(balclr, method = "euclidean"))  

# Bray-curtis dissimilarity 
with_seed(1237845, dist <- vegdist(comp, method = "bray"))  

# can use PCA = vegan::rda with no covariates
# pcoa_ord <- pcoa(dist)  
with_seed(12785, nmds_ord <- metaMDS(dist))

write_rds(nmds_ord, here::here("data", "16S_amplicon", "ordination.rds"))

# permanova ---------------------------------------------------------------

# Results are the same whether you use clr transformation with euclidean 
# distance or proportions with bray-curtis dissimilarity 

# should maybe change to days * pseudomonas_hist * predation 
# to allow for different effect with time - although doesn't seem necessary

with_seed(17853,
  ad.res <- adonis2(
    as.data.frame(comp) ~ days + pseudomonas_hist * predation,
    data = metadf,
    method = "bray"
  )
)

print(ad.res)

# ad.res %>%
#   xtable::xtable()


# permdisp ----------------------------------------------------------------

# Check for equal dispersions among groups. This is important because you 
# can get significant Adonis results if there is no difference in mean 
# centroids between groups but rather a difference in the variances between 
# groups. We need to make sure what we are seeing above is real. 
# Spoiler: it is real. You can read on to see why...

# betadisper is a multivariate analogue of Levene's test for 
# homogeneity of variances. It implements PERMDISP2 for testing of 
# multivariate homogeneity. Apparently you can only check one group at a time...

# First step is to look at all possible groups simultaneously as a combination
# of treatment and day.

# Results show that there is no significant variance differences 
# between treatment-day groupings.

mygroups <- paste(pull(metadf[, 6]), pull(metadf[, 2]), sep = "_")

with_seed(347856, mod0 <-
                     betadisper(vegdist(comp, method = "bray"), mygroups))

# To test if one or more groups is more variable than the others, 
# ANOVA of the distances to group centroids can be performed and parametric 
# theory used to interpret the significance of F. *Not Significant*
print(anova(mod0))
# anova(mod0)%>% xtable::xtable()

# permutest permutes model residuals to generate a permutation distribution of 
# F under the null hypothesis of no difference in dispersion between groups.
# *Not Significant*
print(permutest(mod0, permutations = 999))

# classical comparison of group dispersions calculated using Tukey's Honest 
# Significant Differences between groups, via TukeyHSD.betadisper.
# *No pairwise comparisons Significant*
print(TukeyHSD(mod0))

# visualize results
plot(mod0, ellipse = TRUE, hull = FALSE, conf = 0.90, label = FALSE)
boxplot(mod0)

# permdisp individual -----------------------------------------------------

# Checking the effect of day, SBW25 evolution, and predation individually...
# Read in one of Jari Oksanen's notes that you should check factors individually
# At least here it produces the same results

#### Day ####
with_seed(1324895, mod1 <-
                     betadisper(vegdist(comp, method = "euclidean"),
                                pull(metadf[, 2])))

# ANOVA of the distances to group centroids *Not Significant*
print(anova(mod1))

# permutest *Not Significant*
print(permutest(mod1, permutations = 999))

# Tukey's Honest Significant Differences *No pairwise comparisons Significant*
print(TukeyHSD(mod1))
plot(mod1, ellipse = TRUE, hull = FALSE, conf = 0.90)
boxplot(mod1)


#### SBW25 evolution ####
with_seed(1234895, mod2 <-
                     betadisper(vegdist(comp, method = "euclidean"),
                                pull(metadf[, 3])))

# ANOVA of the distances to group centroids *Not Significant*
print(anova(mod2))

#permutest  *Not Significant*
print(permutest(mod2, permutations = 999))

print(TukeyHSD(mod2))
plot(mod2, ellipse = TRUE, hull = FALSE, conf = 0.90)
boxplot(mod2)


#### Predation ####

with_seed(12785, mod3 <-
                     betadisper(vegdist(comp, method = "euclidean"),
                                pull(metadf[, 4])))

# ANOVA of the distances to group centroids *Not Significant*
print(anova(mod3))

#permutest *Not Significant*
print(permutest(mod3, permutations = 999))
plot(mod3, ellipse = TRUE, hull = FALSE, conf = 0.90)
boxplot(mod3)