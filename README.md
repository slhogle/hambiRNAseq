# hambiRNAseq

<!-- badges: start -->
<!-- badges: end -->

This is project contains data from the HAMBI species gene
expression project. 30 HAMBI bacterial species were grown with
either a low trait-diversity, ancestral Pseudomonas fluorescens SBW25
or a local adapted population of SBW25 that had been co-cultured with 
a *Tetrahymena* ciliate. These two treatments were conducted in the
presence and absence of the ciliate. The experiment duraction was 55 days. No fresh media was added, so the bacteria and ciliate were growing 
in a closed system without new nutrient inputs. RNAseq samples were collectd on days 4 and 45. 16S amplicon samples were collected on days 4, 41, and 45. Ciliate counts, bacterial CFUs, and community ATP concentrations were measured every ~ 4 days.

1. `/sh` contains shell scripts. Mostly from running analysis on the [puhti compute cluster](https://docs.csc.fi/computing/systems-puhti/)
2. `/R` contains R scripts
3. `dataRaw` contains unprocessed data scraped from puhti
4. `/data` contains data that has been processed in some way for later use
5. `/figs` contains figures generated from R scripts
6. `/tables` contains summary tables generated from R scripts
7. `/rendered` contains any markdown or HTML that gets rendered in the process of running the scipts.

# MANUSCRIPT:

""

[Preprint available from bioRxiv]()

Data and code here is provided under the MIT License. Feel free to use or remix as you see fit.

# PROCESSING 16S AMPLICON DATA

## DOWNLOAD
Install the [NCBI SRA Toolkit](https://github.com/ncbi/sra-tools). You can download a prebuilt binary [from here.](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)

Access SRA data following the [instructions here.](https://github.com/ncbi/sra-tools/wiki/HowTo:-Access-SRA-Data)

You will need to [setup your configurations](https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration), but afterwards you could basically do:

```{bash}
prefetch SRR14323812
fasterq-dump SRR14323812
```

You will need to do this for all the SRA accessions associated with BioProject: [PRJNA725120](https://www.ncbi.nlm.nih.gov/bioproject/725120). For example:

```{bash}
cut -f1 00_SRA/sra_records.tsv | while read ID; do
  prefetch $ID
  fasterq-dump $ID
done
```

This will take some time...

# TODO
1. Rename fastq files on puhti for amplicon and RNAseq. Use `16Smapping_newnames.tsv` for name conversions
2. Stage files for upload to NCBI

# ANALYSIS
To run the full analysis (will take some time). Note you must untar `rawData/16SAmplicon/bbmapRPKM.tar.gz` first.

```{bash}
bash sh/runFullAnalysis.sh
```

## SBW25 variant analysis steps

### Run on cluster:
1. `Shovill.sh` - Assemble progenitor reads using shovill
2. `Bwa.sh` - Map reads against Shovill assembly 
3. `BwaVsASM922.sh` - Map reads against NCBI SBW25 assembly ASM922
4. `OctopusIndividual.sh` - call variants on clonal samples
5. `OctopusPolyclone.sh` - call variants on mixed population samples
6. `SNPannCommands.sh` - kind of a hacky pipeline to convert Octopus output to tabular variant format. Uses snEff and snippy-vcf_to_tab from the Snippy software suite.

### Analysis in R

Note this analysis borrows heavily from the [excellent paper by BH Good](https://www.nature.com/articles/nature24287) and the [code that he released publicly](https://github.com/benjaminhgood/LTEE-metagenomic). Most of this analysis I translated from his python code; some I wrote myself.

### Parallelism at nucleotide-level
For each site, we defined them multiplicity, $m_{i}$, as the number of populations with a point mutation detected at that site.  We calculated the multiplicity separately for both the mutator and nonmutator populations, so that the multiplicity could range from 1 to 3. Each mutation was then assigned a multiplicity score according to the site in which it occurred.

So there are 5 sites where the same mutation occurred in 3 independent lines, 11 sites where the same mutation occurred in 2 independent lines, and 632 sites where the mutation was unique in one population.

To put these observations in context, we can compare them to a null model in which mutations are uniformly distributed across the sites in the genome. In this model, the expected fraction of mutations with multiplicity $\geq$ min a sample of size $n_{tot}$ is given by

$S(m) \approx \sum\limits_{n \geq m} \frac{n}{n_{tot}} \times L_{tot} \times \frac{ \left( \frac{n_{tot}}{L_{tot}} \right) ^n}{n!} \times e^{-n_{tot}/L_{tot}}$

### Parallelism at the gene-level

If selection pressures and mutation rates did not vary between genes, the number of mutations in each gene should be proportional to the target size.  While it is difficult to estimate the local target size for beneficial, deleterious, and neutral mutations in any particular gene, we assume that gene length is a primary driver of the target size.  Similar to our nucleotide-level analysis above, we then define a multiplicity for each gene according to

$m_{i} = n_{i} \times \frac{\bar{L}}{L_{i}}$

where $n_{i}$ is  the  number  of  mutations  in  $gene_{i}$ across  all  replicate  populations  (including  indels and structural variants, but excluding synonymous mutations), $L_{i}$ is the total number of nonsynonymous and noncoding sites in $gene_{i}$, and $\bar{L}$ is the average value of $L_{i}$ across all genes in the genome. This definition ensures that under the null hypothesis, all genes have the same expected multiplicity $\bar{m} = \frac{n_{tot}}{n_{genes}}$. 

Using the observed and expected values, we can quantify the net increase of the log-likelihood of the alternative hypothesis relative to the null.

$\Delta l = \sum\limits_{i}n_{i} \log \left( \frac{m_{i}}{\bar{m}} \right)$

where significance is assessed using permutation tests. Because this measure can be sensitive to $n_{tot}$  for comparisons across different strains and treatments we randomly sub-sampled mutations as a multinomial distribution, where the probability of sampling a mutation at gene i was given by $p_{i} = {n_{i}} / {n_{tot}}$  Multinomial sampling was performed 10,000 times with a sub-sampled n tot set to 50.

## Longitudinal Gaussian Process Regression



## Amplicon analysis steps
### Run on cluster:
1. `AmpliconQualityControl.sh`
2. `AmpliconMapping.sh`

### Analysis in R
1. `rpkm2tab.R` -- Prepare amplicon count tables. Count tables saved in `data`
2. `makePhyloseq.R` -- Make phyloseq object
3. `makeFigS2.R` -- Make Fig. S2
4. `shannonDiversity.R` -- Runs DivNet estimate and plots DivNet vs Plugin. Also runs breakawy for testing difference in diversity between samples. Generates results for Table S3
5. `ordination.R` -- Run ordination analysis and PERMANOVA. Generates results for Table S3
6. Run corncob beta binomial regression for differential abundance. Generates results for Table S3. __Fitting corncob models with the bootstrap likelihood ratio test takes about an hour each.__ This was done separately for each model in the scripts below. What each of these scripts does is fit the same full model `~ days + pseudomonas_hist * predation` but with different null models to test the effect of leaving out these different terms by way of parametric bootstrap likelihood ratio tests.
    a. `corncobEvolution.R` 
    b. `corncobPredation.R`
    c. `corncobInteraction.R`
    d. `corncobSBW25.R`
7. `makeFig3.R` -- compile results to make main text Fig. 3

## RNAseq analysis steps

### Run on cluster
1. `process_metatranscriptomes.sh` - quality control RNAseq data
2. `map_count_metaT.sh` - map reads against 30 HAMBI genomes using bbmap

### Analysis in R
1. `01_read_filter_mrna.R` - read the feat counts and format
2. `02_rnaseq_stats.R` - calculate general stats (% noncoding RNA, etc...)
3. `03_species_rna_relative_abundance.R` - make supplementary Fig S3 showing proportion of bacterial species in the RNA dataset (excluding Tetrahymena).
4. `04_prep_deseq_data.R` - prepare all the files necessary to run the rlog transform and the deseq procedure
5. `05_rlog_transform.R` - perform the rlog transform necessary for distatis and to plot the ordination
 - Implements the sum-taxon scaling + species amplicon abundance estimate normalization approach described in [Zhang 2021](https://doi.org/10.1093/bioinformatics/btab327).
 - The idea to modify the internal normalization factors in DESeq2 is [from Mike Love himself](https://support.bioconductor.org/p/99165/). There is additional information about the `normalizationFactors` approach in the [DESeq2 vignette.](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#sample-gene-dependent-normalization-factors). It is important that the normalization matrix has row-wise geometric means of 1, so that the mean of normalized counts for a gene is close to the mean of the unnormalized counts. This is accomplished by dividing out the current row geometric means in the normalization matrix. Modifying the normalization factor matrix replaces the `estimateSizeFactors` step which occurs within the DESeq function. The DESeq function will look for pre-existing normalization factors and use these in the place of size factors (and a message will be printed confirming this). The sizeFactor estimation process is described in very simple terms [here](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html).
 - [This post](http://seqanswers.com/forums/showthread.php?t=39228) introduces the `estimateNormFactors` function from DESeq2 but it is not exported. [See also here.](https://support.bioconductor.org/p/97676/) My implementation is just based on this function. Some additional useful information [here](https://support.bioconductor.org/p/9135244/), [here](https://support.bioconductor.org/p/131438/), and [here](https://support.bioconductor.org/p/97936/)
 - As a final comment... the field is not even clear whether these approaches are necessary? At least there are many 'high profile' papers that don't appear to consider taxon sum scaling and just DESeq2 the same as for single-organism RNAseq experiments. For examples where standard DESeq2 approach is applied to MTX, see [here](https://doi.org/10.1038/s41586-018-0207-y), [here](https://github.com/alexcritschristoph/angelo_biosynthetic_genes_analysis), [here](https://doi.org/10.1038/s41396-020-00820-x), and [here](https://github.com/speeding-up-science-workshops/metatranscriptomics-visualizations)
6. `06_distatis.R` - perform distatis analysis, clustering, plots Fig. 4 from main text, and performs PERMANOVA
7. `07_deseq_tetrahymena.R` - perform simple DeSeq2 analysis and functional enrichment for Tetrahymena using `clusterProfiler`
8. `08_deseq.R` - perform DeSeq2 analysis. Using the same normalization approach from step 5.
9. `09_deseq_contrasts.R` - run contrasts between Evolved and Progenitor Pseudomonas. Requires `apeglm` package to shrink estimated log-fold changes near 0 counts.
10. `10_volcano_plot.R` - Plot Fig. S5
11. `11_fraction_diff_expressed.R` - basic statistics about diff expressed genes.
12. `12_functional_enrichment.R` - perform the `clusterProfiler` functional enrichment analysis

### Notes
As seen from the DiSTATIS ordination analysis there is quite a big difference between days 45
and 4. So much so that there is not much point even plotting them on the same ordination since the 
variation with time swamps out any variation from predation or evolution. Basically, a lot 
happped between those time points including a nutrient addition on day 41. I've found it 
best to analyze the day 4 and day 45 samples separately.

Implements the DiSTATIS method which a 3-way generalization of metric multidimensional scaling (a.k.a. classical MDS or principal coordinate analysis). DiSTATIS takes a set of K distance matrices describing a set of I observations and computes:

1. a set of factor scores that describes the similarity structure of the distance matrices (e.g., what distance matrices describe the observations in the same way, what distance matrices differ from each other)
2. a set of factor scores (called the compromise factor scores) for the observations that best describes the similarity structure of the observations
3. partial factor scores that show how each individual distance matrix "sees" the compromise space. distatis computes the compromise as an optimum linear combination of the cross-product matrices associated to each distance matrix. distatis can also be applied to a set of covariance matrices.

_From: Abdi, H., L. J. Williams, D. Valentin, and M. Bennani-Dosse. 2012. STATIS and DISTATIS: optimum multitable principal component analysis and three way metric multidimensional scaling. WIREs Comp Stat 4: 124–167._

The general idea behind STATIS is to analyze the structure of the individual data sets (i.e., the relation between the individual data sets) and to derive from this structure an optimal set of weights for computing the best common representation of the observations called the compromise, or also sometimes the consensus. To compute this compromise, the elements of each table are multiplied by the optimal weight of this table and the compromise is obtained by the addition of these ‘weighted’ K tables (i.e., the compromise is a linear combination of the tables). These weights are chosen so that the compromise provides the best representation (in a least square sense) of the whole set of tables. The PCA of the compromise decomposes the variance of the compromise into a set of new orthogonal variables called principal components (also often called dimensions, axes, factors, or even latent variables) ordered by the amount of variance that each component explains. The coordinates of the observations on the components are called factor scores and these can be used to plot maps of the observations in which the observations are represented as points such that the distances in the map best reflect the similarities between the observations. The position of the observations ‘as seen by’ each data set can be also represented as points in the compromise. As the components are obtained by combining the original variables, each variable contributes a certain amount to each component. This quantity, called the loading of a variable on a component, reflects the importance of that variable for this component and can also be used to plot maps of the variables that reflect their association. Finally, as a byproduct of the computation of the optimal weights, the data sets can be also represented as points in a multidimensional space.
