<!-- badges: start -->
**The ISME Journal:** [![DOI:10.1038/s41396-023-01361-9](http://img.shields.io/badge/DOI-10.1038/s41396--023--01361--9-4F761F.svg)](https://doi.org/10.1038/s41396-023-01361-9)

**bioRxiv:** [![DOI:10.1101/2022.05.26.493533](http://img.shields.io/badge/DOI-10.1101/2022.05.26.493533-B31B1B.svg)](https://doi.org/10.1101/2022.05.26.493533)

**Code and data archive:** [![DOI](https://zenodo.org/badge/411568623.svg)](https://zenodo.org/badge/latestdoi/411568623)
<!-- badges: end -->

# Publication

"Localized coevolution between microbial predator and prey alters community-wide gene expression and ecosystem function"

- [Preprint available from bioRxiv](https://doi.org/10.1101/2022.05.26.493533)
- [Peer-reviewed version published in The ISME Journal](https://doi.org/10.1038/s41396-023-01361-9)

Data and code here is provided under the MIT License. Feel free to use or remix as you see fit.

# Project description
This is project contains data from the HAMBI species gene expression project. 30 HAMBI bacterial species were grown with either a low trait-diversity, ancestral *Pseudomonas fluorescens* SBW25 or a coevolved population of SBW25 that had been co-cultured with a *Tetrahymena* ciliate. These two treatments were conducted in the presence and absence of the coevolved ciliate. The experiment duration was 55 days. No fresh media was added, so the bacteria and ciliate were growing 
in a closed system without new nutrient inputs. RNAseq samples were collectd on days 4 and 45. 16S amplicon samples were collected on days 4, 41, and 45. Ciliate counts, bacterial CFUs, and community ATP concentrations were measured every ~ 4 days.

# Repository structure

1. `/R` contains R scripts
2. `/data` contains data that has been processed in some way for later use
3. `/data_raw` contains unprocessed data scraped from compute cluster
4. `/figs` contains figures generated from R scripts
5. `/sh` contains shell scripts. Mostly from running analysis on the [puhti compute cluster](https://docs.csc.fi/computing/systems-puhti/)
6. `/tables` contains summary tables generated from R scripts

# Obtaining sequencing data

Sequencing data is available from the NCBI Sequence Read Archive under [Bioproject PRJNA818876](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA818876) which can be viewed on the [NCBI SRA run selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA818876&o=acc_s%3Aa). 

## SRA download
Install the [NCBI SRA Toolkit](https://github.com/ncbi/sra-tools). You can download a prebuilt binary [from here.](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)

Access SRA data following the [instructions here.](https://github.com/ncbi/sra-tools/wiki/HowTo:-Access-SRA-Data)

You will need to [setup your configurations](https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration), but afterwards you can basically do:

```{bash}
prefetch SRR18441242
fasterq-dump SRR18441242
```

You will need to do this for all the SRA accessions associated with BioProject: [PRJNA818876](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA818876). For example:

```{bash}
cut -f1 tables/SraRunTable.tsv | tail -n +2 | while read ID; do
  prefetch $ID
  fasterq-dump $ID
done
```

This sequencing data can be preprocessed and mapped using the scripts in the `/sh/rnaseq`, `/sh/amplicon`, and `/sh/variant` directories.

# SBW25 variant analysis steps

Note this analysis borrows heavily from the [excellent paper by BH Good](https://www.nature.com/articles/nature24287) and the [code that he released publicly](https://github.com/benjaminhgood/LTEE-metagenomic). Most of this analysis I translated from his python code; some I wrote myself.

## Run on cluster:
1. `Shovill.sh` - Assemble progenitor reads using shovill
2. `Bwa.sh` - Map reads against Shovill assembly 
3. `BwaVsASM922.sh` - Map reads against NCBI SBW25 assembly ASM922
4. `OctopusIndividual.sh` - call variants on clonal samples
5. `OctopusPolyclone.sh` - call variants on mixed population samples
6. `SNPannCommands.sh` - kind of a hacky pipeline to convert Octopus output to tabular variant format. Uses snEff and snippy-vcf_to_tab from the Snippy software suite.

## Analysis in R
1. `01_snp_functional_enrichment.R` -- Identify nucleotide sites and genes under parallel evolution, check for functional enrichment of mutated genes, plot Figure S1
2. `02_plot_variant_freq.R` -- Produce Fig. 2 from the main text

# Longitudinal Gaussian Process Regression

## Run on cluster:
1. `submit_lgpr.sh` -- submit steps 2-5 below using this script

## Analysis in R
1. `01_format_data.R` -- Formats data to be used in lgpr
2. `02_puhti_lgpr_ATP.R` -- run lgpr on cluster for ATP
3. `03_puhti_lgpr_ciliate.R` -- run lgpr on cluster for ciliates
4. `04_puhti_lgpr_cfus.R` -- run lgpr on cluster for bacteria colony forming units
5. `05_puhti_lgpr_opd.R` -- run lgpr on cluster for bacteria optical density
6. `06_lgpr_process_plot.R` -- Process output from lgpr and plot Figure 3

# Amplicon analysis steps
## Run on cluster:
1. `AmpliconQualityControl.sh`
2. `AmpliconMapping.sh`

## Analysis in R
1. `01_rpkm2tab.R` -- Prepare amplicon count tables. Count tables saved in `data`
2. `02_make_phyloseq.R` -- Make phyloseq object
3. `03_make_figS3.R` -- Make Fig. S3
4. `04_shannon_diversity.R` -- Runs DivNet estimate and plots DivNet vs Plugin. Also runs breakawy for testing difference in diversity between samples. Generates results for Table S3
5. `05_ordination.R` -- Run ordination analysis and PERMANOVA. Generates results for Table S3
6. Run corncob beta binomial regression for differential abundance. Generates results for Table S3. __Fitting corncob models with the bootstrap likelihood ratio test takes about an hour each.__ This was done separately for each model in the scripts below. What each of these scripts does is fit the same full model `~ days + pseudomonas_hist * predation` but with different null models to test the effect of leaving out these different terms by way of parametric bootstrap likelihood ratio tests.
    a. `06_corncob_SBW25.R`
    b. `07_corncob_evolution.R`
    c. `08_corncob_predation.R`
    d. `09_corncob_interaction.R`
7. `10_make_Fig3.R` -- compile results to make main text Fig. 4

# RNAseq analysis steps

## Run on cluster
1. `process_metatranscriptomes.sh` - quality control RNAseq data
2. `map_count_metaT.sh` - map reads against 30 HAMBI genomes using bbmap

## Analysis in R
1. `01_read_filter_mrna.R` - read the feat counts and format
2. `02_rnaseq_stats.R` - calculate general stats (% noncoding RNA, etc...)
3. `03_species_rna_relative_abundance.R` - make supplementary Fig S3 showing proportion of bacterial species in the RNA dataset (excluding Tetrahymena).
4. `04_prep_deseq_data.R` - prepare all the files necessary to run the rlog transform and the deseq procedure
5. `05_rlog_transform.R` - perform the rlog transform necessary for distatis and to plot the ordination
 - Implements the sum-taxon scaling + species amplicon abundance estimate normalization approach described in [Zhang 2021](https://doi.org/10.1093/bioinformatics/btab327).
 - The idea to modify the internal normalization factors in DESeq2 is [from Mike Love himself](https://support.bioconductor.org/p/99165/). There is additional information about the `normalizationFactors` approach in the [DESeq2 vignette.](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#sample-gene-dependent-normalization-factors). It is important that the normalization matrix has row-wise geometric means of 1, so that the mean of normalized counts for a gene is close to the mean of the unnormalized counts. This is accomplished by dividing out the current row geometric means in the normalization matrix. Modifying the normalization factor matrix replaces the `estimateSizeFactors` step which occurs within the DESeq function. The DESeq function will look for pre-existing normalization factors and use these in the place of size factors (and a message will be printed confirming this). The sizeFactor estimation process is described in very simple terms [here](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html).
 - [This post](http://seqanswers.com/forums/showthread.php?t=39228) introduces the `estimateNormFactors` function from DESeq2 but it is not exported. [See also here.](https://support.bioconductor.org/p/97676/) My implementation is just based on this function. Some additional useful information [here](https://support.bioconductor.org/p/9135244/), [here](https://support.bioconductor.org/p/131438/), and [here](https://support.bioconductor.org/p/97936/)
 - As a final comment... the field is not even clear whether these approaches are necessary? At least there are many 'high profile' papers that don't appear to consider taxon sum scaling and just DESeq2 the same as for single-organism RNAseq experiments. For examples where standard DESeq2 approach is applied to MTX, see [here](https://doi.org/10.1038/s41586-018-0207-y), [here](https://github.com/alexcritschristoph/angelo_biosynthetic_genes_analysis), [here](https://doi.org/10.1038/s41396-020-00820-x), and [here](https://github.com/speeding-up-science-workshops/metatranscriptomics-visualizations)
6. `06_distatis.R` - perform distatis analysis, clustering, plots Fig. 5 from main text, and performs PERMANOVA
7. `07_deseq_tetrahymena.R` - perform simple DeSeq2 analysis and functional enrichment for Tetrahymena using `clusterProfiler`
8. `08_deseq.R` - perform DeSeq2 analysis. Using the same normalization approach from step 5.
9. `09_deseq_contrasts.R` - run contrasts between Evolved and Progenitor Pseudomonas. Requires `apeglm` package to shrink estimated log-fold changes near 0 counts.
10. `10_functional_enrichments.R` - perform functional enrichment analysis for the bacterial community. Produces Figure 6 and supplementary Figures S5 and S6
11. `11_venndiagram.R` - Compare differentially expressed genes from different contrasts in step 9. Produces Figure S7
12. `12_volcano_plot.R` - Plot Fig. S8
13. `13_fraction_diff_expressed.R` - basic statistics about diff expressed genes.
