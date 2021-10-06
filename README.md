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

Note this analysis borrows heavily from the [excellent paper by BH Good](https://www.nature.com/articles/nature24287) and the [code that he released publicly](https://github.com/benjaminhgood/LTEE-metagenomic). Most of this analysis I translated from his python code; some I wrote myself.

### Parallelism at nucleotide-level
For each site, we defined them multiplicity, $m_{i}$, as the number of populations with a point mutation detected at that site.  We calculated the multiplicity separately for both the mutator and nonmutator populations, so that the multiplicity could range from 1 to 3. Each mutation was then assigned a multiplicity score according to the site in which it occurred.

So there are 5 sites where the same mutation occurred in 3 independent lines, 11 sites where the same mutation occurred in 2 independent lines, and 632 sites where the mutation was unique in one population.

To put these observations in context, we can compare them to a null model in which mutations are uniformly distributed across the sites in the genome. In this model, the expected fraction of mutations with multiplicity≥min a sample of size $n_{tot}$ is given by

$S(m) \approx \sum\limits_{n \geq m} \frac{n}{n_{tot}} \times L_{tot} \times \frac{ \left( \frac{n_{tot}}{L_{tot}} \right) ^n}{n!} \times e^{-n_{tot}/L_{tot}}$

### Parallelism at the gene-level

If selection pressures and mutation rates did not vary between genes, the number of mutations in each gene should be proportional to the target size.  While it is difficult to estimate the local target size for beneficial, deleterious, and neutral mutations in any particular gene, we assume that gene length is a primary driver of the target size.  Similar to our nucleotide-level analysis above, we then define a multiplicity for each gene according to

$ m_{i} = n_{i} \times \frac{\bar{L}}{L_{i}} $

where $n_{i}$ is  the  number  of  mutations  in  $gene_{i}$ across  all  replicate  populations  (including  indels and structural variants, but excluding synonymous mutations), $L_{i}$ is the total number of nonsynonymous and noncoding sites in $gene_{i}$, and $\bar{L}$ is the average value of $L_{i}$ across all genes in the genome. This definition ensures that under the null hypothesis, all genes have the same expected multiplicity $\bar{m} = \frac{n_{tot}}{n_{genes}}$. 

Using the observed and expected values, we can quantify the net increase of the log-likelihood of the alternative hypothesis relative to the null.

$ \Delta l = \sum\limits_{i}n_{i} \log \left( \frac{m_{i}}{\bar{m}} \right) $

where significance is assessed using permutation tests. Because this measure can be sensitive to $n_{tot}$  for comparisons across different strains and treatments we randomly sub-sampled mutations as a multinomial distribution, where the probability of sampling a mutation at gene i was given by $ p_{i} = {n_{i}} / {n_{tot}}$  Multinomial sampling was performed 10,000 times with a sub-sampled n tot set to 50.


## Longitudinal Gaussian Process Regression



## Amplicon analysis steps
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


