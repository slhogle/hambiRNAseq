
# DESeq2 Functions --------------------------------------------------------

rlogifySTC <-
  function(data,
           species,
           metadata,
           rnatotals,
           ntop = 500) {
    print(paste0("Working on ", species, "..."))
    
    # needs to be used with imap or equivalent because requires species name in
    # the final formula
    design <- paste("~", species, "+", "summary_cat")
    
    # make summary_cat a factor
    metadata$summary_cat <- factor(metadata$summary_cat)
    
    # remove genes with zero counts in all samples
    data <- data[rowMeans(data) > 0, ]
    
    # create Deseq2 dataset
    dds <- DESeqDataSetFromMatrix(
      countData = data,
      colData = metadata,
      design = as.formula(design),
      tidy = FALSE
    )
    
    # These steps replace estimateSizeFactors which occurs within the DESeq
    # function. The DESeq function will look for pre-existing normalization
    # factors and use these in the place of size factors (and a message will be
    # printed confirming this).
    
    # none of the normalization factors can be zero
    rna_abund <- rnatotals[["coding"]][[species]]
    rna_abund <- replace(rna_abund, rna_abund == 0, NA)
    rna_abund <-
      replace(rna_abund, is.na(rna_abund), min(rna_abund, na.rm = TRUE))
    
    normMatrix  <- t(matrix(rep(rna_abund, nrow(data)), ncol = nrow(data)))
    sizeFactors <- estimateSizeFactorsForMatrix(data / normMatrix)
    normFactors <- t(t(normMatrix) * sizeFactors)
    normFactors <- normFactors / exp(rowMeans(log(normFactors)))
    
    # assign the normalization factors
    normalizationFactors(dds) <- normFactors
    
    # perform regularized log transformation of data
    # important to set blind = FALSE
    rld <- rlog(dds, blind = FALSE)
    
    # Select the top n genes with the highest row variance
    rv <- rowVars(assay(rld))
    select <-
      order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    topN <- t(assay(rld)[select,])
    return(topN)
  }

rlogifyASC <-
  function(data,
           species,
           metadata,
           rnatotals,
           ntop = 500) {
    # needs to be used with imap or equivalent because requires species name in
    # the final formula
    print(paste0("Working on ", species, "..."))
    
    design <- "~ summary_cat"
    
    # make summary_cat a factor
    metadata$summary_cat <- factor(metadata$summary_cat)
    
    # remove genes with zero counts in all samples
    data <- data[rowMeans(data) > 0,]
    
    # create Deseq2 dataset
    dds <- DESeqDataSetFromMatrix(
      countData = data,
      colData = metadata,
      design = as.formula(design),
      tidy = FALSE
    )
    
    # These steps replace estimateSizeFactors which occurs within the DESeq
    # function. The DESeq function will look for pre-existing normalization
    # factors and use these in the place of size factors (and a message will be
    # printed confirming this).
    
    # none of the normalization factors can be zero
    sp_abund <- metadata[[species]]
    sp_abund <- replace(sp_abund, sp_abund == 0, NA)
    sp_abund <-
      replace(sp_abund, is.na(sp_abund), min(sp_abund, na.rm = TRUE))
    
    normMatrix  <-
      t(matrix(rep(sp_abund, nrow(data)), ncol = nrow(data)))
    sizeFactors <- estimateSizeFactorsForMatrix(data / normMatrix)
    normFactors <- t(t(normMatrix) * sizeFactors)
    normFactors <- normFactors / exp(rowMeans(log(normFactors)))
    
    # assign the normalization factors
    normalizationFactors(dds) <- normFactors
    
    # perform regularized log transformation of data
    # important to set blind = FALSE
    rld <- rlog(dds, blind = FALSE)
    
    # Select the top n genes with the highest row variance
    rv <- rowVars(assay(rld))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    topN <- t(assay(rld)[select,])
    return(topN)
  }

deseqify <-
  function(data, metadata, mycontrast, design = "~ summary_cat") {
    # make summary_cat a factor
    metadata$summary_cat <- factor(metadata$summary_cat)
    
    # remove genes with zero counts in all samples
    data <- data[rowMeans(data) > 0, ]
    
    # create Deseq2 dataset
    dds <- DESeqDataSetFromMatrix(
      countData = data,
      colData = metadata,
      design = as.formula(design),
      tidy = FALSE
    )
    # run deseq process
    dds.res <- DESeq(dds)
    # get results as tibble
    dds.res.filt <-
      as_tibble(results(dds.res, contrast = mycontrast), 
                rownames = "Geneid") %>% 
      filter(padj < 0.1)
    return(dds.res.filt)
  }

deseqifySTC <-
  function(data,
           species,
           metadata,
           rnatotals, 
           overall = FALSE) {
    print(paste0("Working on ", species, "..."))
    
    # needs to be used with imap or equivalent because requires species name in
    # the final formula
    
    # + pseudomonas_hist:predation
    if ( overall == TRUE ) {
      design <- paste("~", species, "+", "pseudomonas_hist + predation")
    }  else {
      design <- paste("~", species, "+", "summary_cat")
    }
  
    # make summary_cat a factor
    metadata$summary_cat <- factor(metadata$summary_cat)
    
    # remove genes with zero counts in all samples
    data <- data[rowMeans(data) > 0, ]
    
    # create Deseq2 dataset
    dds <- DESeqDataSetFromMatrix(
      countData = data,
      colData = metadata,
      design = as.formula(design),
      tidy = FALSE
    )
    
    # These steps replace estimateSizeFactors which occurs within the DESeq
    # function. The DESeq function will look for pre-existing normalization
    # factors and use these in the place of size factors (and a message will be
    # printed confirming this).
    
    # none of the normalization factors can be zero
    rna_abund <- rnatotals[["coding"]][[species]]
    rna_abund <- replace(rna_abund, rna_abund == 0, NA)
    rna_abund <-
      replace(rna_abund, is.na(rna_abund), min(rna_abund, na.rm = TRUE))
    
    normMatrix  <- t(matrix(rep(rna_abund, nrow(data)), ncol = nrow(data)))
    sizeFactors <- estimateSizeFactorsForMatrix(data / normMatrix)
    normFactors <- t(t(normMatrix) * sizeFactors)
    normFactors <- normFactors / exp(rowMeans(log(normFactors)))
    
    # assign the normalization factors
    normalizationFactors(dds) <- normFactors
    
    # Run DESEq
    dds <- DESeq(dds)
    
    return(dds)
  }


# Formatting Functions ----------------------------------------------------

# 03_species_rna_relative_aundance.R
# normalizes by genome size and number of noncoding RNAs
comp_lump <- function(.data) {
  .data %>%
    group_by(sample, genome) %>%
    summarize(s = sum(mappedreads)) %>%
    ungroup() %>%
    left_join(., hambi_anno) %>%
    mutate(s = s / count) %>%
    mutate(genome2 = fct_lump(genome, prop = 0.01, w = s)) %>%
    group_by(sample, genome2) %>%
    tally(s) %>%
    mutate(f = n / sum(n) * 100) %>%
    ungroup() %>%
    mutate(genome = as.character(genome2))
}

# 04_prep_deseq_data.R
rnadfformat <- function(.data) {
  .data %>%
    filter(!str_detect(genome, "Tthermo")) %>%
    group_by(sample, rnatype, genome) %>%
    summarize(s = sum(mappedreads)) %>%
    arrange(sample, genome) %>%
    ungroup() %>%
    group_by(rnatype) %>%
    group_split() %>%
    map(.,
        pivot_wider,
        names_from = genome,
        values_from = s) %>%
    map(., select,-rnatype) %>%
    map(., column_to_rownames, var = "sample")
}

split_wide_mat <- function(.data) {
  map(.data, select, Geneid, sample, mappedreads) %>%
    map(., arrange, Geneid, sample) %>%
    map(.,
        pivot_wider,
        names_from = sample,
        values_from = mappedreads) %>%
    map(., column_to_rownames, var = "Geneid") %>%
    map(., as.matrix)
}

# 06_distatis.R
# Get dataframes of partial scores
get_partial <- function(distatis_result, metadata) {
  as.data.frame.table(distatis_result$res4Splus$PartialF) %>%
    as_tibble() %>%
    dplyr::rename(
      sample = Var1,
      factor = Var2,
      table = Var3,
      value = Freq
    ) %>%
    spread(key = "factor", value = "value") %>%
    mutate(sample = factor(sample),
           table = factor(table)) %>%
    left_join(., rownames_to_column(metadata, var = "sample"))
}

# bootstraps for distatis compromise
boot2df <- function(distatis_comp_boot,
                    metadata,
                    axis1 = 1,
                    axis2 = 2) {
  left_join(
    as.data.frame(t(distatis_comp_boot[, axis1,])) %>% pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = 'Factor 1'
    ),
    as.data.frame(t(distatis_comp_boot[, axis2,])) %>% pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = 'Factor 2'
    )
  ) %>% left_join(., rownames_to_column(metadata, var = "sample"))
}

# 09_deseq_counts.R
get_deseq_counts <- function(dds, species) {
  # needs to be used with imap or equivalent because requires species name in
  # the final formula
  print(paste0("Working on ", species, "..."))
  
  as_tibble(counts(dds, normalized=T), rownames = "Geneid") %>%
    pivot_longer(cols=-Geneid, names_to="sample", values_to = "normcounts") %>%
    group_by(Geneid) %>%
    mutate(zscore=(normcounts-mean(normcounts))/sd(normcounts)) %>%
    ungroup() %>%
    mutate(genome = species)
}

# 09_deseq_constrasts.R
get_deseq_results <- function(dds, species, contrast) {
  # needs to be used with imap or equivalent because requires species name in
  # the final formula
  print(paste0("Working on ", species, "..."))
  
  as_tibble(results(dds, contrast = contrast), rownames = "Geneid") %>%
    mutate(genome = species)
}

# 09_deseq_constrasts.R
get_deseq_results_shrunk <-
  function(dds, species, contrast, contrast2coef = FALSE) {
    # needs to be used with imap or equivalent because requires species name in
    # the final formula
    print(paste0("Working on ", species, "..."))
    
    if (contrast2coef == TRUE) {
      dds[[contrast[1]]] <- relevel(factor(dds[[contrast[1]]]), ref = contrast[3])
      coef <- paste(contrast[1], contrast[2], "vs", contrast[3], sep = "_")
      
      # as recommended here
      # http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators
      dds <- nbinomWaldTest(dds)
      
      print(paste0("Using Coefficient ", coef, "..."))
      
      as_tibble(lfcShrink(dds, coef = coef, type = "apeglm"), rownames = "Geneid") %>%
        mutate(genome = species) %>%
        mutate(padj=if_else(is.na(padj), 0.99, padj),
               coefficient = {{ coef }})
    } else {
      coef <- paste(contrast[1], contrast[2], "vs", contrast[3], sep = "_")
      
      print(paste0("Using Coefficient ", coef, "..."))
      
      as_tibble(lfcShrink(dds, coef = coef, type = "apeglm"), rownames = "Geneid") %>%
        mutate(genome = species) %>%
        mutate(padj=if_else(is.na(padj), 0.99, padj),
               coefficient = {{ coef }})
    }
    
  }

# 12_functional_enrichment.R
enrichfun <- function(deseq_res_split, path2gene, path2name) {
  
  genome    <- unique(deseq_res_split$genome)
  coefficient <- unique(deseq_res_split$coefficient)
  day       <- unique(deseq_res_split$days)
  GENES     <- deseq_res_split$Geneid
  T2G       <- path2gene %>% dplyr::filter(str_detect(Geneid, {{ genome }} ))
  
  print(paste("Working on genome", genome, "with condition", coefficient, "on day", day))
  
  x <- try(enricher(
    GENES,
    TERM2GENE = T2G,
    TERM2NAME = path2name,
    pAdjustMethod = "BH"
  ),
  silent = TRUE)
  
  if (!is.null(x)) {
    mydf <- as_tibble(x) %>%
      dplyr::mutate(enrichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio))) %>%
      dplyr::mutate(genome = {{ genome }},
                    coefficient = {{ coefficient }},
                    days = {{ day }} ) %>%
      tidyr::separate_rows(`geneID`, sep = "/") %>%
      dplyr::rename(Geneid = geneID)
    return(mydf)
  }
}

# Plotting Functions ------------------------------------------------------

mybartheme <- function(...){
  theme(
    panel.spacing.x = unit(0,"line"),
    strip.placement = 'outside',
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    #axis.text.x = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    ...)
}

mytheme <- function(...){
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    ...)
}
