# Functions -----------------------------------------------------------------

get_condition_specific_genes <- function(de, day) {
  
  df1 <- de %>%
    mutate(padj = if_else(is.na(padj), 0.99, padj)) %>%
    filter(padj < 0.1 & abs(log2FoldChange) > 2) %>%
    filter(days== {{ day }}) 
  
  ids1 <- df1 %>%
    group_by(days, coefficient) %>%
    group_split() %>%
    map(., pull, Geneid)
  
  # 1 predation_yes_vs_no
  # 2 pseudomonas_hist_coevolved_vs_ancestral
  # 3 summary_cat_COEVO_COEVO_vs_ANC_COEVO
  # 4 summary_cat_COEVO_NONE_vs_ANC_NONE
  
  # treated vs untreated
  # log2FoldChange = log2(treated/untreated).
  
  # unique to predation
  v1 <- setdiff(ids1[[1]], c(ids1[[2]], ids1[[3]], ids1[[4]]))
  
  # unique to evolution
  v2 <- setdiff(ids1[[2]], c(ids1[[1]], ids1[[3]], ids1[[4]]))
  
  # unique to evolution with predation
  v3 <- setdiff(ids1[[3]], c(ids1[[1]], ids1[[2]], ids1[[4]]))
  
  # unique to evolution without predation
  v4 <- setdiff(ids1[[4]], c(ids1[[1]], ids1[[2]], ids1[[3]]))
  
  bind_rows(
    df1 %>% filter(coefficient == "predation_yes_vs_no") %>% filter(Geneid %in% v1),
    df1 %>% filter(coefficient == "pseudomonas_hist_coevolved_vs_ancestral") %>% filter(Geneid %in% v2),
    df1 %>% filter(coefficient == "summary_cat_COEVO_COEVO_vs_ANC_COEVO") %>% filter(Geneid %in% v3),
    df1 %>% filter(coefficient == "summary_cat_COEVO_NONE_vs_ANC_NONE") %>% filter(Geneid %in% v4)
  )
  
}

# order rows hierarchical clustering for plotting
diff_subset <- function(x, day, coeff, genomecutoff=5) {
  x %>%
    mutate(padj = if_else(is.na(padj), 0.99, padj)) %>%
    filter(padj < 0.1 & abs(log2FoldChange) > 2) %>%
    filter(days == {{ day }} & coefficient == {{ coeff }}) %>%
    group_by(genome) %>%
    add_count() %>%
    ungroup() %>%
    filter(n >= {{ genomecutoff }})
}

hyper_tests_subset <- function(x, y, day, coeff, genomes) {
  left_join(x, y, by=c("Geneid")) %>%
    filter(days == {{ day }} & coefficient == {{ coeff }}) %>%
    filter(genome %in% genomes) %>%
    group_by(Geneid, days, coefficient) %>%
    filter(enrichFactor == max(enrichFactor)) %>%
    add_count(Geneid) %>%
    ungroup() %>%
    select(Geneid, Description)
}

gene_counts_subset <- function(x, y, day) {
  x %>%
    filter(days == {{ day }}) %>%
    select(Geneid, sample, zscore) %>%
    pivot_wider(id_cols="Geneid", names_from="sample", values_from="zscore") %>%
    semi_join(., y) %>%
    arrange(Geneid) 
}

geneorder <- function(x) { 
  x %>%
    separate(Geneid, into = c("genome", "num"), sep="_", remove = FALSE) %>%
    group_by(genome) %>%
    group_split() %>%
    map(., select, -genome, -num) %>%
    map(., column_to_rownames, var="Geneid") %>%
    map(., dist) %>%
    map(., hclust) %>%
    map(., function(x) x$labels[x$order]) %>%
    flatten_chr()
}

row_anns <- function(x, y) {
  left_join(x, y) %>%
    separate(Geneid, into = c("genome", "num"), sep="_", remove = FALSE) %>%
    select(Geneid, genome, Description) %>%
    mutate(Description=if_else(is.na(Description), "none", Description)) %>%
    column_to_rownames(var="Geneid")
}

get_palette <- function(row_anns) {
  mylength <- length(unique(row_anns$Description))
  with_seed(1236,
            pal <- unname(createPalette(mylength, c("#F3874AFF", "#FCD125FF"), M=5000)))
  names(pal) <- unique(row_anns$Description)
  pal["none"] <- NA
  pal
}

customheatmap <- function(day, coeff, deseqres, clustprofilres,
                          fun_anns, norm_counts, metadata, 
                          genomecolors, sampleorder, limlo, limhi) {
  
  dexp     <- diff_subset(deseqres, day, coeff, genomecutoff=5)
  genomes  <- pull(distinct(select(dexp, genome)), genome)
  enriched <- hyper_tests_subset(clustprofilres, fun_anns, day, coeff, genomes)
  zscores  <- gene_counts_subset(norm_counts, dexp, day)
  genes    <- geneorder(zscores)
  
  zscores.mat <- zscores %>%
    arrange(match(Geneid, genes)) %>%
    column_to_rownames(var="Geneid") %>%
    as.matrix()
  
  # row annotations
  my_sample_row <- row_anns(zscores, enriched)
  # column annotations
  my_gene_col <- select(metadata, pseudomonas_hist, predation) 
  my_gene_col <- my_gene_col[sampleorder, ]
  
  # colors
  my_color = list(
    predation = c(yes = "#b2df8a", no = "#a6cee3"),
    pseudomonas_hist = c(ancestral = "#cccccc", coevolved = "#969696"),
    genome = genomecolors,
    Description=get_palette(my_sample_row)
  )
  
  pheatmap(zscores.mat[, sampleorder],
           annotation_colors = my_color,
           annotation_col = my_gene_col,
           annotation_row = my_sample_row,
           #cutree_cols = 2,
           cluster_rows=FALSE,
           cluster_cols=FALSE,
           show_rownames=FALSE,
           #cellwidth=1,
           #cellheight=0.1,
           breaks = seq(from=limlo, to=limhi, length.out=100),
           border_color = NA,
           silent=TRUE,
           color=colorRampPalette(diverging_hcl(7, palette = "Blue-Red 3"))(100))
}

plotlfc <- function(a, b, c){
  left_join(a, b) %>%
    filter(Description %in% c) %>%
    mutate(Description=factor(Description, levels=c)) %>%
    ggplot(aes(x=Description, y=log2FoldChange, color=genome, size=baseMean)) +
    geom_jitter(width=0.2, height=0.5) +
    geom_hline(yintercept=0) +
    scale_color_manual(values=mycols19) +
    labs(x="Log Fold Change", y="", size="Normalized mean") + 
    scale_size_continuous(breaks=c(10, 100, 1000, 10000), range=c(1,3)) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2)) + 
    theme_bw() +
    mytheme()
}

