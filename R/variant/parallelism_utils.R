# Sample from multinomial -------------------------------------------------

# sample from multinomial the total number of mutations with prob
# weighted the same as the n_population/n_total
samp_sizes = function(nmutations, population_num, population_probs){
  rsizes = as.list(rmultinom(n = 1, size = nmutations, prob = population_probs))
  names(rsizes) = LETTERS[seq( from = 1, to = population_num )]
  return(rsizes)
}

# GSCORES -----------------------------------------------------------------

G_score = function(gene_df){
  Ls = gene_df$l_i
  ns = gene_df$n_i
  
  Ltot = sum(Ls)
  ntot = sum(ns)
  
  ps = Ls/Ltot
  gs = ns * log(ns/(ntot*ps))
  
  observed_G = sum(gs, na.rm=T)/sum(ns, na.rm=T)
  return(observed_G)
}

# uses m_exp from num total mutations / num genes. Same as
# above just expressed a different way
G_score_1 = function(gene_df){

  L_avg   = mean(gene_df$l_i)
  n_tot   = sum(gene_df$n_i)
  n_genes = nrow(gene_df)
  
  m_exp   = n_tot / n_genes
  m_i = gene_df$n_i * (L_avg / gene_df$l_i)
  g_score = gene_df$n_i * log(m_i / m_exp)
  observed_G = sum(g_score, na.rm=T)/sum(gene_df$n_i, na.rm=T)
  return(observed_G)
}

# SIMULATION --------------------------------------------------------------

# simulate number of mutations per replicate using probabilities observed from data
simmultinom = function(nmutations, population_num, population_probs, gene_probs, gene_names) {
  samp_l = samp_sizes(nmutations, population_num, population_probs)
  cbind(LOCUS_TAG=gene_names, map_df(samp_l, rmultinom, n = 1, prob = gene_probs))
}

# prepare a tibble for input to G_score function
G_prep = function(nmutations, population_num, population_probs, gene_probs, gene_names, gene_df){
  simmultinom(nmutations,
              population_num,
              population_probs,
              gene_probs,
              gene_names) %>%
    pivot_longer(A:last_col(), names_to = "REPLICATE") %>%
    filter(value > 0) %>%
    group_by(LOCUS_TAG) %>%
    summarize(n_i = sum(value), n_populations = n_distinct(REPLICATE)) %>%
    ungroup() %>%
    full_join(., lns, by = "LOCUS_TAG")  %>%
    rename(l_i = length) %>%
    arrange(LOCUS_TAG) %>%
    mutate(across(everything(), ~ replace_na(.x, 0)))
}


# GENE SPECIFIC -----------------------------------------------------------

from_parallelism_statistics = function(gene_df) {
  list(
    "gene_names" = gene_df$LOCUS_TAG,
    "Ls" = gene_df$l_i,
    "ns" = gene_df$n_i,
    "Ltot" = sum(gene_df$l_i),
    "ntot" = sum(gene_df$n_i),
    "ps" = gene_df$l_i / sum(gene_df$l_i),
    "expected_ns" = sum(gene_df$n_i) * (gene_df$l_i / sum(gene_df$l_i))
    
  )
}


calculate_parallelism_qvalues = function(gene_df) {
  fps = from_parallelism_statistics(gene_df)
  # probability of getting ns-0.5 multiplicity in ntot number of trials
  # with the expected sampling probabilities ps
  
  # survival function = 1-CDF
  # using binomial distribution function
  pvalues = (1 - pbinom(fps$ns - 0.5, fps$ntot, fps$ps))
  
  # p.adjust and qvalues are slightly different things: see https://www.biostars.org/p/128931/
  # qvalues = p.adjust(pvalues, method = "BH")
  qvalues = qvalue(p = pvalues)
  qvalues.v = qvalues[['qvalues']]
  names(qvalues.v) = fps$gene_names
  return(qvalues.v)
}

calculate_parallelism_logpvalues = function(gene_df, ethresh=1e-30) {
  fps = from_parallelism_statistics(gene_df)
  
  # survival function = 1-CDF
  survivals = (1 - ppois(fps$ns-0.1, fps$expected_ns))
  # for the case if there are so many mutations in some gene that its survival is 0
  survivals[survivals==0] = min(survivals[survivals>0])/1e10
  logpvalues = rep(0, each = length(survivals))
  logpvalues[survivals > ethresh] = -log(survivals[survivals > ethresh])
  logpvalues[survivals <= ethresh] = -(fps$ns * log(fps$ns / fps$expected_ns +
                                                  (fps$ns == 0)) + fps$ns - fps$expected_ns)[survivals <= ethresh]
  names(logpvalues) = fps$gene_names
  return(logpvalues)
}

calculate_unnormalized_survival_from_vector = function(xs, min_x=NULL, max_x=NULL, min_p=1e-10){
  if (is.null(min_x)){
    min_x = min(xs)-1
  }
  if (is.null(max_x)){
    max_x = max(xs)+1
  }
  
  unique_xs = unique(xs)
  unique_xs = sort(append(unique_xs, c(min_x, max_x)))
  
  num_observations = map_dbl(unique_xs, ~ sum(xs>=.x))
  
  # So that we can plot CDF, SF on log scale
  num_observations[0] = num_observations[0] - min_p
  num_observations[1] = num_observations[1] - min_p
  num_observations[-1] = num_observations[-1] + min_p 
  
  return(list(unique_xs, num_observations))
}

calculate_poisson_log_survival = function(ns, expected_ns, ethresh){
  # survival function = 1-CDF
  survivals = (1 - ppois(ns-0.1, expected_ns))
  # for the case if there are so many mutations in some gene that its survival is 0
  #survivals[survivals==0] = min(survivals[survivals>0])/1e10
  logpvalues = rep(0, each = length(survivals))
  logpvalues[survivals > ethresh] = -log(survivals[survivals > ethresh])
  logpvalues[survivals <= ethresh] = -(ns * log(ns / expected_ns +
                                                  (ns == 0)) + ns - expected_ns)[survivals <= ethresh]
  return(logpvalues)
}

map_pval_and_probs = function(observed_p, logpvalues, probabilities, ns){
  sum(map2_dbl(logpvalues, probabilities, ~ sum((.x>=observed_p)*(ns >= nmin)*.y)))
}

NullGeneLogpSurvivalFunction = function(gene_df, ns = rep(0:399), alpha=0.05){
  fps = from_parallelism_statistics(gene_df)
  logpvalues = map(fps$expected_ns, calculate_poisson_log_survival, ns=ns, ethresh=1e-20)
  logprobabilities = map(fps$expected_ns, ~ ns*log(.x)-lgamma(ns+1)-.x)
  probabilities = map(logprobabilities, exp)
  survivals = map_dbl(observed_ps, map_pval_and_probs, logpvalues, probabilities, ns)
  threshold_idx = which(survivals/observed_pvalue_survival <= alpha)[1]
  pstar = observed_ps[threshold_idx] # lowest value where this is true
  num_significant = observed_pvalue_survival[threshold_idx]
  return(list("survivals"=survivals, "pstar"=pstar, "num_significant"=num_significant))
}


# ggplot theme
mytheme <- function(...){
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank()
  )
}

