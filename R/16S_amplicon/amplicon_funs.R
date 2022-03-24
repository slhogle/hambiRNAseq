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


# ggplot theme
mytheme = function() {
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

# opposite of %in% fuction
`%nin%` = Negate(`%in%`)

# logit transform
logit = function(x){
  log(x/(1-x))
}

minnz = function(V) {
  # Calculates the smallest value of the vector except for 0 (non-zero minumum)
  # Argument: vector
  C <- NULL        # prepare
  k <- length(V)   # count to
  for (i in 1:k) { # check all
    if ((V[i] == 0) == FALSE) (C[i] <- V[i]) else (C[i] <- 9999919) # if V[i] is not 0, add it to C
  }
  m <- min(C)               # minimum of V, not counting 0
  if (max(V) == 1) (m <- 1) # fix for binary vectors (0,1)
  if (m == 9999919) (warning("Error: Minimum calculation failed."))  # warning because of hard-coded replacement
  return(m)
}

quibble95 = function(x, q = c(0.025, 0.5, 0.975)) {
  tibble(x = quantile(x, q), quantile = c("q2.5", "q50", "q97.5"))
}

grouped_bbdml_quantiles95 = function (x, B = 1000, ...) {
  mod = x
  
  M = mod$M
  W = mod$W
  
  ymin = ymax = rep(NA, length(M))
  sims = matrix(NA, nrow = B, ncol = length(W))
  newdat = mod$dat
  
  for (i in 1:B) {
    sim = simulate(mod, nsim = length(W))
    newdat$W = sim
    refit = suppressWarnings(bbdml(mod$formula, phi.formula = mod$phi.formula, 
                                   link = mod$link, phi.link = mod$phi.link, inits = mod$inits, 
                                   data = newdat))
    sims[i, ] = simulate(refit, nsim = length(W))
  }
  
  predall = (t(sims) / M)
  rownames(predall) = rownames(mod$dat)
  
  df = data.frame(cbind(predall, mod$dat)) %>%
    pivot_longer(-W:-last_col()) %>%
    group_by(days, pseudomonas_hist, predation) %>%
    summarize(quibble95(value))
  
  return(df) 
}
