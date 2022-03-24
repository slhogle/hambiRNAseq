### simple modifications of lgpr functions
# 1. plot_posterior_components_sub2
# 2. plot_components
# 3. plot_component

get_posterior_components <- function(fit, PRED = NULL, time_is_xvar = TRUE, n_sd = 2) 
{
  if (class(fit) != "lgpfit") 
    stop("Class of 'fit' must be 'lgpfit'!")
  
  model <- fit@model
  AAA <- lgpr:::PRED_to_arrays(PRED)
  MMM <- AAA$MMM
  SSS <- n_sd * AAA$SSS
  X_test <- PRED$X_test_scaled
  
  GG <- list()
  sum_D <- sum(model@stan_dat$D) + 1
  
  for (d in 1:sum_D) {
    
    gg <- get_component(MMM, SSS, model, d, time_is_xvar, X_test)
    GG[[d]] <- gg
  }
  
  return(GG)
}

get_component <- function(MMM, SSS, model, idx, time_is_xvar, X_test) 
{
  sdat <- model@stan_dat
  
  if (is.null(X_test)) {
    X <- t(sdat$X)
    Xnn <- sdat$X_notnan
  } else {
    X <- X_test
    Xnn <- as.numeric(!is.nan(X_test[, 3]))
  }
  
  D <- sdat$D
  SCL <- model@scalings
  sum_D <- sum(D) + 1
  S <- dim(MMM)[1]
  n <- dim(MMM)[2]
  d <- dim(MMM)[3]
  cpn <- model@info$component_names
  cvn <- model@info$covariate_names
  if (sum_D != d) {
    stop("dim(MMM)[3] must be ", sum_D)
  }
  if (idx < 1 || idx > sum_D) {
    stop("idx must be between 1 and ", sum_D)
  }
  if (idx < sum_D) {
    ctype <- lgpr:::component_index_to_type(D, idx)
    cind <- lgpr:::component_index_to_covariate_index(D, idx)
  } else {
    ctype <- 7
  }
  plot_eb <- !is.null(SSS) && !(ctype %in% c(1, 7))
  MMM_i <- MMM[, , idx]
  f <- as.numeric(t(MMM_i))
  if (plot_eb) {
    UUU_i <- MMM_i + SSS[, , idx]
    LLL_i <- MMM_i - SSS[, , idx]
    ub <- as.numeric(t(UUU_i))
    lb <- as.numeric(t(LLL_i))
  }
  sample <- rep(1:S, each = n)
  id <- rep(X[, 1], S)
  if (ctype == 4) {
    icnt <- idx_to_cont_index(D, idx)
    cscl <- model@scalings$CSCL[[icnt]]
  }
  if (ctype %in% c(1, 2, 5, 6, 7)) {
    time_is_xvar <- TRUE
  }
  if (time_is_xvar) {
    xvar <- rep(SCL$TSCL$fun_inv(X[, 2]), S)
    xvar_name <- model@info$varInfo$time_variable
  }
  else {
    xvar <- rep(X[, cind], S)
    xvar_name <- cvn[cind]
    if (ctype == 4) {
      xvar <- cscl$fun_inv(xvar)
    }
  }
  if(exists("ub")){
    return(data.frame(xvar, f, id, ub, lb))
  }else {
    return(data.frame(xvar, f, id))
  }
}
