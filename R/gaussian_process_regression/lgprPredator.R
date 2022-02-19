.libPaths(c("/projappl/project_2001175/rpackages-r-env-singul-4.0.2", .libPaths()))

library(tidyverse)
library(lgpr)

my_prior.g <- list(
  wrp = log_normal(-0.1, 0.2)
)

lgpinput <- readRDS("cil3") %>% 
  dplyr::select(-count, -countred, -valsqrt, -vallog, -pred, -pred_evo) %>%
  dplyr::rename(y=vallogscale)

fit.g <- lgp(y ~ gp(day) + zs(id)*gp(day) + zs(evo)*gp(day) + gp_vm(recoveryonset),
             lgpinput,
             likelihood = "gaussian",
             prior = my_prior.g,
             verbose = TRUE,
             refresh = 500,
             chains = 4,
             cores = 4,
             control = list(adapt_delta = 0.95),
             iter = 2000)

saveRDS(fit.g, "fit-cil-zs-gs.rds")