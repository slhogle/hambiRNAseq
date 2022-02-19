.libPaths(c("/projappl/project_2001175/rpackages-r-env-singul-4.0.2", .libPaths()))

library(tidyverse)
library(lgpr)

atp3 <- readRDS("atp3")

my_prior.g <- list(
  wrp = log_normal(-0.1, 0.2)
)

lgpinput <- atp3 %>% 
  dplyr::select(-count, -countred) %>%
  dplyr::rename(y=vallogscale)

fit.g <- lgp(y ~ gp(day) + zs(id)*gp(day) + zs(pred)*gp(day) + zs(evo)*gp(day) + 
               zs(pred_evo)*gp(day) + gp_vm(recoveryonset),
             lgpinput, 
             likelihood = "gaussian",
             prior = my_prior.g,
             verbose = TRUE,
             refresh = 500,
             chains = 4,
             cores = 4,
             control = list(adapt_delta = 0.95),
             iter = 2000)

saveRDS(fit.g, "fit-atp-zs-gs.rds")

fit.g1 <- lgp(y ~ gp(day) + zs(id)*gp(day) + categ(pred)*gp(day) + categ(evo)*gp(day) + 
             categ(pred_evo)*gp(day) + gp_vm(recoveryonset),
             lgpinput, 
             likelihood = "gaussian",
             prior = my_prior.g,
             verbose = TRUE,
             refresh = 500,
             chains = 4,
             cores = 4,
             control = list(adapt_delta = 0.95),
             iter = 2000)

saveRDS(fit.g1, "fit-atp-cat-gs.rds")