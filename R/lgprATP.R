.libPaths(c("/projappl/project_2001175/rpackages-r-env-singul-4.0.2", .libPaths()))

library(tidyverse)
library(lgpr)

atp3 <- readRDS("atp3")

my_prior.bb <- list(
  wrp = log_normal(-0.1, 0.2),
  alpha = student_t(20),
  gamma = bet(2, 2)
)

lgpinput <- atp3 %>% 
  dplyr::select(-count, -vallogscale) %>%
  dplyr::rename(y=countred)

fit.bb <- lgp(y ~ gp(day) + zs(id)*gp(day) + zs(pred)*gp(day) + zs(evo)*gp(day) + 
                zs(pred_evo)*gp(day) + gp_vm(recoveryonset),
              lgpinput, 
              likelihood = "bb",
              num_trials = 215,
              prior=my_prior.bb,
              sample_f = TRUE,
              verbose = TRUE,
              refresh = 300,
              chains = 4,
              cores = 4,
              control = list(adapt_delta = 0.95),
              iter = 2000)

saveRDS(fit.bb, "fit.atp.bb")