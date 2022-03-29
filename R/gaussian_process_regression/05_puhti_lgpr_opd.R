.libPaths(c("/projappl/project_2005777/rpackages-r-env-singul-4.0.2", .libPaths()))

library(tidyverse)
library(lgpr)

my_prior <- list(
  wrp = log_normal(-0.1, 0.2)
)

mylist <- read_rds("input/formatted_for_cluster.rds") 

lgpinput <- mylist[["opd"]]

fit <- lgp(y ~ gp(day) + zs(id)*gp(day) + zs(pred)*gp(day) + zs(evo)*gp(day) +
               zs(pred_evo)*gp(day) + gp_vm(recoveryonset),
           lgpinput,
           likelihood = "gaussian",
           prior = my_prior,
           verbose = TRUE,
           refresh = 500,
           chains = 4,
           cores = 4,
           control = list(adapt_delta = 0.95),
           iter = 2000)

# write out with compression since files are large
write_rds(fit, "output/lgpr_fit_opd_xz.rds", "xz", compression = 9L)