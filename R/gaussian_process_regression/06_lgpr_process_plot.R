library(here)
library(tidyverse)
library(patchwork)
library(scales)
library(lgpr)

source(here::here("R", "gaussian_process_regression", ""))


# Read data ---------------------------------------------------------------

lgprinput <- read_rds(here::here("data", "gaussian_process_regression", "formatted_for_cluster.rds"))

fit_atp <- read_rds(here::here("data", "gaussian_process_regression", "lgpr_fit_atp_xz.rds"))
fit_cfu <- read_rds(here::here("data", "gaussian_process_regression", "lgpr_fit_cfu_xz.rds"))
fit_cil <- read_rds(here::here("data", "gaussian_process_regression", "lgpr_fit_cil_xz.rds"))
fit_opd <- read_rds(here::here("data", "gaussian_process_regression", "lgpr_fit_opd_xz.rds"))

fit_list <- list(fit_atp, fit_cfu, fit_cil, fit_opd)
fit_names <- c("atp", "cfu", "cil", "opd")
names(fit_list) <- fit_names

# Quick posterior checks --------------------------------------------------

plot_draws(fit_atp, type = 'dens')
plot_draws(fit_cfu, type = 'dens')
plot_draws(fit_cil, type = 'dens')
plot_draws(fit_opd, type = 'dens')

# Component relevances ----------------------------------------------------

relevance <- imap_dfr(fit_list, ~get_relavance(.x, .y))

# Model posterior prediction: component sum -------------------------------

component_sum <- imap_dfr(fit_list, ~sumfunctiondata(lgprinput, 55, .x, .y)) 

# Model posterior prediction: individual components -----------------------

component_idv <- imap_dfr(fit_list, ~componentdata(lgprinput, 55, .x, .y))

# format for plotting
component_idv_fmt <- component_idv %>%
  group_by(covariate) %>%
  group_split() %>%
  map_df(., componentdata_fmt)

# Plotting ----------------------------------------------------------------

pp2m <- component_sum %>%
  mutate(summary_cat =  case_when((evo == "none" & pred == "none") ~ "ANC_NONE",
                                  (evo == "none" & pred == "predation") ~ "ANC_COEVO",
                                  (evo == "evolution" & pred == "none") ~ "COEVO_NONE",
                                  (evo == "evolution" & pred == "predation") ~ "COEVO_COEVO")) %>%
  mutate(summary_cat = factor(summary_cat, levels=c("ANC_NONE", "COEVO_NONE", "ANC_COEVO", "COEVO_COEVO")),
         measure = factor(measure, levels=c("opd", "cfu", "cil", "atp"))) %>%
  ggplot(aes(x=day, group=id)) + 
  geom_ribbon(aes(ymax=upper, ymin=lower), fill="#C0C0C0", alpha=0.1) +
  geom_jitter(aes(y=obs, color=summary_cat)) + # size=pointsize
  geom_line(aes(y=y, color=summary_cat)) + #, size=linesize
  scale_color_brewer(palette="Paired") +
  facet_wrap( ~ measure, scales="free_y", nrow=1) + 
  labs(x="Day", y="", color="") + 
  pretty_log() + 
  theme_bw() +
  mytheme(legend.position="bottom")

ggsave(plot=pp2m, filename=here::here("figs", "fig2.svg"), device="svg", height = 3, width = 8)


# summed components
p1 <- component_sum %>%
  mutate(summary_cat =  case_when((evo == "none" & pred == "none") ~ "ANC_NONE",
                                  (evo == "none" & pred == "predation") ~ "ANC_COEVO",
                                  (evo == "evolution" & pred == "none") ~ "COEVO_NONE",
                                  (evo == "evolution" & pred == "predation") ~ "COEVO_COEVO")) %>%
  mutate(summary_cat = factor(summary_cat, levels=c("ANC_NONE", "COEVO_NONE", "ANC_COEVO", "COEVO_COEVO")),
         measure = factor(measure, levels=c("opd", "cfu", "cil", "atp"))) %>%
  ggplot(aes(x=day, group=id)) + 
  geom_ribbon(aes(ymax=upper, ymin=lower), fill="#C0C0C0", alpha=0.1) +
  geom_jitter(aes(y=obs, color=summary_cat)) + # size=pointsize
  geom_line(aes(y=y, color=summary_cat)) + #, size=linesize
  scale_color_brewer(palette="Paired") +
  #facet_wrap(~measure, scales="free", ncol=1) + 
  facet_grid(measure ~., scales="free") + 
  labs(x="Day", y="", color="") + 
  pretty_log() + 
  theme_bw() +
  mytheme()

# Relevances
p2 <- relevance %>%
  mutate(var=case_when(str_detect(variable, "id") ~ "ID",
                       str_detect(variable, "recovery") ~ "Recov",
                       str_detect(variable, "noise") ~ "Noise",
                       str_detect(variable, "\\(evo\\)") ~ "Evo",
                       str_detect(variable, "\\(pred\\)") ~ "Pred",
                       str_detect(variable, "pred_evo") ~ "Evo/Pred",
                       TRUE ~ "Day")) %>%
  mutate(var=factor(var, levels=c("Pred", "Evo", "Evo/Pred", "Recov", "Day", "ID", "Noise")),
         measure = factor(measure, levels=c("opd", "cfu", "cil", "atp"))) %>%
  ggplot(aes(x=var, y=relevance, color=selected)) +
  stat_summary(fun.data = "mean_sdl") +
  scale_color_manual(values = c("grey70", "black")) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + #guide_axis(n.dodge = 2)
  facet_grid(measure ~.) + 
  labs(x="", y="") + 
  theme_bw() +
  mytheme()


# individual components
p3 <- component_idv_fmt %>%
  mutate(var=case_when(str_detect(covariate, "id") ~ "ID",
                       str_detect(covariate, "recovery") ~ "Recov",
                       str_detect(covariate, "noise") ~ "Noise",
                       str_detect(covariate, "\\(evo\\)") ~ "Evo",
                       str_detect(covariate, "\\(pred\\)") ~ "Pred",
                       str_detect(covariate, "pred_evo") ~ "Evo/Pred",
                       TRUE ~ "Day")) %>%
  mutate(var=factor(var, levels=c("Pred", "Evo", "Evo/Pred", "Recov", "Day", "ID", "Noise")),
         measure = factor(measure, levels=c("opd", "cfu", "cil", "atp"))) %>%
  ggplot(aes(x=day, y=y, ymax=upper, ymin=lower, fill=condition, color=condition)) + 
  geom_hline(yintercept=0, linetype="dashed") +
  geom_ribbon(alpha=alpha) +
  labs(x="Day", y="") + 
  geom_line(size=linesize) +
  facet_grid(measure~var) +
  theme_bw() +
  mytheme()

# combined plot
figS2 <- p1 + p2 + p3 + 
  plot_layout(guides="collect", widths = c(1, 1, 5)) + 
  plot_annotation(tag_levels = "A") &
  theme(legend.position="bottom")

ggsave(plot=figS2, filename=here::here("figs", "figS2.svg"), device="svg", height = 5, width = 9)

