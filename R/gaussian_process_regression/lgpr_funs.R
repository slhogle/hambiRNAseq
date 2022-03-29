make_id <- function(.data) {
  .data %>%
    group_by(var, pred, evo, replicate) %>%
    mutate(id = dplyr::cur_group_id()) %>%
    ungroup() %>%
    select(y, var, pred, evo, pred_evo, recoveryonset, day, id) %>%
    arrange(id, day)
}

factorize_my_vars <- function(.data){
  .data %>%
    mutate(
      pred = as.factor(pred),
      evo = as.factor(evo),
      pred_evo = as.factor(pred_evo),
      recoveryonset = as.numeric(recoveryonset),
      day = as.numeric(day),
      id = as.factor(id),
      y = as.numeric(y)
    )
}

get_relavance <- function(fit, measure) {
  
  rel <- lgpr:::relevances.default.all(fit) %>% 
    tibble() %>%
    pivot_longer(everything(), names_to = "variable", values_to = "relevance") 
  
  thresh <- select(fit, threshold = 0.95) %>% 
    rownames_to_column(var="variable") %>% 
    tibble() %>%
    rename(selected = 2)
  
  left_join(rel, thresh) %>%
    mutate(measure = {{ measure }})
}

threshold_density <- function(x) {stats::dbeta(x, 20, 2)}


sumfunctiondata <- function(model_input, lastday, model_fit, measure){
  t <- seq(1, lastday, by = 1)
  x_pred <- new_x(model_input[[ {{ measure }} ]], t, x="day", x_ns = "recoveryonset")
  p <- pred(model_fit, x_pred, reduce = mean)
  input <- lgpr:::plot_pred.create_input(model_fit, p, x_pred, draws, reduce, "id", "day")
  # transform data back to original scale
  data <- lgpr:::create_plot_df(model_fit, x="day") %>%
    mutate(y=10^y) %>% dplyr::rename(obs=y)
  # get ribbon for 95% credible interval 
  rib <- lgpr:::plot_pred.create.df_ribbon(model_fit, p, lgpr:::dollar(input, "df_base"), 2) %>%
    mutate(upper=10^upper,
           lower=10^lower)
  # get line from the fit
  lines <- lgpr:::plot_pred.create.df_line(model_fit, p, lgpr:::dollar(input, "df_base")) %>%
    mutate(y=10^y)
  # cobine for final output
  left_join(lines, rib) %>%
    left_join(., data) %>%
    mutate(measure = {{ measure }}) %>%
    left_join(., distinct(dplyr::select(model_input[[ {{ measure }} ]], pred, evo, pred_evo, id))) %>%
    mutate(measure = {{ measure }})
}


ribbon_line <- function(lgprinput, model_fit, model_pred, pred_input, index, covariate, individual=FALSE){
  left_join(lgpr:::plot_f.create.df_ribbon(model_fit, model_pred, lgpr:::dollar(pred_input, "df_base"), comp_idx=index, 2),
            lgpr:::plot_f.create.df_line(model_fit,   model_pred, lgpr:::dollar(pred_input, "df_base"), comp_idx=index)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(covariate = {{ covariate }})
}


componentdata <- function(model_input, lastday, model_fit, measure){
  t <- seq(1, lastday, by = 1)
  x_pred <- new_x(model_input[[ {{ measure }} ]], t, x="day", x_ns = "recoveryonset")
  p <- pred(model_fit, x_pred, reduce = mean)
  input <- lgpr:::plot_pred.create_input(model_fit, p, x_pred, draws, reduce, "id", "day")
  vars2get <- lgpr:::component_names(model_fit)
  imap_dfr(vars2get, ~ribbon_line(model_input[[ {{ measure }} ]], model_fit, p, input, .y, .x)) %>%
    mutate(measure = {{ measure }})
}

componentdata_fmt <- function(.data){
  if ( unique(.data$covariate) == "gp(day)" ) {
    .data %>%
      dplyr::select(-id) %>% 
      dplyr::distinct()
  } else if ( unique(.data$covariate) == "zs(pred)*gp(day)" ) {
    .data %>%
      left_join(., distinct(dplyr::select(lgprinput[[1]], pred, id))) %>%
      dplyr::select(-id) %>% 
      dplyr::distinct() %>%
      dplyr::rename(condition = pred)
  } else if ( unique(.data$covariate) == "zs(evo)*gp(day)" ) {
    .data %>%
      left_join(., distinct(dplyr::select(lgprinput[[1]], evo, id))) %>%
      dplyr::select(-id) %>% 
      dplyr::distinct() %>%
      dplyr::rename(condition = evo)
  } else if ( unique(.data$covariate) == "zs(pred_evo)*gp(day)" ) {
    .data %>%
      left_join(., distinct(dplyr::select(lgprinput[[1]], pred_evo, id))) %>%
      dplyr::select(-id) %>% 
      dplyr::distinct() %>%
      dplyr::rename(condition = pred_evo)
  } else if ( unique(.data$covariate) == "gp_vm(recoveryonset)" ) {
    .data %>%
      dplyr::select(-id) %>% 
      dplyr::distinct()
  }
}


pretty_log <- function(...) {
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10 ^ x),
    labels = trans_format("log10", scales::math_format(10 ^ .x))
  )
}

# ggplot theme
mytheme <- function(...) {
  theme(
    panel.grid = element_blank(),
    #panel.border = element_blank(),
    panel.background = element_blank(),
    strip.placement = 'outside',
    strip.background = element_blank(),
    #axis.line.x = element_line(color = "black"),
    #axis.line.y = element_line(color = "black"),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    ...
  )
}