source("R/default-model/helpers-tables.R")
source("R/helpers-webppl.R")
source("R/helper-functions.R")
library(rwebppl)
library(tidyverse)
fs = .Platform$file.sep
par = expand.grid(alpha=c(1,3,5,10),
                  cost_conditional=c(0, 0.01, 0.1, 0.5), 
                  theta=c(0.9, 0.95, 0.975))
n_iter = nrow(par)
print(paste('# iterations:', n_iter))

params <- configure(c("bias_none", "pl", "targets_paper_config"))
if(!dir.exists(params$target_dir)) dir.create(params$target_dir, recursive = T)
params$target_params <- str_replace(params$target, "results", "params")
params$target <- file.path(params$target_dir, "param-sweep.csv", fsep=fs)

params$utts_path <- file.path(params$target_dir, params$utts_fn, fsep=fs)
##---- Generate/Retrieve tables ----##
if(!"tables_path" %in% names(params)){
  params$tables_path <- paste(params$target_dir, params$tables_fn, sep=fs)
}
tables <- create_tables(params)
# tables = readRDS(params$tables_path)
tbls.map = tables %>% select(bn_id, cn.orig)
params$tables = tables %>% ungroup %>%
  dplyr::select(bn_id, cn, ps, vs, ll) %>% group_by(bn_id)

params$save = FALSE
##---- Generate/Retrieve utterances ----##
utterances <- readRDS(params$utts_path)
print(paste("utterances read from:", params$utts_path))
params$utterances <- utterances

get_vois <- function(data, params){
  df.wide = data %>% ungroup() %>% 
    group_by(bn_id, level) %>% select(-bias) %>% 
    pivot_wider(names_from="cell", values_from="val")

  theta=params$theta
  df = df.wide %>% 
    mutate(p_c=AC+`-AC`, p_a=AC+`A-C`) %>% 
    compute_cond_prob("P(A|C)") %>% rename(p_a_given_c = p) %>% 
    compute_cond_prob("P(-C|-A)") %>% rename(p_nc_given_na = p) %>%
    mutate(uncertainty =
            case_when((p_a<theta & p_a>1-theta) & (p_c<theta & p_c>1-theta) ~ "both")) 
  
  df.uncertain_both = df %>% filter(uncertainty == "both") %>%
    group_by(level) %>% 
    summarize(ev_uncertain_both=sum(prob), .groups="drop_last")
  # expected values
  evs = df %>% group_by(level, cn) %>% 
    transmute(ev_p_rooij = prob*p_rooij, ev_p_delta = prob * p_delta,
              ev_a_given_c = prob * p_a_given_c,
              ev_nc_given_na= prob * p_nc_given_na) %>%
    summarize(ev_p_rooij = sum(ev_p_rooij), ev_p_delta=sum(ev_p_delta),
      ev_a_given_c=sum(ev_a_given_c),
      ev_nc_given_na=sum(ev_nc_given_na), .groups="drop_last") %>%
    pivot_longer(cols=c("ev_p_rooij", "ev_p_delta", "ev_a_given_c", "ev_nc_given_na"),
                 names_to="key", values_to="ev")
  
  evs.cns =  df %>% group_by(level, cn) %>%
    summarize(ev=sum(prob), .groups="drop_last") %>% 
    pivot_wider(names_from="cn", names_prefix="ev_", values_from="ev")
  if(params$level_max == "prior"){
    levels = c("prior")
  } else {
    levels = c("prior", "LL", "PL")
  }
  res = evs %>% filter(level %in% levels) %>% 
    group_by(level, cn) %>% pivot_wider(names_from="key", values_from="ev")
  
  results = left_join(res, evs.cns, by="level") %>% 
    add_column(alpha=params$alpha, cost_conditional=params$cost_conditional)
  results = left_join(results, df.uncertain_both, by="level")
  return(results)
}

# Run Model ---------------------------------------------------------------
results = map_dfr(seq(1, n_iter), function(i){
  params$alpha = par[i,]$alpha
  params$cost_conditional = par[i,]$cost_conditional
  params$theta = par[i,]$theta
  posterior <- run_webppl(params$model_path, params)
  
  # restructure data and save
  data <- posterior %>% structure_listener_data(params, tbls.map)
  vois <- get_vois(data, params)
  return(vois)
});

write_csv(results, params$target)

results %>%
  pivot_longer(cols=starts_with("ev"), names_to="key", values_to="val") %>%
  filter(key=="ev_uncertain_both") %>%
  ggplot(aes(x=alpha, y=val)) +
  geom_point(aes(color=cn), alpha=0.5) +
  ggtitle("ev_uncertain_both") +
  facet_grid(cost_conditional~level)

# results %>%
#   pivot_longer(cols=starts_with("ev"), names_to="key", values_to="val") %>%
#   filter(cn !="A || C" & level == "PL") %>% 
#   ggplot(aes(x=alpha, y=val)) +
#   geom_point(aes(color=cn), alpha=0.5) +
#   ggtitle("dep")
#   facet_grid(cost_conditional~level)

results %>%
  pivot_longer(cols=starts_with("ev"), names_to="key", values_to="val") %>%
  filter(key=="ev_p_rooij") %>% 
  ggplot(aes(x=alpha, y=val)) +
  geom_point(aes(color=cn), alpha=0.5) +
  ggtitle("ev_p_rooij") +
  facet_grid(cost_conditional~level)

results %>%
  pivot_longer(cols=starts_with("ev"), names_to="key", values_to="val") %>%
  filter(key=="ev_A || C") %>% 
  ggplot(aes(x=alpha, y=val)) +
  geom_point(aes(color=cn), alpha=0.5) +
  ggtitle("ev_A || C") +
  facet_grid(cost_conditional~level)

results %>%
  pivot_longer(cols=starts_with("ev"), names_to="key", values_to="val") %>%
  filter(key=="ev_a_given_c") %>% 
  ggplot(aes(x=alpha, y=val)) +
  geom_point(aes(color=cn), alpha=0.5) +
  ggtitle("ev_a_given_c") +
  facet_grid(cost_conditional~level)
