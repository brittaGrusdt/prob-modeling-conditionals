source("R/default-model/helpers-tables.R")
source("R/helpers-webppl.R")
source("R/helper-functions.R")
library(rwebppl)
library(tidyverse)
fs = .Platform$file.sep
par = expand.grid(alpha=c(1,3,5,10), theta=c(0.9, 0.95, 0.975))
n_iter = nrow(par)
print(paste('# iterations:', n_iter))

params <- configure(c("pl", "targets_paper_config"))
if(!dir.exists(params$target_dir)) dir.create(params$target_dir, recursive = T)
params$target_params <- str_replace(params$target, "results", "params")
params$target <- file.path(params$target_dir, "param-sweep.csv", fsep=fs)
params$utts_path <- file.path(params$target_dir, params$utts_fn, fsep=fs)
params$save = FALSE
params$utterances <- readRDS(params$utts_path)
print(paste("utterances read from:", params$utts_path))

# Run Model ---------------------------------------------------------------
# 1. Prior + Listeners
results = map_dfr(seq(1, n_iter), function(i){
  params$alpha = par[i,]$alpha
  params$theta = par[i,]$theta
  posterior <- run_webppl(params$model_path, params)
  
  # restructure data and save
  data <- posterior %>% structure_listener_data(params)
  data.cp <- data_cp_plots(params, data) %>% 
    add_column(iter=i, theta=params$theta, alpha=params$alpha)
  return(data.cp)
});
write_csv(results, params$target)

results %>% plot_cp_probs() + facet_wrap(alpha~theta)
results %>% filter(val=="cns") %>% 
  plot_cp_cns(labels.cp) + facet_wrap(alpha~theta)


# 2. Speaker



