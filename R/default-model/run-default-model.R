source("R/default-model/helpers-tables.R")
source("R/helpers-webppl.R")
source("R/helper-functions.R")
library(rwebppl)
library(tidyverse)

params <- configure(c("bias_none", "pl", "targets_default"))
# params <- configure(c("bias_none", "speaker", "targets_default"))
# params <- configure(c("bias_none", "speaker_literal", "targets_default"))
# params <- configure(c("bias_none", "speaker_p_rooij", "targets_default"))
# params <- configure(c("bias_none", "speaker_uncertain", "targets_default"))
# params <- configure(c("bias_none", "speaker_certain", "targets_default"))
# params <- configure(c("priorN", "targets_default"))
  # Setup -------------------------------------------------------------------
dir.create(params$target_dir, recursive = TRUE)
params$target <- file.path(params$target_dir, params$target_fn, fsep=.Platform$file.sep)
params$target_params <- file.path(params$target_dir, params$target_params, fsep=.Platform$file.sep)
params$utts_path <- file.path(params$target_dir, params$utts_fn, fsep = .Platform$file.sep)

##---- Generate/Retrieve tables ----##
if(!"tables_path" %in% names(params)){
  params$tables_path <- paste(params$target_dir, params$tables_fn, sep=.Platform$file.sep)
}
tables <- create_tables(params)
# tables = readRDS(params$tables_path)
params$tables = tables %>% ungroup %>%
  dplyr::select(bn_id, cn, ps, vs, ll)

##---- Generate/Retrieve utterances ----##
generate_utts <- function(params){
  utterances <- run_webppl("./model/default-model/utterances.wppl", params)
  utterances <- utterances %>% map(function(x){x %>% pull(value)}) %>% unlist()
  utterances %>% save_data(params$utts_path)
  return(utterances)
}
if(!file.exists(params$utts_path)){
  utterances <- generate_utts(params)
} else {
  utterances <- readRDS(params$utts_path)
  print(paste("utterances read from:", params$utts_path))
}
params$utterances <- utterances
 
# Run Model ---------------------------------------------------------------
posterior <- run_webppl(params$model_path, params)

# structure + save data
params$target <- file.path(paste(params$target_dir, params$target_fn, sep=.Platform$file.sep))
params$target_params <- str_replace(params$target, "results", "params")
params$plot_dir = paste(params$target_dir, "figs", sep=.Platform$file.sep)
if(!dir.exists(params$plot_dir)){
  dir.create(params$plot_dir)
}

# restructure data and save
if(params$level_max == "speaker") {
  speaker <- posterior$distributions %>% structure_speaker_data(params)
  save_data(posterior$all_ids %>% rename(stimulus_id=value),
            paste(params$target_dir, .Platform$file.sep,
                  "sample-ids-", params$target_fn, sep=""))
  speaker_avg <- speaker %>% average_speaker(params) %>% arrange(avg)
  speaker_avg
} else if(params$level_max %in% c("priorN")){
    data <- structure_bns(posterior, params)
} else if(params$level_max == "log_likelihood"){
  data <- tibble(id=posterior$id$value, cn=posterior$cn$value,
                 logL=posterior$logL$value)
} else {
  data <- posterior %>% structure_listener_data(params)
  data_voi <- voi_default(data, params)
}
