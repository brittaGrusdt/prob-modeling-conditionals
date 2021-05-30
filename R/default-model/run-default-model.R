source("R/default-model/helpers-tables.R")
source("R/helpers-webppl.R")
source("R/helper-functions.R")
library(rwebppl)
library(tidyverse)
fs = .Platform$file.sep

target="targets_paper_config"
# target="targets_my_config"

# params <- configure(c("bias_none", "pl", target))
# params <- configure(c("bias_none", "speaker", target))
# params <- configure(c("bias_none", "speaker_literal", target))
# params <- configure(c("bias_none", "speaker_p_rooij", target))
##params <- configure(c("bias_none", "speaker_p_rooij_ifac_applicable", target))
# params <- configure(c("bias_none", "speaker_uncertain", target))
params <- configure(c("bias_none", "speaker_certain", target))
# params <- configure(c("priorN", target))

# Setup -------------------------------------------------------------------
if(!dir.exists(params$target_dir)) dir.create(params$target_dir, recursive = TRUE)
params$target <- file.path(
  params$target_dir, paste(params$target_fn, "-", params$fn_suffix, ".rds", sep=""),
  fsep=fs
)
params$target_params <- str_replace(params$target, "results", "params")
params$utts_path <- file.path(params$target_dir, params$utts_fn, fsep=fs)

params$plot_dir = paste(params$target_dir, "figs", sep=fs)
if(!dir.exists(params$plot_dir)) dir.create(params$plot_dir)

##---- Generate/Retrieve tables ----##
if(!"tables_path" %in% names(params)){
  params$tables_path <- paste(params$target_dir, params$tables_fn, sep=fs)
}
# params$seed_tables <- 0907
# params$seed_webppl <- 12345
tables <- create_tables(params)
# tables = readRDS(params$tables_path)
tbls.map = tables %>% select(bn_id, cn.orig)

params$tables = tables %>% ungroup %>%
  dplyr::select(bn_id, cn, ps, vs, ll) %>% group_by(bn_id)

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
if(params$level_max == "speaker") {
  speaker <- posterior$distributions %>%
    structure_speaker_data(params, tbls.map)
  save_data(posterior$all_ids %>% rename(bn_id=value),
            str_replace(params$target, "results", "sample-ids"))
  speaker_avg <- speaker %>% average_speaker(params)
  speaker_avg %>% arrange(desc(avg))
} else if(params$level_max %in% c("priorN", "prior_conditioned")){
    data <- structure_bns(posterior, params)
} else if(params$level_max == "log_likelihood"){
  data <- tibble(id=posterior$id$value, cn=posterior$cn$value,
                 logL=posterior$logL$value)
} else {
  data <- posterior %>% structure_listener_data(params, tbls.map)
  data_voi <- voi_default(data, params)
}
