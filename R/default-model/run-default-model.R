source("R/default-model/helpers-tables.R")
source("R/helpers-webppl.R")
source("R/helper-functions.R")
library(rwebppl)
library(tidyverse)
fs = .Platform$file.sep

target="targets_paper_config"
# target="targets_my_config"
# params <- configure(c("speaker_uncertain_certain", target))
# params <- configure(c("speaker_uncertain", target))
# params <- configure(c("speaker_certain", target))
# params <- configure(c("speaker_p_rooij", target))

# set seed once
seed = as.numeric(Sys.time())

for(i in seq(1,4)){
  if(i == 1){
    params <- configure(c("pl", target))
  } else if(i==2){
    params <- configure(c("speaker", target))
  } else if(i==3){
    params <- configure(c("speaker_literal", target))
  } else if (i==4){
    params <- configure(c("priorN", target))
  }
  params$seed_webppl = seed
  
  # Setup -------------------------------------------------------------------
  if(!dir.exists(params$target_dir)) dir.create(params$target_dir, recursive=T)
  params$target <- file.path(
    params$target_dir, paste("results-", params$fn_suffix, ".rds", sep=""),
    fsep=fs
  )
  params$target_params <- str_replace(params$target, "results", "params")
  params$utts_path <- file.path(params$target_dir, params$utts_fn, fsep=fs)
  
  params$plot_dir = paste(params$target_dir, "figs", sep=fs)
  if(!dir.exists(params$plot_dir)) dir.create(params$plot_dir)
  
  ##---- Generate/Retrieve utterances ----##
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
    speaker <- posterior %>% structure_speaker_data(params)
    save_data(posterior$all_ids %>% rename(bn_id=value),
              str_replace(params$target, "results", "sample-ids"))
    # speaker_avg <- speaker %>% average_speaker(params)
  } else if(params$level_max %in% c("priorN", "prior_conditioned")){
      data <- structure_bns(posterior, params)
  } else {
    data <- posterior %>% structure_listener_data(params)
    data_voi <- voi_default(data, params)
  }
}
