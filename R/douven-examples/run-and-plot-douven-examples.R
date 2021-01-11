library(tidyverse)
library(ggplot2)
library(rwebppl)
source("R/helper-functions.R")
source("R/helpers-webppl.R")

# Run Model ---------------------------------------------------------------
get_douven_data_evs <- function(params){
  params$save <- FALSE
  data <- run_webppl(params$model_path, params) %>%
    structure_listener_data(params) %>% select(-bias) %>%
    group_by(cn, cell, level) %>% summarise(ev=sum(prob*val))
  data %>% save_data(params$target)
  params %>% save_data(params$target_params)
  return(data)
}

listener_beliefs <- function(data, params){
  data <- data %>% filter(level==params$level_max) %>%
    group_by(cn, cell) 
  if(!is.na(params$condition_on)){
  data <- data %>% filter_vars(params$condition_on) %>% filter(keep) %>% 
    # summarise(ev=sum(ev), .groups="drop_last") %>%
    # mutate(ev=ev/sum(ev)) # marginalize
    # group_by(cn) %>%
    ungroup() %>%
    mutate(ev=ev/sum(ev)) %>%  select(-keep)
  }
  return(data %>% mutate(level="trust"))
}


# 1. Skiing ---------------------------------------------------------------
params.skiing <- configure(c("skiing"))
data.skiing <- get_douven_data_evs(params.skiing)
skiing <- bind_rows(data.skiing, listener_beliefs(data.skiing, params.skiing)) 
skiing.marginal = skiing %>% filter_vars(c("E")) %>% filter(keep) %>% 
  select(-keep) %>% group_by(cn, level) %>%
  summarise(ev=sum(ev), .groups="drop_last") %>% 
  add_column(marginal="e", example="skiing")

# 2. Sundowners -----------------------------------------------------------
params.sundowners <- configure(c("sundowners"))
data.sundowners <- get_douven_data_evs(params.sundowners)
sundowners <- bind_rows(data.sundowners, listener_beliefs(data.sundowners, params.sundowners)) 

sundowners.rain = sundowners %>% filter_vars(c("R")) %>% filter(keep) %>% 
  group_by(cn, level) %>% summarise(ev=sum(ev), .groups="drop_last") %>%
  add_column(marginal="r")

sundowners.rain_sundowner = sundowners %>% filter_vars(c("R", "S")) %>%
  filter(keep) %>% 
  group_by(cn, level) %>% summarise(ev=sum(ev), .groups="drop_last") %>%
  add_column(marginal="rs")

sundowners.data <- bind_rows(sundowners.rain, sundowners.rain_sundowner) %>%
  add_column(example="sundowners")

# Joint sundowners+skiing data
douven.data <- bind_rows(sundowners.data, skiing.marginal)

plot_douven_cases <- function(data){
  p <- data %>% mutate(level=as.factor(level)) %>%
    ggplot() +
    geom_bar(mapping = aes(y=level, x=ev, fill=cn), stat="identity", position="stack") +
    facet_wrap(~example) + 
    facet_wrap(~marginal, labeller=
                 as_labeller(c(`r`="Sundowners: P(R)", `rs`="Sundowners: P(R,S)",
                               `e`="Skiing: P(E)"))) +
    scale_y_discrete(
      name = "",
      limits = c("trust", "PL", "LL", "prior"),
      labels = c(
        paste(strwrap("Listener's beliefs", width=20), collapse="\n"),
        paste(strwrap("Pragmatic interpretation", width=20), collapse="\n"),
        paste(strwrap("Literal interpretation", width=20), collapse="\n"),
        "Prior Belief"
      )) +
    scale_fill_discrete(name="causal net",
                        limits=c("R > W > S", "R||S", "E || S>C", "E>S>C"),
                        labels=c("R->W->S", "R,S indep.", "E indep. S, S->C", "E->S->C")) +
    labs(x="Expected value", title="") +
    theme_bw(base_size=25) + theme(legend.position="bottom")
  
  return(p)
}

p <- douven.data %>% plot_douven_cases()
ggsave(here("data", "douven-examples", "douven-cases.png"), p, width=15, height=6)

  

