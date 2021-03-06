library(tidyverse)
library(ggplot2)
library(rwebppl)
library(here)
source("R/helper-functions.R")
source("R/helpers-webppl.R")

# Run Model ---------------------------------------------------------------
get_douven_data_evs <- function(params){
  if(!dir.exists(params$target_dir)) {
    dir.create(params$target_dir, recursive=TRUE)
  }
  data <- run_webppl(params$model_path, params) %>%
    structure_listener_data(params) %>% select(-bias) %>%
    group_by(cn, cell, level) %>% summarise(ev=sum(prob*val))
  return(data)
}

listener_beliefs <- function(data, params){
  data <- data %>% filter(level==params$level_max) %>%
    group_by(cn, cell) 
  if(!is.na(params$condition_on)){
  data <- data %>% filter_vars(params$condition_on) %>%
    # filter(keep) %>% ungroup() %>%
    group_by(keep) %>% 
    mutate(ev=ev/sum(ev), ev=case_when(keep ~ ev, TRUE ~ 0)) %>%
    ungroup() %>% select(-keep)
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

# 2. Garden Party ---------------------------------------------------------
params.gp <- configure(c("garden_party"))
data.gp <- get_douven_data_evs(params.gp)
garden_party <- bind_rows(data.gp, listener_beliefs(data.gp, params.gp)) 
garden_party.marginal = garden_party %>% filter_vars(c("D")) %>%
  filter(keep) %>% 
  select(-keep) %>% group_by(cn, level) %>%
  summarise(ev=sum(ev), .groups="drop_last") %>% 
  add_column(marginal="d", example="gardenParty")

garden_party %>% filter(level=="trust") %>% filter_vars(c("D")) %>%
  filter(keep) %>% 
  select(-keep) %>% group_by(cn, level) %>%
  summarise(ev=sum(ev), .groups="drop_last") %>% 
  add_column(marginal="d", example="gardenParty")

# 3. Sundowners -----------------------------------------------------------
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

# Joint sundowners+skiing+garden party data
douven.data <- bind_rows(sundowners.data, skiing.marginal, garden_party.marginal)

plot_douven_cases <- function(data){
  df = tribble(~limits, ~labels, ~example,
               "R > W > S", "R->W->S", "sundowners",
               "R||S", "R,S indep.", "sundowners", 
               "E || S>C", "E indep. S, S->C", "skiing",
               "E>S>C", "E->S->C", "skiing",
               "D || G>-S", "D indep. G, G->¬S","gardenParty", 
               "D>G>-S", "D->G->¬S", "gardenParty")
  ex.in = data$example %>% unique()
  df = df %>% filter(example %in% ex.in)
  p <- data %>% mutate(level=as.factor(level)) %>%
    ggplot() +
    geom_bar(mapping = aes(y=level, x=ev, fill=cn), stat="identity", position="stack") +
    facet_wrap(~example) + 
    facet_wrap(~marginal, labeller=
                 as_labeller(c(`r`="Sundowners: P(R)", `rs`="Sundowners: P(R,S)",
                               `e`="Skiing: P(E)", `d`="Garden Party: P(D)"))) +
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
                        limits = df$limits,
                        labels = df$labels
                      ) +
    labs(x="Expected value", title="") +
    theme_bw() + theme(legend.position="bottom")
  
  return(p)
}

p <- douven.data %>% filter(example != "gardenParty") %>% 
  plot_douven_cases()
ggsave(here("data", "douven-examples", "skiing-sundowners.png"), p, width=6, height=4)

p <- douven.data %>% filter(example == "gardenParty") %>% 
  plot_douven_cases()
ggsave(here("data", "douven-examples", "gardenParty.png"), p, width=6, height=4)

