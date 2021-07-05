labels.cp = list(
  "A,C indep." = "A,C indep.",
  "A -> C" = expression(A %->% C),
  "A -> ¬C" = expression(A%->%~"¬C"),
  "C -> A" = expression(C %->% A),
  "C -> ¬A" = expression(C%->%~"¬A"),
  # "A,C indep." = "A,C indep.",
  "A implies C" = expression(A %->% C),
  "A implies ¬C" = expression(A%->%~"¬C"),
  "C implies A" = expression(C %->% A),
  "C implies ¬A" = expression(C%->%~"¬A")
)

save_data <- function(data, target_path){
  data %>% write_rds(target_path)
  print(paste("saved to:", target_path))
}

filter_vars <- function(df_long, vars){
  df <- df_long %>% mutate(keep=TRUE)
  for(var in vars){
    if(str_detect(var, "^-")){
      # negative variable, starts with -
      df <- df %>% mutate(keep=case_when(!keep ~ keep, TRUE ~ str_detect(cell, var)))
    }
    else {
      token <- paste("-", var, sep="")
      df <- df %>% mutate(keep=case_when(!keep ~ keep, TRUE ~ !str_detect(cell, token)))
    }
  }
  return(df)
}

# Utterances --------------------------------------------------------------
generate_utts <- function(params){
  utterances <- run_webppl(
    here("model", "default-model", "utterances.wppl", sep=.Platform$file.sep),
    params
  )
  utterances <- utterances %>% map(function(x){x %>% pull(value)}) %>% unlist()
  if(params$save) utterances %>% save_data(params$utts_path)
  return(utterances)
}

# instead of all different utterances, chunk them into categories (for plotting)
chunk_utterances <- function(data, utts_kept=c()){
  levels = c("likely + literal", "conditional", "literal", "conjunction");
  s = paste(utts_kept, collapse="");
  if(str_detect(s, ">") || str_detect(s, "if")){
    levels = c("likely + literal", "other conditional", "literal", "conjunction",
               utts_kept);
  }
  data = data %>% mutate(
    utterance = case_when(
      utterance %in% utts_kept ~ utterance,
      startsWith(utterance, "likely") ~ "likely + literal",
      str_detect(utterance, ">") ~ levels[[2]],
      str_detect(utterance, "and") ~ "conjunction",
      TRUE ~ "literal"
    ),
    utterance = str_replace_all(utterance, "-", "¬"),
    utterance = str_replace(utterance, ">", "->"),
    utterance = factor(utterance, levels=
                         c(map(utts_kept, function(s){
                           s <- str_replace_all(s, "-", "¬")
                           return(str_replace(s, ">", "->"))
                         }),
                         levels)
    )
  );
  return(data)
}

chunk_cns <- function(data) {
  data = data %>% mutate(cn = case_when(cn == "A || C" ~ "A,C independent",
                                        TRUE ~ "A,C dependent"),
                         cn = factor(cn))
  return(data)
}

#@arg tables: with columns AC, A-C, -AC, -A-C
table_to_utts = function(tables, theta){
  tbls = tables %>% add_probs() %>%
    mutate(literal.a=case_when(p_a >= theta ~ "A",
                               p_na >= theta ~ "-A",
                               T ~ ""),
           literal.c = case_when(p_c >= theta ~ "C",
                                 p_nc >= theta ~ "-C",
                                 T~"")) %>%
    pivot_longer(cols=starts_with("literal"), names_to="literal.tmp",
                 values_to="literal") %>%
    
    mutate(conjunction=case_when(AC >= theta ~ "AC",
                                 `A-C` >= theta ~ "A-C",
                                 `-AC` >= theta ~ "-AC", 
                                 `-A-C` >= theta ~ "-A-C", 
                                 T ~ ""),
           likely.a=case_when(p_a >= 0.5 ~ "likely A",
                            p_na >= 0.5 ~ "likely -A",
                            T ~ ""),
           likely.c=case_when(p_c >= 0.5 ~ "likely C",
                              p_nc >= 0.5 ~ "likely -C",
                              T ~ "")
           ) %>%
    dplyr::select(-`literal.tmp`, -p_a, -p_na, -p_c, -p_nc) %>% 
    pivot_longer(cols=starts_with("likely"), names_to="likely.tmp",
                 values_to="likely") %>%
    dplyr::select(-`likely.tmp`) %>% 
    
    mutate(conditional.ca = case_when(p_c_given_a >= theta ~ "A > C", T ~ ""),
           conditional.cna = case_when(p_c_given_na >= theta ~ "-A > C", T ~ ""),
           conditional.nca = case_when(p_nc_given_a >= theta ~ "A > -C", T ~ ""),
           conditional.ncna = case_when(p_nc_given_na >= theta ~ "-A > -C", T ~ ""),
           conditional.ac = case_when(p_a_given_c >= theta ~ "C > A", T ~ ""),
           conditional.anc = case_when(p_a_given_nc >= theta ~ "-C > A", T ~ ""),
           conditional.nac = case_when(p_na_given_c >= theta ~ "C > -A", T ~ ""),
           conditional.nanc = case_when(p_na_given_nc >= theta ~ "-C > -A", T ~ "")
    ) %>% 
    pivot_longer(cols=starts_with("conditional"), names_to="conditional.tmp",
                 values_to="conditional") %>% 
    dplyr::select(bn_id, cn, conjunction, conditional, literal, likely,
                  AC, `A-C`, `-AC`, `-A-C`)
  return(tbls)
}

# Probabilities -----------------------------------------------------------
#@arg vars: list of variables, if more than one, only states where all hold
# are retained
# @arg data: in long format, such that cell is one column and marginals can
# be computed for any cell entries
# @return: in wide format
marginalize <- function(data, vars){
  df <- data %>% filter_vars(vars)
  df <- df %>%  mutate(p=case_when(keep ~ val, TRUE ~ 0)) %>%
          group_by(bn_id, level) %>% mutate(p=sum(p))  %>%
          select(-keep) %>% spread(key=cell, val=val, fill = 0)
        
  return(df)
}

# takes the expected value of column 'p' with probability in column 'prob'
# @args:
#   df_wide; tibble with one bn per row, at least columns: p, prob, level
#   value_str: str describing value, e.g. 'P(A)' for expected val of P(A)
expected_val <- function(df_wide, value_str){
  evs <- df_wide %>% mutate(ev_prod=p * prob)
  evs <- evs %>% group_by(level)
  evs <- evs %>% summarise(ev=sum(ev_prod), .groups="drop") %>% add_column(p=value_str) %>% ungroup()
  
  # fill non-existent levels for plotting
  levels <- evs$level 
  if(is.na(match("prior", levels))){
    evs <- evs %>% add_row(level="prior", ev=0, p=value_str)}
  if(is.na(match("LL", levels))){
    evs <- evs %>% add_row(level="LL", ev=0, p=value_str)}
  if(is.na(match("PL", levels))){
    evs <- evs %>% add_row(level="PL", ev=0, p=value_str)}
  return(evs)
}

compute_cond_prob <- function(distr_wide, prob){
  if(prob=="P(C|A)"){
    distr <- distr_wide %>% mutate(p=`AC`/(`AC`+`A-C`))
  } else if(prob=="P(A|C)"){
    distr <- distr_wide %>% mutate(p=`AC`/(`-AC`+`AC`))
  } else if(prob=="P(C|-A)"){
    distr <- distr_wide %>% mutate(p=`-AC`/(`-AC`+`-A-C`))
  } else if(prob=="P(A|-C)"){
    distr <- distr_wide %>% mutate(p=`A-C`/(`A-C`+`-A-C`))
  
  } else if(prob=="P(-C|A)"){
    distr <- distr_wide %>% mutate(p=`A-C`/(`AC`+`A-C`))
  } else if(prob=="P(-A|C)"){
    distr <- distr_wide %>% mutate(p=`-AC`/(`-AC`+`AC`))
  } else if(prob=="P(-C|-A)"){
    distr <- distr_wide %>% mutate(p=`-A-C`/(`-AC`+`-A-C`))
  } else if(prob=="P(-A|-C)"){
    distr <- distr_wide %>% mutate(p=`-A-C`/(`A-C`+`-A-C`))
  
  }  else{
    stop("not implemented.")
  }
  return(distr)
}
# @arg df: data frame containing columns `AC`, `A-C`, `-AC`
add_probs <- function(df){
  df <- df %>% mutate(p_a = AC + `A-C`, p_c = AC + `-AC`,
                      p_na = `-AC` + `-A-C`, p_nc = `A-C` + `-A-C`) %>%
    mutate(p_c_given_a = AC / p_a,
           p_c_given_na = `-AC` / p_na,
           p_a_given_c = AC / p_c, 
           p_a_given_nc = `A-C` / p_nc, 
           p_nc_given_a = `A-C`/p_a,
           p_nc_given_na = `-A-C`/p_na,
           p_na_given_c = `-AC`/p_c,
           p_na_given_nc = `-A-C`/p_nc,
           p_likely_a = p_a,
           p_likely_na=p_na,
           p_likely_c = p_c,
           p_likely_nc=p_nc
    )
  return(df)
}

# other functions ---------------------------------------------------------

#@arg config_keys: order in config_keys is important since same key values
# are overwritten!
configure <- function(config_keys) {
  key <- config_keys[[1]]
  params <- config::get(config=key)
  # print(key)
  for(key in config_keys[-1]){
    # print(key)
    params2 <- config::get(config=key)
    for (name in names(params2)) {
      params[name] = params2[name]
    }
  }
  return(params)
}

# plotting functions ------------------------------------------------------
plot_speaker_conditions <- function(data) {
  df <- data %>% mutate(p=round(as.numeric(p), 2),
                        utterance=as.character(utterance))
  p <- df %>%
    ggplot(aes(y=utterance, x=p, fill=cn)) +
    guides(fill=guide_legend(title="causal net")) +
    geom_bar(stat="identity",
             position=position_dodge(preserve="single"))  +
    labs(x="proportion", y="best utterance") + theme_minimal() +
    facet_wrap(~speaker_condition, labeller=label_parsed) +
    theme(axis.text.y=element_text(), legend.position="top",
          legend.key.size = unit(0.75,"line")) +
    scale_fill_brewer(palette="Dark2")
  return(p)
}

data_cp_plots <- function(params, data=NA){
  
  if(is.na(data)) data <- read_rds(params$target)
  data <- data %>% ungroup() %>% group_by(bn_id, level) %>%
    select(-p_delta, -p_rooij, -p_diff)
  data.wide <- data %>%
    pivot_wider(names_from = "cell", values_from = "val") %>% 
    compute_cond_prob("P(-C|-A)") %>%  rename(`-C_-A` = p) %>% 
    compute_cond_prob("P(A|C)") %>% rename(`A_C` = p)
  
  # Expected values for P(-C|-A) and P(A|C) and for causal nets
  ev_nc_na = data.wide %>% rename(p=`-C_-A`) %>% expected_val("P(-C|-A)")
  ev_a_c = data.wide %>% rename(p=`A_C`) %>% expected_val("P(A|C)") 
  ev_probs <- bind_rows(ev_a_c, ev_nc_na) %>%
    mutate(level=factor(level, levels=c("PL", "LL", "prior")),
           p=str_replace_all(p, "-", "¬")) %>%
    rename(val_type=p) %>% add_column(val="p")
  
  ev_cns = data.wide %>% ungroup() %>% group_by(level, cn) %>% 
    summarise(ev=sum(prob), .groups="drop_last") %>%
    mutate(level=factor(level, levels=c("PL", "LL", "prior")),
           cn=case_when(cn=="A || C" ~ "A,C indep.", TRUE ~ cn),
           cn=str_replace(cn, "-", "¬"),
           cn=str_replace(cn, "implies", "->")
    ) %>%
    rename(val_type=cn) %>% add_column(val="cns")
  
  cns <- c("A -> ¬C", "C -> ¬A", "C -> A", "A -> C", "A,C indep.")
  data <- bind_rows(ev_probs, ev_cns) %>% 
    mutate(val_type=factor(val_type, levels=c(c("P(A|C)", "P(¬C|¬A)"), cns)))
  
  return(data)  
}

plot_cp_cns <- function(cp.cns, labels){
  p.cns <- cp.cns %>% 
    ggplot(aes(y=level, x=ev, fill=val_type)) + 
    geom_bar(position=position_stack(), stat="identity") +
    scale_fill_brewer(palette="Dark2", name="causal net", labels=labels) +
    labs(x="Degree of belief", y="Interpretation level") +
    theme_minimal() +
    theme(legend.position="top", legend.key.size = unit(0.75,"line")) +
    guides(fill=guide_legend(reverse = TRUE))
  return(p.cns)
}

plot_cp_probs <- function(data.cp){
  p.probs <- data.cp %>% filter(val=="p") %>% 
    ggplot(aes(y=level, x=ev, fill=val_type)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    scale_fill_brewer(palette="Dark2", name="value") +
    labs(x="Degree of belief", y="Interpretation level") +
    theme_minimal() +
    theme(legend.position="top", legend.key.size = unit(0.75,"line")) +
    guides(fill=guide_legend(reverse = TRUE))
  return(p.probs)
}

# Acceptability/Assertability conditions ----------------------------------
# p_rooij: (P(e|i) - P(e|¬i)) / (1-P(e|¬i))
# p_delta: P(e|i) - P(e|¬i)
acceptability_conditions <- function(data_wide){
  df <- data_wide %>% compute_cond_prob("P(C|A)") %>% rename(p_c_given_a=p) %>% 
    compute_cond_prob("P(C|-A)") %>% rename(p_c_given_na=p) %>%
    mutate(p_delta=p_c_given_a - p_c_given_na,
           p_nc_given_na=1-p_c_given_na,
           p_rooij=p_delta/(1-p_c_given_na),
           pc=`AC` + `-AC`,
           p_diff=round(p_c_given_a - pc, 5)) %>%
    select(-p_nc_given_na, -p_c_given_a, -p_c_given_na, -pc)
  return(df)
}

# Douven examples ---------------------------------------------------------
voi_douven <- function(posterior, params, model){
  if(model=="skiing"){voi <- voi_skiing(posterior, params)}
  else if(model=="sundowners"){voi <- voi_sundowners(posterior, params)}
  return(voi)
}

voi_skiing <- function(posterior, params){
  pe <- marginalize(posterior, c("E")) 
  ev_pe <- pe %>% expected_val("E") %>% rename(value=ev, key=p) %>% 
    mutate(alpha=params$alpha, cost=params$cost_conditional, pe=params$prior_pe)
  if(params$save){
    ev_pe %>% save_data(paste(params$target_dir, .Platform$file.sep,
                              params$target, "-voi.rds", sep=""))
  }
  return(ev_pe)
}

voi_sundowners <- function(posterior, params){
  pr <- marginalize(posterior, c("R"))
  ev_pr <- pr %>% expected_val("R") %>% rename(value=ev, key=p)
  
  prs <- marginalize(posterior, c("R", "S"))
  ev_prs <- prs %>% expected_val("R and S") %>% rename(value=ev, key=p)
  
  vois <- bind_rows(ev_prs, ev_pr)  %>% 
    mutate(alpha=params$alpha, cost=params$cost_conditional, 
           pr1=params$prior_pr[1],
           pr2=params$prior_pr[2],
           pr3=params$prior_pr[3]) %>% nest(prior_pr=c(pr1,pr2,pr3))
  if(params$save){
    vois %>%  save_data(paste(params$target_dir, .Platform$file.sep,
                              params$target, "-voi.rds", sep=""))
  }
  return(vois)
}
