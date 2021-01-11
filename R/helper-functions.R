EPSILON = 0.0000001
SEP = .Platform$file.sep

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
sort_utterances <- function(utterances){
  literals <- c("A", "C", "-A", "-C")
  conjs <- c("C and A", "-C and A", "C and -A", "-C and -A")
  likely <- c("likely A", "likely C", "likely -A", "likely -C")
  ifs <- c("A > C", "A > -C", "-A > C", "-A > -C",
           "C > A", "C > -A", "-C > A", "-C > -A")
  return(utterances[order(match(utterances, c(conjs, literals, ifs, likely)))])
}

add_pspeaker_max_conj_lit <- function(df){
  df <- df %>% mutate(pmax_conj_lit =
                        max(A, `-A`, C, `-C`,
                            `C and A`, `C and -A`, `-C and A`,`-C and -A`))
  return(df)
}

generate_utts <- function(params){
  utterances <- run_webppl(
    here("model", "default-model", "utterances.wppl", sep=.Platform$file.sep),
    params
  )
  utterances <- utterances %>% map(function(x){x %>% pull(value)}) %>% unlist()
  utterances %>% save_data(params$utts_path)
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

# model ------------------------------------------------------------------
add_model_params <- function(df, params){
  df <- df %>% mutate(cost=params$cost_conditional,
                      alpha=params$alpha,
                      bias=params$bias,
                      value=as.character(value))
  return(df)
}

filter_by_model_params <- function(df, params){
  df <- df %>% filter(cost==params$cost_conditional &
                      alpha==params$alpha)
  return(df)
}

# ** computes likelihood for independent net and  A->C / C->A for if**
likelihood <- function(df_wide, sigma_indep){
  # prepare
  df <- df_wide %>%
    compute_cond_prob("P(C|A)") %>% rename(p_c_given_a=p) %>% 
    compute_cond_prob("P(C|-A)") %>% rename(p_c_given_na=p) %>% 
    compute_cond_prob("P(A|C)") %>% rename(p_a_given_c=p) %>% 
    compute_cond_prob("P(A|-C)") %>% rename(p_a_given_nc=p) %>%
    mutate(pa=AC+`A-C`, pc=AC+`-AC`,
           ind.lower=case_when(1-(pa+pc) < 0 ~ abs(1-(pa+pc)),
                               TRUE ~ 0),
           ind.upper=pmin(pa, pc))
  
  df <- df %>% 
    mutate(
        p_nc_given_a = 1 - p_c_given_a,
        p_na_given_c = 1 - p_a_given_c,
        p_nc_given_na = 1 - p_c_given_na,
        p_na_given_nc = 1 - p_a_given_nc,
        
        logL_ind=log(dtruncnorm(x=`AC`, a=ind.lower, b=ind.upper, mean=pa*pc, sd=sigma_indep)),
        logL_if_ac = log(dbeta(p_c_given_a, 10, 1))+log(dbeta(p_c_given_na, 1, 10)),
        logL_if_anc = log(dbeta(p_nc_given_a, 10, 1)) + log(dbeta(p_nc_given_na, 1, 10)),
        logL_if_ca = log(dbeta(p_a_given_c, 10, 1)) + log(dbeta(p_a_given_nc, 1, 10)),
        logL_if_cna = log(dbeta(p_na_given_c, 10, 1)) + log(dbeta(p_na_given_nc, 1, 10))
    ) %>% 
    select(-p_c_given_na, -p_c_given_a, -p_a_given_c, -p_a_given_nc, -pa, -pc,
           -p_nc_given_a, -p_na_given_c, -p_nc_given_na, -p_na_given_nc)
  return(df)
}

# other functions ---------------------------------------------------------
hellinger <- function(p, q){
  (1/sqrt(2)) * sqrt(sum((sqrt(p)-sqrt(q))^2))
}

adapt_bn_ids <- function(data_wide){
  # only considers levels PL, LL and prior
  df <- data_wide %>% dplyr::select(-prob, -level, -bn_id, -cn)
  cell_names <- names(df)
  data_wide <- data_wide %>% unite(cells, names(df), sep="__")
  
  # makes sure that bn_ids are identical across levels PL/LL/prior
  prior <- data_wide %>% filter(level=="prior") %>% arrange(cn, cells)
  ll <- data_wide %>% filter(level=="LL") %>% arrange(cn, cells)
  pl <- data_wide %>% filter(level=="PL") %>% arrange(cn, cells)
  df <- bind_rows(ll, prior)
  
  # not all Bayes nets that are in the prior also occur in the literal/pragmatic listener
  # (but there is no diff btw. those in LL/PL)
  idx_dups <- df %>% dplyr::select(-level, -prob, -bn_id) %>% duplicated()
  duplicates_prior <- df[idx_dups, ] %>%  arrange(cn, cells)
  # duplicates_prior$level %>% unique()
  
  pl <- pl %>% mutate(bn_id=duplicates_prior$bn_id)
  ll <- ll %>% mutate(bn_id=duplicates_prior$bn_id)
  
  df <- bind_rows(prior, ll, pl)
  df <- df %>% separate(cells, cell_names, sep="__", convert=TRUE) 
  return(df)
}

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
plot_evs <- function(data){
  p <- data %>% ggplot() +
    geom_bar(mapping = aes(x=level, y=ev, fill=level), stat="identity", position="dodge") +
    labs(x="", y="", title="") +
    coord_flip() +
    theme_classic(base_size = 20) +
    theme(legend.position="none")
  return(p)
}

plot_speaker <- function(data, fn, w, h, plot_dir, legend_pos="none",
                         facets=TRUE, xlab="", ylab="") {
  df <- data %>% mutate(p=round(as.numeric(p), 2))
  if(xlab==""){xlab = TeX("$\\frac{1}{|S|} \\cdot \\sum_{s \\in S} P_S(u|s)$")}
  if(ylab==""){ylab = "utterance"}
  
  if("cn" %in% colnames(df)) {
    p <- df %>%
      ggplot(aes(y=utterance, x=p, fill=cn)) +
      guides(fill=guide_legend(title="causal net"))
  } else if("speaker_condition" %in% colnames(df)) {
    p <-  df %>% ggplot(aes(y=utterance, x=p, fill=speaker_condition))
  } else {
    p <-  df %>% ggplot(aes(y=utterance, x=p))
  }
  p <- p +
    geom_bar(stat="identity", position=position_dodge(preserve = "single"))  +
    labs(x=xlab, y=ylab) + theme_bw(base_size=25)
  if(facets) {p <- p + facet_wrap(~speaker_condition)
  }
  p <- p + theme(axis.text.y=element_text(size=15), legend.position=legend_pos)
  
  ggsave(paste(plot_dir, fn, sep=SEP), p, width=w, height=h)
  return(p)
}

# @arg posterior: in long format, must have columns *cell* and *val*
voi_default <- function(posterior, params){
  df = posterior %>% ungroup() %>% dplyr::select(-starts_with("p_"))
  df.wide = df %>% group_by(bn_id) %>%
    pivot_wider(names_from="cell", values_from="val") %>%
    add_probs() %>% dplyr::select(!starts_with("p_likely")) %>%
    group_by(level)
  
  theta=params$theta
  df = df.wide %>% mutate(uncertainty =
    case_when((p_a<theta & p_a>1-theta) & (p_c<theta & p_c>1-theta) ~ "both",
            (p_a<theta & p_a>1-theta) ~ "only A",
            (p_c<theta & p_c>1-theta) ~ "only C",
            TRUE ~ "none"))
  # bns where certain about both (=uncertain about none) is true
  df.certain_both = df %>% filter(uncertainty == "none") %>%
    summarize(ev=sum(prob), .groups="drop_last") %>% 
    add_column(key="uncertain_none")
  
  df.uncertain_only_a = df %>% filter(uncertainty == "only A") %>%
    summarize(ev=sum(prob), .groups="drop_last") %>%
    add_column(key="uncertain_only_A")
  
  df.uncertain_only_c = df %>% filter(uncertainty == "only C") %>%
    summarize(ev=sum(prob), .groups="drop_last") %>%
    add_column(key="uncertain_only_C")
  
  df.uncertain_both = df %>% filter(uncertainty == "both") %>%
    summarize(ev=sum(prob), .groups="drop_last") %>%
    add_column(key="uncertain_both")
  
  # expected value P(A)
  evs = df.wide %>%
    transmute(ev_a=prob*p_a, ev_c=prob*p_c, ev_a_given_c = prob * p_a_given_c,
              ev_nc_given_na= prob * p_nc_given_na) %>%
    summarize(ev_a=sum(ev_a), ev_c=sum(ev_c), ev_a_given_c=sum(ev_a_given_c),
              ev_nc_given_na=sum(ev_nc_given_na), .groups="keep") %>%
    pivot_longer(cols=c("ev_a", "ev_c", "ev_a_given_c", "ev_nc_given_na"),
                 names_to="key", values_to="ev")
  
  results <- bind_rows(df.uncertain_both, df.uncertain_only_a,
                       df.uncertain_only_c, df.certain_both, evs)
  if(params$level_max == "prior"){
    levels = c("prior")
  } else {
    levels = c("prior", "LL", "PL")
  }
  results = results %>% filter(level %in% levels)
  
  if(params$save){
    results %>%
      save_data(paste(params$target_dir, .Platform$file.sep,
                      str_split(params$target_fn, "\\.")[[1]][1],
                      "-prior-LL-PL-vois.rds", sep=""))
  }
  return(results)
}
# Acceptability/Assertability conditions ----------------------------------
# p_rooij: (P(e|i) - P(e|¬i)) / (1-P(e|¬i))
# p_delta: P(e|i) - P(e|¬i)
acceptability_conditions <- function(data_wide){
  df <- data_wide %>% compute_cond_prob("P(C|A)") %>% rename(p_c_given_a=p) %>% 
    compute_cond_prob("P(C|-A)") %>% rename(p_c_given_na=p) %>%
    mutate(p_delta=round(p_c_given_a - p_c_given_na, 5),
           p_nc_given_na=round(1-p_c_given_na, 5),
           p_rooij=case_when(p_nc_given_na == 0 ~ round(p_delta/0.00001, 5),
                             TRUE ~ round(p_delta/p_nc_given_na, 5)),
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
