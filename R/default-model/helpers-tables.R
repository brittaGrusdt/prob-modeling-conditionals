library(truncnorm)
library(dplyr)
library(here)
source(here("R", "helper-functions.R"))

# Table Generation --------------------------------------------------------
create_dependent_tables <- function(params, cns){
  all_tables <- list()
  idx <- 1
  for(cn in cns){
    theta <- rbeta(params$n_tables, 10, 1)
    beta <- rbeta(params$n_tables, 1, 10)
    p_child_parent <- theta + beta * (1 - theta)
    p_child_neg_parent <- beta
    p_parent <- runif(params$n_tables)
  
    if(cn %in% c("A implies C", "C implies A")){
      probs <- tibble(cond1=p_child_parent, cond2=p_child_neg_parent, marginal=p_parent)
    } else if(cn %in% c("A implies -C", "C implies -A")){
      probs <- tibble(cond1=1-p_child_parent, cond2=1-p_child_neg_parent, marginal=p_parent)
    }  
    # A -> C and -A -> C use the same probabilities (P(C|A), P(C|-A), P(A)/P(-A))
    if(startsWith(cn, "A")){
      probs <- probs %>% mutate(`AC`=cond1 * marginal,
                                `A-C`=(1-cond1) * marginal,
                                `-AC`=cond2 * (1-marginal),
                                `-A-C`=(1-cond2) * (1-marginal))
    } else if(startsWith(cn, "C")){
      # diagonals are switched
      probs <- probs %>% mutate(`AC`=cond1 * marginal,
                                `A-C`=cond2 * (1-marginal),
                                `-AC`=(1-cond1) * marginal,
                                `-A-C`=(1-cond2) * (1-marginal))
    } else {
      stop(paste(cn, "not implemented."))
    }
    tables <- probs %>% dplyr::select(-cond1, -cond2, -marginal) %>%
      rowid_to_column("id")
    tables_long <- tables %>%
      gather(`AC`, `A-C`, `-AC`, `-A-C`, key="cell", val="val") #%>%
      #mutate(val=round(val, 4))
    tables_wide <- tables_long %>% group_by(id) %>%
      summarise(ps = list(val), .groups = 'drop') %>% add_column(cn=(!! cn)) %>%
      mutate(vs=list(c("AC", "A-C", "-AC", "-A-C"))) %>% dplyr::select(-id)
    
    all_tables[[idx]] <- tables_wide
    idx <- idx + 1
  }
  tables <- all_tables %>% bind_rows()
  return(tables)
}

create_independent_tables <- function(params){
  # digits=4
  tables <- tibble(pc=runif(params$n_ind_tables), pa=runif(params$n_ind_tables)) %>%
    rowid_to_column("id") %>%
    mutate(upper_bound = pmin(pa, pc),
           lower_bound = ifelse(1-(pa+pc) <= 0, abs(1-(pa+pc)), 0),
           #noisy samples
           `AC`= rtruncnorm(1, a=lower_bound, b=upper_bound, mean=pa*pc, sd=params$indep_sigma), 
           `-AC`=pc-`AC`,
           `A-C`=pa-`AC`,
           s=`AC` + `-AC` + `A-C`,
           `-A-C`= 1 - s) %>%
    select(-upper_bound, -lower_bound, -pa, -pc, -s)
  tables.mat = tables  %>% select(-id) %>% as.matrix() 
  
  tables = prop.table(tables.mat + epsilon, 1) %>% as_tibble() %>%
    mutate(n=AC + `A-C` + `-AC` + `-A-C`) %>%
    add_column(id=tables$id)
           
  tables_long <- tables %>%
    gather(`AC`, `A-C`, `-AC`, `-A-C`, key="cell", val="val")# %>% 
    # mutate(val=round(val, 4))
  tables_wide <- tables_long %>% group_by(id) %>%
    summarise(ps = list(val), .groups = 'drop') %>% add_column(cn="A || C") %>% 
    mutate(vs=list(c("AC", "A-C", "-AC", "-A-C"))) %>%
    dplyr::select(-id)
  return(tables_wide)
}

create_tables <- function(params){
  set.seed(params$seed_tables)
  cns_dep=params$cns[params$cns != "A || C"]
  tables_all <- list()
  tables_ind <- create_independent_tables(params)
  tables_dep <- create_dependent_tables(params, cns_dep)
  tables <- bind_rows(tables_ind, tables_dep) %>% rowid_to_column("id") %>% 
              mutate(seed=params$seed_tables)
  tables <- tables %>% unnest(c(vs, ps)) %>%
    group_by(id) %>% pivot_wider(names_from="vs", values_from="ps") %>% 
    likelihood(params$indep_sigma) %>% 
    mutate(vs=list(c("AC", "A-C", "-AC", "-A-C")),
           ps=list(c(`AC`, `A-C`, `-AC`, `-A-C`))) %>%
    select(-`AC`, -`A-C`, -`-AC`, -`-A-C`) %>%
    mutate(stimulus_id=case_when(
      cn=="A || C" ~ paste(id, "independent", sep="_"),
      TRUE ~ paste(id, str_replace_all(cn, " ", ""), sep="_"))
      );
  
  tables %>% save_data(params$tables_path)
  return(tables)
}
# 
# filter_tables <- function(tables, params){
#   if(is.na(params$nor_beta)){
#     df <- tables %>% filter(is.na(nor_beta) & is.na(nor_theta))
#   } else{
#     df <- tables %>% filter(nor_beta == params$nor_beta & nor_theta == params$nor_theta)
#   }
#   df <- df %>% filter(n_tables == params$n_tables & indep_sigma == params$indep_sigma)
#   return(df)
# }

unnest_tables <- function(tables){
  tables <- tables %>% rowid_to_column()
  tables_long <- tables %>% unnest(cols=c(vs, ps)) %>% rename(cell=vs, val=ps)
  return(tables_long)
}

# # TODO: compare with same function in helper-functions.R!
# adapt_bn_ids <- function(data_wide){
#   # makes sure that bn_ids are identical across levels PL/LL/prior
#   prior <- data_wide %>% filter(level=="prior") %>% arrange(cn, `AC`, `-AC`, `A-C`, `-A-C`)
#   ll <- data_wide %>% filter(level=="LL") %>% arrange(cn, `AC`, `-AC`, `A-C`, `-A-C`)
#   pl <- data_wide %>% filter(level=="PL") %>% arrange(cn, `AC`, `-AC`, `A-C`, `-A-C`)
#   df <- bind_rows(ll, prior)
#   
#   # not all Bayes nets that are in the prior also occur in the literal/pragmatic listener
#   # (but there is no diff btw. those in LL/PL)
#   idx_dups <- df %>% dplyr::select(-level, -prob, -bn_id) %>% duplicated()
#   duplicates_prior <- df[idx_dups, ] %>%  arrange(cn, `AC`, `AC`, `A-C`, `-AC`, `-A-C`)
#   duplicates_prior$level %>% unique()
#   
#   # pl <- pl %>% mutate(bn_id=duplicates_prior$bn_id)
#   ll <- ll %>% mutate(bn_id=duplicates_prior$bn_id)
#   
#   df <- bind_rows(prior, ll, pl)
#   return(df)
# }

plot_tables <- function(data){
  # data must be in long format with columns *cell* and *val*
  cns <- data$cn %>% as.factor() %>% levels()
  plots <- list(); idx = 1
  for(causal_net in cns){
    if(causal_net == "A || C"){
      cn_title <- "A,C indep."
    } else {
      cn_title = ifelse(endsWith(causal_net, "-A"), TeX("$C\\rightarrow\\neg A$"),
                        ifelse(endsWith(causal_net, "-C"), TeX("$A\\rightarrow\\neg C$"),
                                        parse(text=str_replace(causal_net, " implies ", '%->%'))))
    }
    ylab = ifelse(idx %in% c(1,4), "density", "");
    xlab = ifelse(idx %in% c(3,4,5), "probability", "");
    
    p <- data %>% 
      filter(cn==causal_net) %>%
      ggplot(aes(x=val,  color = cell)) +
      geom_density() +
      facet_wrap(~cell, ncol = 2, scales = "free",
                 labeller = labeller(cell = c(`AC` = "P(A,C)", `A-C` = "P(A,¬C)",
                                              `-AC`= "P(¬A,C)", `-A-C` = "P(¬A,¬C)"))
                 ) +
      labs(title = cn_title, x=xlab, y=ylab) +
      theme_classic(base_size = 20) +
      theme(legend.position = "none", axis.text.x = element_text(size=10))
    plots[[idx]] <- p
    idx <- idx + 1
    print(p)
  }
  return(plots)
}

plot_tables_all_cns <- function(tables_path, plot_dir, w, h){
  tables.wide <- readRDS(tables_path) %>% unnest_tables() %>%
    rename(bn_id=rowid) %>% group_by(bn_id, cn) %>% 
    pivot_wider(names_from = cell, values_from = val) %>% ungroup()
  tables.long <- tables.wide %>% 
    mutate(`-A-C` = case_when(is.na(`-A-C`) ~ rowSums(select(., starts_with("-A-C_"))),
                              TRUE ~ `-A-C`),
           `-AC` = case_when(is.na(`-AC`) ~ rowSums(select(., starts_with("-AC_"))),
                             TRUE ~ `-AC`), 
           `A-C` = case_when(is.na(`A-C`) ~ rowSums(select(., starts_with("A-C_"))),
                             TRUE ~ `A-C`),
           `AC` = case_when(is.na(`AC`) ~ rowSums(select(., starts_with("AC_"))),
                            TRUE ~ `AC`)) %>% 
    group_by(bn_id, cn) %>%
    pivot_longer(cols = c(AC, `A-C`, `-AC`, `-A-C`), names_to = "cell", values_to = "val") %>% 
    ungroup() %>% 
    mutate(cell=factor(cell, levels=c("AC", "A-C", "-AC", "-A-C")),
           cn=case_when(cn=="A || C" ~ "A,C independent",
                        TRUE ~ cn),
           cn=as.factor(cn)) %>% 
    group_by(bn_id, cn)
  
  # tables must be in long format with columns *cell* and *val*
  all_plots = list()
  cns <- list(c("A,C independent"), c("A implies -C", "C implies -A"), c("A implies C", "C implies A"))
  cns.short <- c("indep", "anc-cna", "ac-ca")
  for(i in seq(1,3)) {
    p <- tables.long %>% filter(cn %in% cns[[i]]) %>%
      ggplot(aes(x=val,  fill = cn)) +
      geom_density(alpha=0.5) +
      facet_wrap(~cell, ncol = 2, scales = "free",
                 labeller = labeller(cell = c(`AC` = "P(A,C)", `A-C` = "P(A,¬C)",
                                              `-AC`= "P(¬A,C)", `-A-C` = "P(¬A,¬C)"))
      ) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
      labs(x="probability", y="density") +
      theme_classic(base_size = 20) +
      theme(legend.position = "bottom")
    all_plots[[i]] = p
    
    save_to = paste(plot_dir, paste("tables-", cns.short[[i]], ".png", sep=""), sep=SEP)
    ggsave(save_to, p, width=w, height=h)
    print(paste('saved to', save_to))
  }
  return(all_plots)
}


makeTables <- function(){
  params <- configure(c("bias_none", "tables"))
  tables <- create_tables(params)
  return(tables)  
}

# params = configure(c("bias_none", "tables"))
# TABLES <- readRDS("./data/default-model/tables-default.rds") %>%
#   select(id, cn, vs, ps) %>% unnest(c(ps, vs)) %>% group_by(id)
# TABLES.wide <-  TABLES %>% pivot_wider(names_from = vs, values_from = ps)  








# tables_data <- marginalize(tables_long %>% rename(level=cn), c("A")) %>%
#   mutate(p=case_when(p==0 ~ 0.00001, TRUE~p), `P(C|A)`=AC/p)
# tables_data %>% group_by(level) %>% summarise(m=mean(`P(C|A)`))

# df <- tables %>% spread(key=cell, val=val) %>% mutate(pca=AC/(AC+`A-C`), pac=AC/(AC+`-AC`))
# df %>% group_by(cn) %>% summarise(mean_pca=mean(pca, na.rm = TRUE), mean_pac=mean(pac, na.rm = TRUE))