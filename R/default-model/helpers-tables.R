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
      gather(`AC`, `A-C`, `-AC`, `-A-C`, key="cell", val="val") %>%
      filter(val != 0)
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
  tables <- tibble(pc=runif(params$n_ind_tables), pa=runif(params$n_ind_tables)) %>%
    rowid_to_column("id") %>%
    mutate(upper_bound = pmin(pa, pc),
           lower_bound = ifelse(1-(pa+pc) < 0, abs(1-(pa+pc)), 0),
           #noisy samples
           `AC`= rtruncnorm(1, a=lower_bound, b=upper_bound, mean=pa*pc, sd=params$indep_sigma), 
           `-AC`=pc-`AC`, `A-C`=pa-`AC`, s=`AC` + `-AC` + `A-C`, `-A-C`= 1 - s) %>%
    select(-upper_bound, -lower_bound, -pa, -pc, -s)
  tables.mat = tables  %>% select(-id) %>% as.matrix() 
  
  tables = prop.table(tables.mat, 1) %>% as_tibble() %>%
    mutate(n=AC + `A-C` + `-AC` + `-A-C`) %>%
    add_column(id=tables$id)
           
  tables_long <- tables %>%
    gather(`AC`, `A-C`, `-AC`, `-A-C`, key="cell", val="val") %>%
    filter(val != 0)
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
  tables <- bind_rows(tables_ind, tables_dep) %>% rowid_to_column("table_id") %>% 
              mutate(seed=params$seed_tables)
  tables <- tables %>% unnest(c(vs, ps)) %>%
    group_by(table_id) %>% pivot_wider(names_from="vs", values_from="ps") %>% 
    likelihood(params$indep_sigma) %>% 
    mutate(vs=list(c("AC", "A-C", "-AC", "-A-C")),
           ps=list(c(`AC`, `A-C`, `-AC`, `-A-C`))) %>%
    select(-`AC`, -`A-C`, -`-AC`, -`-A-C`) %>%
    rename(cn.orig=cn)
  tables = tables_to_bns(tables, params)  %>%
    mutate(bn_id=case_when(
      cn=="A || C" ~ paste(table_id, "independent", sep="_"),
      TRUE ~ paste(table_id, str_replace_all(cn, " ", ""), sep="_")
    )) 
  tables %>% save_data(params$tables_path)
  return(tables)
}

tables_to_bns = function(tables, params) {
  tables.ll = tables %>% group_by(table_id) %>%
    pivot_longer(cols=starts_with("logL_"), names_to="ll_cn", values_to="ll")
  # filter out the n worst causal nets, if there are several cns with the same
  # ll, make sure that at least the cn with the best ll is kept!
  # (e.g. for 0.25-0.25-0.25-0.25 ll for all dependent nets identical, but keep ind!)
  tbls = tables.ll  %>%
    mutate(cn=case_when(ll_cn=="logL_ind" ~ "A || C",
                        ll_cn=="logL_if_ac" ~ "A implies C",
                        ll_cn=="logL_if_ca" ~ "C implies A",
                        ll_cn=="logL_if_anc" ~ "A implies -C",
                        ll_cn=="logL_if_cna" ~ "C implies -A")) %>%
    dplyr::select(-ll_cn, -ind.lower, -ind.upper) %>% 
    group_by(table_id) %>% mutate(best.cn=ll==max(ll)) %>%
    arrange(desc(ll)) %>% mutate(rank=seq(1:n()))
  
  bns = tbls %>%
    mutate(bn_id=case_when(
      cn=="A || C" ~ paste(table_id, "independent", sep="_"),
      TRUE ~ paste(table_id, str_replace_all(cn, " ", ""), sep="_")
    )) %>% group_by(bn_id) %>% ungroup() %>%
    filter(rank <= params$n_best_cns) %>%
    dplyr::select(bn_id, table_id, ps, vs, ll, cn, cn.orig)
  return(bns)
}

unnest_tables <- function(tables){
  tables <- tables %>% rowid_to_column()
  tables_long <- tables %>% unnest(cols=c(vs, ps)) %>% rename(cell=vs, val=ps)
  return(tables_long)
}

# data must be in long format with columns 'cell', 'val'
plot_tables <- function(data){
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
                 labeller = labeller(cell =
                                       c(`-A-C` = "P(¬A,¬C)", `-AC`= "P(¬A,C)",
                                         `A-C` = "P(A,¬C)", `AC` = "P(A,C)"))
                 ) +
      labs(title = cn_title, x=xlab, y=ylab) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(size=10)) +
      scale_color_brewer(palette="Dark2")
    plots[[idx]] <- p
    idx <- idx + 1
  }
  return(plots)
}

# plot densities of generated tables for different causal nets
# @arg tables: long format
plot_tables_cns <- function(tables_path, plot_dir, w, h, tables=NA){
  if(is.na(tables)){
    tables <- readRDS(tables_path)
  }
  tables.wide <- tables %>% unnest_tables() %>% ungroup() %>%
    select(-table_id) %>% group_by(rowid) %>% 
    pivot_wider(names_from = cell, values_from = val)
  tables.long <- tables.wide %>% 
    pivot_longer(cols = c(AC, `A-C`, `-AC`, `-A-C`), names_to = "cell",
                 values_to = "val") %>% 
    group_by(bn_id, cn) %>% 
    mutate(cell=factor(cell, levels=c("AC", "A-C", "-AC", "-A-C")),
           cn=case_when(cn=="A || C" ~ "A,C independent",
                        cn=="A implies -C" ~ "A implies ¬C",
                        TRUE ~ cn),
           cn=as.factor(cn))
  all_plots = list()
  cns <- list(c("A,C independent"), c("A implies ¬C"), c("A implies C"))
  cns.expr <- list(c("A,C independent"), c(expression(paste(A %->%"", "¬C"))),
                   c(expression(A %->% C)))
  cns.short <- c("indep", "anc", "ac")
  for(i in seq(1,3)) {
    p <- tables.long %>% filter(cn %in% cns[[i]]) %>%
      ggplot(aes(x=val,  fill = cn)) +
      geom_density() +
      facet_wrap(~cell, ncol = 2, scales = "free",
                 labeller = labeller(cell = c(`AC` = "P(A,C)", `A-C` = "P(A,¬C)",
                                              `-AC`= "P(¬A,C)", `-A-C` = "P(¬A,¬C)"))
      ) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
      scale_fill_brewer(palette = "Dark2") +
      labs(x="probability", y="density") +
      theme_minimal() +
      theme(legend.position = "none") +
      ggtitle(cns.expr[[i]])
    all_plots[[i]] = p
    
    save_to = paste(plot_dir, paste("tables-", cns.short[[i]], ".png", sep=""),
                    sep=SEP)
    ggsave(save_to, p, width=w, height=h)
    print(paste('saved to', save_to))
  }
  return(all_plots)
}

# Analyze generated Tables ------------------------------------------------
analyze_tables <- function(params){
  theta = params$theta
  prior = readRDS(params$target) %>% filter(level=="prior") %>%
    pivot_wider(names_from="cell", values_from="val") %>% ungroup() %>% 
    dplyr::select(prob, bn_id)
  tables <- read_rds(params$tables_path) %>% ungroup() %>% 
      select(bn_id, cn.orig, cn, vs, ps) %>% group_by(bn_id, cn) %>% unnest(c(ps, vs))
  tables = left_join(tables, prior, by=c("bn_id"))
  tables.wide <-  tables %>% pivot_wider(names_from = vs, values_from = ps)  
  
  results = tibble()
  n.wide = nrow(tables.wide)
  n.long = nrow(tables)
  # conjunctions
  df.conj <- tables %>% mutate(conj=case_when(ps >= theta ~ TRUE,
                                              TRUE ~ FALSE))
  count_conj=function(df.conj, cell){
    df = df.conj %>% filter(conj & vs == (!! cell)) 
    prob = sum(df$prob)
    n = df %>% nrow()
    ratio = round(n / nrow(df.conj), 5)
    return(tibble(n=n, ratio=ratio, key=cell))
  }
  
  conj = pmap_dfr(tibble(cell=c("AC", "A-C", "-AC", "-A-C")), function(cell){
    return(count_conj(df.conj, cell))
  })
  
  conditionals <- tables.wide %>%
    add_column(pca=compute_cond_prob(tables.wide, "P(C|A)") %>% pull(p) > theta,
               pac=compute_cond_prob(tables.wide, "P(A|C)") %>% pull(p) > theta,
               pcna = compute_cond_prob(tables.wide, "P(C|-A)") %>% pull(p) > theta,
               panc = compute_cond_prob(tables.wide, "P(A|-C)") %>% pull(p) > theta
    )
  ifs = tribble(~n, ~key,
          conditionals %>% filter(pca) %>% nrow(), "if_ac",
          conditionals %>% filter(pac) %>% nrow(), "if_ca",
          conditionals %>% filter(pcna) %>% nrow(), "if_nac",
          conditionals %>% filter(panc) %>% nrow(), "if_ncna") %>% 
    mutate(ratio = round(n/n.wide, 2))

  literals <- tables.wide %>%
    mutate(a=`AC` + `A-C` > 0.5,
           c=`AC` + `-AC` > 0.5,
           na=`-AC` + `-A-C` > 0.5,
           nc=`A-C` + `-A-C` > 0.5)
  n.likely_a = literals %>% filter(a) %>% nrow
  n.likely_c = literals %>% filter(c) %>% nrow
  n.likely_na = literals %>% filter(na) %>% nrow
  n.likely_nc = literals %>% filter(nc) %>% nrow
  
  likely = tribble(~n, ~key,
                    n.likely_a, "likely_a",
                    n.likely_c, "likely_c",
                    n.likely_na, "likely_na",
                    n.likely_nc, "likely_nc") %>% 
    mutate(ratio=round(n/n.wide, 2))
  
  literals <- tables.wide %>%
    mutate(a=`AC` + `A-C` > theta,
           c=`AC` + `-AC` > theta,
           na=`-AC` + `-A-C` > theta,
           nc=`A-C` + `-A-C` > theta)
  n.a = literals %>% filter(a) %>% nrow
  n.c = literals %>% filter(c) %>% nrow
  n.na = literals %>% filter(na) %>% nrow
  n.nc = literals %>% filter(nc) %>% nrow

  lits = tribble(~n, ~key,
                    n.a, "A",
                    n.c, "C",
                    n.na, "-A",
                    n.nc, "-C") %>% 
    mutate(ratio=round(n/n.wide, 2))
  return(bind_rows(conj, lits, likely, ifs) %>% arrange(ratio))
}

analyze_table_likelihoods = function(params, w=5, h=5){
  tables <- read_rds(params$tables_path) %>% unnest(c(ps, vs)) %>% group_by(id)
  tables.wide <-  tables %>% pivot_wider(names_from = "vs", values_from = "ps") %>%
    mutate(cn.orig=case_when(cn.orig=="A || C" ~ "A,C independent",
                             cn.orig=="A implies C" ~ "A implies C",
                             cn.orig=="A implies -C" ~ "A implies ¬C",
                             cn.orig=="C implies A" ~ "C implies A",
                             cn.orig=="C implies -A" ~ "C implies ¬A"))
  tbls.best.cns = tables.wide %>% mutate(best_cn= ll==max(ll)) %>%
    filter(best_cn) %>% dplyr::select(-best_cn)
  tables.long = tbls.best.cns %>%
    pivot_longer(c(AC, `A-C`, `-AC`, `-A-C`), names_to="cell", values_to="val") %>%
    mutate(cell=factor(cell, levels=c("AC", "A-C", "-AC", "-A-C")),
           cn.orig=as.factor(cn.orig)) %>% 
    group_by(id, cn.orig)
  
  all_plots = list()
  cns <- list(c("A,C independent"), c("A implies ¬C"), c("A implies C"))
  cns.short <- c("indep", "anc", "ac")
  for(i in seq(1,3)) {
    p <- tables.long %>% filter(cn.orig %in% cns[[i]]) %>%
      ggplot(aes(x=val,  fill = cn.orig)) +
      geom_density() +
      facet_wrap(~cell, ncol = 2, scales = "free",
                 labeller = labeller(cell = c(`AC` = "P(A,C)", `A-C` = "P(A,¬C)",
                                              `-AC`= "P(¬A,C)", `-A-C` = "P(¬A,¬C)"))
      ) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
      labs(x="probability", y="density") +
      theme_classic(base_size = 20) +
      theme(legend.position = "none") +
      ggtitle(cns[[i]])
    all_plots[[i]] = p
    
    save_to = paste(params$plot_dir, paste("tables-", cns.short[[i]], ".png", sep=""), sep=SEP)
    ggsave(save_to, p, width=w, height=h)
    print(paste('saved to', save_to))
  }
}

plot_sigma_ind_tables = function(){
  x=seq(0,0.3, length=1000)
  x_labels <-seq(0, 0.3, by=0.05)
  x_breaks <- seq(0, 0.3, by=0.05)
  sigmas = c(0.001, 0.005, 0.01, 0.1)
  for(sigma in sigmas){
    y=dtruncnorm(x, a=0, b=0.3, mean=0.12, sd=sigma)
    p = ggplot() + 
      geom_point(aes(x=x, y=y), size=0.5) +
      scale_x_continuous(breaks=x_breaks, labels=x_labels) +
      theme_classic() +
      theme(axis.text.x=element_text(angle=30, vjust=-.01)) +
      labs(x="") 
    
    ggsave(here("figs", paste("truncnorm_ind_sigma_", sigma, ".png", sep="")), p)
  }
}

