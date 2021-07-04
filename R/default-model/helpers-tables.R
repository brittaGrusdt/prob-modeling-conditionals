library(truncnorm)
library(dplyr)
library(here)
source(here("R", "helper-functions.R"))

unnest_tables <- function(tables){
  tables <- tables %>% rowid_to_column()
  tables_long <- tables %>% unnest(cols=c(vs, ps)) %>% rename(cell=vs, val=ps)
  return(tables_long)
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
                    sep=.Platform$file.sep)
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

