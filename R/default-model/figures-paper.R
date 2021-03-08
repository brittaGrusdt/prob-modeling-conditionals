library(tidyverse)
library(config)
library(ggplot2)
library(latex2exp)
library(cowplot)
library(grid)
source("R/helper-functions.R")
source("R/default-model/helpers-tables.R")

SEP = .Platform$file.sep
data_dir = here("data", "default-model", "paper-config")
# data_dir = here("data", "default-model", "my-config")

plot_dir = paste(data_dir, "figs", sep=SEP)
if(!dir.exists(plot_dir)) dir.create(plot_dir)
params <- read_rds(paste(data_dir, "params-none.rds", sep=SEP))
theta=params$theta

params.speaker <- read_rds(paste(data_dir, "params-speaker.rds", sep=SEP))
UTTERANCES <- read_rds(paste(params.speaker$target_dir, params.speaker$utts_fn, sep=SEP))

# Figure 2 ----------------------------------------------------------------
# for plotting densities sample more tables 
params.tables = list(
  seed_tables=params$seed_tables, cns=params$cns, indep_sigma=params$indep_sigma,
  n_best_cns=params$n_best_cns,
  n_tables=10000, n_ind_tables=10000,
  tables_path=paste(params$target_dir, "tables-default-samples-for-plots.rds",
                    sep=SEP)
)
tables.plot = create_tables(params.tables) %>% select(-cn, -ll) %>%
  rename(cn=cn.orig)
save_data(tables.plot, params.tables$tables_path)
plot_tables_cns(params.tables$tables_path, params$plot_dir, w=5, h=5)

# plot tables for samples from prior
tbls.prior = readRDS(here("data", "default-model", "paper-config",
                   "results-none-priorN.rds"))
plot_tables(tbls.prior)

# analyze_table_likelihoods(params)

# Figure 3 ----------------------------------------------------------------
# pragmatic/literal interpretations of conditional If A, C
# probability assigned to states where speaker is uncertain about A and C
dat.none.voi <- read_rds(
  paste(params$target_dir, "results-none-prior-LL-PL-vois.rds", sep=SEP)
  ) %>%
  filter((key == "uncertain_both") & level !="prior") %>% 
  mutate(ev=round(ev, 2))

p <- dat.none.voi %>% 
  ggplot(aes(y=ev, x=level)) + 
  geom_bar(stat="identity")  +
  geom_text(aes(label = ev, x = level,  y = ev), hjust=-0.1, size=6) + 
  scale_x_discrete(limits = c("LL", "PL"),               
                   labels=c(paste(strwrap("Literal interpretation", width=20), collapse="\n"),
                            paste(strwrap("Pragmatic interpretation", width=20), collapse="\n"))
  ) + 
  labs(y= TeX("$\\sum_{s\\in Uncertain_s(A) \\bigcap Uncertain_s(C) \\}
            Pr(s|u=A\\rightarrow C)$"), x="") +
  theme_classic() +
  theme(legend.position = "none") +
  coord_flip()
ggsave(paste(params$plot_dir, "ignorance-inferences.png", sep=SEP), p,
       width=13.5, height=4)


# Figure 4 ----------------------------------------------------------------
data.speaker <- read_rds(params.speaker$target) %>% select(-level, -bias) %>%
  select(-p_delta, -p_diff)
bn_ids = read_rds(paste(params.speaker$target_dir, SEP, "sample-ids-",
                        params.speaker$target_fn, sep="")) %>% 
  group_by_all() %>% summarize(n_sampled=n(), .groups="drop_last")

data.speaker.best <- data.speaker %>% group_by(bn_id) %>%
  mutate(p_best=max(probs), u_best=list(utterance[probs == max(probs)])) %>%
  unnest(c(u_best)) %>% select(-p_best)

data.speaker.best = left_join(data.speaker.best, bn_ids, by=c("stimulus_id"))

dat <- data.speaker.best %>%
  mutate(pa=`AC` + `A-C`, pc=`AC` + `-AC`,
         certainA = pa >= theta | pa <= 1-theta,
         certainC = pc >= theta | pc <= 1-theta,
         uncA = pa > 1 - theta & pa < theta,
         uncC = pc > 1-theta & pc < theta,
         certain=certainA & certainC,
         uncertain=uncA & uncC, 
         both=!certain & !uncertain,
         speaker_condition = case_when(certain ~ "A,C certain",
                                       uncertain ~ "A,C uncertain",
                                       both ~ "A (C) certain, C (A) uncertain")) %>% 
  select(-uncA, -uncC, -certainA, -certainC, -certain, -uncertain) %>%
  select(-utterance) %>%
  rename(utterance = u_best) %>% distinct(bn_id, .keep_all = TRUE) %>%
  chunk_utterances()

df <- dat %>%
  group_by(speaker_condition, cn, utterance) %>%
  mutate(p_ifac=(AC/pa) >= 0.9) %>% 
  summarise(p=sum(n_sampled), n.ifac=sum(p_ifac*n_sampled), .groups = "drop_last") %>%
  arrange(p) %>% 
  mutate(N=sum(p), ratio=p/N) %>% rename(n=p, p=ratio) %>%
  mutate(speaker_condition = factor(speaker_condition,
                                    levels=c("A,C certain", "A,C uncertain",
                                             "A (C) certain, C (A) uncertain"))) %>%
  chunk_cns()

p = plot_speaker(df, "bottom", facets=TRUE, xlab="proportion", ylab="best utterance")
ggsave(paste(params.speaker$plot_dir, "speaker_freq_best_un_certain_other.png",
             sep=SEP), p)

# Figure 5 ----------------------------------------------------------------
plot_evs_cp <- function(data, facets){
  df = data %>%
    mutate(kind=case_when(val == "p" ~ "CP-propositions",
                          val_type == "A,C indep." ~ "cn.ind",
                          val_type %in% c("A -> C", "C -> A") ~ "cns.dep.pos", 
                          val_type %in% c("A -> ¬C", "C -> ¬A") ~ "cns.dep.neg"),
           val=recode(val, "cns"="causal nets", "p"="CP-probabilities"))
  
  p <- df %>% ggplot(aes(y=level, x=ev, fill=val_type)) + 
    geom_bar(position=position_stack(), stat="identity") +
    # scale_y_discrete(name=paste(strwrap("values", width=20),
    #                             collapse="\n")) +
    scale_fill_discrete(name="value") +
    #                     breaks=c("prior", "LL", "PL"),
    #                     labels=c("A priori", "Literal", "Pragmatic")) +
    labs(x="Degree of belief", y="Interpretation level") +
    theme_bw() +
    theme(legend.position="top") +
    facet_wrap(~val, labeller=as_labeller(facets))
    # theme(legend.position = c(.95, .95), legend.justification=c("right", "top"))
  return(p)
}

data_cp_plots <- function(params){
  data <- read_rds(params$target) %>%
    ungroup() %>% group_by(bn_id, level) %>% select(-p_delta, -p_rooij, -p_diff, -bias)
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
  
data.cp.none = data_cp_plots(params)
facets = list("cns"="causal nets", "p"="CP-probabilities")
p <- plot_evs_cp(data.cp.none, facets) 
  
ggsave(paste(params$plot_dir, "none-evs-cp.png", sep=SEP), p, width=6.7, height=3)


# Figure 6 ----------------------------------------------------------------
plot_accept_conditions <- function(dat, fn){
  x_labels <- c(c(-100000, -1000, -50, -10, -1),
                seq(from=-0.75, by=0.25, to=0.7),
                seq(from=0.75, by=0.05, to=1))
  x_breaks <- seq(from=0.5, by=1, length.out=length(x_labels))
  level_labels = as_labeller(
    c(`literal-speaker` = paste(strwrap("Literal speaker", width=15), collapse="\n"),
      `pragmatic-speaker` = paste(strwrap("Pragmatic speaker", width=15), collapse="\n"),
      `prior` = "Prior")
  )
  condition_labels = as_labeller(c(`p_rooij`= "△*P", `p_delta`="△P", `p_diff`="P(C|A)-P(C)"))
  cn_labels = as_labeller(dat$cn %>% unique())
  # x-axis binned into different non-linear intervals (groups)
  getGroup = function(val) {
    i<-1
    x <- x_labels[[1]]
    while(val > x){
      i <- i + 1
      x <- x_labels[[i]]
    } 
    return(i-1)
  }
  
  data <- dat %>% group_by(bn_id, level, condition) %>%
    mutate(group = getGroup(val)) %>%
    group_by(level, condition, group, cn) %>%
    mutate(n_sampled=as.numeric(n_sampled)) %>% 
    mutate(n_sampled=case_when(is.na(n_sampled) ~ 1,
                               TRUE ~ n_sampled))
  data.sum <- data %>% 
    summarise(count=sum(n_sampled), .groups="drop_last") %>%
    group_by(level, condition) %>%
    mutate(ratio = count/sum(count))
  
  p <- data.sum %>% ggplot(aes(x=group, fill=cn)) +
    facet_grid(level~condition, scales="free", labeller=
                 labeller(level=level_labels, condition=condition_labels))
  
  p <- p + geom_bar(aes(y=ratio), stat="identity", position="dodge") +
    scale_x_continuous(breaks=x_breaks, labels=x_labels) +
    labs(x=paste(strwrap("accept/assert condition value intervals", width=25), collapse="\n"),
         y="ratio", fill="causal net") +
    theme_bw() +
    theme(legend.position="bottom", axis.text.x=element_text(angle=0, size=11))
  
  ggsave(paste(params$plot_dir, fn, sep=SEP), p, width=18, height=10)
  return(p)
}

params.prior <- read_rds(paste(data_dir, "params-none-priorN.rds", sep=SEP))
prior <-  read_rds(params.prior$target) %>%
  pivot_wider(names_from = "cell", values_from = "val") %>% 
  mutate(level = "prior") %>% select(-bias)
params.sp_literal <- read_rds(paste(data_dir, "params-speaker-literal.rds", sep=SEP))

formatSpeaker <- function(params, cat) {
  bn_ids = read_rds(paste(params$target_dir, .Platform$file.sep, "sample-ids-",
                          params$target_fn, sep="")) %>% 
    group_by_all() %>% summarize(n_sampled=n(), .groups="drop_last")
  speaker <- read_rds(file.path(params$target_dir, params$target_fn)) %>%
    select(-bias) %>% group_by(bn_id, cn) %>%
    mutate(utterance = paste("utt", utterance, sep="_")) %>% 
    pivot_wider(names_from = utterance, values_from=probs) %>%
    mutate(level = cat)
  sp = left_join(speaker, bn_ids, by=c("stimulus_id"))
  return(sp)
}

# get data from implemented literal speaker (conditioned s.t. A > C is true)
speaker.literal <- formatSpeaker(params.sp_literal, "literal-speaker")
speaker.literal.best = speaker.literal %>% 
  pivot_longer(cols = starts_with("utt_"), names_to = "utterance",
             values_to = "probs") %>% 
  group_by(bn_id, level, cn) %>% 
  mutate(p_best=max(probs), u_best=list(utterance[probs == max(probs)])) %>%
  unnest(u_best)

speaker.pragmatic <- speaker.literal.best %>%
  filter(u_best == "utt_A > C") %>%
  mutate(level = "pragmatic-speaker") %>%
  select(-u_best, -p_best) %>%
  pivot_wider(names_from = "utterance", values_from = "probs") 

speakers = bind_rows(speaker.literal, speaker.pragmatic)

df <- bind_rows(prior, speaker.literal, speaker.pragmatic) %>%
  group_by(bn_id, cn, level) %>%
  pivot_longer(cols=c(p_delta, p_rooij, p_diff), names_to="condition", values_to = "val") %>% 
  select(-starts_with("utt")) %>%
  mutate(
    level=factor(level, levels = c("prior", "literal-speaker", "pragmatic-speaker")), 
    cn = case_when(cn == "A || C" ~ "A,C indep.",
                   TRUE ~str_replace(cn, "-", "¬")
    ),
    condition = factor(condition, levels = c("p_rooij", "p_delta", "p_diff"))
  ) %>% 
  group_by(bn_id, level, condition)

plot_accept_conditions(df %>% filter(condition != "p_diff"), "accept-conditions.png")

# -------------- additional checks on tables of pragmatic speaker condition
speaker.pragmatic.long = speaker.pragmatic %>% 
   pivot_longer(cols=starts_with("utt_"), names_to="utterance", values_to="probs")
 
prag.ind= speaker.pragmatic.long %>% filter(cn=="A || C" & probs>0) 
# analyze prior given tables where no literal no conjunction is true
prior = readRDS(params$target) %>% filter(level=="prior") %>%
  pivot_wider(names_from="cell", values_from="val")

prior.unc = prior %>%
  mutate(conj=AC>theta | `A-C`>theta | `-AC`>theta | `-A-C`>theta,
         lit=(AC+`A-C` > theta) | (AC+`-AC`>theta) |
           (AC+`A-C` < 1-theta) | (AC+`-AC`<1-theta)) %>%
  filter(!conj &!lit) %>% ungroup() %>%
  filter(AC/(AC+`A-C`) > theta)

ids = prior.unc %>% filter(cn=="A || C") %>% pull(bn.id)
ratio = nrow(prag.ind) / nrow(prag.ind %>% filter(stimulus_id %in% ids))
print(paste(ratio, 'cant say conjunction or literal'))

speaker.pragmatic %>% filter(p_rooij>0.9) %>% pull(n_sampled) %>% sum() /
  speaker.pragmatic %>% pull(n_sampled) %>% sum()
speaker.literal %>% filter(p_rooij<0) %>% pull(n_sampled) %>% sum() /
  speaker.literal %>% pull(n_sampled) %>% sum()

prior %>% filter(p_rooij < 0) %>% ungroup() %>% nrow() / prior  %>% nrow()
prior %>% filter(p_rooij > 0.9) %>% ungroup() %>% nrow() / prior  %>% nrow()


# Check almost true states Appendix
df= speaker.pragmatic %>% filter(cn=="A || C") %>% select(-starts_with("utt_")) %>% 
  ungroup()
df.lit = df %>% group_by(bn_id) %>% 
  mutate(pa=`AC`+`A-C`, pc=`AC`+`-AC`, pna=`-AC`+`-A-C`, pnc=`A-C`+`-A-C`) %>%
  select(bn_id, pa, pc, pna, pnc) %>% 
  pivot_longer(cols=c(-bn_id), names_to="p", values_to="val") %>%
  arrange(desc(val)) %>% 
  distinct_at(vars(c("bn_id")), .keep_all = TRUE)
df.almost_true = df.lit %>% filter(val>=0.85) 


# Figure 7 ----------------------------------------------------------------

# speaker p_rooij filtered from literal speaker condition
speaker.p_rooij.best.from_lit = speaker.literal %>% filter(p_rooij>0.9) %>% 
  group_by(stimulus_id, cn) %>%
  pivot_longer(cols=starts_with("utt_"), names_prefix="utt_", names_to="utterance", values_to="probs") %>%
  mutate(p_best=max(probs), u_best=list(utterance[probs == max(probs)])) %>%
  unnest(u_best) %>% select(-p_best)

# frequency utterance type is the speaker's best utterance
# GIVEN A > C is NOT the speaker's best utterance
df <- speaker.p_rooij.best.from_lit %>%
  filter(utterance == u_best) %>% select(-u_best) %>%
  filter(utterance != "A > C") %>%
  chunk_utterances(c("-C > -A", "-A > -C", "C > A"))

df.sum <- df %>% group_by(utterance, cn) %>%
  summarise(p=sum(n_sampled), .groups = "drop_last") %>% arrange(p) %>% 
  mutate(N=sum(p), p=p/N)

p = plot_speaker(df.sum, "bottom", facets=FALSE, xlab="proportion",
                 ylab="best utterance")

ggsave(paste(params.speaker$plot_dir,
             "literal-speaker_prooij_large_freq_best_not_ac.png",
             sep=SEP), p, width=13.5, height=4.5)

