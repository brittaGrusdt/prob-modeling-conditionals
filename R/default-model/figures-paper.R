library(tidyverse)
library(config)
library(ggplot2)
library(latex2exp)
library(cowplot)
library(grid)
source("R/helper-functions.R")
source("R/default-model/helpers-tables.R")

# target <- "my-config"
target <- "paper-config"

# seed <- "seed-1625405121"
# seed <- "seed-1625403074"
seed <- ""
  
# setup -------------------------------------------------------------------
target_seed_dir = function(params, seed){
  subdirs = str_split(params$target, SEP)[[1]]
  if(seed != "") params$target <- paste(params$target_dir, seed,
                                        subdirs[[length(subdirs)]], sep=SEP)
  return(params)
}

# load parameters
SEP = .Platform$file.sep
data_dir = here("data", "default-model", target)
if(seed != "") data_dir <- paste(data_dir, SEP, seed, sep="")
plot_dir = paste(data_dir, "figs", sep=SEP)
if(!dir.exists(plot_dir)) dir.create(plot_dir)

params <- read_rds(paste(data_dir, "params-prior-ll-pl.rds", sep=SEP))
params <- target_seed_dir(params, seed)
theta=params$theta
data <- read_rds(params$target)
                                     
params.speaker <- read_rds(paste(data_dir, "params-speaker.rds", sep=SEP))
params.speaker <- target_seed_dir(params.speaker, seed)

params.prior <- read_rds(paste(data_dir, "params-priorN.rds", sep=SEP))
params.prior <- target_seed_dir(params.prior, seed)

# get data from implemented literal speaker (conditioned s.t. A > C is true)
params.sp_literal <- read_rds(paste(data_dir, "params-speaker-literal.rds",sep=SEP))
params.sp_literal <- target_seed_dir(params.sp_literal, seed)

# Figure 2 ----------------------------------------------------------------
tables = data %>% filter(level == "prior") %>% 
  ungroup() %>% dplyr::select(cn, bn_id, cell, val) %>% 
  group_by(bn_id) %>% pivot_wider(names_from="cell", values_from="val") %>% 
  mutate(vs=list(c("AC", "A-C", "-AC", "-A-C")),
         ps=list(c(`AC`, `A-C`, `-AC`, `-A-C`))) %>%
  select(-`AC`, -`A-C`, -`-AC`, -`-A-C`) %>% 
  mutate(bn_id.tmp=bn_id) %>% 
  separate(bn_id, into=c("cn", "t0", "t1", "t2", "t3"), sep="_") %>% 
  unite("table_id", "t0", "t1", "t2", "t3", sep="_") %>% 
  rename(bn_id=bn_id.tmp)
plot_tables_cns("", plot_dir, w=3, h=3.5, tables)

# Figure 3 ----------------------------------------------------------------
data.speaker <- read_rds(params.speaker$target) %>%
  select(-level, -p_delta, -p_diff) %>% ungroup() %>% 
  distinct_at(vars(c(bn_id, utterance)), .keep_all = T)
# bn_ids of sampled states with nb of occurrence in overall set of samples
bn_ids = read_rds(str_replace(params.speaker$target, "results", "sample-ids")) %>% 
  group_by_all() %>% summarize(n_sampled=n(), .groups="drop_last")

data.speaker.best <- data.speaker %>% group_by(bn_id) %>%
  mutate(p_best=max(probs), u_best=list(utterance[probs == max(probs)])) %>%
  unnest(c(u_best)) %>% select(-p_best) %>% 
  filter(utterance == u_best)
data.speaker.best = left_join(data.speaker.best, bn_ids, by=c("bn_id"))

sp.best.conditions <- data.speaker.best %>%
  mutate(pa=`AC` + `A-C`, pc=`AC` + `-AC`,
         certainA = pa >= theta | pa <= 1-theta,
         certainC = pc >= theta | pc <= 1-theta,
         certain=certainA & certainC,
         uncertain=!certainA & !certainC, 
         speaker_condition = case_when(
           certain ~ "certain",
           uncertain ~ "uncertain",
           T ~ "xor"
         ),
         sp_condition=speaker_condition) %>% 
  select(-certainA, -certainC, -certain, -uncertain) %>%
  chunk_utterances()

df <- sp.best.conditions %>%
  chunk_cns() %>% 
  group_by(speaker_condition, cn, utterance) %>%
  dplyr::select(cn, utterance, speaker_condition, n_sampled) %>%
  summarize(n_sampled=sum(n_sampled), .groups="drop_last") %>%
  mutate(N=sum(n_sampled), p=n_sampled/N) %>%
  mutate(cn = factor(cn, levels=c("A,C independent", "A,C dependent")),
         speaker_condition=factor(speaker_condition,
                                  levels=c("certain", "uncertain", "xor")))
levels(df$speaker_condition) <-
  c("certain" = expression(atop("A,C certain")),
    "uncertain" = expression(atop("A,C uncertain")),
    "xor" = expression(atop("A xor C certain")))

p = plot_speaker_conditions(df)
ggsave(paste(plot_dir, "speaker_freq_best_un_certain_other.png", sep=SEP),
       width=7, height=2.5)
# some checks
dat.utts = table_to_utts(sp.best.conditions, params$theta)
# 1. check A,C uncertain + independent + best utterance is conditional
check = sp.best.conditions %>% 
  filter(sp_condition=="uncertain" & utterance == "conditional" & cn=="A || C") %>%
  ungroup() %>% dplyr::select(bn_id) %>% distinct() %>% pull(bn_id)
dat.utts %>% filter(bn_id %in% check)

# Figure 4+5 ----------------------------------------------------------------
data.cp = data_cp_plots(params) 
cp.cns = data.cp %>% filter(val=="cns") %>% 
  mutate(val_type=factor(val_type, levels=c("A,C indep.", "A -> C",
                                            "A -> ¬C", "C -> A", "C -> ¬A")))
p.probs <- data.cp %>% filter(val=="p") %>% plot_cp_probs()
p.cns <- cp.cns %>% plot_cp_cns(labels.cp)
ggsave(paste(plot_dir, "cp-evs-probs.png", sep=SEP), p.probs,width=7, height=2.5)
ggsave(paste(plot_dir, "cp-evs-cns.png", sep=SEP), p.cns, width=7, height=2.5)

# Figure 6 ----------------------------------------------------------------
prior <-  read_rds(params.prior$target) %>%
  pivot_wider(names_from = "cell", values_from = "val") %>% 
  mutate(level = "prior") %>%
  pivot_longer(cols=c(p_delta, p_rooij, p_diff),names_to="condition", values_to = "val")

# speaker literal: only total of 18 utterances used (2 not applicable) 
speaker.literal <- read_rds(params.sp_literal$target) %>%
  mutate(utterance = paste("utt", utterance, sep="_"),
         level="literal-speaker") %>% 
  group_by(rowid) %>% 
  mutate(p_best=max(probs), u_best = probs == p_best) %>%
  pivot_longer(cols=c(p_delta, p_rooij, p_diff),names_to="condition",
               values_to = "val")
# pragmatic speaker condition: best utterance is A > C
speaker.literal.best = speaker.literal %>%
  filter(u_best) %>% dplyr::select(-p_best, -u_best)
speaker.pragmatic <- speaker.literal.best %>%
  filter(utterance == "utt_A > C") %>% 
  mutate(level = "pragmatic-speaker")

speakers = bind_rows(speaker.literal, speaker.pragmatic) %>% group_by(rowid,level)
df <- bind_rows(prior, speakers) %>%
  mutate(
    level=factor(level, levels = c("prior", "literal-speaker", "pragmatic-speaker")), 
    cn = case_when(cn == "A || C" ~ "A,C indep.",
                   TRUE ~str_replace(cn, "-", "¬")
    ),
    condition = factor(condition, levels = c("p_rooij", "p_delta", "p_diff"))
  ) %>% 
  group_by(rowid, level, condition)

plot_accept_conditions <- function(dat){
  x_labels <- c(c(-100000, -1000, -50, -10, -1),
                seq(from=-0.75, by=0.25, to=0.7),
                seq(from=0.75, by=0.05, to=1))
  x_breaks <- seq(from=0.5, by=1, length.out=length(x_labels))
  level_labels = as_labeller(
    c(`literal-speaker` = paste(strwrap("Literal speaker", width=25), collapse="\n"),
      `pragmatic-speaker` = paste(strwrap("Pragmatic speaker", width=25), collapse="\n"),
      `prior` = "Prior")
  )
  condition_labels = as_labeller(c(`p_rooij`= "△*P",
                                   `p_delta`="△P",
                                   `p_diff`="P(C|A)-P(C)"))
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
  data <- dat %>% rowid_to_column("idx") %>% group_by(idx) %>% 
    mutate(group = getGroup(val))
  # group_by(rowid, level, condition) %>%
  data.sum <- data %>% group_by(level, condition, group, cn) %>% 
    summarise(count=n(), .groups="drop_last") %>%
    group_by(level, condition) %>%
    mutate(ratio = count/sum(count))
  
  p <- data.sum %>% ggplot(aes(x=group, fill=cn)) +
    facet_wrap(~level, scales="free_x", labeller=
                 labeller(level=level_labels)) +
    geom_bar(aes(y=ratio), stat="identity", position="dodge") +
    scale_x_continuous(breaks=x_breaks, labels=x_labels) +
    labs(x=TeX("$\\Delta^{*}  P$"),
         y="relative frequencey", fill="causal net") +
    scale_fill_brewer(palette="Dark2", labels=labels.cp) +
    theme_minimal() +
    theme(legend.position="top",
          axis.text.x=element_text(angle=40, size=6),
          legend.key.size = unit(0.75,"line"))
  return(p)
}
p = plot_accept_conditions(df %>% filter(condition == "p_rooij"))
ggsave(paste(plot_dir, "accept-conditions.png", sep=SEP), p, width=7, height=3)

# -------------- additional checks on tables of pragmatic speaker condition
prag.ind= speaker.pragmatic %>% filter(cn=="A || C" & probs>0) 
# analyze prior given tables where no literal no conjunction is true
prior = data %>% filter(level=="prior") %>%
  pivot_wider(names_from="cell", values_from="val")

# A->C is applicable and neither conj nor literal true
prior.unc = prior %>%
  mutate(conj=AC>theta | `A-C`>theta | `-AC`>theta | `-A-C`>theta,
         lit=(AC+`A-C` > theta) | (AC+`-AC`>theta) |
           (AC+`A-C` < 1-theta) | (AC+`-AC`<1-theta)) %>%
  filter(!conj &!lit) %>% ungroup() %>%
  filter(AC/(AC+`A-C`) > theta)

ids = prior.unc %>% filter(cn=="A || C") %>% pull(bn_id)
ratio = nrow(prag.ind %>% filter(bn_id %in% ids)) / nrow(prag.ind)
print(paste("in", ratio*100, '% of cases, neither conjunction nor literal applicable'))

# some checks
speaker.pragmatic %>% filter(condition=="p_rooij" & val >0.9) %>% nrow() /
  (speaker.pragmatic %>% filter(condition=="p_rooij") %>% nrow())
speaker.literal %>% filter(condition=="p_rooij" & val<0) %>% nrow() /
  (speaker.literal %>% filter(condition=="p_rooij") %>% nrow())

prior %>% filter(p_rooij < 0) %>% ungroup() %>% nrow() / (prior  %>% nrow())
prior %>% filter(p_rooij > 0.9) %>% ungroup() %>% nrow() / (prior  %>% nrow())

# Check almost true states Appendix
df= speaker.pragmatic %>% filter(cn=="A || C") %>% select(-utterance) %>% 
  ungroup()
df.lit = df %>% group_by(bn_id) %>% 
  mutate(pa=`AC`+`A-C`, pc=`AC`+`-AC`, pna=`-AC`+`-A-C`, pnc=`A-C`+`-A-C`) %>%
  select(bn_id, pa, pc, pna, pnc) %>% 
  pivot_longer(cols=c(-bn_id), names_to="p", values_to="val") %>%
  arrange(desc(val)) %>% 
  distinct_at(vars(c("bn_id")), .keep_all = TRUE)
df.almost_true = df.lit %>% filter(val>=0.85) 

nrow(df.almost_true) / df.lit %>% nrow()

# Figure 7 ----------------------------------------------------------------
# V1: speaker results (not literal speaker condition, just samples from prior)
# df.sp <- read_rds(params.speaker$target) %>%
#   select(-level, -p_delta, -p_diff) %>% ungroup() %>% 
#   mutate(utterance = paste("utt", utterance, sep="_")) %>% 
#   filter(AC/(AC+`A-C`)>=0.8) %>% 
#   group_by(rowid) %>% 
#   mutate(p_best=max(probs), u_best = probs == p_best) %>% 
#   filter(u_best) %>% 
#   filter(p_rooij >= 0.8 & utterance != "utt_A > C") %>% 
#   select(-p_best, -u_best) %>%
#   mutate(utterance = str_replace(utterance, "utt_", "")) %>% 
#   chunk_utterances(c("-C > -A", "-A > -C", "C > A"))

# V2: best or second best utterance if best is not A>C
# df.sp = speaker.literal %>%
#   filter(condition=="p_rooij" & val >= 0.9 & utterance != "utt_A > C") %>%
#   group_by(rowid) %>%
#   mutate(p_best=max(probs), u_best = probs == p_best) %>%
#   filter(u_best) %>%
#   select(-p_best, -u_best) %>%
#   mutate(utterance = str_replace(utterance, "utt_", "")) %>%
#   chunk_utterances(c("-C > -A", "-A > -C", "C > A"))

# speaker results for states from literal speaker condition filtered s.t.
# p_rooij>0.9 and the speaker's best utterance is NOT A->C
df.sp = speaker.literal %>%
  filter(condition=="p_rooij" & u_best) %>%
  filter(val >= 0.9 & utterance != "utt_A > C") %>%
  select(-p_best, -u_best) %>%
  mutate(utterance = str_replace(utterance, "utt_", "")) %>%
  chunk_utterances()

# frequency of each utterance type for these states
df.sum <- df.sp %>% group_by(utterance, cn) %>%
  summarise(count=n(), .groups = "drop_last") %>%
  mutate(freq=count/sum(count), count.utt=sum(count)) %>% ungroup() %>%
  mutate(N=sum(count), freq.utt=count.utt/sum(count))
p <- df.sum %>%
  mutate(utterance = as.character(utterance), 
         utterance = case_when(utterance == "conditional" ~ "other conditional",
                               T ~ utterance),
         utterance = as.factor(utterance),
         cn=case_when(cn=="A || C" ~ "A,C indep.", T ~ cn),
         cn=factor(cn, levels=c("A,C indep.", "A implies C", "C implies A"))) %>% 
  ggplot(aes(y=utterance, x=count, fill=cn)) +
  geom_bar(stat="identity", position=position_stack())  +
  labs(x="count", y="best utterance") + theme_minimal() +
  theme(axis.text.y=element_text(), legend.position="top",
        legend.key.size = unit(0.75,"line")) +
  scale_fill_brewer(palette="Dark2", name="causal net", labels=labels.cp,
                    guide = guide_legend(reverse = TRUE))

ggsave(paste(plot_dir, "literal-speaker_prooij_large_freq_best_not_ac.png",
             sep=SEP), p, width=7, height=2.5)

# only independent tables sampled -----------------------------------------
freqs.best.utts.ind = data.speaker.best %>% 
  filter(utterance==u_best & cn=="A || C") %>%
  group_by(utterance) %>% 
  summarize(n=sum(n_sampled), .groups="drop")  %>% 
  mutate(N=sum(n), ratio=n/N, utt=utterance) %>% 
  chunk_utterances() 

# from independent states sampled from prior, ratio for each
# utterance type to be hyperrational speakers choice
freqs.best.utts.ind %>% group_by(utterance) %>% 
  summarize(p=sum(ratio))

# for specific utterance A > C
freqs.best.utts.ind %>% filter(utt == "A > C")
