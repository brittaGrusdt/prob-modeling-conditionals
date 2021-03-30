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
params <- read_rds(paste(data_dir, "params-none-prior-ll-pl.rds", sep=SEP))
theta=params$theta

params.speaker <- read_rds(paste(data_dir, "params-none-speaker.rds", sep=SEP))
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
plot_tables_cns(params.tables$tables_path, params$plot_dir, w=3, h=3)

#-- plot tables for samples from prior --#
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
  theme_minimal() +
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

data.speaker.best = left_join(data.speaker.best, bn_ids, by=c("bn_id"))

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

p = plot_speaker_conditions(df)
ggsave(paste(params.speaker$plot_dir, "speaker_freq_best_un_certain_other.png",
             sep=SEP), width=7, height=2.5)

# Figure 5 ----------------------------------------------------------------
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
  
data.cp = data_cp_plots(params) 
cp.cns = data.cp %>% filter(val=="cns") %>% 
  mutate(val_type=factor(val_type, levels=c("A,C indep.", "A -> C",
                                               "A -> ¬C", "C -> A", "C -> ¬A")))
labels.cp = tribble(~ val, ~label,
                    "A,C indep.", "A,C indep.",
                    "A -> C", expression(A %->% C),
                    "A -> ¬C", expression(A%->%~"¬C"),
                    "C -> A",expression(C %->% A),
                    "C -> ¬A", expression(C%->%~"¬A"))

p.probs <- data.cp %>% filter(val=="p") %>% 
  ggplot(aes(y=level, x=ev, fill=val_type)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  scale_fill_brewer(palette="Dark2", name="value") +
  labs(x="Degree of belief", y="Interpretation level") +
  theme_minimal() + theme(legend.position="top", legend.key.size = unit(0.75,"line"))

p.cns <- cp.cns %>% 
  ggplot(aes(y=level, x=ev, fill=val_type)) + 
  geom_bar(position=position_stack(), stat="identity") +
  scale_fill_brewer(palette="Dark2", name="value", labels=labels.cp$label) +
  labs(x="Degree of belief", y="Interpretation level") +
  theme_minimal() + theme(legend.position="top", legend.key.size = unit(0.75,"line"))

ggsave(paste(params$plot_dir, "cp-evs-probs.png", sep=SEP), p.probs, width=5, height=2)
ggsave(paste(params$plot_dir, "cp-evs-cns.png", sep=SEP), p.cns, width=5, height=2)


# Figure 6 ----------------------------------------------------------------
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
    facet_wrap(~level, scales="free", labeller=
                 labeller(level=level_labels)) +
    geom_bar(aes(y=ratio), stat="identity", position="dodge") +
    scale_x_continuous(breaks=x_breaks, labels=x_labels) +
    labs(x=paste(strwrap("accept/assert condition value intervals", width=25),
                 collapse="\n"),
         y="ratio", fill="causal net") +
    scale_fill_brewer(palette="Dark2") +
    theme_minimal() +
    theme(legend.position="top",
          axis.text.x=element_text(angle=40, size=8))
  return(p)
}

params.prior <- read_rds(paste(data_dir, "params-none-priorN.rds", sep=SEP))
prior <-  read_rds(params.prior$target) %>%
  pivot_wider(names_from = "cell", values_from = "val") %>% 
  mutate(level = "prior") %>% select(-bias) %>%
  pivot_longer(cols=c(p_delta, p_rooij, p_diff),names_to="condition", values_to = "val")

# get data from implemented literal speaker (conditioned s.t. A > C is true)
params.sp_literal <- read_rds(paste(data_dir, "params-speaker-literal.rds", sep=SEP))

speaker.literal <- read_rds(file.path(params$target_dir, params.sp_literal$target_fn)) %>%
  select(-bias) %>%
  mutate(utterance = paste("utt", utterance, sep="_"),
         level="literal-speaker") %>% 
  group_by(rowid) %>% 
  mutate(p_best=max(probs), u_best = probs == p_best) %>%
  pivot_longer(cols=c(p_delta, p_rooij, p_diff),names_to="condition", values_to = "val")

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

p = plot_accept_conditions(df %>% filter(condition == "p_rooij"))
ggsave(paste(params$plot_dir, "accept-conditions.png", sep=SEP), p, width=10, height=4)

# -------------- additional checks on tables of pragmatic speaker condition
prag.ind= speaker.pragmatic %>% filter(cn=="A || C" & probs>0) 
# analyze prior given tables where no literal no conjunction is true
prior = readRDS(params$target) %>% filter(level=="prior") %>%
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


# Figure 7 ----------------------------------------------------------------
# speaker results for states from literal speaker condition filtered s.t.
# p_rooij>0.9 and the speaker's best utterance is NOT A->C
df.sp = speaker.literal %>%
  filter(condition=="p_rooij" & val > 0.9 & u_best) %>%
  filter(utterance != "utt_A > C") %>%
  select(-p_best, -u_best) %>%
  mutate(utterance = str_replace(utterance, "utt_", "")) %>% 
  chunk_utterances(c("-C > -A", "-A > -C", "C > A"))

# frequency of each utterance type for these states
df.sum <- df.sp %>% group_by(utterance, cn) %>%
  summarise(count=n(), .groups = "drop_last") %>%
  mutate(freq=count/sum(count), count.utt=sum(count)) %>% ungroup() %>%
  mutate(N=sum(count), freq.utt=count.utt/sum(count))
p <- df.sum %>%
    mutate(cn=case_when(cn=="A || C" ~ "A,C indep.", T ~ cn),
           cn=factor(cn, levels=c("A,C indep.", "A implies C", "C implies A"))) %>% 
    ggplot(aes(y=utterance, x=count, fill=cn)) +
    geom_bar(stat="identity", position=position_stack())  +
    labs(x="count", y="best utterance") + theme_minimal() +
    theme(axis.text.y=element_text(), legend.position="top") +
    scale_fill_brewer(palette="Dark2")

ggsave(paste(params.speaker$plot_dir,
             "literal-speaker_prooij_large_freq_best_not_ac.png",
             sep=SEP), p, width=5, height=2)


# Fig listener surprisal ml -----------------------------------------------
data.speaker <- read_rds(params.speaker$target) %>% select(-level, -bias) %>%
  select(-starts_with("p_"))
sp.ind = data.speaker %>% filter(cn == "A || C")

freqs.best.utts = data.speaker %>% group_by(rowid) %>%
  mutate(best=probs==max(probs)) %>%
  filter(best) %>%  group_by(utterance, cn) %>% 
  chunk_utterances() %>% mutate(cn=case_when(cn=="A || C" ~ "A,C indep.",
                                             T ~ "A,C dep.")) %>% 
  summarize(freq=n()) %>% 
  arrange(freq) 

p = freqs.best.utts %>% ggplot(aes(x=utterance, y=freq, fill=cn)) +
  geom_bar(stat="identity") + theme_minimal() +
  theme(legend.position="top") +
  labs(y="frequency", x="best utterance") +
  scale_fill_brewer(name="causal net", palette="Dark2")

ggsave(paste(params.speaker$plot_dir,
             "listener-surprisal-speakers-utt-choice.png",
             sep=SEP), p, width=5, height=2)



