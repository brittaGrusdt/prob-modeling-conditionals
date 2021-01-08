library(tidyverse)
library(config)
library(ggplot2)
library(latex2exp)
library(cowplot)
library(grid)
source("R/helper-functions.R")
source("R/default-model/helpers-tables.R")

data_dir = here("data", "default-model")
plot_dir = here("figs")
params <- read_rds(paste(data_dir, "params-none.rds", sep=SEP))
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

params.speaker <- read_rds(paste(data_dir, "params-speaker.rds", sep=SEP))
UTTERANCES <- read_rds(paste(params.speaker$target_dir, params.speaker$utts_fn, sep=SEP))
THETA = 0.9

# Figure 2 ----------------------------------------------------------------
# for plotting densities sample more tables 
params.tables = list(
  seed_tables=params$seed_tables, cns=params$cns, indep_sigma=params$indep_sigma,
  n_tables=10000, n_ind_tables=10000,
  tables_path=paste(params$target_dir, "tables-default-samples-for-plots.rds",
                    sep=.Platform$file.sep)
)
tables.plot = create_tables(params.tables)
plot_tables_cns(params.tables$tables_path, plot_dir, w=5, h=5)


# Figure 3 ----------------------------------------------------------------
# pragmatic/literal interpretations of conditional If A, C
# probability assigned to states where speaker is uncertain about A and C
dat.none.voi <- read_rds(
  paste(params$target_dir, "results-none-prior-LL-PL-vois.rds",
        sep=.Platform$file.sep)) %>%
  filter((key == "uncertain_both") & level!="prior") %>% 
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
  theme_classic(base_size=25) +
  theme(legend.position = "none") +
  coord_flip()
p
ggsave(paste(plot_dir, "ignorance-inferences.png", sep=SEP), p,
       width=13.5, height=4)



# Figure 4 ----------------------------------------------------------------
data.speaker <- read_rds(params.speaker$target) %>% select(-level, -bias) %>%
  select(-p_delta, -p_diff)
data.speaker.best <- data.speaker %>% group_by(bn_id) %>%
  mutate(p_best=max(probs), u_best=list(utterance[probs == max(probs)])) %>%
  unnest(c(u_best)) %>% select(-p_best)

dat <- data.speaker.best %>%
  mutate(pa=`AC` + `A-C`, pc=`AC` + `-AC`,
         certainA = pa >= THETA | pa <= 1-THETA,
         certainC = pc >= THETA | pc <= 1-THETA,
         uncA = pa > 1 - THETA & pa < THETA,
         uncC = pc > 1-THETA & pc < THETA,
         certain=certainA & certainC,
         uncertain=uncA & uncC, 
         both=!certain & !uncertain,
         speaker_condition = case_when(certain ~ "certain",
                                       uncertain ~ "uncertain",
                                       both ~ "both")) %>% 
  select(-uncA, -uncC, -certainA, -certainC, -certain, -uncertain) %>%
  select(-utterance) %>%
  rename(utterance = u_best) %>% distinct(bn_id, .keep_all = TRUE) %>%
  chunk_utterances()

df <- dat %>% group_by(speaker_condition, cn, utterance) %>%
  summarise(p=n(), .groups = "drop_last") %>% arrange(p) %>% 
  mutate(N=sum(p), ratio=p/N) %>% rename(n=p, p=ratio) %>%
  mutate(speaker_condition = factor(speaker_condition,
                                    levels=c("certain", "uncertain", "both"))) %>%
  chunk_cns()


plot_speaker(df, "speaker_freq_best_un_certain_other.png", w=13.5, h=5, plot_dir,
             "bottom", TRUE, "proportion", "best utterance")


# Figure 5 ----------------------------------------------------------------
plot_evs_cp <- function(data){
  lim = ifelse(max(data$ev) < 1-0.1, max(data$ev)+0.1, 1)
  p <- data %>% ggplot(aes(y=val_type, x=ev, fill=level)) + 
    geom_bar(position=position_dodge2(preserve = "single"), stat="identity") +
    scale_x_continuous(limits=c(0, lim)) +
    scale_y_discrete(name=paste(strwrap("causal nets / conditional probs", width=20),
                                collapse="\n")) +
    scale_fill_discrete(name="Interpretation Level",
                        breaks=c("prior", "LL", "PL"),
                        labels=c("A priori", "Literal", "Pragmatic")) +
    labs(x="Degree of belief") +
    theme_bw(base_size=25) +
    theme(legend.position = c(.95, .95), legend.justification=c("right", "top"))
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
p <- plot_evs_cp(data.cp.none)
ggsave(paste(plot_dir, "none-evs-cp.png", sep=SEP), p, width=16, height=8)


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
    group_by(level, condition, group, cn)
  data.sum <- data %>% 
    summarise(count=n(), .groups="drop_last") %>%
    group_by(level, condition) %>%
    mutate(ratio = count/sum(count))
  
  p <- data.sum %>% ggplot(aes(x=group, fill=cn)) +
    facet_grid(level~condition, scales="free", labeller=
                 labeller(level=level_labels, condition=condition_labels))
  
  p <- p + geom_bar(aes(y=ratio), stat="identity", position="dodge") +
    scale_x_continuous(breaks=x_breaks, labels=x_labels) +
    labs(x=paste(strwrap("accept/assert condition value intervals", width=25), collapse="\n"),
         y="ratio", fill="causal net") +
    theme_bw(base_size=25) +
    theme(legend.position="bottom", axis.text.x=element_text(angle=0, size=11))
  
  ggsave(paste(plot_dir, fn, sep=SEP), p, width=18, height=10)
  return(p)
}

params.prior <- read_rds(paste(data_dir, "params-none-priorN.rds", sep=SEP))
prior <-  read_rds(params.prior$target) %>%
  pivot_wider(names_from = "cell", values_from = "val") %>% 
  mutate(level = "prior") %>% select(-bias)
params.sp_literal <- read_rds(paste(data_dir, "params-speaker-literal.rds", sep=SEP))

formatSpeaker <- function(params, cat) {
  speaker <- read_rds(file.path(params$target_dir, params$target_fn)) %>%
    select(-bias) %>% group_by(bn_id, cn) %>%
    mutate(utterance = paste("utt", utterance, sep="_")) %>% 
    pivot_wider(names_from = utterance, values_from=probs) %>%
    mutate(level = cat)
  return(speaker)
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



# Figure 7 ----------------------------------------------------------------
params.speaker.p_rooij = read_rds(
  paste(data_dir, "params-speaker-p_rooij-large.rds", sep=SEP)
)
speaker.p_rooij <- read_rds(params.speaker.p_rooij$target) %>%
  select(-level, -bias, -p_delta, -p_diff)

speaker.p_rooij.best <- speaker.p_rooij %>% group_by(bn_id) %>%
  mutate(p_best=max(probs), u_best=list(utterance[probs == max(probs)])) %>%
  unnest(u_best) %>% select(-p_best)

# frequency utterance type is the speaker's best utterance
# GIVEN A > C is NOT the speaker's best utterance
df <- speaker.p_rooij.best %>%
  filter(utterance == u_best) %>% 
  select(-u_best) %>%
  distinct(bn_id, .keep_all = TRUE) %>%
  filter(utterance != "A > C" & utterance != "C > A") %>%
  chunk_utterances()

df.sum <- df %>% group_by(utterance) %>%
  summarise(p=n(), .groups = "drop_last") %>% arrange(p) %>% 
  mutate(N=sum(p), p=p/N)

plot_speaker(df.sum, "speaker_prooij_large_freq_best_not_ac.png", w=13.5, h=4,
             plot_dir, "bottom", FALSE, "proportion", "best utterance")









