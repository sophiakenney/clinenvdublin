# Check Read Loss

# ---- Load Packages ----
library(tidyverse)
library(tidylog)
library(reshape2)
library(ggpubr)

# Set working directory 
setwd("./R/")

# Load tables - 
raw <- read.delim("readqc/rawstats.txt", sep = "\t")
trim <- read.delim("readqc/trimstats.txt", sep = "\t")
filt <- read.delim("readqc/k2filtenterostats.txt", sep = "\t") # this is the only that was changed 9/20

# ----- Combine and create summary table ----
all <- raw %>%
  mutate(ID = str_split_fixed(file, "\\_", 2)[,1]) %>%
  mutate(read = case_when(
    str_detect(file, "_1.fq.gz") ~ "forward",
    str_detect(file, "_2.fq.gz") ~ "reverse"
  )) %>%
  mutate(raw_sumlen = sum_len) %>%
  mutate(raw_numseqs = num_seqs) %>%
  select(ID, read, raw_numseqs, raw_sumlen) %>%
  inner_join(
    trim %>%
      mutate(ID = str_split_fixed(file, "\\_", 2)[,1]) %>%
      filter(str_detect(file, "P")) %>%
      mutate(read = case_when(
        str_detect(file, "_1P.fq.gz") ~ "forward",
        str_detect(file, "_2P.fq.gz") ~ "reverse"
      )) %>%
      mutate(trim_sumlen = sum_len) %>%
      mutate(trim_numseqs = num_seqs) %>%
      select(ID, read, trim_numseqs, trim_sumlen), by = c("ID", "read")
  ) %>%
  inner_join(
    filt %>%
      mutate(ID = str_split_fixed(file, "\\.", 2)[,1]) %>%
      mutate(read = case_when(
        str_detect(file, "_1P.fq") ~ "forward",
        str_detect(file, "_2P.fq") ~ "reverse"
      )) %>%
      mutate(filt_sumlen = sum_len) %>%
      mutate(filt_numseqs = num_seqs) %>%
      select(ID, read, filt_numseqs, filt_sumlen), by = c("ID", "read")
  )

#saved as "readqc/readloss_throughk2filtentero.txt"

all %>%
  mutate(rawtotrim = (raw_numseqs-trim_numseqs)/raw_numseqs) %>%
  mutate(trimtofilt = (trim_numseqs-filt_numseqs)/trim_numseqs) %>%
  mutate(rawtofilt = (raw_numseqs-filt_numseqs)/raw_numseqs)

# ---- Plot ---- 

all.m <- reshape2::melt(all)


ggarrange(
  ggplot(all.m %>%
           filter(str_detect(variable, "numseqs")), aes(x=variable, y=value))+
    theme_classic() +
    geom_boxplot() +
    geom_jitter(width = 0.1, aes(colour = read)) +
    theme(axis.title= element_blank(),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_blank()), 
  ggplot(all.m %>%
           filter(str_detect(variable, "sumlen")), aes(x=variable, y=value))+
    theme_classic() +
    geom_boxplot() +
    geom_jitter(width = 0.1, aes(colour = read)) +
    theme(axis.title= element_blank(),
          axis.text = element_text(size = 12)), common.legend = TRUE)

# ---- Check Post Kraken2 Coverage and Depth ----- 
cov1 <- filt %>%
  mutate(ID = str_split_fixed(file, "\\.", 2)[,1]) %>%
  mutate(read = case_when(
    str_detect(file, "_1P.fq") ~ "forward",
    str_detect(file, "_2P.fq") ~ "reverse"
  )) %>%
  mutate(totalseq = avg_len*num_seqs) %>%
  group_by(ID) %>%
  mutate(totalseq_bothreads = sum(totalseq)) %>%
  ungroup() %>%
  select(ID, totalseq_bothreads) %>%
  unique() %>%
  mutate(coverage = totalseq_bothreads/4857000) # taken from LT2 ATCC reference length

# saved as "readqc/enterodepthall.txt"

# how many are under 30x coverage 
nrow(cov1 %>%
       filter(coverage <=30)) # 1 

# how many are under 50x coverage
nrow(cov1 %>%
       filter(coverage <50)) # 1

# which are under 30X
sub30 <- cov1 %>%
  filter(coverage <=30)

covno16 <- cov1 %>%
  filter(ID != "ARS16")

summary(covno16$coverage)

# ---- Descriptive stats ----- 

no16 <- all %>%
  filter(ID != "ARS16")

summarydf.no16 <- as.data.frame(t(summary(no16[,sapply(no16, is.numeric)])))

summarydf.no16 <- summarydf.no16 %>%
  mutate(measure = str_split_fixed(Freq, ":", 2)[,1]) %>%
  mutate(value = round(as.numeric(str_split_fixed(Freq, ":", 2)[,2]), 2)) %>%
  select(-c(Var2, Freq))

summarydf.no16 <- summarydf.no16 %>%
  pivot_wider(names_from = measure)

colnames(summarydf.no16)[1] <- "Parameter"
summarydf.no16 <- summarydf.no16[,1:7]

# saved as "readqc/readsummarystats_noARS16.txt"

