# Check Read Loss - All

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
filt <- read.delim("readqc/k2filtenterostats.txt", sep = "\t")
pb <- read.delim("readqc/pbrawstats.txt")

# ----- Combine and create summary table ----
#updated to incorporate long reads
all <- raw %>%
  mutate(ID = str_split_fixed(file, "\\_", 2)[,1]) %>%
  mutate(read = case_when(
    str_detect(file, "_1.fq.gz") ~ "forward",
    str_detect(file, "_2.fq.gz") ~ "reverse"
  )) %>%
  group_by(ID) %>%
  mutate(raw_sumlen = sum(sum_len)) %>%
  mutate(raw_numseqs = sum(num_seqs)) %>%
  mutate(raw_minlen = min(min_len)) %>%
  mutate(raw_avglen = mean(avg_len)) %>%
  mutate(raw_maxlen = max(max_len)) %>%
  select(ID, read, raw_numseqs, raw_sumlen, raw_minlen, raw_avglen, raw_maxlen) %>%
  inner_join(
    trim %>%
      mutate(ID = str_split_fixed(file, "\\_", 2)[,1]) %>%
      filter(str_detect(file, "P")) %>%
      mutate(read = case_when(
        str_detect(file, "_1P.fq.gz") ~ "forward",
        str_detect(file, "_2P.fq.gz") ~ "reverse"
      )) %>%
      group_by(ID) %>%
      mutate(trim_sumlen = sum(sum_len)) %>%
      mutate(trim_numseqs = sum(num_seqs)) %>%
      mutate(trim_minlen = min(min_len)) %>%
      mutate(trim_avglen = mean(avg_len)) %>%
      mutate(trim_maxlen = max(max_len)) %>%
      select(ID, read, trim_numseqs, trim_sumlen, trim_minlen, trim_avglen, trim_maxlen), by = c("ID", "read")
  ) %>%
  inner_join(
    filt %>%
      mutate(ID = str_split_fixed(file, "\\.", 2)[,1]) %>%
      mutate(read = case_when(
        str_detect(file, "_1P.fq") ~ "forward",
        str_detect(file, "_2P.fq") ~ "reverse"
      )) %>%
      group_by(ID) %>%
      mutate(filt_sumlen = sum(sum_len)) %>%
      mutate(filt_numseqs = sum(num_seqs)) %>%
      mutate(filt_minlen = min(min_len)) %>%
      mutate(filt_avglen = mean(avg_len)) %>%
      mutate(filt_maxlen = max(max_len)) %>%
      select(ID, read, filt_numseqs, filt_sumlen, filt_minlen, filt_avglen, filt_maxlen), by = c("ID", "read")
  ) %>%
  ungroup() %>%
  select(-read) %>%
  unique() %>%
  inner_join(pb %>%
               mutate(ID = str_remove_all(file, "\\.fastq.gz")) %>%
               mutate(pb_sumlen = sum(sum_len)) %>%
               mutate(pb_numseqs = sum(num_seqs)) %>%
               mutate(pb_minlen = min(min_len)) %>%
               mutate(pb_avglen = mean(avg_len)) %>%
               mutate(pb_maxlen = max(max_len)) %>%
               select(ID, pb_numseqs, pb_sumlen, pb_minlen, pb_avglen, pb_maxlen), by = c("ID")) %>%
  mutate(totalassembseqs = filt_numseqs + pb_numseqs) %>%
  mutate(avgcoverage = ((filt_numseqs*filt_avglen) + (pb_numseqs*pb_avglen))/4857000) %>%
  mutate(loss_rawtotrim = (raw_numseqs-trim_numseqs)/raw_numseqs) %>%
  mutate(loss_trimtofilt = (trim_numseqs-filt_numseqs)/trim_numseqs) %>%
  mutate(loss_rawtofilt = (raw_numseqs-filt_numseqs)/raw_numseqs)
               

#write.table(all, "readqc/readloss_throughk2filtentero_withpb.txt")

# ---- Plot ---- 

all.m <- reshape2::melt(all)


ggarrange(
  ggplot(all.m %>%
           filter(str_detect(variable, "numseqs")), aes(x=variable, y=value))+
    theme_classic() +
    geom_boxplot() +
    #geom_jitter(width = 0.1, aes(colour = read)) +
    theme(axis.title= element_blank(),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_blank()), 
  ggplot(all.m %>%
           filter(str_detect(variable, "sumlen")), aes(x=variable, y=value))+
    theme_classic() +
    geom_boxplot() +
    #geom_jitter(width = 0.1, aes(colour = read)) +
    theme(axis.title= element_blank(),
          axis.text = element_text(size = 12)), common.legend = TRUE)

# ---- Check Post Kraken2 Coverage and Depth ----- 
#completed before with short read only - no longer needed
#cov1 <- filt %>%
#  mutate(ID = str_split_fixed(file, "\\.", 2)[,1]) %>%
# mutate(read = case_when(
#   str_detect(file, "_1P.fq") ~ "forward",
#   str_detect(file, "_2P.fq") ~ "reverse"
# )) %>%
# mutate(totalseq = avg_len*num_seqs) %>%
# group_by(ID) %>%
# mutate(totalseq_bothreads = sum(totalseq)) %>%
# ungroup() %>%
# select(ID, totalseq_bothreads) %>%
# unique() %>%
# mutate(coverage = totalseq_bothreads/4857000) # taken from LT2 ATCC reference length

# saved as "readqc/enterodepthall.txt"

# how many are under 30x coverage 
#nrow(cov1 %>%
#      filter(coverage <=30)) # 1 

# how many are under 50x coverage
#nrow(cov1 %>%
#      filter(coverage <50)) # 1

# which are under 30X
#sub30 <- cov1 %>%
# filter(coverage <=30)

#covno16 <- cov1 %>%
# filter(ID != "ARS16")

#summary(covno16$coverage)

# ---- Descriptive stats ----- 
#nolonger needed
#no16 <- all %>%
#filter(ID != "ARS16")

#update to reflect keeping ARS16
summarydf <- as.data.frame(t(summary(all[,sapply(all, is.numeric)])))

summarydf <- summarydf %>%
  mutate(measure = str_split_fixed(Freq, ":", 2)[,1]) %>%
  mutate(value = round(as.numeric(str_split_fixed(Freq, ":", 2)[,2]), 2)) %>%
  select(-c(Var2, Freq))

summarydf <- summarydf %>%
  pivot_wider(names_from = measure)

colnames(summarydf)[1] <- "Parameter"
summarydf <- summarydf[,1:7]

#write.table(summarydf, "readqc/readsummarystats_withpb.txt")

