# Genome Assembly QC

# ---- Load Packages ----
library(tidylog)
library(tidyverse)
library(data.table)


# ---- Load Tables ----
quast <- read.delim("assembqc/report.tsv", sep = "\t")
checkm <- read.delim("assembqc/genphenentero_checkm.tsv", sep = "\t")
prk <- read.delim("assembqc/enteroprokka_summary.txt", sep = "\t")
ss <- read.delim("assembqc/enteroseqsero_summary.tsv", sep = "\t")
sr <- read.delim("assembqc/enterosistr_output.txt", sep = "\t")
st <- read.delim("assembqc/enteromlst.tsv", sep = "\t", header = FALSE)


# ---- Format to combine ----
# CheckM
colnames(checkm)

checkm2 <- checkm %>%
  mutate(ID=str_split_fixed(Bin.Id, "\\_", 2)[,1]) %>% # cut _assembly off name
  select(ID, Completeness, Contamination, Strain.heterogeneity) # select most relevant columns

# Prokka 
colnames(prk)

prk2 <- prk %>%
  mutate(ID=str_split_fixed(Filename, "\\.", 2)[,1]) %>% # cut .txt off name
  select(ID, contigs,bases,CDS,rRNA, tRNA, repeat_region, tmRNA) # drop the file name and organism columns

# QUAST 
colnames(quast)
quast2<- data.table::transpose(quast)
quast2 <- as.data.frame(apply(quast2, 2, as.numeric))
colnames(quast2) <- quast$Assembly
quast2$ID <- colnames(quast)
quast2$ID <- str_remove_all(quast2$ID, "_assembly")
rownames(quast2) <- quast2$ID
quast2 <- quast2[-1,]

# SISTR
colnames(sr)

sr2 <- sr %>%
  mutate(ID=str_split_fixed(genome, "\\_", 2)[,1]) %>% # cut _assembly off name
  select(ID, o_antigen, h1, serogroup, serovar)

# SeqSero2 
colnames(ss)

ss2 <- ss %>%
  mutate(ID=Sample.name) %>%
  select(ID, Predicted.antigenic.profile, Predicted.serotype)

# MLST
colnames(st)

st2 <- st %>%
  mutate(ID=str_split_fixed(V1, "\\_", 2)[,1]) %>%
  mutate(ST=V3) %>%
  select(ID, ST)

# Combine

all <- checkm2 %>%
  inner_join(prk2, by = "ID") %>%
  inner_join(quast2, by = "ID") %>%
  inner_join(sr2, by = "ID") %>%
  inner_join(ss2, by = "ID") %>%
  inner_join(st2, by = "ID")

# ----- Summary Stats ----

summarydf.all <- as.data.frame(t(summary(all[,sapply(all, is.numeric)])))

summarydf.all <- summarydf.all %>%
  mutate(measure = str_split_fixed(Freq, ":", 2)[,1]) %>%
  mutate(value = round(as.numeric(str_split_fixed(Freq, ":", 2)[,2]), 2)) %>%
  select(-c(Var2, Freq))

summarydf.all <- summarydf.all %>%
  pivot_wider(names_from = measure)

colnames(summarydf.all)[1] <- "Parameter"
summarydf.all <- summarydf.all[,1:7]

# repeat without ARS16
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

# saved as "allcombinedentero_summarystats_noars16.txt"