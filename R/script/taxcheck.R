# Check Taxonomic Classification

# ---- Load Packages and Tables ----

library(tidyverse)
library(tidylog)
library(reshape2)

# load tables 
tax <- read.delim("readqc/taxassignmenttab.tsv" ,sep = "\t")
classified <- read.csv("readqc/classificationsummary.csv")

tab <- reshape2::melt(tax %>%
                        select(-c(taxID,Max)))

tab2 <- merge(tab %>%
                mutate(Name = str_split_fixed(variable, "\\.", 2)[,1]) %>%
                select(-c(variable, lineage)),
              classified %>%
                select(Name, Classified.reads) %>%
                mutate(classified = as.numeric(str_remove_all(Classified.reads, ","))) %>%
                select(-c(Classified.reads)), by = "Name")

# subset by strain
bystrain <- list()

for (i in unique(tab2$Name)) {
  df <- tab2[tab2$Name == i, ]
  bystrain[[i]] <- df
}

bystrain[[i]] # check

# ---- Top Family ----

top_ppns_list <- list() # make list

# loop to get the top genera in each by ppn and combine into one table
for (i in unique(tab2$Name)){
  result <- bystrain[[i]] %>%
    filter(taxRank == "F") %>% # Family
    mutate(ppn = value / classified) %>%
    filter(ppn > 0.6) %>%
    unique()
  
  # Add the result to the list
  top_ppns_list[[i]] <- result
}

# Combine all results into a single data frame
top.ppns <- do.call(rbind, top_ppns_list)

summary(top.ppns$ppn)

no16 <- top.ppns %>%
  filter(Name != "ARS16")

summary(no16$ppn)
