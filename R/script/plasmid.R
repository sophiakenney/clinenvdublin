# Plasmid Profiles 

# ---- Load Packages and Tables ----
library(tidyverse)

pld <- read.delim("R/abricate/pld_long.txt", sep = "\t")

# ----- Plot -----

sourcepal = c("Clinical" = "#1E5A46FF",
              "Environmental" = "#75884BFF")

ggplot(pld %>%
         filter(ID != "ARS16") %>%
         mutate(Source = case_when(
           str_detect(ID, "ARS") ~ "Environmental",
           str_detect(ID, "NDSU") ~ "Clinical")) %>%
         filter(value != 0) %>%
         group_by(Source, variable) %>%
         mutate(ppn = case_when(
           Source == "Clinical" ~ sum(value)/26,
           Source == "Environmental" ~ sum(value)/17
         )) %>%
         ungroup() %>%
         select(Source, variable, ppn) %>%
         unique() %>%
         mutate(var2 = case_when(
           variable == "Col440I_1" ~ "Col440I",
           variable == "IncA.C2_1" ~ "IncA/C2",
           variable == "IncFIA.HI1._1_HI1" ~ "IncFIA/HI1",
           variable == "IncFII.S._1" ~ "IncFII/S",
           variable == "IncHI1A_1" ~ "IncHI1A",
           variable == "IncHI1B.R27._1_R27" ~ "IncHI1B/R27",
           variable == "IncN_1" ~ "IncN",
           variable == "IncX1_1" ~ "IncX1",
         )), aes(x=fct_infreq(var2), y= round(ppn*100, digits = 2), fill=Source)) +
  geom_bar(position="dodge", stat="identity") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  theme_classic() +
  xlab("Plasmid") +
  ylab("Proportion of Strains") +
  scale_fill_manual(values = sourcepal)+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 20, colour = "black"))


# ----- Statistics: Is MDR IncA/C2 Plasmid more prevalent in clinical strains? ----

prop.test(c(pld %>%
              filter(variable == "IncA.C2_1") %>%
              filter(value == 1) %>%
              mutate(source = ifelse(str_detect(ID, "ARS"), "Environmental", "Clinical"
              )) %>%
              filter(source == "Clinical") %>%
              nrow(),
            
            pld %>%
              filter(variable == "IncA.C2_1") %>%
              filter(value == 1) %>%
              mutate(source = ifelse(str_detect(ID, "ARS"), "Environmental", "Clinical"
              )) %>%
              filter(source == "Environmental") %>%
              nrow()), c(26,17), alternative = "greater") # p-vale 0.01112



