# Phenotype Analysis

# ---- Load Packages and Tables -----
library(tidyverse)
library(purrr)
library(tidylog)

# Load Table 
ast <- readxl::read_xlsx("R/ast/ast_results.xlsx", sheet = "cutoffs")

# ---- Correct to NARMS Interps ----

ast2 <- rbind(
  # NARMS interps
  ast %>%
  select(-c(`ABX_Concentrations(mg/L)`)) %>%
  filter(!ABX %in% c("Amikacin", "Gentamicin", "Rifampin", "Cefazolin",
                     "Erythromycin", "Oxacillin_2pNaCL", "Penicillin", "Ticarcillin", "Timentin (Ticarcillin/clauvanic acid)"))%>%
  mutate(across(everything(), ~ if_else(ABX == "Ampicillin", str_replace_all(., ">32", "R"), .)))  %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ampicillin", str_replace_all(., "2", "S"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ampicillin", str_replace_all(., "1", "S"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ampicillin", str_replace_all(., "0.5", "S"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Azithromycin", str_replace_all(., "2", "S"), .))) %>% #NARMS -
  mutate(across(everything(), ~ if_else(ABX == "Azithromycin", str_replace_all(., "4", "R"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ceftazidime", str_replace_all(., "32", "R"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ceftazidime", str_replace_all(., "64", "R"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ceftazidime", str_replace_all(., "16", "R"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ceftazidime", str_replace_all(., "8", "I"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ceftazidime", str_replace_all(., "4", "S"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ceftazidime", str_replace_all(., "2", "S"), .)))%>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ceftazidime", str_replace_all(., "<=1", "S"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ceftazidime", str_replace_all(., "<=1", "S"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ceftiofur", str_replace_all(., ">4", "R"), .))) %>% # NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ceftiofur", str_replace_all(., "4", "I"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ceftiofur", str_replace_all(., "0.5", "S"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ceftiofur", str_replace_all(., "1", "S"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Ceftiofur", str_replace_all(., "<=0.25", "S"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Chloramphenicol", str_replace_all(., ">32", "R"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Chloramphenicol", str_replace_all(., "<=4", "S"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Chloramphenicol", str_replace_all(., "32", "R"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Chloramphenicol", str_replace_all(., "16", "I"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Clarithromycin", str_replace_all(., ">8", "No interpretation"), .))) %>% # No Interpretation
  mutate(across(everything(), ~ if_else(ABX == "Doxycycline", str_replace_all(., ">16", "R"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Doxycycline", str_replace_all(., "16", "R"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Doxycycline", str_replace_all(., "8", "I"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Doxycycline", str_replace_all(., "<=2", "S"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Enrofloxacin", str_replace_all(., "<=0.25", "No interpretation"), .))) %>% # use NARMS once found
  mutate(across(everything(), ~ if_else(ABX == "Enrofloxacin", str_replace_all(., "0.5", "No interpretation"), .))) %>% # use NARMS once found
  mutate(across(everything(), ~ if_else(ABX == "Imipenem", str_replace_all(., "<=1", "S"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Tetracycline", str_replace_all(., ">8", "R"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Tetracycline", str_replace_all(., "<=2", "S"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Trim/Sulfa", str_replace_all(., ">4", "R"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Trim/Sulfa", str_replace_all(., "4", "R"), .))) %>% #NARMS
  mutate(across(everything(), ~ if_else(ABX == "Trim/Sulfa", str_replace_all(., "<=0.5", "S"), .))), #NARMS
  
  # non-NARMS interps
  ast %>%
    select(-c(`ABX_Concentrations(mg/L)`)) %>%
    filter(ABX %in% c("Amikacin", "Gentamicin", "Rifampin", "Cefazolin", "Clarithromycin",
                      "Erythromycin", "Oxacillin_2pNaCL", "Penicillin", "Ticarcillin", "Timentin (Ticarcillin/clauvanic acid)")) %>%
    mutate(across(everything(), ~ if_else(ABX == "Amikacin", str_replace_all(., "<=4", "S"), .))) %>% 
    mutate(across(everything(), ~ if_else(ABX == "Cefazolin", str_replace_all(., ">16", "R"), .))) %>% # default to CLSI M100 interpretation - no NARMS
    mutate(across(everything(), ~ if_else(ABX == "Cefazolin", str_replace_all(., "<=4", "S"), .))) %>% # default to CLSI M100 interpretation - no NARMS
    mutate(across(everything(), ~ if_else(ABX == "Cefazolin", str_replace_all(., "8", "R"), .))) %>%  # ddefault to CLSI M100 interpretation - no NARMS
    mutate(across(everything(), ~ if_else(ABX == "Clarithromycin", str_replace_all(., ">8", "R"), .))) %>% # default to adl interpretation - no NARMS/CLSI
    mutate(across(everything(), ~ if_else(ABX == "Erythromycin", str_replace_all(., ">8", "R"), .))) %>% # default to adl interpretation - no NARMS/CLSI
    mutate(across(everything(), ~ if_else(ABX == "Gentamicin", str_replace_all(., "<=1", "S"), .))) %>% # default to CLSI M100 Enterobacterales since full NARMS range wasnt tested
    mutate(across(everything(), ~ if_else(ABX == "Gentamicin", str_replace_all(., ">8", "R"), .))) %>% # default to CLSI M100 /Enterobacterales since full NARMS range wasnt tested
    mutate(across(everything(), ~ if_else(ABX == "Oxacillin_2pNaCL", str_replace_all(., ">4", "R"), .))) %>% # default to adl interpretation - no NARMS/CLSI
    mutate(across(everything(), ~ if_else(ABX == "Penicillin", str_replace_all(., ">8", "No interpretation"), .))) %>% # default to adl interpretation - no NARMS/CLSI
    mutate(across(everything(), ~ if_else(ABX == "Penicillin", str_replace_all(., "4", "No interpretation"), .))) %>% # default to adl interpretation - no NARMS/CLSI
    mutate(across(everything(), ~ if_else(ABX == "Penicillin", str_replace_all(., "8", "No interpretation"), .))) %>% # default to adl interpretation - no NARMS/CLSI
    mutate(across(everything(), ~ if_else(ABX == "Penicillin", str_replace_all(., "2", "No interpretation"), .))) %>% # default to adl interpretation - no NARMS/CLSI
    mutate(across(everything(), ~ if_else(ABX == "Rifampin", str_replace_all(., ">4", "R"), .))) %>% # default to adl interpretation - no NARMS/CLSI
    mutate(across(everything(), ~ if_else(ABX == "Rifampin", str_replace_all(., "4", "R"), .))) %>% # default to adl interpretation - no NARMS/CLSI
    mutate(across(everything(), ~ if_else(ABX == "Ticarcillin", str_replace_all(., ">64", "R"), .))) %>% # default to adl interpretation - no NARMS/CLSI
    mutate(across(everything(), ~ if_else(ABX == "Ticarcillin", str_replace_all(., "64", "I"), .))) %>% # default to adl interpretation - no NARMS/CLSI
    mutate(across(everything(), ~ if_else(ABX == "Ticarcillin", str_replace_all(., "32", "I"), .))) %>% # default to adl interpretation - no NARMS/CLSI
    mutate(across(everything(), ~ if_else(ABX == "Ticarcillin", str_replace_all(., "16", "S"), .))) %>% # default to adl interpretation - n - no NARMS/CLSIo NARMS/CLSI
    mutate(across(everything(), ~ if_else(ABX == "Ticarcillin", str_replace_all(., "<=8", "S"), .))) %>% # default to adl interpretation
    mutate(across(everything(), ~ if_else(ABX == "Timentin (Ticarcillin/clauvanic acid)", str_replace_all(., ">64", "R"), .))) %>% # default to CLSI M100 interpretation - no NARMS
    mutate(across(everything(), ~ if_else(ABX == "Timentin (Ticarcillin/clauvanic acid)", str_replace_all(., "64", "I"), .))) %>% # default to CLSI M100 interpretation - no NARMS
    mutate(across(everything(), ~ if_else(ABX == "Timentin (Ticarcillin/clauvanic acid)", str_replace_all(., "32", "I"), .))) %>% # default to CLSI M100 interpretation - no NARMS
    mutate(across(everything(), ~ if_else(ABX == "Timentin (Ticarcillin/clauvanic acid)", str_replace_all(., "16", "S"), .))) %>% # default to CLSI M100 interpretation - no NARMS
    mutate(across(everything(), ~ if_else(ABX == "Timentin (Ticarcillin/clauvanic acid)", str_replace_all(., "<=8", "S"), .)))) # default to CLSI M100 interpretation - no NARMS
  
  
# ---- Reformat to Long ----

# change interps to numeric
ast3 <- cbind(
  ast2 %>% 
  select(-c(ABX, ABX_CLASS)) %>%
  mutate(across(everything(), ~ case_when(
    . == "R" ~ 1,
    . == "S" ~ -1,
    . == "I" ~ 0,
    . == "No interpretation" ~ 5,
    TRUE ~ as.numeric(.) # Retain original value if it doesn't match
  ))),
  ast2 %>%
    select(ABX, ABX_CLASS))

  
# reformat 
astm <- reshape2::melt(ast3)

# add interp column
astm <- astm %>%
  mutate(SIR = case_when(
    value == 0 ~ "I",
    value == 1 ~ "R", 
    value == -1 ~ "S",
    value == 5 ~ "No interpretation"
  ))

#write.table(astm, "R/ast/ast_longfmt.txt", sep = "\t")

# ---- Plot AST Interps ----

# Plot all for supplemental figure 
ggplot(astm %>%
         filter(variable != "ARS16") %>% # exclude ARS16 since excluded from the study 
         filter(SIR != "No interpretation") %>%
         mutate(ABX = if_else(str_detect(ABX, "Oxa"), "Oxacillin (2% NaCl)", ABX)) %>%
         mutate(ABX = if_else(str_detect(ABX, "Timentin"), "Timentin", ABX))%>%
         mutate(ABX = if_else(str_detect(ABX, "Trim/Sulfa"), "Trimethoprim Sulfa", ABX)), aes(y = fct_infreq(ABX), x = fct_infreq(variable), fill = SIR)) +
  geom_tile(color = "white") +
  coord_fixed() + 
  #geom_vline(xintercept = 26.5, color = "red") +
  theme_classic() + 
  scale_y_discrete(limits = rev) + 
  scale_fill_manual(values = c(
    "I" = "grey", 
    "R" = "black", 
    "S" = "white"
  )) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        axis.ticks = element_blank()) 


SIRpal <- c("S" = "#EDEDEF",
            "R" = "#232C43FF",
            "I" = "#4C6C94FF",
            "No interpretation" = "#A4ABB0FF")

# Plot Figure 1

ggplot(astm %>%
         #filter(SIR != "No interpretation") %>%
         filter(SIR != "S")%>%
         filter(!str_detect(variable, "ARS16")) %>% 
         filter(ABX %in% c("Chloramphenicol", "Doxycycline", "Tetracycline", "Ampicillin", "Ceftiofur", "Ceftazidime", "Azithromycin", "Trim/Sulfa")) %>%
         mutate(ABX = if_else(str_detect(ABX, "Trim/Sulfa"), "Trimethoprim/Sulfa", ABX)) %>%
         mutate(variable = factor(variable, levels = rev(c(paste0("ARS", 1:18), paste0("NDSU", 1:30))))),
       aes( y = fct_infreq(ABX), x = fct_infreq(variable), fill = SIR)) +
  geom_tile(color = "white") +
  coord_fixed() + 
  #geom_vline(xintercept = 26.5, color = "red") +
  theme_classic() + 
  scale_y_discrete(limits = rev) + 
  scale_fill_manual(values = SIRpal[1:3], limits = c("R", "I", "S")) +
  guides(fill = guide_legend(nrow = 1, title = "Interpretation"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14, color= "black"),
        axis.text.y = element_text(size = 14, color= "black"),
        axis.title = element_blank(),
        legend.position = "bottom", 
        text = element_text(size = 14))


# ---- Plot Comparisons AMR/MDR between Clinical and Environmental ----

sourcepal = c("Clinical" = "#1E5A46FF",
              "Environmental" = "#75884BFF")


# Figure 2A
 ggplot(astm %>%
  filter(ABX %in% c("Chloramphenicol", "Doxycycline", "Tetracycline", "Ampicillin", "Ceftiofur", "Ceftazidime", "Azithromycin", "Trim/Sulfa")) %>%
  group_by(ABX_CLASS, variable) %>%
  mutate(classtotal = sum(SIR == "R") + sum(SIR == "I")) %>%
  ungroup() %>%
  mutate(classPA = ifelse(
    classtotal > 0, 1, 0
  )) %>%
  select(variable, ABX_CLASS, classtotal, classPA) %>%
  unique() %>%
  group_by(variable) %>%
  mutate(straintotal = sum(classPA == 1)) %>%
  ungroup() %>%
  select(variable, straintotal) %>%
  unique() %>%
  filter(variable != "ARS16") %>%
  mutate(Clinical = ifelse(str_detect(variable, "NDSU"), "Clinical", "Environmental")) %>%
  filter(Clinical == "Clinical") %>%
  group_by(straintotal) %>%
  count() %>%
  ungroup() %>%
  mutate(Clinical = n) %>%
  mutate(Nclasses = straintotal) %>%
  select(Nclasses, Clinical) %>%
  full_join(
    astm %>%
      filter(ABX %in% c("Chloramphenicol", "Doxycycline", "Tetracycline", "Ampicillin", "Ceftiofur", "Ceftazidime", "Azithromycin", "Trim/Sulfa")) %>%
      group_by(ABX_CLASS, variable) %>%
      mutate(classtotal = sum(SIR == "R") + sum(SIR == "I")) %>%
      ungroup() %>%
      mutate(classPA = ifelse(
        classtotal > 0, 1, 0
      )) %>%
      select(variable, ABX_CLASS, classtotal, classPA) %>%
      unique() %>%
      group_by(variable) %>%
      mutate(straintotal = sum(classPA == 1)) %>%
      ungroup() %>%
      select(variable, straintotal) %>%
      unique() %>%
      filter(variable != "ARS16") %>%
      mutate(Environmental = ifelse(str_detect(variable, "NDSU"), "Clinical", "Environmental")) %>%
      filter(Environmental == "Environmental") %>%
      group_by(straintotal) %>%
      count() %>%
      ungroup() %>%
      mutate(Environmental = n) %>%
      mutate(Nclasses = straintotal) %>%
      select(Nclasses, Environmental), by = "Nclasses") %>%
  arrange(Nclasses) %>%
  mutate(across(everything(), ~ replace_na(., 0))) %>%
  pivot_longer(!Nclasses, names_to="Source", values_to="Count") %>%
  mutate(ppn = ifelse(Source == "Clinical", Count/26, Count/17)) %>%
   mutate(ppn = round(ppn*100, digits = 1)) %>%
    mutate(Pattern = ifelse(Nclasses == 0, "Pansusceptible", paste0(Nclasses, " Classes"))) %>%
   mutate(Pattern = ifelse(Nclasses == 1, "1 Class", Pattern)), aes(x=Pattern, y=ppn, fill = Source)) +
   geom_bar(position="dodge", stat="identity") +
   theme_classic() +
   scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
   scale_x_discrete(limits=c("Pansusceptible", "1 Class", "2 Classes", "3 Classes", "4 Classes", "5 Classes")) +
   xlab("Resistance Pattern") +
   ylab("Proportion of Strains") +
   scale_fill_manual(values = sourcepal)+
   theme(axis.text = element_text(size = 20, color = "black"),
         axis.title = element_text(size = 20, color = "black"),
         legend.position = "bottom",
         legend.title = element_text(size = 20, colour = "black"),
         legend.text = element_text(size = 20, colour = "black"))


# Figure 2B
 
 mdr <- astm %>%
   filter(ABX %in% c("Chloramphenicol", "Doxycycline", "Tetracycline", "Ampicillin", "Ceftiofur", "Ceftazidime", "Azithromycin", "Trim/Sulfa")) %>%
   group_by(ABX_CLASS, variable) %>%
   mutate(classtotal = sum(SIR == "R") + sum(SIR == "I")) %>%
   ungroup() %>%
   mutate(classPA = ifelse(
     classtotal > 0, 1, 0
   )) %>%
   select(variable, ABX_CLASS, classtotal, classPA) %>%
   unique() %>%
   group_by(variable) %>%
   mutate(straintotal = sum(classPA == 1)) %>%
   ungroup() %>%
   select(variable, straintotal) %>%
   unique() %>%
   filter(variable != "ARS16") %>%
   mutate(Clinical = ifelse(str_detect(variable, "NDSU"), "Clinical", "Environmental")) %>%
   filter(Clinical == "Clinical") %>%
   group_by(straintotal) %>%
   count() %>%
   ungroup() %>%
   mutate(Clinical = n) %>%
   mutate(Nclasses = straintotal) %>%
   select(Nclasses, Clinical) %>%
   full_join(
     astm %>%
       filter(ABX %in% c("Chloramphenicol", "Doxycycline", "Tetracycline", "Ampicillin", "Ceftiofur", "Ceftazidime", "Azithromycin", "Trim/Sulfa")) %>%
       group_by(ABX_CLASS, variable) %>%
       mutate(classtotal = sum(SIR == "R") + sum(SIR == "I")) %>%
       ungroup() %>%
       mutate(classPA = ifelse(
         classtotal > 0, 1, 0
       )) %>%
       select(variable, ABX_CLASS, classtotal, classPA) %>%
       unique() %>%
       group_by(variable) %>%
       mutate(straintotal = sum(classPA == 1)) %>%
       ungroup() %>%
       select(variable, straintotal) %>%
       unique() %>%
       filter(variable != "ARS16") %>%
       mutate(Environmental = ifelse(str_detect(variable, "NDSU"), "Clinical", "Environmental")) %>%
       filter(Environmental == "Environmental") %>%
       group_by(straintotal) %>%
       count() %>%
       ungroup() %>%
       mutate(Environmental = n) %>%
       mutate(Nclasses = straintotal) %>%
       select(Nclasses, Environmental), by = "Nclasses") %>%
   arrange(Nclasses) %>%
   mutate(across(everything(), ~ replace_na(., 0))) %>%
   pivot_longer(!Nclasses, names_to="Source", values_to="Count") %>%
   mutate(ppn = ifelse(Source == "Clinical", Count/26, Count/17)) %>%
   mutate(ppn = round(ppn*100, digits = 1)) %>%
   mutate(Pattern = ifelse(Nclasses == 0, "Pansusceptible", paste0(Nclasses, " Classes"))) %>%
   mutate(Pattern = ifelse(Nclasses == 1, "1 Class", Pattern)) %>%
   mutate(MDR = case_when(
     Pattern %in% c("Pansusceptible", "1 Class", "2 Classes") ~ "Not MDR",
     Pattern %in% c("3 Classes", "4 Classes", "5 Classes") ~ "MDR"
   )) %>%
   group_by(MDR, Source) %>%
   mutate(PPN_MDR = sum(ppn)) %>%
   ungroup()
 
 ggplot(mdr, aes(x=MDR, y=PPN_MDR, fill=Source)) +
   geom_bar(stat="identity", position="dodge") +
   theme_classic() +
   scale_fill_manual(values = sourcepal)+
   scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
   scale_x_discrete(limits = c("Not MDR", "MDR"))+
   xlab("Resistance Pattern") +
   ylab("Proportion of Strains") +
   theme(axis.text = element_text(size = 20, color = "black"),
         axis.title = element_text(size = 20, color = "black"),
         legend.position = "bottom",
         legend.title = element_text(size = 20, colour = "black"),
         legend.text = element_text(size = 20, colour = "black"))

  

# Figure 2C
ggplot(astm %>%
  filter(!str_detect(variable, "ARS16")) %>% 
  filter(ABX %in% c("Chloramphenicol", "Doxycycline", "Tetracycline", "Ampicillin", "Ceftiofur", "Ceftazidime", "Azithromycin", "Trim/Sulfa")) %>%
  filter(SIR != "No interpretation") %>%
  mutate(Source = ifelse(str_detect(variable, "NDSU"), "Clinical", "Environmental")) %>%
  mutate(SIR = ifelse(SIR == "I", "R", SIR)) %>%
  filter(SIR == "R") %>%
  select(variable, Source, ABX, SIR) %>%
  unique() %>%
  group_by(ABX, Source) %>%
  mutate(N_strains = sum(SIR=="R")) %>%
  select(-c(variable)) %>%
  unique() %>%
  ungroup() %>%
  mutate(ppn = case_when(
    Source == "Clinical" ~ round((N_strains/26)*100, digits = 1),
    Source == "Environmental" ~ round((N_strains/17)*100, digits = 1)
  )) %>%
  bind_rows(data.frame(Source = "Environmental",
                       ABX = "Trim/Sulfa",
                       SIR = "R",
                       N_strains = 0,
                       ppn = 0)) %>%
  mutate(ABX = if_else(str_detect(ABX, "Trim/Sulfa"), "Trimethoprim/Sulfa", ABX)) ,
  aes(x=ABX, y=ppn, fill=Source)) +
  geom_bar(position="dodge", stat="identity") +
  theme_classic() +
  #scale_y_continuous(expand = c(0,0), limits = c(0,60)) +
  #scale_x_discrete(limits=c("Pansusceptible", "1 Class", "2 Classes", "3 Classes", "4 Classes", "5 Classes")) +
  xlab("Antimicrobial") +
  ylab("Proportion of Strains") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  scale_fill_manual(values = sourcepal)+
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 20, colour = "black"))
  

# ---- Statistical Comparisons ---- 
 
# Stats for MDR: Are a higher ppn of Clinical strains MDR than Environmental?
mdr %>%
  filter(Nclasses >=3) %>%
  filter(Source == "Clinical") %>%
  mutate(total = sum(Count)) %>%
  select(total) %>%
  unique() # 23

mdr %>%
  filter(Nclasses >=3) %>%
  filter(Source == "Environmental") %>%
  mutate(total = sum(Count)) %>%
  select(total) %>%
  unique() # 10


prop.test(c(23, 10), c(26, 17), alternative = "greater") # p val 0.03005


# Stats for each drug: Are a higher ppn of Clinical strains resistance to certain drugs than Environmental? 
ppnRbyclass <- merge(astm %>%
        filter(str_detect(variable, "NDSU")) %>%
        filter(SIR != "No interpretation") %>%
        filter(ABX %in% c("Chloramphenicol", "Doxycycline", "Tetracycline", "Ampicillin", "Ceftiofur", "Ceftazidime", "Azithromycin", "Trim/Sulfa")) %>%
        group_by(ABX) %>%
        mutate(NDSU_R = sum(SIR %in% c("I", "R"))/sum(SIR %in% c("S", "I", "R"))) %>%
        mutate(NDSU_total = 26) %>%
        select(ABX, ABX_CLASS, NDSU_R, NDSU_total) %>%
        unique(),
      
      astm %>%
        filter(str_detect(variable, "ARS")) %>% 
        filter(SIR != "No interpretation") %>%
        filter(variable != "ARS16") %>%
        filter(ABX %in% c("Chloramphenicol", "Doxycycline", "Tetracycline", "Ampicillin", "Ceftiofur", "Ceftazidime", "Azithromycin", "Trim/Sulfa")) %>%
        group_by(ABX) %>%
        mutate(ARS_R = sum(SIR %in% c("I", "R"))/sum(SIR %in% c("S", "I", "R"))) %>%
        mutate(ARS_total = 17) %>%
        select(ABX, ABX_CLASS, ARS_R, ARS_total) %>%
        unique(), by = c("ABX", "ABX_CLASS")) %>%
  mutate(NDSU_R = ifelse(NDSU_R == 1, 0.9999999, NDSU_R)) %>%
  mutate(ARS_R = ifelse(ARS_R == 1, 0.9999999, ARS_R)) %>%
  mutate(NDSU_R = ifelse(NDSU_R == 0, 0.0000001, NDSU_R)) %>%
  mutate(ARS_R = ifelse(ARS_R == 0, 0.0000001, ARS_R))

 
    
# Apply prop.test() across the dataframe
results <- apply(ppnRbyclass, 1, function(row) {
  # Convert proportions to counts of resistant cases
  resistant1 <- as.numeric(row["NDSU_R"]) * as.numeric(row["NDSU_total"])
  resistant2 <- as.numeric(row["ARS_R"]) * as.numeric(row["ARS_total"])
  
  # Perform the two-proportion z-test
  test_result <- prop.test(c(resistant1, resistant2), c(as.numeric(row["NDSU_total"]), as.numeric(row["ARS_total"])), 
                            alternative = "greater")
  
  # Return the p-value (or any other part of the test result)
  return(test_result$p.value)
})

# View the results (p-values)
ppnRbyclass$pval <- results
ppnRbyclass$BH_adjusted <- p.adjust(results, method = "BH")  

  
  