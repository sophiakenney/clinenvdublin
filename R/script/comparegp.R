# Compare Genotype Phenotype 

# ---- Load Packages ----
library(tidyverse)
library(tidylog)
library(ggpubr)
library(irr)
library(gt)

# ---- Load Tables and Format -----
arg <- read.delim("R/amrfinder/amrfinder_entero90.txt", sep = "\t")
arp <- read.delim("R/ast/ast_longfmt.txt", sep = "\t")


# ---- Merge ARG and AST Tables ----- 
# rename abx classes so they match
table(arp$ABX_CLASS)
table(arg$Class)

arg <- arg %>%
  mutate(Class= str_to_sentence(Class)) %>%
  mutate(Subclass = str_to_sentence(Subclass)) %>%
  mutate(Class = if_else(Subclass == "Sulfonamide", "Folate Pathway Inhibitor", Class)) 

arp <- arp %>%
  filter(ABX %in% c("Chloramphenicol", "Doxycycline", "Tetracycline", "Ampicillin", "Ceftiofur", "Ceftazidime", "Azithromycin", "Trim/Sulfa")) %>%
  mutate(Class = case_when(
    str_detect(ABX_CLASS, "Fluoroquinolone") ~ "Quinolone",
    str_detect(ABX_CLASS, "ß-lactam") ~ "Beta-lactam", 
    str_detect(ABX_CLASS, "Trimethoprim/Sulfonamide") ~ "Folate Pathway Inhibitor", 
    .default = ABX_CLASS))

# make sure all geneXdrug combos exist 
allgenes<- unique(arg %>%
                    filter(Class %in% arp$Class) %>% select(Gene.symbol))
allstrains <- unique(arp$variable)

allcombos <- expand.grid(ID = allstrains, Gene.symbol = allgenes$Gene.symbol)

combined <- arg %>%
  filter(Class %in% arp$Class) %>%
  mutate(SIR_gen = "R") %>%
  full_join(allcombos, by = c("ID", "Gene.symbol"), relationship = "one-to-many")

combined$SIR_gen[is.na(combined$SIR_gen)] <- "S" # fix NAs

# add the corresponding gene info based on other columns
combined <- combined %>%
  group_by(Gene.symbol) %>%
  mutate(Element.type = ifelse(is.na(Element.type), first(na.omit(Element.type)), Element.type)) %>%
  mutate(Element.subtype = ifelse(is.na(Element.subtype), first(na.omit(Element.subtype)), Element.subtype)) %>%
  mutate(Class = ifelse(is.na(Class), first(na.omit(Class)), Class)) %>%
  mutate(Subclass = ifelse(is.na(Subclass), first(na.omit(Subclass)), Subclass)) %>%
  ungroup()


combined <- arp %>%
  filter(Class %in% arg$Class) %>%
  mutate(ID=variable) %>%
  mutate(SIR_phen = SIR) %>%
  select(ID, Class, ABX, SIR_phen) %>%
  unique() %>%
  full_join(combined, by = c("ID", "Class"), relationship = "many-to-many")


# where intermediate is considered resistant
combined2 <- combined %>%
  filter(ID != "ARS16") %>%
  filter(SIR_phen != "No interpretation") %>%
  mutate(SIR_phen = if_else(SIR_phen == "I", "R", SIR_phen)) %>%
  select(ID, Class, ABX, Gene.symbol, SIR_gen, SIR_phen, Element.type, Element.subtype, Scope)

# where intermediate is retained but still filtering out no interpretation and ARS16
combined3 <- combined %>%
  filter(ID != "ARS16") %>%
  #filter(SIR_phen != "No interpretation") %>%
  select(ID, Class, ABX, Gene.symbol, SIR_gen, SIR_phen, Element.type, Element.subtype, Scope)

# Group genes with redundant functions
combined4 <- combined3 %>%
  filter(ABX != "Enrofloxacin") %>%
  mutate(Gene.symbol2 = case_when(
    str_detect(Gene.symbol, "CMY") ~ "ampC_blaCMY",
    str_detect(Gene.symbol, "TEM") ~ "ESBL_blaTEM",
    str_detect(Gene.symbol, "tet") ~ "tetA_B",
    str_detect(Gene.symbol, "sul") ~ "sul1_2",
    Gene.symbol %in% c("floR", "catA1") ~ "floR_catA1"
  )) %>%
  group_by(ID, Gene.symbol2) %>%
  mutate(
    SIR_gen2 = case_when(
      any(SIR_gen == "R") ~ "R",    # If any row in the group has "R", set "R"
      all(SIR_gen == "S") ~ "S",    # If all are "S", set "S"
      TRUE ~ SIR_gen                # Otherwise, retain the original value
    )) %>%
  select(-c(Gene.symbol, SIR_gen)) %>%
  unique()

combined4$ABX <- factor(combined4$ABX, levels = c("Ampicillin", "Imipenem", "Ceftazidime", "Ceftiofur", "Chloramphenicol", 
                                                  "Doxycycline", "Tetracycline", "Trim/Sulfa"))

# ---- Plot Comparison -----

classpal = c("Beta-lactam" = "#AF4F2FFF",
             "Phenicol" = "#D8B847FF",
             "Folate Pathway Inhibitor" = "#59385CFF",
             "Tetracycline" = "#732F30FF", 
             "No ARG" = "white")

# Figure 4
ggplot(combined4 %>%
         mutate(Class = ifelse(SIR_gen2 == "S", "No ARG", Class)) %>%
         mutate(Gene.symbol2 = case_when(
           Gene.symbol2 == "ampC_blaCMY" ~ "AmpC (blaCMY)", 
           Gene.symbol2 == "ESBL_blaTEM" ~ "ESBL (blaTEM)", 
           Gene.symbol2 == "floR_catA1" ~ "floR/catA1", 
           Gene.symbol2 == "tetA_B" ~ "tetA/B", 
           Gene.symbol2 == "sul1_2" ~ "sul1/2")) %>%
         #mutate(Gene.symbol = ifelse(Gene.symbol == "Not present", "No ARG", Gene.symbol)) %>%
         mutate(ABX = if_else(str_detect(ABX, "Trim/Sulfa"), "Trimethoprim Sulfa", ABX),
                ABX=factor(ABX, levels = c("Ampicillin", "Imipenem", "Ceftazidime", "Ceftiofur", "Chloramphenicol",
                                           "Doxycycline", "Tetracycline", "Trimethoprim Sulfa", "Plasmid"))), aes(y=fct_infreq(ID), x=Gene.symbol2, fill = Class)) + 
  geom_tile(color="black")+
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_manual(values = classpal)+
  #coord_fixed()+
  theme_classic() +
  facet_grid(~ABX, scales = "free_x", space = "free_x")+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black", fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(angle = 90, hjust = 0,  size = 18, color = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black")) +
  geom_text(data = subset(combined4 %>%
                            mutate(Gene.symbol2 = case_when(
                              Gene.symbol2 == "ampC_blaCMY" ~ "AmpC (blaCMY)", 
                              Gene.symbol2 == "ESBL_blaTEM" ~ "ESBL (blaTEM)", 
                              Gene.symbol2 == "floR_catA1" ~ "floR/catA1", 
                              Gene.symbol2 == "tetA_B" ~ "tetA/B", 
                              Gene.symbol2 == "sul1_2" ~ "sul1/2")) %>%
                            mutate(ABX = if_else(str_detect(ABX, "Trim/Sulfa"), "Trimethoprim Sulfa", ABX),
                                   ABX=factor(ABX, levels = c("Ampicillin", "Imipenem", "Ceftazidime", "Ceftiofur", "Chloramphenicol", 
                                                              "Doxycycline", "Tetracycline", "Trimethoprim Sulfa", "Plasmid"))),
                          SIR_phen == "R"), aes(label = "R"), color = "black", fontface = "bold", size = 6) +
  geom_text(data = subset(combined4 %>%
                            mutate(Gene.symbol2 = case_when(
                              Gene.symbol2 == "ampC_blaCMY" ~ "AmpC (blaCMY)", 
                              Gene.symbol2 == "ESBL_blaTEM" ~ "ESBL (blaTEM)", 
                              Gene.symbol2 == "floR_catA1" ~ "floR/catA1", 
                              Gene.symbol2 == "tetA_B" ~ "tetA/B", 
                              Gene.symbol2 == "sul1_2" ~ "sul1/2")) %>%
                            mutate(ABX = if_else(str_detect(ABX, "Trim/Sulfa"), "Trimethoprim Sulfa", ABX),
                                   ABX=factor(ABX, levels = c("Ampicillin", "Imipenem", "Ceftazidime", "Ceftiofur", "Chloramphenicol", 
                                                              "Doxycycline", "Tetracycline", "Trimethoprim Sulfa", "Plasmid"))), 
                          SIR_phen == "I"), aes(label = "I"), color = "black", fontface = "bold", size = 6) +
  geom_text(data = subset(combined4 %>%
                            mutate(Gene.symbol2 = case_when(
                              Gene.symbol2 == "ampC_blaCMY" ~ "AmpC (blaCMY)", 
                              Gene.symbol2 == "ESBL_blaTEM" ~ "ESBL (blaTEM)", 
                              Gene.symbol2 == "floR_catA1" ~ "floR/catA1", 
                              Gene.symbol2 == "tetA_B" ~ "tetA/B", 
                              Gene.symbol2 == "sul1_2" ~ "sul1/2")) %>%
                            mutate(ABX = if_else(str_detect(ABX, "Trim/Sulfa"), "Trimethoprim Sulfa", ABX),
                                   ABX=factor(ABX, levels = c("Ampicillin", "Imipenem", "Ceftazidime", "Ceftiofur", "Chloramphenicol", 
                                                              "Doxycycline", "Tetracycline", "Trimethoprim Sulfa", "Plasmid"))), 
                          SIR_phen == "No interpretation"), aes(label = "NI"), color = "black", fontface = "bold", size = 6)



# ---- Concordance - Statistical Analysis

cohens <- data.frame()
# now add for by ABX
for (i in unique(combined4$ABX)) {
  cohens <- rbind(cohens, data.frame(Class = paste(i), 
                                     kappa = kappa2(combined4 %>%
                                                      ungroup() %>%
                                                      filter(SIR_phen != "No interpretation") %>% 
                                                      mutate(SIR_phen = ifelse(SIR_phen == "I", "R", SIR_phen)) %>%
                                                      filter(ABX == i) %>%
                                                      select(SIR_gen2, SIR_phen), weight = "unweighted")$value, 
                                     p.value = kappa2(combined4 %>%
                                                        ungroup () %>%
                                                        filter(SIR_phen != "No interpretation") %>% 
                                                        mutate(SIR_phen = ifelse(SIR_phen == "I", "R", SIR_phen)) %>%
                                                        filter(ABX == i) %>%
                                                        select(SIR_gen2, SIR_phen), weight = "unweighted")$p.value,
                                     level = "ABX"))
  
}

senspectab <- combined4 %>%
  #group_by(Gene.symbol) %>%
  mutate(hitR_R = ifelse(
    SIR_phen == "R" & SIR_gen2 == "R", T, F)) %>%
  mutate(hitS_S = ifelse(
    SIR_phen == "S" & SIR_gen2 == "S", T, F)) %>%
  mutate(missR_S = ifelse(
    SIR_phen == "R" & SIR_gen2 == "S", T, F)) %>%
  mutate(missS_R = ifelse(
    SIR_phen == "S" & SIR_gen2 == "R", T, F)) %>%
  group_by(ABX) %>%
  mutate(GR_PR = sum(hitR_R == T)) %>%
  mutate(GR_PS = sum(missS_R  == T)) %>%
  mutate(GS_PR = sum(missR_S == T)) %>%
  mutate(GS_PS = sum(hitS_S == T)) %>%
  select(ABX, GR_PR, GS_PR, GR_PS, GS_PS) %>%
  unique() %>%
  mutate(sensitivity = round(GR_PR / (GR_PR + GR_PS), digits = 2)) %>%
  mutate(specificity = round(GS_PS / (GS_PR + GS_PS), digits = 2)) %>%
  mutate(PPV = round(GR_PR  / (GR_PR  + GS_PR), digits = 2)) %>%
  mutate(NPV = round(GR_PS / (GS_PS + GR_PS), digits = 2)) %>%
  mutate(accuracy = round((GS_PS+GR_PR)  / (GS_PS+GR_PR+GS_PR+GR_PS), digits = 2))

ss_cohens <- senspectab %>%
  inner_join(cohens %>%
               mutate(ABX=Class) %>%
               filter(level == "ABX") %>%
               select(ABX, kappa), by = "ABX")


# Table 3
gt(ss_cohens %>% ungroup())

