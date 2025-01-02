# Genotype Analysis

# ----- Load Packages and Data ---- 
library(tidyverse)
library(tidylog)
library(ggpubr)
library(reshape2)
library(gt)

arg <- read.delim("amrfinder/amrfinder_entero90.txt", sep = "\t")

# Format 

# wide format 
argmat <- arg %>%
  mutate(value = 1) %>%
  pivot_wider(id_cols = "ID", names_from= "Gene.symbol", values_from = "value") %>%
  filter(!str_detect(ID, "ARS16")) # filter ARS16 from dataset

# correct rownames and table issues from reformatting
argmat <- as.data.frame(lapply(argmat, function(x) str_replace_all(x, "NULL", "0")))
rownames(argmat) <- argmat$ID
argmat2 <- as.data.frame(lapply(argmat, function(x) as.numeric(x)))
argmat2$ID <- argmat$ID


# ---- Visualize Strain Genotypes ----

# Perform hierarchical clustering to order axis by strain similarity

# create distance matrix (presence/absence)
dmat <- dist(argmat[,-1], method = "binary")

# perform hierarchical clustering
hc <- hclust(dmat)

# reorder the rows based on clustering
argmat2$ID <- factor(argmat2$ID , levels = argmat2$ID [hc$order])


# Reformat to long for plotting
arg2 <- melt(argmat2, id.vars = "ID")

arg2 <- merge(arg2 %>% mutate(Gene.symbol = variable) %>% select(-c(variable)), # rename column
              # merge with gene metadata
              arg %>% filter(!str_detect(ID, "ARS16")) %>% # except for this strain
                select(ID, Gene.symbol, Scope, Element.type, 
                       Element.subtype, Class, Subclass, X..Identity.to.reference.sequence) %>%
                unique() %>%
                group_by(Element.subtype) %>%
                arrange(Gene.symbol, .by_group = TRUE) %>%
                mutate(Gene.symbol = factor(Gene.symbol, levels = unique(Gene.symbol))) %>%
                mutate(Element.subtype = str_to_sentence(Element.subtype)) %>% # format naming
                mutate(Element.subtype = ifelse(Element.subtype == "Amr", "AMR", Element.subtype)), by=c("ID", "Gene.symbol")) # format naming


# Plot Dissertation Figure 5.3
ggplot(arg %>%
         mutate(ID = factor(arg$ID , levels = argmat2$ID [hc$order])) %>%
         filter(ID != "ARS16") %>%
         mutate(Element.subtype = str_to_sentence(Element.subtype)) %>%
         mutate(Element.subtype = ifelse(Element.subtype == "Amr", "AMR", Element.subtype))%>%
         filter(Element.subtype %in% c("AMR","Point")) %>%
         mutate(Class = str_to_sentence(Class)) %>%
         mutate(Class = if_else(str_detect(Class, "uino"), "Quinolone", Class)),
       aes(x=ID, y= Gene.symbol, shape = Element.subtype, color = Class)) +
  geom_point(size = 6)+
  coord_fixed() + 
  scale_color_viridis_d(option = "A", end = 0.8) + 
  theme_classic() +
  guides(shape = guide_legend(ncol = 1, title = "ARG Type"),
         color = guide_legend(ncol = 1, title = "Drug Class"))+
  theme(axis.text.y = element_text(face = "italic", color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        text=element_text(size = 18),
        axis.title = element_blank(),
        legend.position = "right") 


# ---- Plot Dissertation Table 5.2 ----

efflux <- rbind(arg %>%
                  filter(!str_detect(ID, "ARS16")) %>% 
                  filter(Element.type == "AMR") %>%
                  select(ID, Gene.symbol, Element.type, Class) %>%
                  unique(),
                #add dataframe for efflux genes identified by prokka
                expand.grid(ID = unique(arg$ID), Gene.symbol = c("emrA", "emrB", "mdtA", "mdtB", "acrA", "acrB", "acrD", "acrE", "acrF")) %>%
                  full_join(data.frame(Class = c(rep("EFFLUX", 9)),
                                       Element.type = c(rep("AMR ", )),
                                       Gene.symbol = c("emrA", "emrB", "mdtA", "mdtB", "acrA", "acrB", "acrD", "acrE", "acrF")
                  ), by = "Gene.symbol"))


gt(efflux %>%
     group_by(Gene.symbol, Class) %>%
     filter(ID != "ARS16") %>%
     mutate(N_strains = n()) %>%
     mutate(ppn = round((N_strains/43)*100, digits = 1)) %>%
     mutate(ppn = paste0(ppn, "%")) %>%
     ungroup() %>%
     mutate(Class = str_to_sentence(Class)) %>%
     mutate(Class = ifelse(str_detect(Class, "Quin"), "Quinolone", Class)) %>%
     select(Class, Gene.symbol, N_strains, ppn) %>%
     group_by(Class) %>%
     arrange(.by_group = TRUE) %>%
     unique()) %>%
  cols_label(
    N_strains = "Number of isolates") %>%
  cols_label(
    ppn = "Positive rate (%)") %>%
  cols_label(
    Gene.symbol = "Resistance Gene") %>%
  tab_style(
    style = cell_text(align = "left", weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style = cell_text(align = "center", weight = "bold"),
    locations = cells_column_labels(columns = c(Gene.symbol, Class))
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_body(columns = c(N_strains, ppn))
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = c(Gene.symbol))
  ) %>%
  tab_style(
    style = cell_text(font = "Times New Roman"),
    locations = cells_body()
  ) %>%
  tab_style(
    style = cell_text(font = "Times New Roman"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(font = "Times New Roman"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style = cell_borders(sides = c("left", "right"), color = "lightgrey", weight = px(1)),
    locations = cells_body(columns = everything())
  ) %>%
  tab_style(
    style = cell_borders(sides = c("left", "right"), color = "lightgrey", weight = px(1)),
    locations = cells_column_labels(columns = everything())
  ) %>%
  tab_style(
    style = cell_borders(sides= c( "top", "bottom"), color = "lightgrey", weight = px(1)),
    locations = cells_body(columns = everything())
  ) %>%
  opt_vertical_padding(scale = 0) %>%
  tab_footnote(
    footnote = "Genes detected by AMRFinderPlus and/or Prokka",
    locations = cells_column_labels(
      columns =Gene.symbol)
  )



