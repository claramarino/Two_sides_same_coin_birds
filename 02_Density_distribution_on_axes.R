# Visualisation of impact and mecha groups along PCA axis

rm(list=ls())

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(dummies)
library(Hmisc)
source("R/R_rainclouds.R")

#### Load data ####

# PCA from pcoa function
coordinates <- readRDS( "Output/01_PCA_coord_10D_all_birds_pcoa")
# bird groups and impact information
eicat_iast <- readRDS("Data/00_Alien_IAST_birds_impact_group")

# Impact magnitude EICAT
eicat <- read.csv2("Data/EICAT_birds_2016.csv") %>%
  mutate(binomial = paste(Gender, Species, sep=" ")) %>%
  dplyr::select(binomial, Impact_cate) %>% distinct()
# high impact ias-t
high_impact_ias <- readRDS("Output/01_IAST_high_impact_list") %>%
  pull(scientificName) %>% unique()

# taxonomic correspondence avo_traits x eicat/ias-t
avo_iast <- readRDS("Output/Data_clean/03_traits_AVONET_iast")
avo_eicat <- readRDS("Output/Data_clean/05_traits_AVONET_eicat")
taxo_corresp <- bind_rows(avo_eicat %>% mutate(group = "EICAT"), 
                     avo_iast %>%  mutate(bino_eicat = Species1, group = "IAS-T")) %>%
  select(Species1, bino_eicat) %>% distinct()


#### Set up species groups ####

all_groups <- left_join(eicat_iast, eicat, by="binomial") %>%
  # classify in 3 impact categories: Impact (MO, MR, MV), no-impact (MC, MN), DD
  mutate(impact = case_when(
    Impact_cate %in% c("MC", "MN") ~ "no_imp",
    Impact_cate=="DD" ~ "DD",
    Impact_cate %in% c("MO", "MV", "MR") ~ "imp")) %>%
  mutate(impact = if_else(group=="IAS-T",'NA', impact)) %>%
  mutate(group2 = paste(group, impact, sep="_")) %>%
  dplyr::select(-c(Impact_cate, impact)) %>%
  mutate(group2 = if_else(group2=="IAS-T_NA", 
                          if_else(binomial %in% high_impact_ias,
                                  "IAS-T_imp","IAS-T_no_imp"), group2))

str(all_groups)
summary(all_groups %>% mutate_if(is.character, as.factor))

#### Plot mechanisms along pca axis ####

# attribute taxonomic corresp to all_groups
syno_iucn_eicat <- readRDS("Output/Synonyms/05_synonyms_iucn_eicat")
syno_iucn_eicat <- syno_iucn_eicat %>%
  distinct(accepted_name, synonym) %>%
  rename(binomial=synonym)
# add synonyms founds manually outside of iucn
manual_syno = data.frame(
  binomial = c("Streptopelia risoria", "Nesoenas picturata", 
               "Streptopelia chinensis", "Porphyrio poliocephalus",
               "Sturnus burmannicus", "Copsychus malabaricus",
               "Aratinga holochlora",
               "Collocalia vanikorensis", "Coturnix chinensis",
               "Dendragapus canadensis", "Sturnus melanopterus",
               "Sturnus contra", "Urocissa erythrorhyncha",
               "Estrilda caerulescens", "Parus varius", "Bowdleria punctata"),
  accepted_name = c("Streptopelia roseogrisea", "Nesoenas picturatus",
                    "Spilopelia chinensis", "Porphyrio porphyrio",
                    "Acridotheres burmannicus", "Kittacincla malabarica", 
                    "Psittacara holochlorus",
                    "Aerodramus vanikorensis", "Excalfactoria chinensis",
                    "Falcipennis canadensis", "Acridotheres melanopterus",
                    "Gracupica contra", "Urocissa erythroryncha", 
                    "Glaucestrilda caerulescens" ,
                    "Sittiparus varius","Poodytes punctatus")
)

syno_all <- bind_rows(syno_iucn_eicat, manual_syno) %>% distinct()

for(i in 1:nrow(taxo_corresp)){
  for (j in 1:nrow(syno_all)){
    if(taxo_corresp$Species1[[i]]==syno_all$accepted_name[[j]]){
      taxo_corresp$bino_eicat[[i]]=syno_all$binomial[[j]]
    }
    if(taxo_corresp$Species1[[i]]==syno_all$binomial[[j]]){
      taxo_corresp$bino_eicat[[i]]=syno_all$accepted_name[[j]]
    }
  }
}

all_groups_tax <- left_join(all_groups, 
                            taxo_corresp %>% rename(binomial = bino_eicat),
                            by = "binomial")
colSums(is.na(all_groups_tax))


# create group db without dd eicat
ei_no_dd <- left_join(all_groups_tax %>% filter(combi!="0_0_0"),
                      coordinates %>% select(binomial, V1:V4) %>%
                        rename(Species1=binomial),
                      by= 'Species1') %>%
  mutate(combi2 = case_when(
    combi=="0_0_1" ~ "Only indirect",
    combi=="0_1_1" ~ "Direct_Indirect",
    combi=="1_1_1" ~ "All mecha",
    combi=="0_1_0" ~ "Only direct",
    combi=="1_0_0" ~ "Only ecosystem",
    combi=="1_1_0" ~ "Ecosystem_Direct",
    combi=="1_0_1" ~ "Ecosystem_Indirect")) %>%
  mutate(group2 = case_when(
    group2=="EICAT_imp" ~ "Alien birds with impact",
    group2=="EICAT_no_imp" ~ "Alien birds without impact",
    group2=="IAS-T_imp" ~ "Native birds impacted by IAS",
    group2=="IAS-T_no_imp" ~ "Native birds not impacted by IAS"
  )) %>%
  mutate_if(is.character, as.factor) %>%
  # remove sp with no match in traits (7 sp)
  filter(!is.na(V1))

colSums(is.na(ei_no_dd))

# remove all groups with less than ndim+1 species ?
table(ei_no_dd$combi, ei_no_dd$group)
ei_no_dd_m <- ei_no_dd %>% filter(!(group=="EICAT" & combi =="1_1_0")) %>%
  filter(!(group=="IAS-T" & combi =="1_0_1")) %>%
  filter(combi!="1_1_1") # remove cate all because not usefull
table(ei_no_dd_m$combi, ei_no_dd_m$group)

group_color = c("Direct_Indirect" = "darkorange1", "Ecosystem_Direct" = 
                  "darkorchid3", "Ecosystem_Indirect" = "chartreuse3",
                "Only direct" = "firebrick2","Only ecosystem" = "dodgerblue3",
                "Only indirect" = "gold")

# keep only alone mecha
ei_no_dd_alone <- ei_no_dd_m %>%
  filter(combi2 %in% c("Only direct", "Only ecosystem", "Only indirect"))
table(ei_no_dd_alone$combi2, ei_no_dd_alone$group)


#### Graphs ####

p1 <- ggplot(ei_no_dd_alone, aes(x = group, y = V1, fill = combi2)) +
  geom_flat_violin(aes(fill = combi2),
                   position = position_nudge(x = .15, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_boxplot(aes(x = group, y = V1, fill = combi2),
               outlier.shape = NA, alpha = .8, width = .3, colour = "grey30")+
  scale_fill_manual(values = group_color)+
  ylab("PC1") +
  labs(fill="", colour="") +
  theme_classic()+
  coord_flip()


p2 <- ggplot(ei_no_dd_alone, aes(x = group, y = V2, fill = combi2)) +
  geom_flat_violin(aes(fill = combi2),
                   position = position_nudge(x = .15, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_boxplot(aes(x = group, y = V2, fill = combi2),
               outlier.shape = NA, alpha = .8, width = .3, colour = "grey30")+
  scale_fill_manual(values = group_color)+ 
  ylab("PC2") +
  labs(fill="", colour="") +
  theme_classic()+
  coord_flip()

p3 <- ggplot(ei_no_dd_alone, aes(x = group, y = V3, fill = combi2)) +
  geom_flat_violin(position = position_nudge(x = 0.15, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_boxplot(aes(x = group, y = V3, fill = combi2),
               outlier.shape = NA, alpha = .8, width = .3, colour = "grey30")+
  scale_fill_manual(values = group_color)+
  ylab("PC3") +
  labs(fill="") +
  theme_classic()+
  coord_flip()

p4 <- ggplot(ei_no_dd_alone, aes(x = group, y = V4, fill = combi2)) +
  geom_flat_violin(position = position_nudge(x = 0.15, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_boxplot(aes(x = group, y = V4, fill = combi2),
               outlier.shape = NA, alpha = .8, width = .3, colour = "grey30")+
  scale_fill_manual(values = group_color)+
  ylab("PC4") +
  labs(fill="") +
  theme_classic()+
  coord_flip()

p1
p2
p3

b = ggarrange(p1, p2, p3, p4, ncol = 4, legend = "top", common.legend = T)

# save plot
pdf(file = "Fig/10_bird_mecha_alone_groups_pca.pdf", 10, 6)
print(b)
dev.off()

#### Statistical tests ####

axis = c("V1", "V2", "V3", "V4")
groups = c("EICAT","IAS-T")
duo = combn(c("1_0_0", "0_1_0", "0_0_1"), 2)

ks_df <- data.frame(
  Axis = rep(axis, each = 6),
  Group = rep(groups, times = 4, each = 3),
  a = character(24),
  b = character(24),
  D = numeric(24),
  p_val = numeric(24)
)
for (i in 1:ncol(duo)){
  ks_df[seq(i,24, by = 3),"a"] <- duo[1,i]
  ks_df[seq(i,24, by = 3),"b"] <- duo[2,i]
}

for(i in 1:nrow(ks_df)){
  gpa = pull(ei_no_dd %>% 
               filter(group==ks_df$Group[i] & combi==ks_df$a[i]), 
             ks_df$Axis[i])
  gpb = pull(ei_no_dd %>% 
               filter(group==ks_df$Group[i] & combi==ks_df$b[i]), 
             ks_df$Axis[i])
  k = ks.test(gpa, gpb)
  ks_df$D[i] = k$statistic
  ks_df$p_val[i] = k$p.value
}

ks_df$p_val_round = round(ks_df$p_val, 4)

write.table(ks_df, "Output/10_KS_test_mecha_alone.csv", sep = ";", row.names = F)

#### Correlation traits axis ####
df_ei_imp_final <- readRDS("Data/00_Complete_trait_data_birds")

eicat_iast_ft <- df_ei_imp_final %>%
  mutate_if(is.character, as.factor) %>%
  column_to_rownames("binomial") %>%
  # convert Mass, Beak & Clutch with log 
  mutate(ln.Mass = log(Mass),
         ln.Clutch = log(Clutch),
         ln.Beak.Depth = log(Beak.Depth),
         ln.Beak.Length_Nares = log(Beak.Length_Nares)) %>%
  # remove converted var + trophic niche
  select(-c(Mass, Clutch, Beak.Depth, Beak.Length_Nares, Trophic.Niche))

ei_ft_num <- bind_cols(eicat_iast_ft %>%
                         dplyr::select_if(is.numeric),
                       as.data.frame(dummy(eicat_iast_ft$Trophic.Level, sep="_")), 
                       as.data.frame(dummy(eicat_iast_ft$Primary.Lifestyle, sep="_")))

colnames(ei_ft_num)
# fix df names
names(ei_ft_num) <- c("Hand.Wing.Index", "Tail.Length", "hab_sum",
                      "insular_endemic", "volant", "ln.Mass",
                      "ln.Clutch", "ln.Beak.Depth", "ln.Beak.Length_Nares",
                      paste0("Diet.",levels(eicat_iast_ft$Trophic.Level)), 
                      paste0("ForN.",levels(eicat_iast_ft$Primary.Lifestyle)))
# change to matrix
data_ft_num <- as.matrix(ei_ft_num)

# Compute correlations between species coordinates and trait values
cor_mat<-rcorr(as.matrix(cbind(data_ft_num, coordinates %>% dplyr::select(V1:V4))), 
                          type = 'spearman')
# Matrix of correlation coeficients
R <- as.data.frame(round(cor_mat$r, 2)) %>% 
  select(V1:V4) %>%
  rownames_to_column("Trait") %>%
  filter(!(Trait %in% c("V1", "V2", "V3", "V4")))
# Matrix of p-value 
p <- as.data.frame(round(cor_mat$P,3))  %>% 
  select(V1:V4) %>%
  rownames_to_column("Trait") %>%
  filter(!(Trait %in% c("V1", "V2", "V3", "V4")))


cor_tab <- left_join(R, p, by="Trait")
str(cor_tab)

# write.table(cor_tab, "Output/10_Corr_traits_axis_pcoa.csv", sep=";", row.names = F)


#### Plot impacts along PCA axis ####

# IAS-T-imp, IAS-T-no-imp, EICAT-imp, EICAT-no-imp

p1i <- ggplot(ei_no_dd, aes(x= 1, y = V1, fill = group2)) +
  geom_flat_violin(aes(fill = group2),
                   position = position_nudge(x = .1, y = 0), 
                   adjust = 1, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = 1-.15, y = V1, colour = group2),
             position = position_jitter(width = .05), size = .25, shape = 20)+
  geom_boxplot(aes(x=1, y = V1, fill = group2),
               outlier.shape = NA, alpha = .5, width = .2, colour = "grey30")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab("PC1") +
  labs(fill="", colour="") +
  coord_flip()

p2i <- ggplot(ei_no_dd, aes(x= 1, y = V2, fill = group2)) +
  geom_flat_violin(aes(fill = group2),
                   position = position_nudge(x = .1, y = 0), 
                   adjust = 1, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = 1-.15, y = V2, colour = group2),
             position = position_jitter(width = .05), size = .25, shape = 20)+
  geom_boxplot(aes(x=1, y = V2, fill = group2),
               outlier.shape = NA, alpha = .5, width = .2, colour = "grey30")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ylab("PC2") +
  labs(fill="", colour="") +
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_flip()

p3i <- ggplot(ei_no_dd, aes(x= 1, y = V3, fill = group2)) +
  geom_flat_violin(aes(fill = group2),
                   position = position_nudge(x = .1, y = 0), 
                   adjust = 1, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = 1-.15, y = V3, colour = group2),
             position = position_jitter(width = .05), size = .25, shape = 20)+
  geom_boxplot(aes(x=1, y = V3, fill = group2),
               outlier.shape = NA, alpha = .5, width = .2, colour = "grey30")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ylab("PC3") +
  labs(fill="", colour="") +
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_flip()

p4i <- ggplot(ei_no_dd, aes(x= 1, y = V4, fill = group2)) +
  geom_flat_violin(aes(fill = group2),
                   position = position_nudge(x = .1, y = 0), 
                   adjust = 1, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = 1-.15, y = V4, colour = group2),
             position = position_jitter(width = .05), size = .25, shape = 20)+
  geom_boxplot(aes(x=1, y = V4, fill = group2),
               outlier.shape = NA, alpha = .5, width = .2, colour = "grey30")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ylab("PC4") +
  labs(fill="", colour="") +
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_flip()

a = ggarrange(p1i, p2i, p3i, p4i, ncol = 4, legend = "top", common.legend = T)


# save plot
pdf(file = "Fig/10_bird_impact_groups_pca.pdf", 8, 4)
print(a)
dev.off()

?ggarrange

p1i
p2i
p3i


#### Statistical test for comparing distributions ####

axis = c("V1", "V2", "V3", "V4")
duo_imp = combn(c("Alien birds with impact", "Alien birds without impact",
              "Native birds impacted by IAS", "Native birds not impacted by IAS"), 2)

ks_df_imp <- data.frame(
  Axis = rep(axis, each = 6),
  a = character(24),
  b = character(24),
  D = numeric(24),
  p_val = numeric(24)
)
for (i in 1:ncol(duo_imp)){
  ks_df_imp[seq(i,24, by = 6),"a"] <- duo_imp[1,i]
  ks_df_imp[seq(i,24, by = 6),"b"] <- duo_imp[2,i]
}

for(i in 1:nrow(ks_df_imp)){
  gpa = pull(ei_no_dd %>% 
               filter(group2 == ks_df_imp$a[i]), 
             ks_df_imp$Axis[i])
  gpb = pull(ei_no_dd %>% 
               filter(group2 == ks_df_imp$b[i]), 
             ks_df_imp$Axis[i])
  k = ks.test(gpa, gpb)
  ks_df_imp$D[i] = k$statistic
  ks_df_imp$p_val[i] = k$p.value
}

ks_df_imp$p_val_round = round(ks_df_imp$p_val, 4)

write.table(ks_df_imp, "Output/10_KS_test_impact.csv", sep = ";", row.names = F)




#### Plot ALIEN/IAS-T along PCA axis ####

# IAS-T all, EICAT all

p1 <- ggplot(ei_no_dd, aes(x= 1, y = V1, fill = group)) +
  geom_flat_violin(aes(fill = group),
                   position = position_nudge(x = .1, y = 0), 
                   adjust = 1, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = 1-.15, y = V1, colour = group),
             position = position_jitter(width = .05), size = .25, shape = 20)+
  geom_boxplot(aes(x=1, y = V1, fill = group),
               outlier.shape = NA, alpha = .5, width = .2, colour = "grey30")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab("PC1") +
  labs(fill="", colour="") +
  coord_flip()

p2 <- ggplot(ei_no_dd, aes(x= 1, y = V2, fill = group)) +
  geom_flat_violin(aes(fill = group),
                   position = position_nudge(x = .1, y = 0), 
                   adjust = 1, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = 1-.15, y = V2, colour = group),
             position = position_jitter(width = .05), size = .25, shape = 20)+
  geom_boxplot(aes(x=1, y = V2, fill = group),
               outlier.shape = NA, alpha = .5, width = .2, colour = "grey30")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ylab("PC2") +
  labs(fill="", colour="") +
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_flip()

p3 <- ggplot(ei_no_dd, aes(x= 1, y = V3, fill = group)) +
  geom_flat_violin(aes(fill = group),
                   position = position_nudge(x = .1, y = 0), 
                   adjust = 1, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = 1-.15, y = V3, colour = group),
             position = position_jitter(width = .05), size = .25, shape = 20)+
  geom_boxplot(aes(x=1, y = V3, fill = group),
               outlier.shape = NA, alpha = .5, width = .2, colour = "grey30")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ylab("PC3") +
  labs(fill="", colour="") +
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_flip()

p4 <- ggplot(ei_no_dd, aes(x= 1, y = V4, fill = group)) +
  geom_flat_violin(aes(fill = group),
                   position = position_nudge(x = .1, y = 0), 
                   adjust = 1, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = 1-.15, y = V4, colour = group),
             position = position_jitter(width = .05), size = .25, shape = 20)+
  geom_boxplot(aes(x=1, y = V4, fill = group),
               outlier.shape = NA, alpha = .5, width = .2, colour = "grey30")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ylab("PC4") +
  labs(fill="", colour="") +
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_flip()

a = ggarrange(p1, p2, p3, p4, ncol = 4, legend = "top", common.legend = T)


# save plot
pdf(file = "Fig/10_bird_role_groups_pca.pdf", 8, 4)
print(a)
dev.off()

#### Statistical test for comparing distributions ####

axis = c("V1", "V2", "V3", "V4")

ks_df_imp <- data.frame(
  Axis = axis,
  D = numeric(4),
  p_val = numeric(4)
)

for(i in 1:length(axis)){
  eicat = pull(ei_no_dd %>% filter(group == "EICAT"), 
             ks_df_imp$Axis[i])
  iast = pull(ei_no_dd %>% filter(group == "IAS-T"), 
             ks_df_imp$Axis[i])
  k = ks.test(eicat, iast)
  ks_df_imp$D[i] = k$statistic
  ks_df_imp$p_val[i] = k$p.value
}

ks_df_imp$p_val_round = round(ks_df_imp$p_val, 4)


#### Mecha && Magnitude of impact #######

# how many species in each group ?
table(ei_no_dd_alone$combi2, ei_no_dd_alone$group2)

# remove groups with less than 4 sp
# Native imp - indirect + alien imp - ecosys

ei_no_dd_alone_mi <- ei_no_dd_alone %>% 
  filter(!(group2=="Native birds impacted by IAS" & combi2 =="Only indirect")) %>%
  filter(!(group2=="Alien birds with impact" & combi2 =="Only ecosystem"))
table(ei_no_dd_alone_mi$combi2, ei_no_dd_alone_mi$group2)



pmi1 <- ggplot(ei_no_dd_alone_mi, aes(x = group2, y = V1, fill = combi2)) +
  geom_flat_violin(aes(fill = combi2),
                   position = position_nudge(x = .15, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_boxplot(aes(x = group2, y = V1, fill = combi2),
               outlier.shape = NA, alpha = .8, width = .3, colour = "grey30")+
  scale_fill_manual(values = group_color)+
  ylab("PC1 (22.5%)") + xlab("") +
  labs(fill="", colour="") +
  theme_classic()+
  theme(axis.text.y = element_blank())+
  coord_flip()


pmi2 <- ggplot(ei_no_dd_alone_mi, aes(x = group2, y = V2, fill = combi2)) +
  geom_flat_violin(aes(fill = combi2),
                   position = position_nudge(x = .15, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_boxplot(aes(x = group2, y = V2, fill = combi2),
               outlier.shape = NA, alpha = .8, width = .3, colour = "grey30")+
  scale_fill_manual(values = group_color)+ 
  ylab("PC2 (17.3%)") + xlab("") +
  labs(fill="", colour="") +
  theme_classic()+
  theme(axis.text.y = element_blank())+
  coord_flip()

pmi3 <- ggplot(ei_no_dd_alone_mi, aes(x = group2, y = V3, fill = combi2)) +
  geom_flat_violin(position = position_nudge(x = 0.15, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_boxplot(aes(x = group2, y = V3, fill = combi2),
               outlier.shape = NA, alpha = .8, width = .3, colour = "grey30")+
  scale_fill_manual(values = group_color)+
  ylab("PC3 (13.3%)") + xlab("") +
  labs(fill="") +
  theme_classic()+ 
  theme(axis.text.y = element_blank())+
  coord_flip()

pmi1
pmi2
pmi3

c = ggarrange(pmi1, pmi2, pmi3, ncol = 3, legend = "top", common.legend = T)

pdf(file = "Fig/15_bird_impact_meca_pca.pdf", 8, 6)
print(c)
dev.off()



######## Save data for traits within groups #######

saveRDS(left_join(all_groups_tax,
          coordinates %>%
            rename(Species1=binomial),
          by= 'Species1') %>%
  mutate(combi2 = case_when(
    combi=="0_0_1" ~ "Only indirect",
    combi=="0_1_1" ~ "Direct_Indirect",
    combi=="1_1_1" ~ "All mecha",
    combi=="0_1_0" ~ "Only direct",
    combi=="1_0_0" ~ "Only ecosystem",
    combi=="1_1_0" ~ "Ecosystem_Direct",
    combi=="1_0_1" ~ "Ecosystem_Indirect")),
  "Output/Data_clean/10_data_for_traits_groups")
  

