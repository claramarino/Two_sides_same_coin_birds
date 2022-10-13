# Visualisation of impact and mecha groups along PCA axis

rm(list=ls())

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
#library(dummies)
library(Hmisc)
source("R/R_rainclouds.R")

#### Load data ####

# PCA from pcoa function
coordinates <- readRDS( "Output/01_PCA_coord_10D_all_birds_pcoa")
# bird groups and impact information
alien_iast <- readRDS("Data/00_Alien_IAST_birds_impact_group")


#### Bind tables ####

# join bird groups and bird traits
# create group database without dd eicat
ei_no_dd <- left_join(alien_iast %>% filter(combi!="0_0_0"),
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



#### Fig 1 -- Impact groups along PCA axis ####

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

a = ggarrange(p1i, p2i, p3i, ncol = 3, legend = "top", common.legend = T)


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





#### Graphs - Figure 4 ####

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

b = ggarrange(p1, p2, p3, ncol = 3, legend = "top", common.legend = T)

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




#### Table 1 -- Correlation traits axis ####

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



######## Save data for traits within groups #######

saveRDS(left_join(alien_iast,
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
  "Output/02_data_for_traits_groups")
  

