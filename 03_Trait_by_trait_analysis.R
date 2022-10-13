# assess traits of mecha & impact magni groups

rm(list=ls())

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(car)
library(emmeans)


# load group data 
tg_combi <- readRDS("Output/02_data_for_traits_groups")
tg <- tg_combi %>% filter(!is.na(hab_sum)) %>%
  filter(combi!="0_0_0") %>%
  select(group2:ln.Beak.Length_Nares) %>%
  rename(binomial = Species1)

# load trait data
all_birds <- readRDS("Data/00_Complete_trait_data_birds") %>%
  mutate_if(is.character, as.factor) %>%
  # convert Mass, Beak & Clutch with log 
  mutate(ln.Mass = log(Mass),
         ln.Clutch = log(Clutch),
         ln.Beak.Depth = log(Beak.Depth),
         ln.Beak.Length_Nares = log(Beak.Length_Nares)) %>%
  # remove converted var + trophic niche
  select(-c(Mass, Clutch, Beak.Depth, Beak.Length_Nares, Trophic.Niche)) %>%
  mutate(group2 = "All_birds")

tg_all <- bind_rows(tg, all_birds)

str(tg_all)

length(tg_all$binomial[duplicated(tg_all$binomial)])

str(tg_all)
group_order <- c("All_birds", "IAS-T_imp", "IAS-T_no_imp",
                 "EICAT_no_imp", "EICAT_imp")
tg_all <- tg_all %>%
  mutate(group2 = as.factor(group2)) %>%
  mutate(group2 = factor(group2, level = group_order))

#####################################################################

# Impact versus no impact

#---------------------- numeric var

# HWI => ns 
a <- ggplot(tg_all, aes(x= group2, y= Hand.Wing.Index)) +
  geom_boxplot(alpha = 0.5)+
  geom_point(data=tg_all %>% filter(group2 != "All_birds"),
             position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) +
  geom_boxplot(alpha = 0.5) +
  stat_summary(fun = "mean", color="violetred3", shape=18, size = 1)+ 
  #xlab("") +
  theme_classic2()

# Tail Length => signif EICAT vs IAS-T
b <- ggplot(tg_all, aes(group2, Tail.Length)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(data=tg_all %>% filter(group2 != "All_birds"),
             position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) + 
  geom_boxplot(alpha = 0.5) +
  stat_summary(fun = "mean", color="violetred3", shape=18, size = 1)+ 
  xlab("") +
  theme_classic2()

lm_var <- lm(Tail.Length ~ group2, data=tg)
lsmeans(lm_var, pairwise ~ group2, adjust = "bonferroni")

# Correct by body mass 
ggplot(tg_all, aes(group2, Tail.Length/ln.Mass)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(data=tg_all %>% filter(group2 != "All_birds"),
             position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) + 
  geom_boxplot(alpha = 0.5) +
  xlab("") +
  theme_classic2()
lm_var <- lm(Tail.Length/ln.Mass ~ group2, data=tg)
lsmeans(lm_var, pairwise ~ group2, adjust = "bonferroni")

# Hab sum => signif eicat vs IAS-T
ggplot(tg_all, aes(group2, hab_sum)) +
  geom_boxplot()+
  geom_point(position = position_jitter(width = .1)) + 
  xlab("")

lm_var <- lm(hab_sum ~ group2, data=tg)
lsmeans(lm_var, pairwise ~ group2, adjust = "bonferroni")

# Mass => signif EICAT > IAS-T imp
e <- ggplot(tg_all, aes(group2, ln.Mass)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(data=tg_all %>% filter(group2 != "All_birds"),
             position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) + 
  geom_boxplot(alpha = 0.5) +
  stat_summary(fun = "mean", color="violetred3", shape=18, size = 1)+ 
  xlab("") +
  theme_classic2()

lm_var <- lm(ln.Mass ~ group2, data=tg )
lsmeans(lm_var, pairwise ~ group2, adjust = "bonferroni")

# Clutch => signif EICAT > IAS-T
f <- ggplot(tg_all, aes(group2, ln.Clutch)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(data=tg_all %>% filter(group2 != "All_birds"),
             position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) + 
  geom_boxplot(alpha = 0.5) +
  stat_summary(fun = "mean", color="violetred3", shape=18, size = 1)+ 
  xlab("") +
  theme_classic2()

lm_var <- lm(ln.Clutch ~ group2, data=tg)
Anova(lm_var)
lsmeans(lm_var, pairwise ~ group2, adjust = "bonferroni")

# Beak Depth => EICAT no imp > IAS-T
c <- ggplot(tg_all, aes(group2, ln.Beak.Depth)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(data=tg_all %>% filter(group2 != "All_birds"),
             position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) + 
  geom_boxplot(alpha = 0.5) +
  stat_summary(fun = "mean", color="violetred3", shape=18, size = 1)+ 
  xlab("") +
  theme_classic2()
lm_var <- lm(ln.Beak.Depth ~ group2, data=tg)
lsmeans(lm_var, pairwise ~ group2, adjust = "bonferroni")

# correct by bodymass
# Eicat no imp > ias-t no imp
# eicat no imp > eicat imp
ggplot(tg_all, aes(group2, ln.Beak.Depth/ln.Mass)) +
  geom_point(data=tg_all %>% filter(group2 != "All_birds"),
             position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) + 
  geom_boxplot(alpha = 0.5) +
  xlab("") +
  theme_classic2()
lm_var <- lm(ln.Beak.Depth/ln.Mass ~ group2, data=tg)
lsmeans(lm_var, pairwise ~ group2, adjust = "bonferroni")

# Beak Length => ns
d <- ggplot(tg_all, aes(group2, ln.Beak.Length_Nares))+
  geom_boxplot(alpha = 0.5) +
  geom_point(data=tg_all %>% filter(group2 != "All_birds"),
             position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) +
  geom_boxplot(alpha = 0.5) +
  stat_summary(fun = "mean", color="violetred3", shape=18, size = 1)+ 
  xlab("") +
  theme_classic2()

# correct by body mass
#eicat < ias-t imp
ggplot(tg_all, aes(group2, ln.Beak.Length_Nares/ln.Mass))+
  geom_point(data=tg_all %>% filter(group2 != "All_birds"),
             position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) +
  geom_boxplot(alpha = 0.5) +
  stat_summary(fun = "mean", color="violetred3", shape=18, size = 1)+ 
  xlab("") +
  theme_classic2()
lm_var <- lm(ln.Beak.Length_Nares/ln.Mass ~ group2, data=tg)
lsmeans(lm_var, pairwise ~ group2, adjust = "bonferroni")

pdf("Fig/11_Single_traits_num_impact_w.pdf", 10, 7)
ggarrange(a, b, e, c, d, f, ncol=3, nrow = 2 )
dev.off()

# test mean and sd

tg_all_mean <- tg_all %>%
  group_by(group2) %>%
  summarise(Hand.Wing.Index = mean(Hand.Wing.Index),
            Tail.Length = mean(Tail.Length),
            ln.Mass = mean(ln.Mass),
            ln.Clutch = mean(ln.Clutch),
            ln.Beak.Depth = mean(ln.Beak.Depth),
            ln.Beak.Length_Nares = mean(ln.Beak.Length_Nares)) %>%
  pivot_longer(!group2, names_to = "metric", values_to = "mean")

tg_all_sd <- tg_all %>%
  group_by(group2) %>%
  summarise(Hand.Wing.Index = sd(Hand.Wing.Index),
            Tail.Length = sd(Tail.Length),
            ln.Mass = sd(ln.Mass),
            ln.Clutch = sd(ln.Clutch),
            ln.Beak.Depth = sd(ln.Beak.Depth),
            ln.Beak.Length_Nares = sd(ln.Beak.Length_Nares)) %>%
  pivot_longer(!group2, names_to = "metric", values_to = "sd")

tg_all_summ <- left_join(tg_all_mean, tg_all_sd, by = c("group2", "metric"))
  
ggplot(tg_all_summ, aes(group2, mean)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  facet_wrap(~metric, scales = "free_y", nrow=2)


#---------------------- Categorial var


# trophic level
a1 <- ggplot(tg_all, aes(group2, fill=Trophic.Level)) +
  geom_bar(aes(y=..count..), position = "fill") +
  xlab("") +ylab("")+
  theme_classic()+
  theme(legend.position = "top")+
  scale_fill_viridis_d(option = "viridis", direction = -1) +
  scale_color_viridis_d(option = "viridis", direction = -1)
# ggplot(tg, aes(group2, fill=Trophic.Niche)) +
#   geom_bar(aes(y=..count..), position = "fill")

# habitat
b1 <- ggplot(tg_all %>%
         mutate(hab_cate = as.character(hab_sum)) %>%
         mutate(hab_cate = if_else(hab_sum>4, "5+", hab_cate)), 
       aes(group2, fill=hab_cate)) +
  geom_bar(aes(y=..count..), position = "fill") +
  xlab("") +ylab("")+
  theme_classic()+
  theme(legend.position = "top")+
  scale_fill_viridis_d(option = "viridis", direction = -1) +
  scale_color_viridis_d(option = "viridis", direction = -1)

# primary lifestyle
c1 <- ggplot(tg_all,
       aes(group2, fill=Primary.Lifestyle)) +
  geom_bar(aes(y=..count..), position = "fill") +
  xlab("") +ylab("")+
  theme_classic()+
  theme(legend.position = "top")+
  scale_fill_viridis_d(option = "viridis", direction = -1) +
  scale_color_viridis_d(option = "viridis", direction = -1)

pdf("Fig/11_Single_traits_cate_impact.pdf", 7, 3)
ggarrange(a1, b1, c1, ncol=3, nrow = 1)
dev.off()
#---------------------- binary var

# volant
table(tg$volant,tg$group2)
vol <- ggplot(tg_all, aes(group2, fill=as.character(volant))) +
  geom_bar(aes(y=..count..), position = "fill") +
  xlab("") +ylab("")+
  theme_classic()+
  theme(legend.position = "top")+
  scale_fill_viridis_d(option = "viridis", direction = -1) +
  scale_color_viridis_d(option = "viridis", direction = -1)

# insular endemic
table(tg$insular_endemic, tg$group2)
ins <- ggplot(tg_all, aes(group2, fill=as.character(insular_endemic))) +
  geom_bar(aes(y=..count..), position = "fill") +
  xlab("") + ylab("")+
  theme_classic()+
  theme(legend.position = "top")+
  scale_fill_viridis_d(option = "viridis", direction = -1) +
  scale_color_viridis_d(option = "viridis", direction = -1)

pdf("Fig/11_Single_traits_bin_impact.pdf", 6, 4)
ggarrange(vol, ins, ncol=2, nrow = 1)
dev.off()

# test for IAS-T imp / no-imp (insular endemic)
x <- matrix(c(57, 131, 111, 157), ncol = 2)
chisq.test(x)
x


pdf("Fig/11_Single_traits_cate_bin_impact.pdf", 7, 6)
ggarrange(a1, b1, c1, vol, ins, ncol=3, nrow = 2)
dev.off()
##################################################################

# Memes analyses avec les mécanismes
# alone et combinés (attention au nb de points par groupes)

tg_m <- tg_combi %>%
  filter(complete.cases(.)) %>%
  filter(! combi %in% c("0_0_0","1_1_1")) %>%
  mutate(mecha = paste(group, combi, sep="_")) %>%
  filter(! mecha %in% c("EICAT_1_1_0","IAS-T_1_0_1", "EICAT_1_0_1",
                        "EICAT_0_1_1", "IAS-T_0_1_1", "IAS-T_1_1_0")) %>%
  mutate(hab_cate = as.character(hab_sum)) %>%
  mutate(hab_cate = if_else(hab_sum>4, "5+", hab_cate))

table(tg_m$mecha)

#---------------------- numeric var

# HWI => ns sauf IAST indir > EICAT ecosys
ggplot(tg_m, aes(mecha, Hand.Wing.Index)) +
  geom_boxplot()+
  geom_point(position = position_jitter(width = .1))
lm_var <- lm(Hand.Wing.Index ~ mecha, data=tg_m)
lsmeans(lm_var, pairwise ~ mecha, adjust = "bonferroni")

# Tail Length => ns sauf EICAT vs IAS-T (but nothing on mecha)
ggplot(tg_m, aes(mecha, Tail.Length)) +
  geom_boxplot()+
  geom_point(position = position_jitter(width = .1))
lm_var <- lm(Tail.Length ~ mecha, data=tg_m)
lsmeans(lm_var, pairwise ~ mecha, adjust = "bonferroni")

# Hab sum => signif eicat dir vs eicat indir
ggplot(tg_m, aes(mecha, hab_sum)) +
  geom_boxplot()+
  geom_point(position = position_jitter(width = .1))
lm_var <- lm(hab_sum ~ mecha, data=tg_m)
lsmeans(lm_var, pairwise ~ mecha, adjust = "bonferroni")

# Mass => IAS-T dir > IAS-T ecosys :O + eicat vs iast
ggplot(tg_m, aes(mecha, ln.Mass)) +
  geom_boxplot()+
  geom_point(position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) + 
  xlab("") +
  theme_classic2()
lm_var <- lm(ln.Mass ~ mecha, data=tg_m)
lsmeans(lm_var, pairwise ~ mecha, adjust = "bonferroni")

# Clutch => ns sauf pour les EICAT vs IAS-T
ggplot(tg_m, aes(mecha, ln.Clutch)) +
  geom_boxplot()+
  geom_point(position = position_jitter(width = .1))
lm_var <- lm(ln.Clutch ~ mecha, data=tg_m)
Anova(lm_var)
lsmeans(lm_var, pairwise ~ mecha, adjust = "bonferroni")

# Beak Depth => IAS-T indir > IAS-T ecosys
ggplot(tg_m, aes(mecha, ln.Beak.Depth)) +
  geom_boxplot()+
  geom_point(position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) + 
  xlab("") +
  theme_classic2()
lm_var <- lm(ln.Beak.Depth ~ mecha, data=tg_m %>% filter(combi!="0_0_0"))
lsmeans(lm_var, pairwise ~ mecha, adjust = "bonferroni")

# Beak Length => IAS-T dir > IAS-T ecosys
ggplot(tg_m, aes(mecha, ln.Beak.Length_Nares)) +
  geom_boxplot()+
  geom_point(position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) + 
  xlab("") +
  theme_classic2()
lm_var <- lm(ln.Beak.Length_Nares ~ mecha, data=tg_m %>% filter(combi!="0_0_0"))
lsmeans(lm_var, pairwise ~ mecha, adjust = "bonferroni")

# coreected by body mass => IAS-T dir < IAS-T ecosys !
ggplot(tg_m, aes(mecha, ln.Beak.Length_Nares/ln.Mass)) +
  geom_boxplot()+
  geom_point(position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) + 
  xlab("") +
  theme_classic2()
lm_var <- lm(ln.Beak.Length_Nares/ln.Mass ~ mecha, data=tg_m %>% filter(combi!="0_0_0"))
lsmeans(lm_var, pairwise ~ mecha, adjust = "bonferroni")

#---------------------- Categorial var

# reprendre ici avec test de khi2

# trophic level
ggplot(tg_m, aes(mecha, fill=Trophic.Level)) +
  geom_bar(aes(y=..count..), position = "fill")

chisq.test(table(tg_m$mecha, tg_m$Trophic.Level))

# habitat
ggplot(tg_m, aes(mecha, fill=hab_cate)) +
  geom_bar(aes(y=..count..), position = "fill")
chisq.test(table(tg_m$mecha, tg_m$hab_cate))


# primary lifestyle
ggplot(tg_m, aes(mecha, fill=Primary.Lifestyle)) +
  geom_bar(aes(y=..count..), position = "fill")
chisq.test(table(tg_m$mecha, tg_m$Primary.Lifestyle))


#---------------------- binary var

# volant
table(tg_m$volant,tg_m$mecha)
ggplot(tg_m, aes(mecha, fill=as.character(volant))) +
  geom_bar(aes(y=..count..), position = "fill")
chisq.test(table(tg_m$mecha, as.character(tg_m$volant)))

# insular endemic
table(tg_m$insular_endemic, tg_m$mecha)
ggplot(tg_m, aes(mecha, fill=as.character(insular_endemic))) +
  geom_bar(aes(y=..count..), position = "fill")
chisq.test(table(tg_m$mecha, as.character(tg_m$insular_endemic)))
