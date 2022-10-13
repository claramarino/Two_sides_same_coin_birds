# Functional space in n dim using pca

rm(list=ls())

#__________________________________________________

# load df_eicat_traits & ias_birds_mecha_traits
# bind in one single df with traits only
# compute distance matrix between species
# compute PCA on distance matrix
# keep 10 first axis
# compute correlation between traits and PCA 10 first axis
# save outputs to open in 11_PCA_ecol_surf_plots
#__________________________________________________

library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(cluster)
library(tibble)
library(dummies)

#### Load data ####

# load data with complete traits for IAS-T and alien species 
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

str(eicat_iast_ft)

#### Run PCoA on traits ####

# 1. compute the dissimilarity matrix between species 
# (Gower distance for non numeric var) & compute PCoA on dis matrix
# 2. use a pca function which can include mixed data

# Mouillot et al 2020
mat_dis <- cluster::daisy(eicat_iast_ft, metric = "gower")
mat_pcoa <- ape::pcoa(mat_dis)

nbdim = 10
nbdim<-min(nbdim,ncol(mat_pcoa$vectors) )
# keeping species coordinates on the 'nbdim' axes
mat_coord<-mat_pcoa$vectors[,1:nbdim]
row.names(mat_coord)<-row.names(eicat_iast_ft)
colnames(mat_coord)<-paste("PC",1:nbdim,sep="")

sum(mat_pcoa$values$Eigenvalues)
mat_pcoa$values$Eigenvalues[1]/sum(mat_pcoa$values$Eigenvalues)



#### Add trait correlation with ecological space ####

# select coordinates from pcoa
coordinates = as.data.frame(mat_coord) 
colnames(coordinates) <- paste0("V", 1:ncol(coordinates))

# Shape traits as numeric or binary variables (with dummy)
# to compute correlation with species
str(eicat_iast_ft)

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
sc_ar <- 0.5 # scale arrows
corr_sp_traits = as.data.frame(
  cor(data_ft_num, coordinates, method = 'spearman')*sc_ar) %>%
  rownames_to_column("trait")
# method = spearman for non parametric test

# Calculate p.values for each correlation
p_values_spe <- sapply(data.frame(coordinates),
                       function(x) Map(function(a,b) cor.test(a,b, method='spearman')$p.value,
                                       list(x),as.data.frame(data_ft_num)))
p_values <- as.data.frame(p_values_spe)
colnames(p_values) = colnames(corr_sp_traits %>% select(-trait))
rownames(p_values) = corr_sp_traits$trait
p_values[p_values>0.05]<-'NA'


# save PCA data & outputs
pc_to_save <- bind_cols(eicat_iast_ft, coordinates) %>%
  rownames_to_column("binomial")

saveRDS(pc_to_save, "Output/01_PCA_coord_10D_all_birds_pcoa")

