# DD prediction

rm(list=ls())

library(dplyr)
library(tidyr)
library(tibble)
library(cluster)
library(ggplot2)
library(dbscan)

#### Load data ####

df_ei_imp_final <- readRDS("Data/00_Complete_trait_data_birds")

tg <- readRDS("Output/02_data_for_traits_groups")
all_groups <- tg %>% filter(!is.na(hab_sum))

#### Set up species groups ####

ei_no_dd <- all_groups %>% 
  filter(combi!="0_0_0") %>%
  mutate_if(is.character, as.factor)

#### Compute distance matrix ####

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


ft_from_1_4 <- all_groups %>% 
  select(Species1, V1,V2,V3,V4) %>% 
  distinct() %>%
  column_to_rownames("Species1")

# distance matrix based on traits 
dis.gower <- daisy(ft_from_1_4,metric='gower')

summary(dis.gower)
matrix = as.matrix(dis.gower)
dim(matrix)
matrix[1,2]

# select only DD species in rows, IAS-T & EICAT no dd sp in columns
no_dd <- unique(as.character(pull(ei_no_dd, Species1)))
dd <- unique(as.character(pull(all_groups %>% filter(combi=="0_0_0"), Species1)))

df_dist <- as.data.frame(matrix) %>%
  select(all_of(no_dd)) %>%
  rownames_to_column("DD_sp") %>%
  filter(DD_sp %in% dd) %>%
  column_to_rownames("DD_sp")

# create species to group correspondence
corresp <- ei_no_dd %>% distinct(Species1, group2)


#### Calculate metrics ####

#---------------------- First closest sp

closest <- data.frame(dd_sp = row.names(df_dist), closest_sp = "")
for(i in 1:length(dd)){
  closest$closest_sp[i] <- names(df_dist)[which.min(df_dist[i,])]
}
# add group to which closest species belongs
closest_gp <- left_join(closest, ei_no_dd %>% select(Species1, group2) %>%
                             rename(closest_sp = Species1), 
                           by = "closest_sp")

table(closest_gp$group2)/length(dd)

#---------------------- find k closest species

k_closest <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(k_closest) <- paste("closest_", 1:10, sep="")

for(i in 1: length(dd)){
  col1 <- as.data.frame(t(df_dist[i,])) %>% rownames_to_column("sp")
  sorted <- col1[order(col1[,2]),]
  for (k in 1:10){
    k_closest[i,k]<- sorted$sp[k]
  } 
}

row.names(k_closest) <- row.names(df_dist)

k_closest_lg <- k_closest %>%
  rownames_to_column("DD_sp") %>%
  pivot_longer(!DD_sp, names_to = "rank", values_to = "sp_2")
  
k_closest_gps <- left_join(k_closest_lg, ei_no_dd %>% 
                             select(Species1, group2) %>%
                             rename(sp_2 = Species1), 
                           by = "sp_2") %>% 
  group_by(DD_sp, group2) %>% summarise(n_10_closest = n())


#----------------------- closest species in a buffer

# Calibrate buffer diameter
min(as.matrix(df_dist))
mean(as.matrix(df_dist))
median(as.matrix(df_dist))
max(as.matrix(df_dist))
# histogram of distance with mean distance
hist(as.matrix(df_dist), breaks = 100)
abline(v=mean(as.matrix(df_dist)), col='red')
abline(v=0.21, col='blue')
# cumulative curve of distance
# a=sort(as.matrix(df_dist))
# plot(a)

# Initialisation: set up d_buffer
d_buffer <- c("d1" = 0.21, 
              "d2" = mean(as.matrix(df_dist))) 

#source function for calculating sp in buffer, mean dist & prop in buffer
source("R/find_sp_in_buffer.R")

buffer_d1 <- find_sp_in_buffer_gp(df_dist, dd, d_buffer[1])

buffer_d2 <- find_sp_in_buffer_gp(df_dist, dd, d_buffer[2])

for(i in 1:length(buffer_d1)){
  buffer_d1[[i]] <- buffer_d1[[i]] %>% 
    mutate(DD_sp = names(buffer_d1[i]))
  buffer_d2[[i]] <- buffer_d2[[i]] %>% 
    mutate(DD_sp = names(buffer_d2[i]))
}

buffer_d1_df <- do.call(rbind.data.frame, buffer_d1) %>%
  rename(mean_b1 = mean_buffer,
         n_b1 = n_b,
         prop_b1 = prop_b)
buffer_d2_df <- do.call(rbind.data.frame, buffer_d2) %>%
  rename(mean_b2 = mean_buffer,
         n_b2 = n_b,
         prop_b2 = prop_b)

# How many sp in buffers?
# mean(unlist(lapply(sp_in_buffer, nrow)))
# min(unlist(lapply(sp_in_buffer, nrow)))
# max(unlist(lapply(sp_in_buffer, nrow)))
# median(unlist(lapply(sp_in_buffer, nrow)))


#----------- Mean & min distance with groups

long_df <- df_dist %>% rownames_to_column("DD_sp") %>%
  pivot_longer(!DD_sp, names_to = "sp2", values_to = "dist")

long_df_gp <- left_join(long_df, ei_no_dd %>% select(Species1, group2) %>%
                          rename(sp2 = Species1), 
                        by = "sp2") %>%
  group_by(DD_sp, group2) %>%
  summarise(mean_dist = mean(dist), sd_dist = sd(dist), min_dist = min(dist))

head(long_df_gp)


#-------------- Final table with all methods for DD prediction

predict_dd <- left_join(
  # Add k closest neighbours
  left_join(long_df_gp %>% select(-sd_dist), k_closest_gps, by=c("DD_sp","group2")),
  # Add buffer metrics
  # join with d2 first because no NA
  left_join(buffer_d2_df, buffer_d1_df, by=c("DD_sp","group2"))) %>%
  mutate_at(c("n_10_closest","n_b1","prop_b1"),  ~ ifelse(is.na(.), 0, .))

saveRDS(predict_dd, "Output/04_metrics_predict_dd_pc1to4")



#### Attribute a status to each method ####

# for min_dist & mean_dist, attribute status of minimal distance group
# idem for mean_b1 & mean_b2: minimal mean distance in each buffer

# for proportion of closest species:
# take group status for which the proportion is maximal when removing null_prop

# load metrics file
predict_dd <- readRDS("Output/04_metrics_predict_dd_pc1to4") %>%
  mutate(group2 = as.character(group2))

# calculate null proportions of group representation
null_prop <- c(table(ei_no_dd$group2)/nrow(ei_no_dd))

predict_dd <- left_join(predict_dd,
                        as.data.frame(null_prop) %>% rownames_to_column("group2"),
                        by = "group2")

predict_dd <- predict_dd %>%
  mutate(prop_10_closest = n_10_closest/10) %>%
  mutate(prop_b1_corr = prop_b1 - null_prop,
         prop_b2_corr = prop_b2 - null_prop,
         prop_10_closest_corr = prop_10_closest - null_prop) %>%
  mutate(prop_b1_corr = if_else(is.na(mean_b1), mean_b1, prop_b1_corr),
         prop_b2_corr = if_else(is.na(mean_b2), mean_b2, prop_b2_corr))

# initialize empty df to fill with predicted status
status_all_meth <- data.frame(matrix(ncol = 8, nrow = length(dd)))
names(status_all_meth) <- c("DD_sp","min_dist","mean_dist","mean_b1", 
                            "prop_b1_corr", "mean_b2", "prop_b2_corr",
                            "prop_10_closest_corr")
status_all_meth$DD_sp <- dd

# attribute status for each dd sp & each method
for(sp in dd){
  predict_dd_sp <- predict_dd %>% filter(DD_sp == sp)
  
  # vars for which to take minimal value
  vars_min <- c("min_dist","mean_dist","mean_b1","mean_b2")
  
  for(col in vars_min){
    if(sum(is.na(pull(predict_dd_sp,col)))!=3){
      status_all_meth[which(status_all_meth$DD_sp==sp),col] <-
        predict_dd_sp$group2[which.min(pull(predict_dd_sp,col))]
    } else {
      status_all_meth[which(status_all_meth$DD_sp==sp),col] <- NA
    }
  }
  # vars for which to take maximal value after correction by null_prop
  vars_max <- c("prop_b1_corr", "prop_b2_corr", "prop_10_closest_corr")
  for(col2 in vars_max){
    if(sum(is.na(pull(predict_dd_sp,col2)))!=3){
      status_all_meth[which(status_all_meth$DD_sp==sp),col2] <-
        predict_dd_sp$group2[which.max(pull(predict_dd_sp,col2))]
    } else {
      status_all_meth[which(status_all_meth$DD_sp==sp),col2] <- NA
    }
  }

}

write.table(status_all_meth, "Output/04_Pred_status_all_methods_pc1to4.csv",
            sep=";", row.names = F)
saveRDS(status_all_meth, "Output/04_Pred_status_all_methods_pc1to4")


#### Predict a final status ####

# Select best metrics (following 31_Evaluate_DD_prediction)

best_metrics <- status_all_meth %>% select(DD_sp, min_dist, mean_b2)

best_metrics$IAS_T_no_imp <- rowSums(best_metrics == "IAS-T_no_imp", na.rm = T)
best_metrics$IAS_T_imp <- rowSums(best_metrics == "IAS-T_imp", na.rm = T)
best_metrics$EICAT_imp <- rowSums(best_metrics == "EICAT_imp", na.rm = T)
best_metrics$EICAT_no_imp  <- rowSums(best_metrics == "EICAT_no_imp", na.rm = T)

# Absolute majority (all methods vote for one status)
status_vote <- c("IAS_T_imp", "IAS_T_no_imp", "EICAT_imp", "EICAT_no_imp")

for(i in 1:nrow(best_metrics)){
  if(max(best_metrics[i,status_vote])>1){
    best_metrics$majo_abs_vote[i] <-
      status_vote[which.max(best_metrics[i,status_vote])]
  } else {
    best_metrics$majo_abs_vote[i] <- "undefined"
  }
}
table(best_metrics$majo_abs_vote)

# add majority for 6 x 6 + absolute vote 
best_metrics2 <- status_all_meth %>% 
  select(DD_sp, min_dist, mean_b1, mean_b2,
         prop_b1_corr, prop_b2_corr, prop_10_closest_corr)

best_metrics2$IAS_T_no_imp <- rowSums(best_metrics2 == "IAS-T_no_imp", na.rm = T)
best_metrics2$IAS_T_imp <- rowSums(best_metrics2 == "IAS-T_imp", na.rm = T)
best_metrics2$EICAT_imp <- rowSums(best_metrics2 == "EICAT_imp", na.rm = T)
best_metrics2$EICAT_no_imp  <- rowSums(best_metrics2 == "EICAT_no_imp", na.rm = T)

# Absolute majority (all methods vote for one status)
status_vote <- c("IAS_T_imp", "IAS_T_no_imp", "EICAT_imp", "EICAT_no_imp")
for(i in 1:nrow(best_metrics2)){
  if(max(best_metrics2[i,status_vote])>4){
    best_metrics2$majo_abs_vote[i] <-
      status_vote[which.max(best_metrics2[i,status_vote])]
  } else {
    best_metrics2$majo_abs_vote[i] <- "undefined"
  }
}
table(best_metrics2$majo_abs_vote)

for(i in 1:nrow(best_metrics2)){
  best_metrics2$demo_vote[i] <- 
    status_vote[which.max(best_metrics2[i,status_vote])]
}
table(best_metrics2$demo_vote)



final_pred <- left_join(best_metrics %>% 
                          select(DD_sp, majo_abs_vote) %>%
                          rename(majo_abs_2var = majo_abs_vote), 
                        best_metrics2 %>% 
                          select(DD_sp, majo_abs_vote, demo_vote)%>%
                          rename(majo_abs_6var = majo_abs_vote,
                                 demo_6var = demo_vote),
                        by="DD_sp")

write.table(final_pred, "Output/04_Predicted_status_all_dd_pc1to4.csv", sep=";",
            row.names = F)




#### Vizualise Mean distance with EICAT imp/no_imp & IAS-T groups ####

# vizualise mean distance with each group for all eicat_dd species
ggplot(long_df_gp, aes(x = mean_dist)) +
  geom_histogram() +
  facet_wrap(~group2)

# check if there is a group of species close to EICAT_imp but not to 
# other groups (IAS-T and EICAT_no_imp)
wide_df_gp  <- long_df_gp %>% select(-sd_dist) %>%
  pivot_wider(names_from = group2, values_from = mean_dist)

ggplot(wide_df_gp) +
  geom_point(aes(x=EICAT_imp, y=`IAS-T`, color = EICAT_no_imp )) +
  geom_abline(slope=1, intercept = 0)
ggplot(wide_df_gp) +
  geom_point(aes(x=EICAT_imp, y=EICAT_no_imp)) +
  geom_abline(slope=1, intercept = 0)

ggplot(wide_df_gp) +
  geom_point(aes(x=EICAT_imp, y=EICAT_no_imp, color= `IAS-T`)) +
  geom_abline(slope=1, intercept = 0)

