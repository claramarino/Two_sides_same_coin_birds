# Evaluation of status prediction for DD species

rm(list=ls())

library(dplyr)
options(dplyr.summarise.inform = FALSE)
library(tidyr)
library(cluster)
library(tibble)
library(ggplot2)

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

# distance matrix based on traits 
dis.gower <- daisy(eicat_iast_ft,metric='gower')
summary(dis.gower)
matrix = as.matrix(dis.gower)
dim(matrix)
matrix[1,2]

# select only no_DD species in rows & in columns
# we will evaluate the prediction using only species with a known status

no_dd <- unique(as.character(pull(ei_no_dd, Species1)))

df_dist <- as.data.frame(matrix)  %>%
  select(all_of(no_dd)) %>%
  rownames_to_column("no_DD_sp") %>%
  filter(no_DD_sp %in% no_dd) %>%
  column_to_rownames("no_DD_sp")

# create species to group correspondence
corresp <- ei_no_dd %>% distinct(Species1, group)

#### Method 'Out of the bag' ####

# fill the  table for each of the 574 species (to remove in the predictors)

#------------- Initialize tables for metrics

# for min and mean distance
long_df_all_sp <- data.frame(
  no_DD_sp = character(0), group = character(0),
  mean_dist = numeric(0), min_dist = numeric(0))
# 10 closest neighbors
k_closest <- k_closest_sp <- data.frame(matrix(ncol = 10, nrow = 0)) %>%
  mutate_all(as.character)
colnames(k_closest) <- colnames(k_closest_sp) <- paste("closest_", 1:10, sep="")
# species in buffers 
d_buffer <- c("d1" = 0.21, 
              "d2" = mean(as.matrix(df_dist))) 
buffer_df_all <- data.frame(
  group = character(0), mean_buffer = numeric(0), n_b = numeric(0),
  prop_b = numeric(0), no_DD_sp = character(0),d_buffer=character(0))
#source function for calculating sp in buffer, mean dist & prop in buffer
source("R/find_sp_in_buffer.R")

#------------ Loop on all species

for (sp in no_dd){
  
  # test iteration
  #sp = no_dd[6]
  
  # select focus species & remove it from the predictors
  df_dist_out <- df_dist %>%
    select(-all_of(sp))%>%
    rownames_to_column("no_DD_sp") %>%
    filter(no_DD_sp == sp) %>%
    column_to_rownames("no_DD_sp")
  
  #___________________ Min & mean distance with groups
  
  long_df <- df_dist_out %>% rownames_to_column("no_DD_sp") %>%
    pivot_longer(!no_DD_sp, names_to = "sp2", values_to = "dist")
  
  long_df_sp <- left_join(long_df, corresp %>% rename(sp2 = Species1), 
                          by = "sp2") %>%
    group_by(no_DD_sp, group) %>%
    summarise(mean_dist = mean(dist), min_dist = min(dist))
  
  # implement in final long df with all no_dd species
  long_df_all_sp <- bind_rows(long_df_all_sp, long_df_sp)
  
  #____________________ 10 closest neighbours
  
  col1 <- as.data.frame(t(df_dist_out)) %>% rownames_to_column("sp")
  sorted <- col1[order(col1[,2]),]
  for (k in 1:10){
      k_closest_sp[1,k]<- sorted$sp[k]
  } 
  row.names(k_closest_sp) <- sp
  
  # implement in final k_closest dataset
  k_closest <- bind_rows(k_closest, k_closest_sp)
  
  #____________________ species in buffers
  
  buffer_d1 <- find_sp_in_buffer_role(df_dist_out, sp, d_buffer[1])
  buffer_d2 <- find_sp_in_buffer_role(df_dist_out, sp, d_buffer[2])
  
  buffer_d1[[1]] <- buffer_d1[[1]] %>% 
    mutate(no_DD_sp = sp, d_buffer = "d1")
  buffer_d2[[1]] <- buffer_d2[[1]] %>% 
    mutate(no_DD_sp = sp, d_buffer = "d2")
  
  buffer_df_all <- bind_rows(buffer_df_all, 
                             do.call(rbind.data.frame, buffer_d1),
                             do.call(rbind.data.frame, buffer_d2))
  
}

#___________________ Final metrics dataset

k_closest_lg <- k_closest %>%
  rownames_to_column("no_DD_sp") %>%
  pivot_longer(!no_DD_sp, names_to = "rank", values_to = "sp_2")

k_closest_gps <- left_join(k_closest_lg, corresp %>%
                             rename(sp_2 = Species1), 
                           by = "sp_2") %>% 
  group_by(no_DD_sp, group) %>% summarise(n_10_closest = n())

buffer_d1_df <- buffer_df_all %>%
  filter(d_buffer == "d1") %>%
  rename(mean_b1 = mean_buffer,
         n_b1 = n_b,
         prop_b1 = prop_b) %>%
  select(-d_buffer)
  
buffer_d2_df <- buffer_df_all %>%
  filter(d_buffer == "d2") %>%
  rename(mean_b2 = mean_buffer,
         n_b2 = n_b,
         prop_b2 = prop_b) %>%
  select(-d_buffer)
  
evaluate_no_dd <- full_join(
  # Add k closest neighbours
  full_join(long_df_all_sp, k_closest_gps, by=c("no_DD_sp","group")),
  # Add buffer metrics
  # join with d2 first because no NA
  full_join(buffer_d2_df, buffer_d1_df, by=c("no_DD_sp","group"))) %>%
  # replace NA by 0 for proportions (but nor for means)
  mutate_at(c("n_10_closest","n_b1","prop_b1", "n_b2","prop_b2"), 
            ~ ifelse(is.na(.), 0, .))

saveRDS(evaluate_no_dd, "Output/04_metrics_evaluate_no_dd_ft_role")


#### Attribute a status for each metric ####

# for min_dist & mean_dist, attribute status of minimal distance group
# idem for mean_b1 & mean_b2: minimal mean distance in each buffer

# for proportion of closest species:
# take group status for which the proportion is maximal when removing null_prop

# load metrics file
evaluate_no_dd <- readRDS("Output/04_metrics_evaluate_no_dd_ft_role") %>%
  mutate(group = as.character(group))

# calculate null proportions of group representation
null_prop <- c(table(ei_no_dd$group)/nrow(ei_no_dd))

evaluate_no_dd <- left_join(evaluate_no_dd,
                            as.data.frame(null_prop) %>% rownames_to_column("group"),
                            by = "group")

evaluate_no_dd <- evaluate_no_dd %>%
  mutate(prop_10_closest = n_10_closest/10) %>%
  mutate(prop_b1_corr = prop_b1 - null_prop,
         prop_b2_corr = prop_b2 - null_prop,
         prop_10_closest_corr = prop_10_closest - null_prop) %>%
  mutate(prop_b1_corr = if_else(is.na(mean_b1), mean_b1, prop_b1_corr),
         prop_b2_corr = if_else(is.na(mean_b2), mean_b2, prop_b2_corr))

# initialize empty df to fill with predicted status
status_all_meth <- data.frame(matrix(ncol = 8, nrow = length(no_dd)))
names(status_all_meth) <- c("no_DD_sp","min_dist","mean_dist","mean_b1", 
                            "prop_b1_corr", "mean_b2", "prop_b2_corr",
                            "prop_10_closest_corr")
status_all_meth$no_DD_sp <- no_dd

# attribute status for each dd sp & each method
for(sp in no_dd){
  
  evaluate_no_dd_sp <- evaluate_no_dd %>% filter(no_DD_sp == sp)
  
  # vars for which to take minimal value
  vars_min <- c("min_dist","mean_dist","mean_b1","mean_b2")
  
  for(col in vars_min){
    # predict only if at least 1 group is rpzted in the metric
    if(sum(is.na(pull(evaluate_no_dd_sp,col)))!=3){
      status_all_meth[which(status_all_meth$no_DD_sp==sp),col] <-
        evaluate_no_dd_sp$group[which.min(pull(evaluate_no_dd_sp,col))]
    } else {
      status_all_meth[which(status_all_meth$no_DD_sp==sp),col] <- NA
    }
  }
  # vars for which to take maximal value after correction by null_prop
  vars_max <- c("prop_b1_corr", "prop_b2_corr", "prop_10_closest_corr")
  for(col2 in vars_max){
    if(sum(is.na(pull(evaluate_no_dd_sp,col2)))!=3){
      status_all_meth[which(status_all_meth$no_DD_sp==sp),col2] <-
        evaluate_no_dd_sp$group[which.max(pull(evaluate_no_dd_sp,col2))]
    } else {
      status_all_meth[which(status_all_meth$no_DD_sp==sp),col2] <- NA
    }
  }
  
}


#### Votes & Contingency tables ####


# Democracy vote
status_all_meth$IAS_T <- rowSums(status_all_meth == "IAS-T", na.rm = T)
status_all_meth$EICAT <- rowSums(status_all_meth == "EICAT", na.rm = T)

status_vote <- c("IAS_T", "EICAT")
for(i in 1:nrow(status_all_meth)){
    status_all_meth$demo_vote[i] <- 
      status_vote[which.max(status_all_meth[i,status_vote])]
}

# Absolute majority (all method vote for one status)
for(i in 1:nrow(status_all_meth)){
  if(max(status_all_meth[i,status_vote])>5){
    status_all_meth$majo_abs_vote[i] <-
      status_vote[which.max(status_all_meth[i,status_vote])]
  } else {
    status_all_meth$majo_abs_vote[i] <- "undefined"
  }
}


# Compare democracy vote with observed status
compare <- full_join(status_all_meth %>%
                       rename(Species1 = no_DD_sp),
                     corresp, by = "Species1")

table(compare$group, compare$demo_vote)

sum(diag(as.matrix(table(compare$group, compare$demo_vote))))/length(no_dd)

table(compare$group, compare$majo_abs_vote)

sum(diag(as.matrix(table(compare$group, compare$majo_abs_vote))[1:2,1:2]))/
  sum(as.matrix(table(compare$group, compare$majo_abs_vote))[1:2,1:2])

sum(as.matrix(table(compare$group, compare$majo_abs_vote))[1:2,3])/length(no_dd)


# Make a contingency table per metric
metric_names <- c("min_dist","mean_dist","mean_b1", "prop_b1_corr",
                  "mean_b2","prop_b2_corr", "prop_10_closest_corr")
ct_list <- vector(mode="list", length(metric_names))
names(ct_list) <- metric_names
  
for(col in metric_names){
  ct_list[[col]] <- table(compare$group, pull(compare, col))
}

ct_list


# Evaluate which metrics are the best predictors

# percent of well predicted values (matrix diagonal)
lapply(ct_list, function(x){
  sum(diag(as.matrix(x)))/length(no_dd)
})

# Focus on EICAT_imp --> TP, TN, FP, FN => specificity, sensibility
calculate_eicat_imp_success <- function(mat){
  TP <- mat[1,1]
  TN <- sum(mat[2,2])
  FP <- sum(mat[2,1])
  FN <- sum(mat[1,2])
  index <- c("TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN,
             "sensitivity" = TP/(TP+FN), "specificity" = TN/(TN+FP))
  return(index)
}

lapply(ct_list, function(x){
  mat <- as.matrix(x)
  index <- calculate_eicat_imp_success(mat)
  return(index)
  })

# New vote with best metrics?

# Majority vote seems to be the best one
# Remove metrics one by one and see which combination gives the best 
# prediction rate ?

best_status <- status_all_meth %>%
  select(no_DD_sp:prop_10_closest_corr)

col_to_combi <- colnames(best_status %>% select(-no_DD_sp))

#initialize
line_count = 0 # reach 120 for k in 2:7
df_best_combi <- data.frame(matrix(ncol = 17, nrow = 0))
names(df_best_combi) <- c("k","j","success_demo","TP_demo", "TN_demo", 
                          "FP_demo", "FN_demo", "sensitivity_demo", "specificity_demo",
                          "success_abs","undef_rate", "TP_abs", "TN_abs", "FP_abs", 
                          "FN_abs","sensitivity_abs", "specificity_abs")

for(k in 2:7){
  
  k_combi <- combn(col_to_combi,k)
  
  for (j in 1:ncol(k_combi)){
    line_count <- line_count+1
    
    best_status_to_test <- best_status %>% select(no_DD_sp,k_combi[,j])
    
    # Democracy vote
    best_status_to_test$IAS_T_imp <- rowSums(best_status_to_test == "IAS-T_imp", na.rm = T)
    best_status_to_test$IAS_T_no_imp <- rowSums(best_status_to_test == "IAS-T_no_imp", na.rm = T)
    best_status_to_test$EICAT_imp <- rowSums(best_status_to_test == "EICAT_imp", na.rm = T)
    best_status_to_test$EICAT_no_imp  <- rowSums(best_status_to_test == "EICAT_no_imp", na.rm = T)
    
    status_vote <- c("IAS_T_imp", "IAS_T_no_imp", "EICAT_imp", "EICAT_no_imp")
    for(i in 1:nrow(best_status_to_test)){
      best_status_to_test$demo_vote[i] <- 
        status_vote[which.max(best_status_to_test[i,status_vote])]
    }
    
    # Absolute majority (all method vote for one status)
    if(k>2){
      for(i in 1:nrow(best_status_to_test)){
      if(max(best_status_to_test[i,status_vote])>(length(k_combi[,j])-2)){
        best_status_to_test$majo_abs_vote[i] <-
          status_vote[which.max(best_status_to_test[i,status_vote])]
      } else {
        best_status_to_test$majo_abs_vote[i] <- "undefined" 
      }
        }
    } else {
      for(i in 1:nrow(best_status_to_test)){
        if(max(best_status_to_test[i,status_vote])>1){
          best_status_to_test$majo_abs_vote[i] <-
            status_vote[which.max(best_status_to_test[i,status_vote])]
        } else {
          best_status_to_test$majo_abs_vote[i] <- "undefined" }
      }
    }
    
    # Compare democracy vote with observed status
    compare_to_test <- full_join(best_status_to_test %>%
                                   rename(Species1 = no_DD_sp),
                                 corresp, by = "Species1")
    
    mat_demo <- as.matrix(table(compare_to_test$group,compare_to_test$demo_vote))
    mat_abs <- as.matrix(table(compare_to_test$group, compare_to_test$majo_abs_vote))
    
    success_rate_demo <- sum(diag(mat_demo))/length(no_dd)
    
    success_rate_abs <- sum(diag(mat_abs[1:4,1:4]))/sum(mat_abs[1:4,1:4])
    
    undef_rate_abs <- sum(mat_abs[1:4,5])/length(no_dd)
    
    index_abs <- calculate_eicat_imp_success(mat_abs)
    index_demo <- calculate_eicat_imp_success(mat_demo)
    names(index_abs) <- paste0(names(index_abs), "_abs")
    names(index_demo) <- paste0(names(index_demo), "_demo")
    
    df_best_combi <- bind_rows(df_best_combi, 
              c("k" = k , "j" = j, "success_demo" = success_rate_demo, index_demo,
                "success_abs" = success_rate_abs, "undef_rate" = undef_rate_abs,
                index_abs))
   print(line_count)
   }
  
}


plot(df_best_combi$success_abs, df_best_combi$undef_rate)
plot(df_best_combi$success_abs, df_best_combi$success_demo) 
abline(a=0, b=1)

ggplot(df_best_combi)+
  geom_point(aes(success_abs, success_demo)) +
  geom_abline(slope = 1)


df_best_combi %>% 
  filter(success_demo>0.75 & success_abs > 0.83 &
           undef_rate < 0.3)

df_best_combi %>% 
  filter(success_demo>0.75 & success_abs > 0.83)


combn(col_to_combi,2)[,4]
k=2
j=4
combn(col_to_combi,6)[,6]


ggplot(df_best_combi)+
  geom_point(aes(FP_abs, TP_abs)) +
  geom_abline(slope = 1)

ggplot(df_best_combi)+
  geom_point(aes(FP_demo, TP_demo)) +
  geom_abline(slope = 1)

