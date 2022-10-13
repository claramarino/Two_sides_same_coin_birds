# find species in buffer and calculate group proportions


find_sp_in_buffer_gp <- function(df_dist, dd, d_buffer){
  
  # empty list that will contain sp in buffer for each dd
  sp_in_buffer <- vector(mode="list", length = length(dd))
  names(sp_in_buffer) <- dd
  
  # find all sp in buffer for each dd
  for(i in 1: length(dd)){
    col1 <- as.data.frame(t(df_dist[i,])) %>% rownames_to_column("sp")
    names(col1) <- c("Species1","dist")
    sp_in_buffer[[i]] <- col1 %>% filter(dist<d_buffer)
  }
  
  # calculate proportion and mean for each buffer type
  sp_in_buffer_gp <- lapply(sp_in_buffer, function(x){
    left_join(x, corresp, by="Species1") %>%
      group_by(group2) %>%
      summarise(mean_buffer = mean(dist),
                n_b = n()) %>%
      mutate(prop_b = n_b / sum(n_b))})
  
  return(sp_in_buffer_gp)
}

# adapt function for role in biol inv (alien vs IAS-T)

find_sp_in_buffer_role <- function(df_dist, dd, d_buffer){
  
  # empty list that will contain sp in buffer for each dd
  sp_in_buffer <- vector(mode="list", length = length(dd))
  names(sp_in_buffer) <- dd
  
  # find all sp in buffer for each dd
  for(i in 1: length(dd)){
    col1 <- as.data.frame(t(df_dist[i,])) %>% rownames_to_column("sp")
    names(col1) <- c("Species1","dist")
    sp_in_buffer[[i]] <- col1 %>% filter(dist<d_buffer)
  }
  
  # calculate proportion and mean for each buffer type
  sp_in_buffer_gp <- lapply(sp_in_buffer, function(x){
    left_join(x, corresp, by="Species1") %>%
      group_by(group) %>%
      summarise(mean_buffer = mean(dist),
                n_b = n()) %>%
      mutate(prop_b = n_b / sum(n_b))})
  
  return(sp_in_buffer_gp)
}