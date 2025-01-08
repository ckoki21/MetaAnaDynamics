get_df_sisters_more =  function(jobset_str, Mean_anaphase=NULL, max_anaphase_frames = 45, p_missing_final = 0.2, prop_missing_initial = 0.3){
   library(caTools)
  # Finds the relevant folder/csv with the csv file containing the data for a particular cell
  # Creates the data_single_pair (removing the trajectories (sister pair kts) with high number of NaNs) 
  # Takes the dataframe with the timings for the anaphase onset (already calculated)
  # Sister distance (X1 - X2) for each time point (or frame if you like)
  # Calculates the radius for each time point (sqrt(y^2+z^2))
  # Calculates the twist angle for each time point
  # Joints the anaphase onset for each sister pair.
  #Mean_anaphase is a data frame!!!
  
  K=Inf
  #data_single_pair <- process_jobset(jobset_str,K=K,max_missing=0.20) %>%
  #  filter(!is.na(SisterID)) #omit unpaired KTs
  data_single_pair <-
    filter_cell_data(jobset_str, max_anaphase_frames = 45, p_missing_final = p_missing_final, prop_missing_initial = prop_missing_initial)%>%
    filter(!is.na(SisterID))    
  
  
  data_1 = data_single_pair %>% filter(SisterID == 1) 
  data_2 = data_single_pair %>% filter(SisterID == 2) 
  
  
  P1 = data_1 %>% dplyr::select(-c(movieID, proportionNaN, Amplitude_1, Amplitude_2, Amplitude_3))%>% dplyr::rename(X1 = Position_1,Y1 = Position_2, Z1 = Position_3)
  P2 = data_2 %>% dplyr::select(-c(movieID, proportionNaN, Amplitude_1, Amplitude_2, Amplitude_3)) %>% dplyr::rename(X2 = Position_1,Y2 = Position_2, Z2 = Position_3)
  
  
  df_sisters_more = cbind(P1[c("SisterPairID", "Frame", "Time", "X1", "Y1", "Z1")], P2[c("X2", "Y2", "Z2")]) 
  df_sisters_more = df_sisters_more %>% dplyr::select(X1,X2)%>% mutate(df_sisters_more, Mean_X_sisters = (X1+X2)/2)
  df_sisters_more = df_sisters_more %>%mutate(Radius_mean_sisters = sqrt(((Y1+Y2)/2)^2 + ((Z1+Z2)/2)^2), .keep = "all")
  df_sisters_more = df_sisters_more %>%mutate(Radius_sister_1 = sqrt((Y1)^2 + (Z1)^2), .keep = "all")
  df_sisters_more = df_sisters_more %>%mutate(Radius_sister_2 = sqrt((Y2)^2 + (Z2)^2), .keep = "all")
  
  
  pairIDs <- unique(data_single_pair$SisterPairID) 
  cos_theta_all <- c()
  phi_all <- c()
  phi_mean_all <- c()
  kk_distance_all <- c()
  for( i in pairIDs){
    trajectory_index <- which(pairIDs==i)
    cos_theta <- get_cos_theta(data_single_pair,pairIDs[trajectory_index])
    spherical_phi <- get_spherical_phi (data_single_pair,pairIDs[trajectory_index])
    spherical_phi_mean_sisters <- get_spherical_phi_mean_sisters(data_single_pair,pairIDs[trajectory_index])
    kk_distance <- get_kk_distance(data_single_pair,pairIDs[trajectory_index])
    
    kk_distance_all =  append(kk_distance_all, kk_distance)
    cos_theta_all =  append(cos_theta_all, cos_theta)
    phi_all = append(phi_all, spherical_phi)
    phi_mean_all = append(phi_mean_all, spherical_phi_mean_sisters)
  }   
  
  df_sisters_more$Cos_angle = cos_theta_all 
  df_sisters_more$Spherical_phi_rad = phi_all
  df_sisters_more$Spherical_phi_mean_sisters_rad = phi_mean_all
  df_sisters_more$KK_distance = kk_distance_all
  
  window_size = 20
  df_sisters_more = df_sisters_more %>% 
    mutate(amplitude_x1 = 0.5*(caTools::runmax(X1,window_size) - 
                                 caTools::runmin(X1,window_size)))
  A = df_sisters_more %>% group_by(SisterPairID) %>% summarise(#period = get_period(Position_1,dt),
      amplitude_median_x1 = median(amplitude_x1,na.rm=TRUE))
  df_sisters_more = left_join(df_sisters_more, A)
  
  df_sisters_more = df_sisters_more %>% 
    mutate(amplitude_x2 = 0.5*(caTools::runmax(X2,window_size) - 
                                 caTools::runmin(X2,window_size)), .keep = "all")
    B=  df_sisters_more %>% group_by(SisterPairID) %>% summarise(#period = get_period(Position_1,dt),
      amplitude_median_x2 = median(amplitude_x2,na.rm=TRUE))
  
    df_sisters_more = left_join(df_sisters_more, B)
  #FOR Twist angle in degrees see here:: Kinetochore life histories reveal an Aurora-B-dependent error correction mechanism in anaphase (Sen et al.2021
  #https://doi.org/10.1016/j.devcel.2021.10.007)
  
    # df_sisters_more = df_sisters_more %>% select(Cos_angle) %>% mutate(df_sisters_more, Twist_angle =  (180/pi)*acos(Cos_angle))
    # df_sisters_more$Median_Twist = (180/pi)*median((df_sisters_more %>% filter(Cos_angle!=0))$Cos_angle)
  
   df_sisters_more = df_sisters_more %>% left_join(Mean_anaphase, by = "SisterPairID")
  
  df_sisters_more = df_sisters_more %>% relocate(c("SisterPairID", "Frame", "Time", "X1", "X2", "Y1", "Y2", "Z1", "Z2", "Mean_X_sisters", "Radius_mean_sisters", "Cos_angle", "Spherical_phi_rad", "KK_distance"))
  return(df_sisters_more)
}