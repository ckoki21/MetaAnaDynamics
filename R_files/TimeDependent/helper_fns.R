
interpolate_missing_data <- function(Position,Time){
  Y <- zoo::zoo(cbind(Time,Position)) #see https://stackoverflow.com/questions/7188807/interpolate-na-values
  index(Y) <- Y[,1]
  Y_approx <- na.approx(Y)
  return(as.numeric(Y_approx[,2]))
}

extract_long_tracks <- function(Data,K=Inf,T0=0,max_missing=0){
  #take str giving reference to a jobset converted to csv and return all tracks beyond a certain length
  #defaults to using all the available data
  how_much_missing <- Data %>%
    dplyr::filter(Frame<=(T0+K), Frame>T0) %>% #assess only on first K frames, starting at frame T0 + 1
    dplyr::group_by(SisterPairID,SisterID) %>% #assess each track individually
    dplyr::summarise(proportionNaN = sum(is.na(Position_1))/length(Position_1)) %>%
    dplyr::group_by(SisterPairID) %>% #combine to consider pairs together
    dplyr::summarise(proportionNaN = max(proportionNaN))
  print(how_much_missing)
  print(dim(how_much_missing))
  Data <- dplyr::left_join(Data,how_much_missing) %>%
    dplyr::group_by(SisterPairID,SisterID) %>%
    dplyr::filter(proportionNaN <= max_missing,
                 Frame<=(T0+K), Frame>T0) %>%
    dplyr::ungroup()
  
  return(Data)
}

extract_long_tracks_check_missing_up_to_K <- function(Data,K=Inf,T0=0,max_missing=0){
  #take str giving reference to a jobset converted to csv and return all tracks beyond a certain length
  #defaults to using all the available data
  how_much_missing <- Data %>%
    dplyr::filter(Frame<=(T0+K), Frame>T0) %>% #assess only on first K frames, starting at frame T0 + 1
    dplyr::group_by(SisterPairID,SisterID) %>% #assess each track individually
    dplyr::summarise(proportionNaN = sum(is.na(Position_1))/length(Position_1)) %>%
    dplyr::group_by(SisterPairID) %>% #combine to consider pairs together
    dplyr::summarise(proportionNaN = max(proportionNaN))
  print(how_much_missing)
  print(dim(how_much_missing))
  Data <- dplyr::left_join(Data,how_much_missing) %>%
    dplyr::group_by(SisterPairID,SisterID) %>%
    dplyr::filter(proportionNaN <= max_missing,
                  Frame<=(T0+K), Frame>T0) %>%
#    Frame<=(T0+K), Frame>T0) %>%
  dplyr::ungroup()

return(Data)
}
reorder_sisters <- function(Data){
  #ensures that sister 1 has on average the larger x (Position_1) coordinate
  Data %>%
    dplyr::group_by(SisterPairID,SisterID) %>%
    dplyr::summarise(x = mean(Position_1,na.rm=T)) %>%
    dplyr::group_by(SisterPairID) %>%
    dplyr::summarise(need_to_switch = (first(x) - last(x))<0) %>%
    dplyr::right_join(Data) %>%
    dplyr::mutate(SisterID = as.integer(dplyr::case_when(
      need_to_switch & (SisterID==1) ~ 2,
      need_to_switch & (SisterID==2) ~ 1,
      !need_to_switch ~ as.double(SisterID),
      TRUE ~ as.double(NA)))) %>%
    dplyr::select(-need_to_switch)
}

process_jobset <- function(jobset_str,K = Inf,max_missing=0,start_from=0,plot_opt=0){
  #this function makes it easier to read tracking output of tracked kinetochores
Data <- read.csv(jobset_str,header=TRUE) %>%
  dplyr::group_by(SisterPairID,SisterID) %>%
  dplyr::mutate(Position_1=interpolate_missing_data(Position_1,Time)) %>%
  dplyr::mutate(Position_2=interpolate_missing_data(Position_2,Time)) %>%
  dplyr::mutate(Position_3=interpolate_missing_data(Position_3,Time)) %>%
    dplyr::ungroup() %>%
    extract_long_tracks(K,start_from,max_missing) %>%
    reorder_sisters()
  if (plot_opt){
    g <- ggplot(Data, aes(x=Time, y=Position_1,color=factor(SisterID))) +
      geom_line() +
      facet_wrap(.~SisterPairID) +
      theme_bw() +
      theme(legend.position = "None",
      strip.text.x = element_text(size = 8),
      axis.text.x = element_text(angle=90)) +
      labs(x="Time (s)",y="Position (um)")
    #print(g)
    ggsave(stringr::str_replace(jobset_str,".csv",
           paste("_processed_tracks_length",K,".eps",sep="")))
  }
return(Data)
}

process_jobset_keep_missing <- function(jobset_str, K= Inf,max_missing=0,start_from=0,plot_opt=0){
  #this function makes it easier to read tracking output of tracked kinetochores
  #the frame up to where we want to check the missing data
  Data <- read.csv(jobset_str,header=TRUE) %>%
    dplyr::group_by(SisterPairID,SisterID) %>%
#    dplyr::mutate(Position_1=interpolate_missing_data(Position_1,Time)) %>%
#    dplyr::mutate(Position_2=interpolate_missing_data(Position_2,Time)) %>%
#    dplyr::mutate(Position_3=interpolate_missing_data(Position_3,Time)) %>%
#    dplyr::ungroup() %>%
    extract_long_tracks_check_missing_up_to_K(K,start_from,max_missing) %>% 
    reorder_sisters()
  if (plot_opt){
    g <- ggplot(Data, aes(x=Time, y=Position_1,color=factor(SisterID))) +
      geom_line() +
      facet_wrap(.~SisterPairID) +
      theme_bw() +
      theme(legend.position = "None",
            strip.text.x = element_text(size = 8),
            axis.text.x = element_text(angle=90)) +
      labs(x="Time (s)",y="Position (um)")
    print(g)
    ggsave(stringr::str_replace(jobset_str,".csv",
                                paste("_processed_tracks_length",K,".eps",sep="")))
  }
  return(Data)
}



#
# df = data_files_interpolation[[1]]
# input_df =  df %>% annotate_anaphase_onset_for_cell(method="manual",t_ana_frame=278)


rough_anaphase_frame = function(jobset_str = NULL, df_1 = NULL, prop_missing){
  library("brms") #https://bookdown.org/content/3686/jags-brms.html
  library("mcp") #https://lindeloev.github.io/mcp/
  library(bayesplot)
  
  if(is.null(df_1)&&is.null(jobset_str)){
    stop("Give: either the dataset or the directory of data to import")
  }else{
    if(is.null(df_1)){
      df_1 <- process_jobset(jobset_str,K=Inf, start_from = 0,max_missing = prop_missing) %>%
        filter(!is.na(SisterID)) #omit unpaired KTsD))
      cell_name=str_split(jobset_str,"kittracking")[[1]][2]
      cell_name = str_split(cell_name,".ome")[[1]][1]
      cell_name = str_split(cell_name,"pair")[[1]][1]
    }
  }
  df = df_1  %>% mutate(filename = paste0("kittracking",cell_name))
  
  df = df %>%dplyr::rename(kittracking_file_str = filename) %>% mutate(kts = length(unique(SisterPairID)))
  df_m= df %>%group_by(Frame,SisterID,kittracking_file_str) %>%
    mutate(med=median(Position_1, na.rm=TRUE),
           spread=mad(Position_1,  constant = 1, na.rm=TRUE)) %>% ungroup()
  
  cluster_positions_df_m = df %>%
    group_by(Frame,SisterID,kittracking_file_str) %>%
    summarise(med=median(Position_1,na.rm=TRUE),
              spread=mad(Position_1, constant = 1, na.rm=TRUE),
              s_deviation = sd(Position_1,na.rm=TRUE)) %>% ungroup()
  
  median_positions <- cluster_positions_df_m %>%
    dplyr::select(-c("s_deviation"))%>%
    pivot_wider(names_from = c("SisterID"), values_from = c("med", "spread")) %>%
    mutate(diff_median = med_1-med_2)
  
  model = list(
    diff_median ~ Frame+ar(1),  # plateau (int_1)
    ~ rel(1) + (Frame-rel(1)) /(1-Frame-rel(1))
  )
  fit = mcp(model, median_positions)
  
  summary_fit = summary(fit)
  if(summary_fit %>%filter(name == "cp_1") %>%select(Rhat)>1.1){
    print("Anaphase time not reliable. Taking max recorded frame")
    min_laziness_frame = max(df$Frame)
  }else{
  draws <- as_draws_df(fit$mcmc_post)
  draws %>%  mutate(chain = .chain) %>% mcmc_dens_overlay(pars = vars(cp_1))
  draws %>% mutate(chain = .chain) %>%mcmc_dens_overlay(pars = vars(cp_1))
  posterior_summary(draws, robust=T,probs = c(0.025, 0.5,0.975) )["cp_1",]
  min_laziness_frame = posterior_summary(draws, robust=T,probs = c(0.025, 0.5,0.975))["cp_1","Estimate"]
  }  
  
  return(min_laziness_frame)
}

filter_cell_data <- function(jobset_str = NULL,df_1 = NULL, max_anaphase_frames =30, p_missing_final =0.2, prop_missing_initial =0.3){
  if(is.null(df_1)&&is.null(jobset_str)){
    stop("Give: either the dataset or the directory of data to import")
  }
  ana_frame_r = rough_anaphase_frame(jobset_str,df_1,prop_missing_initial)
  data_single_pair = process_jobset_keep_missing (jobset_str, K= ana_frame_r + max_anaphase_frames,max_missing=p_missing_final,start_from=0,plot_opt=0)
  return(data_single_pair)  
}


# ind_to_binary <- function(ind,nStates=5){
#  aux = rep(0,nStates)
#  aux[ind] = 1.0
#  return(aux)
# }


get_single_pair <- function(Data,id) {
  Y1 <- Data %>% dplyr::filter(SisterPairID==id,
                        SisterID==1) %>%
    dplyr::pull(Position_1)
  Y2 <- Data %>% dplyr::filter(SisterPairID==id,
                        SisterID==2) %>%
    dplyr::pull(Position_1)
  Y = cbind(Y1,Y2)
  return(Y)
}
# 
prepare_for_stan_format <- function(Data,IDs=NA){
  if (is.na(IDs)){
    #use default of all available sisters
    IDs <- unique(Data$SisterPairID)
  }
  Y_list <- purrr::map(IDs,function(id) get_single_pair(Data,id))
  return(Y_list)
}

get_start_end_of_nonmissing_data <- function(y_missing){
#convert a binary vector of whether data is missing into two integers T0 and T1
#which indicate the final missing point at the start, and the final existing point at the end
#T0>=0, T1<=K
######################
y_not_missing_ind <- which(y_missing<=0)
T0 = y_not_missing_ind[1] - 1
T1 = y_not_missing_ind[length(y_not_missing_ind)]
return(list(T0=T0,T1=T1))
}

compute_angle_theta <- function(Position_1,Position_2,Position_3){
  #computes at a single time
  #assumes Position_1 is a length 2 vector with positions for each sister
  #returns cos of the angle from kk-axis to normal of metaphase plate
  ##########
  stopifnot(length(Position_1)==2)
  if (any(is.na(c(Position_1,Position_2,Position_3)))) {return(NA)}
  #When we run get_cos_phi the data are ordered as sister 1 sister 2 sister 1 sister 2
  #for T =0 0 2 2 4 4 etc. The function diff calculates the difference
  #  x[(1+lag):n] - x[1:(n-lag)]. given that we have grouped by time it will take only
  #the difference between the same times lets say  = 0 . Hence as diff(Position_1) = X2_t0 - X1_t0
  #but we want X1_t0- X2_t0. hence -diff
  inter_kt_vec = c(-diff(Position_1),-diff(Position_2),-diff(Position_3))
  cos_theta = inter_kt_vec[1]/norm(inter_kt_vec,type="2")
  return(cos_theta)
}

#I added line 209  in Jonathan's code and grouped by Time AND SisterID because when we reorder sister 2 is on the top leading to positive cosphi
# get_cos_phi <- function(Data,id){
#   cos_phi <- Data %>%
#     dplyr::filter(SisterPairID==id) %>% group_by(Time) %>%
#     summarise(cos_phi = compute_angle_phi(Position_1,Position_2,Position_3)) 
#     pull(cos_phi)
# }


#I added line 209  in Jonathan's code and grouped by Time AND SisterID because when we reorder sister 2 is on the top leading to positive cosphi
get_cos_theta <- function(Data,id){
  cos_theta <- Data %>%
    dplyr::filter(SisterPairID==id) %>% group_by(Time,SisterID) %>%
    summarise(Position_1, Position_2, Position_3) %>% dplyr::select(Time, Position_1, Position_2, Position_3) %>%
    summarise(cos_theta = compute_angle_theta(Position_1,Position_2,Position_3)) %>%
    pull(cos_theta)
}



compute_spherical_angle_phi <- function(Position_1,Position_2,Position_3){
  #computes at a single time
  ##########
  stopifnot(length(Position_1)==2)
  if (any(is.na(c(Position_1,Position_2,Position_3)))) {return(NA)}
  #When we run get_cos_phi the data are ordered as sister 1 sister 2 sister 1 sister 2
  #for T =0 0 2 2 4 4 etc. The function diff calculates the difference
  #  x[(1+lag):n] - x[1:(n-lag)]. given that we have grouped by time it will take only
  #the difference between the same times lets say  = 0 . Hence as diff(Position_1) = X2_t0 - X1_t0
  #but we want X1_t0- X2_t0. hence -diff
  inter_kt_vec = c(-diff(Position_1),-diff(Position_2),-diff(Position_3))
  cos_spherical_phi = sign(inter_kt_vec[2])*acos(inter_kt_vec[3]/sqrt(inter_kt_vec[2]^2+inter_kt_vec[3]^2))
  return(cos_spherical_phi)
}

#I added line 209  in Jonathan's code and grouped by Time AND SisterID because when we reorder sister 2 is on the top leading to positive cosphi
get_spherical_phi <- function(Data,id){
  cos_phi <- Data %>%
    dplyr::filter(SisterPairID==id) %>% group_by(Time,SisterID) %>%
    summarise(Position_1, Position_2, Position_3) %>% dplyr::select(Time, Position_1, Position_2, Position_3) %>%
    summarise(cos_spherical_phi = compute_spherical_angle_phi(Position_1,Position_2,Position_3)) %>%
    pull(cos_spherical_phi)
}


compute_spherical_angle_phi_mean_sisters <- function(Position_1,Position_2,Position_3){
  #computes at a single time
  ##########
  stopifnot(length(Position_1)==2)
  if (any(is.na(c(Position_1,Position_2,Position_3)))) {return(NA)}
  #When we run get_cos_phi the data are ordered as sister 1 sister 2 sister 1 sister 2
  #for T =0 0 2 2 4 4 etc. The function diff calculates the difference
  #  x[(1+lag):n] - x[1:(n-lag)]. given that we have grouped by time it will take only
  #the difference between the same times lets say  = 0 . Hence as diff(Position_1) = X2_t0 - X1_t0
  #but we want X1_t0- X2_t0. hence -diff
  inter_kt_vec = c(mean(Position_1),mean(Position_2),mean(Position_3))
  cos_spherical_phi = sign(inter_kt_vec[2])*acos(inter_kt_vec[3]/sqrt(inter_kt_vec[2]^2+inter_kt_vec[3]^2))
  return(cos_spherical_phi)
}

#I added line 209  in Jonathan's code and grouped by Time AND SisterID because when we reorder sister 2 is on the top leading to positive cosphi
get_spherical_phi_mean_sisters <- function(Data,id){
  cos_phi <- Data %>%
    dplyr::filter(SisterPairID==id) %>% group_by(Time,SisterID) %>%
    summarise(Position_1, Position_2, Position_3) %>% dplyr::select(Time, Position_1, Position_2, Position_3) %>%
    summarise(cos_spherical_phi_mean = compute_spherical_angle_phi_mean_sisters(Position_1,Position_2,Position_3)) %>%
    pull(cos_spherical_phi_mean)
}


get_kk_distance <- function(Data,id){
  kk_distance <- Data %>%
    dplyr::filter(SisterPairID==id) %>% group_by(Time,SisterID) %>%
    summarise(Position_1, Position_2, Position_3) %>% dplyr::select(Time, Position_1, Position_2, Position_3) %>%
    summarise(kk_distance =ifelse( any(is.na(c(Position_1,Position_2,Position_3))),NA,norm(c(diff(Position_1),diff(Position_2),diff(Position_3)), type = "2"))) %>%
    pull(kk_distance)
}

how_much_missing <- function(jobset_str,K=Inf,T0=0, max_mis_prop = 0.25){
  #take str giving reference to a jobset converted to csv and return all tracks beyond a certain length
  #defaults to using all the available data
  Data <- read.csv(jobset_str,header=TRUE)
  where_missing <-  Data %>%
    dplyr::filter(Frame<=(T0+K), Frame>T0) %>% #assess only on first K frames, starting at frame T0 + 1
    dplyr::group_by(SisterPairID,SisterID) %>% #assess each track individually
    dplyr::summarise(where_na = which(is.na(Position_1) ==TRUE, arr.ind = TRUE))  
  
  how_much_missing <- Data %>%
    dplyr::filter(Frame<=(T0+K), Frame>T0) %>% #assess only on first K frames, starting at frame T0 + 1
    dplyr::group_by(SisterPairID,SisterID) %>% #assess each track individually
    dplyr::summarise(proportionNaN = sum(is.na(Position_1))/length(Position_1)) %>%
    dplyr::group_by(SisterPairID) %>% #combine to consider pairs together
    dplyr::summarise(proportionNaN = max(proportionNaN))
  
  
  how_much_missing_data_single_pair = how_much_missing %>% filter(proportionNaN <= max_mis_prop)
  ds_missing_vector = where_missing %>% group_by(SisterPairID, SisterID) %>% select(SisterPairID,SisterID, where_na) %>% summarise(which_missing = list(where_na))
  missing_data_single_pair = how_much_missing_data_single_pair %>% left_join(ds_missing_vector)
  
  
  return(missing_data_single_pair)
}


# for the cells with less than q = 25% of missing data (or any other predefiined percentage, maybe q = 5%) 
#what i want is to know 
#1. how many missing obesrvations
#2. the length of consecutive missing numbers and how many times
#3. where are the missing numbers? middle, start, end ?
#note by start we mean occures at the first 10% of the sequence
#by end we mean after the anaphase time! So we need to run the missing_data
#functions after the change point model when we keep the metaphase period

get_statistics_missing_data  = function(jobset_str, full_cell_name,max_mis_prop, dt = 2.05,T_0 = 0,mean_anaphase_path = NA){
  #Mean anaphase_time is the dataframe with the anaphase times as derived by the change-point model
  #we will only use the Average(mean_anaphase_time$Mean_t_ana)
  if(is.na(mean_anaphase_path)){
    K = Inf
    print(K) 
  }else{
    mean_t_ana_df = read.csv(mean_anaphase_path, header =TRUE)
    K = (mean(mean_t_ana_df[,1])-30)/dt
    print(K)
    #Equivalent
    #    K = colMeans(mean_t_ana_df)[1]
    #    K = colMeans(mean_t_ana_df)["Mean_t_ana"]
  }
  #CAREFUL: T0 = 0 ALWAYS 
  df_missing = how_much_missing(jobset_str, K  , T0=T_0,max_mis_prop)
  #empty df
  df_vector_missing = df_missing %>% mutate(new = as.character(which_missing), .keep = "none")
  
  Missing_statistics<- data.frame(matrix(ncol = 16, nrow = 0))
  #provide column names
  colnames(Missing_statistics) <- c('Cell_name','SisterPairID','SisterID',"Missing values",'proportion','total', 'length_1', 'length_2', 'length_3', 'length_4', 
                                    'length_5', 'length_6', 'length_7', 'length_8', 'length_9', 'length_10+')
  
  # length of missing numbers array (convert this from the list for each sister and each kt-pair)
  for(j in 1:dim(df_missing)[1]){
    Missing_statistics[j,"Cell_name"] = full_cell_name 
    Missing_statistics[j,"SisterPairID"] = df_missing$SisterPairID[j]
    Missing_statistics[j,"SisterID"] = df_missing$SisterID[j]
    Missing_statistics[j,"proportion"] = df_missing$proportionNaN[j]
    Missing_statistics[j,"Missing values"] = df_vector_missing[j,1]
    
    
    ls_which = df_missing$which_missing[[j]]
    if(!is.null(ls_which)){
      Missing_statistics[j,"first_missing"] = min(ls_which)
      Missing_statistics[j,"total"] = length(ls_which)
      #if values ==1 then there are consecutive numbers . if values>1 then single missing data

     ls_which_aug = c(-10, ls_which, 1000)
      A = rle(diff(ls_which_aug))
      A$values[A$values != 1] = 0
      Single_miss = A$values
      Missing_statistics[j,"length 1"] = sum((rle(Single_miss)$lengths[rle(Single_miss)$values ==0 ]-1))
      for(i in 1:9){ 
        if(i<9){
          mis_lengths = rle(diff(ls_which_aug))$lengths[rle(diff(ls_which_aug))$values==1]
          Missing_statistics[j, 7+i] = sum(1*(mis_lengths == i))
        }else{
          Missing_statistics[j, 7+i] = sum(1*(mis_lengths >= i))}
      }
    }
  }
  return(Missing_statistics)
}
