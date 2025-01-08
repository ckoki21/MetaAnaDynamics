TimeDependentModelFunction3 <- function(jobset_str, warm_up_iter = 8000, total_iter = 12000, 
                                       warm_up_iter_anaphase = 3000, total_iter_anaphase = 5000,
                                       dt0 = 2.05, K=Inf,p_missing_final = 0.2, prop_missing_initial =0.3){


  proportional_missing = 1
  whole_name=str_split(jobset_str,"per_frame/kittracking0")
  #wl = str_split(whole_name[[1]][2],"Untreated")[[1]][1]
  split_name=str_split(whole_name[[1]][2],".ome.csv")
  cell_name_old=paste(split_name[[1]][1], sep="")
  split_name= str_split(cell_name_old,"-kitjobset")
  cell_name_old=paste(split_name[[1]][1], sep="")
  
  whole_name2 = str_split(jobset_str,"per_frame/")[[1]][2]
  whole_name2 = str_split(whole_name2,".ome.csv")[[1]][1]
  cell_name = whole_name2
  print(whole_name2)
  
p_missing_interpolate = 0.25

data_single_pair = filter_cell_data(jobset_str, max_anaphase_frames = 45, p_missing_final = p_missing_final, prop_missing = prop_missing_initial)

# Identify mephase window
pairs_to_include = unique(data_single_pair$SisterPairID)
data_single_pair_interpolate <- process_jobset(jobset_str,K=K,max_missing=1) %>%
  filter(!is.na(SisterID)) 
#omit unpaired KTs 
data_single_pair_interpolate <- data_single_pair_interpolate %>%filter(SisterPairID %in% (pairs_to_include)) #keep sisters that have pass the criterion of filter_cell_data

pairIDs <- unique(data_single_pair_interpolate$SisterPairID)
K=max(data_single_pair_interpolate$Frame)
# #Remember this is 2 times the frames because we have data for sister 1 and 2
#2K=dim(data_single_pair %>% filter(SisterPairID==pairIDs[2]))[1] 
nTracks = length(pairIDs)
y_interpolate= prepare_for_stan_format(data_single_pair_interpolate)
y_interpolate_missing = purrr::map(y_interpolate, is.na) %>% purrr::map(function(x) x[,1] | x[,2])
for (i in 1:length(y_interpolate)) {
  y_interpolate[[i]][y_interpolate_missing[[i]],] <- 0
} #replace missing values with zero for stan compatibility (stan does not like NA values)
y_interpolate_missing <- purrr::map(y_interpolate_missing,as.integer)
min_tana<-1000
sister_t_ana<-0
anaphase_times<-matrix(0,2,length(pairIDs))

j<-1
changept_estimate_m <-
  stan_model(file='~/Documents/GitHub/MetaAnaDynamics/STAN_files/TimeDependent/anaphase_changepoint.stan', model_name = "anaphase_changepoint_model")

for(i in pairIDs){
  stan_input = list(dt=dt0, T=K,
                    nTracks = 1,
                    y = y_interpolate[[which(pairIDs==i)]],
                    y_missing = y_interpolate_missing[[which(pairIDs==i)]]
  )
  # file_name <- 'anaphase_changepoint.stan'
  # changept_estimate_m <- stan_model(file_name)
  changept_estimate <- sampling(changept_estimate_m,
                                data=stan_input,
                                seed = 42,
                                chains =4,
                                warmup = warm_up_iter_anaphase,
                                iter = total_iter_anaphase,
                                control=list(adapt_delta=0.95,max_treedepth=12))
  posterior <- as.array(changept_estimate)
  ch_es<-as.data.frame(changept_estimate)
  min2tana<-mean(ch_es$t_ana)
  anaphase_times[,j]<-c(min2tana,i)
  print(min2tana)
  print(i)
  if(j==1){
    min_tana = min2tana
    df_t_ana<- as.data.frame(ch_es$t_ana)
    names(df_t_ana) = c(i)
  }else{
    min_tana = min_tana+min2tana
    df_t_ana <- cbind(df_t_ana, c(ch_es$t_ana))
    names(df_t_ana)[j] = c(i)
  }
  j<-j+1
}

min_tana = min_tana/j   #This is the mean anaphase time from all sisters (or we could take the median. We have noticed that there are some sisters
#with outliers anaphase times and hence minimum is not the right choice. )

#Another way to go is to replace the outlier with the median time and then find the minimum anaphase time. 
Mean_anaphase = data.frame("Mean_t_ana"= colMeans(df_t_ana),"SisterPairID" = as.numeric(names(df_t_ana))) #Dataframe with mean anaphase times for all sister pairs.
#Remove outliers. I define them as the values which are < than (median(anaphase_times)/2)
Mean_anaphase = Mean_anaphase %>% 
  select(Mean_t_ana)%>%
  mutate(Mean_anaphase, Mean_t_ana = ifelse((Mean_t_ana< median(Mean_anaphase$Mean_t_ana)-100), median(Mean_anaphase$Mean_t_ana)-30,Mean_t_ana ))

#write_xlsx(Mean_anaphase,here::here(paste0("Mean_anaphase_times_",whole_name2,"_all.xlsx")))
write.csv(Mean_anaphase,here::here(paste0("Mean_anaphase_times_",whole_name2,".csv")),row.names=FALSE)

#Maybe more accurate to go with median instead of minimum to reduce the effect of a possible outlier/false change point.

min_tana_new =  median(Mean_anaphase$Mean_t_ana)

#!!!!!!!!!!these data are the ones with missing data, without interpolation!!!!!!!!# 
data_single_pair <- data_single_pair[(data_single_pair$Time<floor(min_tana_new-30))==TRUE,]
rm(posterior)



pairIDs <- unique(data_single_pair$SisterPairID)
y = prepare_for_stan_format(data_single_pair)
y_missing = purrr::map(y, is.na) %>% purrr::map(function(x) x[,1] | x[,2])
for (i in 1:length(y)) {
  y[[i]][y_missing[[i]],] <- 10^9
} #replace missing values high values--- see also in .stan file to treat missing (stan does not like NA values)
y_missing <- purrr::map(y_missing,as.integer)
K=max(data_single_pair$Frame)
########################Bridge Sampling############################################
Marginal_likelihoods_bridge <- data.frame(SisterPairID=integer(length(pairIDs)),
                                          constant=double(length(pairIDs)), 
                                          alpha=double(length(pairIDs)), 
                                          kappa=double(length(pairIDs)),
                                          v_minus=double(length(pairIDs)), 
                                          v_plus=double(length(pairIDs)),
                                          v_minus_plus=double(length(pairIDs)),
                                          tau=double(length(pairIDs)),
                                          v_minus_tau=double(length(pairIDs)),
                                          v_plus_tau=double(length(pairIDs)),
                                          v_minus_plus_tau=double(length(pairIDs)))

constant_m <- stan_model('~/Documents/GitHub/MetaAnaDynamics/STAN_files/TimeDependent/metaphase_symmetric_full_likelihood_missing.stan', model_name = "constant_model")
alpha_m <- stan_model('~/Documents/GitHub/MetaAnaDynamics/STAN_files/TimeDependent/time_dependent_alpha.stan', model_name = "alpha_model")
kappa_m <- stan_model('~/Documents/GitHub/MetaAnaDynamics/STAN_files/TimeDependent/time_dependent_kappa.stan', model_name = "kappa_model")
v_minus_m <- stan_model('~/Documents/GitHub/MetaAnaDynamics/STAN_files/TimeDependent/time_dependent_v_minus.stan', model_name = "v_minus_model")
v_plus_m <- stan_model('~/Documents/GitHub/MetaAnaDynamics/STAN_files/TimeDependent/time_dependent_v_plus.stan', model_name = "v_plus_model")
v_minus_plus_m <- stan_model('~/Documents/GitHub/MetaAnaDynamics/STAN_files/TimeDependent/time_dependent_v_both.stan', model_name = "v_minus_plus_model")
tau_m <- stan_model('~/Documents/GitHub/MetaAnaDynamics/STAN_files/TimeDependent/time_dependent_tau.stan', model_name = "tau_model")
v_minus_tau_m <- stan_model('~/Documents/GitHub/MetaAnaDynamics/STAN_files/TimeDependent/time_dependent_v_minus_tau.stan', model_name = "v_minus_tau_model")
v_plus_tau_m <- stan_model('~/Documents/GitHub/MetaAnaDynamics/STAN_files/TimeDependent/time_dependent_v_plus_tau.stan', model_name = "v_plus_tau_model")
v_minus_plus_tau_m <- stan_model('~/Documents/GitHub/MetaAnaDynamics/STAN_files/TimeDependent/time_dependent_v_both_tau.stan', model_name = "v_minus_plus_tau_model")

j=1
for(i in pairIDs){
  Marginal_likelihoods_bridge[j,1]=i
  trajectory_index <- which(pairIDs==i)
  start_end <- purrr::map_df(y_missing,get_start_end_of_nonmissing_data,.id="traj_index") %>%
    mutate(traj_index=as.integer(traj_index)) #use this info to only fit to the existing tracks
  cos_theta <- get_cos_theta(data_single_pair,pairIDs[trajectory_index])
  cos_theta[is.na(cos_theta)] <- 0 #replace NAs as stan does not deal with these, should not come into the fit due to T0 and T1
  T0 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T0);
  T1 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T1);
  Frames=K
  sigma0 = c(0.0,0.5,0.5,0.0)
  stan_input = list(dt=dt0, Frames=K, nStates = 4,
                    y = y[[trajectory_index]],
                    y_missing = y_missing[[trajectory_index]],
                    sigma0 =sigma0, #initial state probabilities
                    T0 =T0, #shall i say T0=T0 based on lines 205,206?
                    T1 = T1,cos_theta = cos_theta )
  
  estimate_constant<- sampling(constant_m,
                             data=stan_input,
                             seed = 42,
                             chains = 4,
                             warmup = warm_up_iter,
                             iter = total_iter,
                             init_r=0.1,
                             pars = c("f","eta","xi","auxStates","P","aux"),
                             include=FALSE, #avoid saving the params listed above
                             control=list(adapt_delta=0.95,
                                          max_treedepth=12))
  name2save <- paste("Constant_", whole_name2,'_',"pair_",i,".rds",sep="")       
  saveRDS(estimate_constant,name2save)       
  bridge1 = bridge_sampler(estimate_constant)
  print(bridge1$logml)
  #bridge1=bridge_sampler(estimate_asym1,method = "warp3")
  Marginal_likelihoods_bridge[j,2] = bridge1$logml
  
 ### # time dependent alpha # ####
  estimate_alpha <- sampling(alpha_m,data=stan_input,
                             seed = 42,
                             chains = 4,
                             warmup = warm_up_iter,
                             iter = total_iter,init_r=0.1,
                             pars = c("eta","f","xi","auxStates","P","aux"),
                             include=FALSE, #avoid saving the params listed above
                             control=list(adapt_delta=0.95,max_treedepth=12))
  name2save <- paste("Temporal_alpha_", whole_name2,'_',"pair_",i,".rds",sep="")       
  saveRDS(estimate_alpha,name2save)  
  bridge2 = bridge_sampler(estimate_alpha)
  print(bridge2$logml)
  Marginal_likelihoods_bridge[j,3] = bridge2$logml
  

###  # time dependent kappa # ###
  estimate_kappa<- sampling(kappa_m,
                             data=stan_input,
                             seed = 42,
                             chains = 4,
                             warmup = warm_up_iter,
                             iter = total_iter,init_r=0.1,
                             pars = c("eta","f","xi","auxStates","P","aux"),
                             include=FALSE, #avoid saving the params listed above
                             control=list(adapt_delta=0.95,
                                          max_treedepth=12))
  name2save <- paste("Temporal_kappa_", whole_name2,'_',"pair_",i,".rds",sep="")       
  saveRDS(estimate_kappa,name2save)  
  bridge3 = bridge_sampler(estimate_kappa)
  print(bridge3$logml)  
  Marginal_likelihoods_bridge[j,4] = bridge3$logml
  
  # time dependent v_minus
  estimate_vminus <- sampling(v_minus_m,
                             data=stan_input,
                             seed = 42,
                             chains = 4,
                             warmup = warm_up_iter,
                             iter = total_iter,init_r=0.1,
                             pars = c("eta","f","xi","auxStates","P","aux"),
                             include=FALSE, #avoid saving the params listed above
                             control=list(adapt_delta=0.95,max_treedepth=12))
  name2save <- paste("Temporal_vminus_", whole_name2,'_',"pair_",i,".rds",sep="")       
  saveRDS(estimate_vminus,name2save)
  bridge4 = bridge_sampler(estimate_vminus)
  print(bridge4$logml)
  
  Marginal_likelihoods_bridge[j,5] = bridge4$logml
  
  # time dependent v_plus
  estimate_vplus <- sampling(v_plus_m,data=stan_input,
                        seed = 42,
                        chains = 4,
                        warmup = warm_up_iter,
                        iter = total_iter,init_r=0.1,
                        pars = c("eta","f","xi","auxStates","P","aux"),
                        include=FALSE, #avoid saving the params listed above
                        control=list(adapt_delta=0.95,max_treedepth=12))
  name2save <- paste("Temporal_vplus_", whole_name2,'_',"pair_",i,".rds",sep="")       
  saveRDS(estimate_vplus,name2save)  
  bridge5 = bridge_sampler(estimate_vplus)
  print(bridge5$logml)
  Marginal_likelihoods_bridge[j,6] = bridge5$logml 
  
  # time dependent v_minus_plus
  estimate_vminus_vplus <- sampling(v_minus_plus_m,data=stan_input,
                        seed = 42,
                        chains = 4,
                        warmup = warm_up_iter,
                        iter = total_iter,init_r=0.1,
                        pars = c("eta","f","xi","auxStates","P","aux"),
                        include=FALSE, #avoid saving the params listed above
                        control=list(adapt_delta=0.95,max_treedepth=12))
  name2save <- paste("Temporal_vminus_vplus_", whole_name2,'_',"pair_",i,".rds",sep="")       
  saveRDS(estimate_vminus_vplus,name2save)  
  bridge6 = bridge_sampler(estimate_vminus_vplus)
  print(bridge6$logml)
  Marginal_likelihoods_bridge[j,7] = bridge6$logml
  
  # time dependent tau
  estimate_tau <- sampling(tau_m,data=stan_input,
                        seed = 42,
                        chains = 4,
                        warmup = warm_up_iter,
                        iter = total_iter,init_r=0.1,
                        pars = c("eta","f","xi","auxStates","P","aux"),
                        include=FALSE, #avoid saving the params listed above
                        control=list(adapt_delta=0.95,max_treedepth=12))
  name2save <- paste("Temporal_tau_", whole_name2,'_',"pair_",i,".rds",sep="")       
  saveRDS(estimate_tau,name2save)  
  bridge7 = bridge_sampler(estimate_tau)
  print(bridge7$logml)
  Marginal_likelihoods_bridge[j,8] = bridge7$logml 
  
  # time dependent vminus and tau
  estimate_vminus_tau<- sampling(v_minus_tau_m,
                             data=stan_input,
                             seed = 42,
                             chains = 4,
                             warmup = warm_up_iter,
                             iter = total_iter,
                             init_r=0.1,
                             pars = c("f","eta","xi","auxStates","P","aux"),
                             include=FALSE, #avoid saving the params listed above
                             control=list(adapt_delta=0.95,
                                          max_treedepth=12))
  name2save <- paste("Temporal_vminus_tau_", whole_name2,'_',"pair_",i,".rds",sep="")       
  saveRDS(estimate_vminus_tau,name2save)  
  bridge8 = bridge_sampler(estimate_vminus_tau)
  print(bridge8$logml)
  Marginal_likelihoods_bridge[j,9] = bridge8$logml
  
  # time dependent vplus and tau
  estimate_vplus_tau <- sampling(v_plus_tau_m,data=stan_input,
                             seed = 42,
                             chains = 4,
                             warmup = warm_up_iter,
                             iter = total_iter,init_r=0.1,
                             pars = c("eta","f","xi","auxStates","P","aux"),
                             include=FALSE, #avoid saving the params listed above
                             control=list(adapt_delta=0.95,max_treedepth=12))
  name2save <- paste("Temporal_vplus_tau_", whole_name2,'_',"pair_",i,".rds",sep="")       
  saveRDS(estimate_vplus_tau,name2save)  
  bridge9  = bridge_sampler(estimate_vplus_tau)
  print(bridge9$logml)
  Marginal_likelihoods_bridge[j,10]=bridge9$logml
  

  # v_minus, v_plus and tau
  estimate_vminus_vplus_tau <- sampling(v_minus_plus_tau_m,
                             data=stan_input,
                             seed = 42,
                             chains = 4,
                             warmup = warm_up_iter,
                             iter = total_iter,init_r=0.1,
                             pars = c("eta","f","xi","auxStates","P","aux"),
                             include=FALSE, #avoid saving the params listed above
                             control=list(adapt_delta=0.95,
                                          max_treedepth=12))
  name2save <- paste("Temporal_vminus_vplus_tau_", whole_name2,'_',"pair_",i,".rds",sep="")       
  saveRDS(estimate_vminus_vplus_tau,name2save)     
  bridge10=bridge_sampler(estimate_vminus_vplus_tau)
  print(bridge10$logml)
  
  Marginal_likelihoods_bridge[j,11]=bridge10$logml
  
  j = j+1
  print(j)
  write_xlsx(Marginal_likelihoods_bridge,paste0("Marginal_likelihoods_bridge_",cell_name,".xlsx"))
  write.csv( Marginal_likelihoods_bridge,here::here(paste0("Marginal_likelihoods_bridge_",whole_name2,".csv")), row.names = FALSE)
}

write_xlsx(Marginal_likelihoods_bridge,paste0("Marginal_likelihoods_bridge_",cell_name,".xlsx"))
  write.csv( Marginal_likelihoods_bridge,here::here(paste0("Marginal_likelihoods_bridge_",whole_name2,".csv")), row.names = FALSE)
}