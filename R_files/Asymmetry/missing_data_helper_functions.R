

single_gap_missing_data_run <- function(data_full, prop_missing, burn_in = 5000,iterations = 15000,dt = 2.05){
  
  d_length = dim(data_full)[1]
  n_gaps = d_length*prop_missing
  yA_to_change = data_full
  
  for(i in 1:(n_gaps)){ 
    k = rdunif(1,(1+(i-1)*floor(d_length/n_gaps)), i*floor(d_length/n_gaps) )
    yA_to_change[k,] =NA    
  }
  
  data_single_pair= data.frame(SisterPairID=double(2*nFrames),Position_1=double(2*nFrames),
                               Frame=integer(2*nFrames),Time=integer(2*nFrames),
                               SisterID=integer(2*nFrames))
  data_single_pair[,1:5]=c(matrix(1,nrow=2*nFrames,ncol=1),cbind(c(yA_to_change[,1],yA_to_change[,2])),rep(c(1:nFrames),times=2),
                           rep(seq(from=0,to=((nFrames-1)*dt),by=dt),times=2),
                           cbind(c(rep(1,times=nFrames),rep(2,times=nFrames))))
  pairIDs <- unique(data_single_pair$SisterPairID)
  y = prepare_for_stan_format(data_single_pair)
  y_missing = purrr::map(y, is.na) %>% purrr::map(function(x) x[,1] | x[,2])
  for (i in 1:length(y)) {
    y[[i]][y_missing[[i]],] <- 10^3
  }
  y_missing <- purrr::map(y_missing,as.integer)
  K=max(data_single_pair$Frame)
  T0=0
  T1=nFrames
  y_missing=rep(0,times=nFrames)
  sigma_0 = c(0.0,0.5,0.5,0.0)
  stan_input = list(dt=dt, Frames=K, nStates = 4,
                    y = y[[1]],
                    y_missing = y_missing,
                    sigma0 =sigma_0, #assume each state equally likely initially
                    T0 =T0, #shall i say T0=T0 based on lines 205,206?
                    T1 = T1,cos_phi = cos_phi)
  
  ##Asymmetric metaphase model
  stan_file= here::here('src/stan_reparametrise_missing/metaphase_asymmetric_full_likelihood_reparametrised_missing.stan')
  m <- stan_model(stan_file) #slightly different way to call a stan model, this compiles the model and the next line does the sampling
  estimate_asym <- sampling(m,
                            data=stan_input,
                            seed = 25,
                            chains = 4,
                            warmup = burn_in,
                            iter = iterations,
                            pars = c("eta","xi","f","auxStates","P","aux"),
                            include=FALSE, #avoid saving the params listed above
                            control=list(adapt_delta=0.95,
                                         max_treedepth=12))
  bridge1<-bridge_sampler(estimate_asym)
  print(bridge1$logml)
  estimate_single_missing_summary = summary(estimate_asym) 
  post_asym_1=as.matrix( estimate_asym)
  states_pos=which(1*(colnames(post_asym_1)=="sigma_sim[1]")==1)
  nstates= 4
  smo_prob_asymmetric_1=matrix(NA,nrow=nFrames,ncol=nstates)
  for(frames in 1:nFrames){
    for(st in 1:nstates){
      smo_prob_asymmetric_1[frames,st]=sum(1*(post_asym_1[,states_pos+frames-1]==st))/(dim(post_asym_1)[1])
    }
  }
  Mp_states=max.col(as.data.frame(smo_prob_asymmetric_1))
  #  par(mfrow=c(2,1))
  #  plot(h_states,col="blue",type="l",pch=3)
  #  plot(Mp_states,col="red",type="l",pch=10)
  
  #  par(mfrow=c(1,1))
  #  plot(h_states,ylim=c(-0.5,4.5),ylab="")
  #  par(new=TRUE)
  #  plot((post_asym_1[3000,c(states_pos:(states_pos+nFrames-1))]-0.2),col="red",ylim=c(-0.5,4.5),ylab="")
  #  par(new=TRUE)
  #  plot((post_asym_1[2900,c(states_pos:(states_pos+nFrames-1))]-0.4),col="blue",ylim=c(-0.5,4.5),ylab="")
  
  return(list(estimate_asym, bridge1))
  ############################ single Missing points #####################################################################
  ##################################################################################################################
}


double_gap_missing_data_run <- function(data_full, prop_missing, burn_in = 5000,iterations = 15000,dt = 2.05){
  
  d_length = dim(data_full)[1]
  n_gaps = floor(d_length*prop_missing/2)
  yA_to_change = data_full
  i=1
  while(i<= n_gaps){ 
    k = rdunif(1,(1+(i-1)*floor(d_length/n_gaps)), i*floor(d_length/n_gaps) )
    print(k)
    if(k<=(d_length-1)){ 
      yA_to_change[k:(k+1),] =NA  
      i=i+1
      } 
    }
  
  data_single_pair= data.frame(SisterPairID=double(2*nFrames),Position_1=double(2*nFrames),
                               Frame=integer(2*nFrames),Time=integer(2*nFrames),
                               SisterID=integer(2*nFrames))
  data_single_pair[,1:5]=c(matrix(1,nrow=2*nFrames,ncol=1),cbind(c(yA_to_change[,1],yA_to_change[,2])),rep(c(1:nFrames),times=2),
                           rep(seq(from=0,to=((nFrames-1)*dt),by=dt),times=2),
                           cbind(c(rep(1,times=nFrames),rep(2,times=nFrames))))
  pairIDs <- unique(data_single_pair$SisterPairID)
  y = prepare_for_stan_format(data_single_pair)
  y_missing = purrr::map(y, is.na) %>% purrr::map(function(x) x[,1] | x[,2])
  for (i in 1:length(y)) {
    y[[i]][y_missing[[i]],] <- 10^3
  }
  y_missing <- purrr::map(y_missing,as.integer)
  K=max(data_single_pair$Frame)
  T0=0
  T1=nFrames
  y_missing=rep(0,times=nFrames)
  sigma_0 = c(0.0,0.5,0.5,0.0)
  stan_input = list(dt=dt, Frames=K, nStates = 4,
                    y = y[[1]],
                    y_missing = y_missing,
                    sigma0 =sigma_0, #assume each state equally likely initially
                    T0 =T0, #shall i say T0=T0 based on lines 205,206?
                    T1 = T1,cos_phi = cos_phi)
  
  ##Asymmetric metaphase model
  stan_file= here::here('src/stan_reparametrise_missing/metaphase_asymmetric_full_likelihood_reparametrised_missing.stan')
  m <- stan_model(stan_file) #slightly different way to call a stan model, this compiles the model and the next line does the sampling
  estimate_asym <- sampling(m,
                            data=stan_input,
                            seed = 25,
                            chains = 4,
                            warmup = burn_in,
                            iter = iterations,
                            pars = c("eta","xi","f","auxStates","P","aux"),
                            include=FALSE, #avoid saving the params listed above
                            control=list(adapt_delta=0.95,
                                         max_treedepth=12))
  bridge1<-bridge_sampler(estimate_asym)
  print(bridge1$logml)
  estimate_single_missing_summary = summary(estimate_asym) 
  post_asym_1=as.matrix( estimate_asym)
  states_pos=which(1*(colnames(post_asym_1)=="sigma_sim[1]")==1)
  nstates= 4
  smo_prob_asymmetric_1=matrix(NA,nrow=nFrames,ncol=nstates)
  for(frames in 1:nFrames){
    for(st in 1:nstates){
      smo_prob_asymmetric_1[frames,st]=sum(1*(post_asym_1[,states_pos+frames-1]==st))/(dim(post_asym_1)[1])
    }
  }
  Mp_states=max.col(as.data.frame(smo_prob_asymmetric_1))
  #  par(mfrow=c(2,1))
  #  plot(h_states,col="blue",type="l",pch=3)
  #  plot(Mp_states,col="red",type="l",pch=10)
  
  #  par(mfrow=c(1,1))
  #  plot(h_states,ylim=c(-0.5,4.5),ylab="")
  #  par(new=TRUE)
  #  plot((post_asym_1[3000,c(states_pos:(states_pos+nFrames-1))]-0.2),col="red",ylim=c(-0.5,4.5),ylab="")
  #  par(new=TRUE)
  #  plot((post_asym_1[2900,c(states_pos:(states_pos+nFrames-1))]-0.4),col="blue",ylim=c(-0.5,4.5),ylab="")
  
  return(list(estimate_asym, bridge1))
  ########################### double Missing point #####################################################################
  ##################################################################################################################
}




triple_gap_missing_data_run <- function(data_full, prop_missing, burn_in = 5000,iterations = 15000,dt = 2.05){
  
  d_length = dim(data_full)[1]
  n_gaps = floor(d_length*prop_missing/3)
  yA_to_change = data_full
  i=1
  while(i<= n_gaps){ 
    k = rdunif(1,(1+(i-1)*floor(d_length/n_gaps)), i*floor(d_length/n_gaps) )
    print(k)
    if(k<=(d_length-2)){ 
      yA_to_change[k:(k+2),] =NA  
      i=i+1
    } 
  }
  data_single_pair= data.frame(SisterPairID=double(2*nFrames),Position_1=double(2*nFrames),
                               Frame=integer(2*nFrames),Time=integer(2*nFrames),
                               SisterID=integer(2*nFrames))
  data_single_pair[,1:5]=c(matrix(1,nrow=2*nFrames,ncol=1),cbind(c(yA_to_change[,1],yA_to_change[,2])),rep(c(1:nFrames),times=2),
                           rep(seq(from=0,to=((nFrames-1)*dt),by=dt),times=2),
                           cbind(c(rep(1,times=nFrames),rep(2,times=nFrames))))
  pairIDs <- unique(data_single_pair$SisterPairID)
  y = prepare_for_stan_format(data_single_pair)
  y_missing = purrr::map(y, is.na) %>% purrr::map(function(x) x[,1] | x[,2])
  for (i in 1:length(y)) {
    y[[i]][y_missing[[i]],] <- 10^3
  }
  y_missing <- purrr::map(y_missing,as.integer)
  K=max(data_single_pair$Frame)
  T0=0
  T1=nFrames
  y_missing=rep(0,times=nFrames)
  sigma_0 = c(0.0,0.5,0.5,0.0)
  stan_input = list(dt=dt, Frames=K, nStates = 4,
                    y = y[[1]],
                    y_missing = y_missing,
                    sigma0 =sigma_0, #assume each state equally likely initially
                    T0 =T0, #shall i say T0=T0 based on lines 205,206?
                    T1 = T1,cos_phi = cos_phi)
  
  ##Asymmetric metaphase model
  stan_file= here::here('src/stan_reparametrise_missing/metaphase_asymmetric_full_likelihood_reparametrised_missing.stan')
  m <- stan_model(stan_file) #slightly different way to call a stan model, this compiles the model and the next line does the sampling
  estimate_asym <- sampling(m,
                            data=stan_input,
                            seed = 25,
                            chains = 4,
                            warmup = burn_in,
                            iter = iterations,
                            pars = c("eta","xi","f","auxStates","P","aux"),
                            include=FALSE, #avoid saving the params listed above
                            control=list(adapt_delta=0.95,
                                         max_treedepth=12))
  bridge1<-bridge_sampler(estimate_asym)
  print(bridge1$logml)
  estimate_single_missing_summary = summary(estimate_asym) 
  post_asym_1=as.matrix( estimate_asym)
  states_pos=which(1*(colnames(post_asym_1)=="sigma_sim[1]")==1)
  nstates= 4
  smo_prob_asymmetric_1=matrix(NA,nrow=nFrames,ncol=nstates)
  for(frames in 1:nFrames){
    for(st in 1:nstates){
      smo_prob_asymmetric_1[frames,st]=sum(1*(post_asym_1[,states_pos+frames-1]==st))/(dim(post_asym_1)[1])
    }
  }
  Mp_states=max.col(as.data.frame(smo_prob_asymmetric_1))

  return(list(estimate_asym, bridge1))
  ############################ triple Missing point #####################################################################
  ##################################################################################################################
}



all_in_a_row_gap_missing_data_run <- function(data_full, prop_missing, burn_in = 5000,iterations = 15000,dt = 2.05){
  
  d_length = dim(data_full)[1]
  n_missing = floor(d_length*prop_missing)
  yA_to_change = data_full
  
    k = rdunif(1,1, d_length-2*n_missing)
     yA_to_change[k:(k+n_missing),] =NA    
  
  
  data_single_pair= data.frame(SisterPairID=double(2*nFrames),Position_1=double(2*nFrames),
                               Frame=integer(2*nFrames),Time=integer(2*nFrames),
                               SisterID=integer(2*nFrames))
  data_single_pair[,1:5]=c(matrix(1,nrow=2*nFrames,ncol=1),cbind(c(yA_to_change[,1],yA_to_change[,2])),rep(c(1:nFrames),times=2),
                           rep(seq(from=0,to=((nFrames-1)*dt),by=dt),times=2),
                           cbind(c(rep(1,times=nFrames),rep(2,times=nFrames))))
  pairIDs <- unique(data_single_pair$SisterPairID)
  y = prepare_for_stan_format(data_single_pair)
  y_missing = purrr::map(y, is.na) %>% purrr::map(function(x) x[,1] | x[,2])
  for (i in 1:length(y)) {
    y[[i]][y_missing[[i]],] <- 10^3
  }
  y_missing <- purrr::map(y_missing,as.integer)
  K=max(data_single_pair$Frame)
  T0=0
  T1=nFrames
  y_missing=rep(0,times=nFrames)
  sigma_0 = c(0.0,0.5,0.5,0.0)
  stan_input = list(dt=dt, Frames=K, nStates = 4,
                    y = y[[1]],
                    y_missing = y_missing,
                    sigma0 =sigma_0, #assume each state equally likely initially
                    T0 =T0, #shall i say T0=T0 based on lines 205,206?
                    T1 = T1,cos_phi = cos_phi)
  
  ##Asymmetric metaphase model
  stan_file= here::here('src/stan_reparametrise_missing/metaphase_asymmetric_full_likelihood_reparametrised_missing.stan')
  m <- stan_model(stan_file) #slightly different way to call a stan model, this compiles the model and the next line does the sampling
  estimate_asym <- sampling(m,
                            data=stan_input,
                            seed = 25,
                            chains = 4,
                            warmup = burn_in,
                            iter = iterations,
                            pars = c("eta","xi","f","auxStates","P","aux"),
                            include=FALSE, #avoid saving the params listed above
                            control=list(adapt_delta=0.95,
                                         max_treedepth=12))
  bridge1<-bridge_sampler(estimate_asym)
  print(bridge1$logml)
  estimate_single_missing_summary = summary(estimate_asym) 
  post_asym_1=as.matrix( estimate_asym)
  states_pos=which(1*(colnames(post_asym_1)=="sigma_sim[1]")==1)
  nstates= 4
  smo_prob_asymmetric_1=matrix(NA,nrow=nFrames,ncol=nstates)
  for(frames in 1:nFrames){
    for(st in 1:nstates){
      smo_prob_asymmetric_1[frames,st]=sum(1*(post_asym_1[,states_pos+frames-1]==st))/(dim(post_asym_1)[1])
    }
  }
  Mp_states=max.col(as.data.frame(smo_prob_asymmetric_1))
  #  par(mfrow=c(2,1))
  #  plot(h_states,col="blue",type="l",pch=3)
  #  plot(Mp_states,col="red",type="l",pch=10)
  
  #  par(mfrow=c(1,1))
  #  plot(h_states,ylim=c(-0.5,4.5),ylab="")
  #  par(new=TRUE)
  #  plot((post_asym_1[3000,c(states_pos:(states_pos+nFrames-1))]-0.2),col="red",ylim=c(-0.5,4.5),ylab="")
  #  par(new=TRUE)
  #  plot((post_asym_1[2900,c(states_pos:(states_pos+nFrames-1))]-0.4),col="blue",ylim=c(-0.5,4.5),ylab="")
  
  return(list(estimate_asym, bridge1))
  ############################ all in a row Missing point #####################################################################
  ##################################################################################################################
}


get_missing_estimations_dfs_mean_se = function(list_runs){
#get the dataframe with all the mean estimations and se for the mean estiamtion
#from the full asymmetric model
#takes as input the list with the stan (S4) files for all the experiments

no_experiments = length(list_runs)
summary_runs = sapply(list_runs,"[",1) %>% purrr::map(summary) %>% sapply("[",1)
names(summary_runs) = seq(1,no_experiments,1)
df_tau1 = summary_runs %>% map_df(., function(x) c(x["theta[1,1]",c(1,3,4,8)]))%>%rename(c("tau1" ="mean","tau1_sd"="sd", "tau1_2.5" = "2.5%", "tau1_97.5"= "97.5%" ))
df_tau2 = summary_runs %>% map_df(., function(x) c(x["theta[1,2]",c(1,3,4,8)]))%>%rename(c("tau2" ="mean","tau2_sd"="sd", "tau2_2.5" = "2.5%", "tau2_97.5"= "97.5%" ))
df_alpha = summary_runs %>% map_df(., function(x) c(x["theta[2,1]",c(1,3,4,8)]))%>%rename(c("alpha" ="mean","alpha_sd"="sd", "alpha_2.5" = "2.5%", "alpha_97.5"= "97.5%" ))
df_kappa = summary_runs %>% map_df(., function(x) c(x["theta[3,1]",c(1,3,4,8)]))%>%rename(c("kappa" ="mean","kappa_sd"="sd", "kappa_2.5" = "2.5%", "kappa_97.5"= "97.5%" ))
df_L = summary_runs %>% map_df(., function(x) c(x["theta[8,1]",c(1,3,4,8)]))%>%rename(c("L" ="mean","L_sd"="sd", "L_2.5" = "2.5%", "L_97.5"= "97.5%" ))
df_picoh = summary_runs %>% map_df(., function(x) c(x["theta[6,1]",c(1,3,4,8)]))%>%rename(c("picoh" ="mean","picoh_sd"="sd", "picoh_2.5" = "2.5%", "picoh_97.5"= "97.5%" ))
df_pcoh = summary_runs %>% map_df(., function(x) c(x["theta[7,1]",c(1,3,4,8)]))%>%rename(c("pcoh" ="mean","pcoh_sd"="sd", "pcoh_2.5" = "2.5%", "pcoh_97.5"= "97.5%" ))
df_vminus1 = summary_runs %>% map_df(., function(x) c(x["theta[4,1]",c(1,3,4,8)]))%>%rename(c("vminus1" ="mean","vminus1_sd"="sd", "vminus1_2.5" = "2.5%", "vminus1_97.5"= "97.5%" ))
df_vminus2 = summary_runs %>% map_df(., function(x) c(x["theta[4,2]",c(1,3,4,8)]))%>%rename(c("vminus2" ="mean","vminus2_sd"="sd", "vminus2_2.5" = "2.5%", "vminus2_97.5"= "97.5%" ))
df_vplus1 = summary_runs %>% map_df(., function(x) c(x["theta[5,1]",c(1,3,4,8)]))%>%rename(c("vplus1" ="mean","vplus1_sd"="sd", "vplus1_2.5" = "2.5%", "vplus1_97.5"= "97.5%" ))
df_vplus2 = summary_runs %>% map_df(., function(x) c(x["theta[5,2]",c(1,3,4,8)]))%>%rename(c("vplus2" ="mean","vplus2_sd"="sd", "vplus2_2.5" = "2.5%", "vplus2_97.5"= "97.5%" ))


df = cbind(df_tau1, df_tau2, df_alpha, df_kappa, df_L, df_picoh, df_pcoh, df_vminus1,df_vminus2, df_vplus1,df_vplus2)

return(df)

}


 
get_missing_estimations_dfs_mean_sd = function(df, var1 = "tau1", scale_plt = 0.2, actual_value){
  no_experiments = dim(df)[1]
  if(var1 =="tau1"){ylabel = expression(tau[1])
  }else if(var1 == "tau2"){ylabel = expression(tau[2])
  }else if(var1 =="vminus1"){ylabel = expression("v"["-"]^1)
  }else if(var1 =="vminus2"){ylabel = expression("v"["-"]^{2})
  }else if(var1 =="vplus1"){ylabel = expression("v"["+"]^{1})
  }else if(var1 =="vplus2"){ylabel = expression("v"["+"]^{2})
  }else if(var1 =="picoh"){ylabel = expression("p"["icoh"])
  }else if(var1 =="pcoh"){ylabel = expression("p"["coh"])
  }else if(var1 =="L"){ylabel = "L"
  }else if(var1 =="alpha"){ylabel = expression(alpha)
  }else{ylabel = expression(kappa)}
  var_se = paste0(var1,"_sd")
  var2 = paste0(var1,"_sd/", "scale_plt")
  df = df  %>% mutate(!! rlang::sym(var_se)/scale_plt, .keep = "all")
  h = ggplot(data = df,aes(x = seq(1,no_experiments,1),y = .data[[var1]])) + geom_point(aes(color = "mean"), size = 3)+
    geom_point(aes(y = .data[[var2]],color= "sd"), shape =2)+      
    geom_line(aes(x=seq(1,no_experiments,1), y = actual_value, color = "actual"), alpha = 0.3)+
    scale_y_continuous(sec.axis = sec_axis(~.*scale_plt, name="standard deviation"))+ 
    labs(x = "experiment", y = ylabel)+scale_color_manual(values = c("darkorange4", "orange2","gray30")) +
    theme_bw()+theme(text = element_text(size = 9),
                     legend.justification = c("right", "top"),
                     legend.margin = margin(5, 5, 5, 5),legend.direction = "vertical",legend.title=element_blank())
  
  return(h)
}



get_missing_estimations_dfs_mean_confidence = function(df, var1 = "tau1", actual_value){
  no_experiments = dim(df)[1]
  var_sd = paste0(var1,"_sd")
  if(var1 =="tau1"){ylabel = expression(tau[1])
  }else if(var1 == "tau2"){ylabel = expression(tau[2])
  }else if(var1 =="vminus1"){ylabel = expression("v"["-"]^1)
  }else if(var1 =="vminus2"){ylabel = expression("v"["-"]^{2})
  }else if(var1 =="vplus1"){ylabel = expression("v"["+"]^{1})
  }else if(var1 =="vplus2"){ylabel = expression("v"["+"]^{2})
  }else if(var1 =="picoh"){ylabel = expression("p"["icoh"])
  }else if(var1 =="pcoh"){ylabel = expression("p"["coh"])
  }else if(var1 =="L"){ylabel = "L"
  }else if(var1 =="alpha"){ylabel = expression(alpha)
  }else{ylabel = expression(kappa)}
#  df = df  %>% mutate(!! rlang::sym(var_sd)/scale_plt, .keep = "all")
  h = ggplot(data = df,aes(x = seq(1,no_experiments,1),y = .data[[var1]])) + geom_point(aes(color = "mean"))+
    geom_errorbar(aes(ymin = .data[[var1]]-.data[[var_sd]],ymax = .data[[var1]]+.data[[var_sd]], color = "+/- sd"),width=0.2)+      
    geom_line(aes(x=seq(1,no_experiments,1), y = actual_value, color = "actual"), alpha = 0.3)+
    scale_x_continuous(limits=c(0,10), breaks=seq(1,10,1))+
    labs(x = "experiment", y = ylabel)+scale_color_manual(values = c("darkorange4", "orange2","gray30"))+
    theme_bw()+theme(text = element_text(size = 9),  legend.position = c(0.95,0.95),
     legend.justification = c("right", "top"),
     legend.margin = margin(5, 5, 5, 5),legend.direction = "horizontal",legend.title=element_blank()) 
  return(h)
}


level_order <- factor(iris$Species, level = c('virginica', 'versicolor', 'setosa'))


get_missing_plots_different_missing_setting_c_frames = function(df, missing_settings = c("no_missing", "singles","doubles","triples", "all_in_a_row"),var1 = "tau1", actual_value){
  library(forcats)
  var_sd = paste0(var1,"_sd")
  if(var1 =="tau1"){ylabel = expression(tau[1])
  }else if(var1 == "tau2"){ylabel = expression(tau[2])
  }else if(var1 =="vminus1"){ylabel = expression("v"["-"]^1)
  }else if(var1 =="vminus2"){ylabel = expression("v"["-"]^{2})
  }else if(var1 =="vplus1"){ylabel = expression("v"["+"]^{1})
  }else if(var1 =="vplus2"){ylabel = expression("v"["+"]^{2})
  }else if(var1 =="picoh"){ylabel = expression("p"["icoh"])
  }else if(var1 =="pcoh"){ylabel = expression("p"["coh"])
  }else if(var1 =="L"){ylabel = "L"
  }else if(var1 =="alpha"){ylabel = expression(alpha)
  }else{ylabel = expression(kappa)}
  
  #  df = df  %>% mutate(!! rlang::sym(var_sd)/scale_plt, .keep = "all")
  h = ggplot(data = df,aes(x = fct_inorder(missing_settings),y = .data[[var1]])) + geom_point(aes(color = "mean")) +
    geom_errorbar(aes(ymin = .data[[var1]]-.data[[var_sd]],ymax = .data[[var1]]+.data[[var_sd]], color = "+/- sd"),width=0.2)+      
    geom_line(aes(x= seq(1,5,1), y = actual_value, group = 1,color = "actual"), alpha = 0.3)+
    labs(x = "missing setting", y = ylabel)+scale_color_manual(values = c("darkorange4","orange", "gray30"))+
    theme_bw()+theme(axis.text.x = element_text(angle = 0,face ="bold"),text = element_text(size = 9),legend.title=element_blank(), legend.position ="none")
  #legend.position = c(0.90,0.85),
  # legend.justification = c("right", "top"),
  # legend.margin = margin(5, 5, 5, 5),legend.direction = "vertical",legend.title=element_blank()) 
  return(h)
}

get_missing_plots_different_missing_setting_sd = function(df, missing_settings = c("no_missing", "singles","doubles","triples", "all_in_a_row"),var1 = "tau1"){
  library(forcats)
  var_sd = paste0(var1,"_sd")
  if(var1 =="tau1"){ylabel = expression(tau[1])
  }else if(var1 == "tau2"){ylabel = expression(tau[2])
  }else if(var1 =="vminus1"){ylabel = expression("v"["-"]^1)
  }else if(var1 =="vminus2"){ylabel = expression("v"["-"]^{2})
  }else if(var1 =="vplus1"){ylabel = expression("v"["+"]^{1})
  }else if(var1 =="vplus2"){ylabel = expression("v"["+"]^{2})
  }else if(var1 =="picoh"){ylabel = expression("p"["icoh"])
  }else if(var1 =="pcoh"){ylabel = expression("p"["coh"])
  }else if(var1 =="L"){ylabel = "L"
  }else if(var1 =="alpha"){ylabel = expression(alpha)
  }else{ylabel = expression(kappa)}
  h = ggplot(data = df,aes(x = fct_inorder(missing_settings),y = .data[[var_sd]])) + geom_line() 
 return(h) 
}

get_missing_plots_different_missing_setting_c_frames_quantiles = function(df, missing_settings = c("no_missing", "singles","doubles","triples", "all_in_a_row"),var1 = "tau1", actual_value){
  library(forcats)
  q025 = paste0(var1,"_2.5")
  q975 = paste0(var1,"_97.5")
  if(var1 =="tau1"){ylabel = expression(tau[1])
  }else if(var1 == "tau2"){ylabel = expression(tau[2])
  }else if(var1 =="vminus1"){ylabel = expression("v"["-"]^1)
  }else if(var1 =="vminus2"){ylabel = expression("v"["-"]^{2})
  }else if(var1 =="vplus1"){ylabel = expression("v"["+"]^{1})
  }else if(var1 =="vplus2"){ylabel = expression("v"["+"]^{2})
  }else if(var1 =="picoh"){ylabel = expression("p"["icoh"])
  }else if(var1 =="pcoh"){ylabel = expression("p"["coh"])
  }else if(var1 =="L"){ylabel = "L"
  }else if(var1 =="alpha"){ylabel = expression(alpha)
  }else{ylabel = expression(kappa)}
  #  df = df  %>% mutate(!! rlang::sym(var_sd)/scale_plt, .keep = "all")
  h = ggplot(data = df,aes(x = fct_inorder(missing_settings),y = .data[[var1]])) + geom_point(aes(color = "mean")) +
    geom_errorbar(aes(ymin = .data[[q025]],ymax = .data[[q975]], color = "quantiles"),width=0.2)+      
    geom_line(aes(x= seq(1,5,1), y = actual_value, group = 1,color = "actual"), alpha = 0.3)+
    labs(x = "missing setting", y = ylabel)+scale_color_manual(values = c("darkorange3","orange", "gray30"))+
    theme_bw()+theme(axis.text.x = element_text(angle = 0,face ="bold"),text = element_text(size = 9),legend.title=element_blank(), legend.position ="none")
  #legend.position = c(0.90,0.85),
  # legend.justification = c("right", "top"),
  # legend.margin = margin(5, 5, 5, 5),legend.direction = "vertical",legend.title=element_blank()) 
  return(h)
}

get_missing_plots_2df_different_ms_comparison = function(df1, df2, missing_settings = c("no_missing", "singles","doubles","triples", "all_in_a_row"),var1 = "tau1", actual_value){
  library(forcats)
  var_sd = paste0(var1,"_sd")
  if(var1 =="tau1"){ylabel = expression(tau[1])
  }else if(var1 == "tau2"){ylabel = expression(tau[2])
  }else if(var1 =="vminus1"){ylabel = expression("v"["-"]^1)
  }else if(var1 =="vminus2"){ylabel = expression("v"["-"]^{2})
  }else if(var1 =="vplus1"){ylabel = expression("v"["+"]^{1})
  }else if(var1 =="vplus2"){ylabel = expression("v"["+"]^{2})
  }else if(var1 =="picoh"){ylabel = expression("p"["icoh"])
  }else if(var1 =="pcoh"){ylabel = expression("p"["coh"])
  }else if(var1 =="L"){ylabel = "L"
  }else if(var1 =="alpha"){ylabel = expression(alpha)
  }else{ylabel = expression(kappa)}
  h = ggplot(data = df1,aes(x = fct_inorder(missing_settings),y = .data[[var1]])) + geom_point(aes(color = "mean")) +  
    geom_errorbar(aes(ymin = .data[[var1]]-.data[[var_sd]],ymax = .data[[var1]]+.data[[var_sd]], color = "+/- sd"),width=0.2) + 
    geom_point(data = df2,aes(color = "mean_inter"), shape = 17)+
    geom_errorbar(data = df2,aes(ymin = .data[[var1]]-.data[[var_sd]],ymax = .data[[var1]]+.data[[var_sd]], color = "+/- sd_inter"),width=0.2, linetype = "dashed")+
    geom_line(aes(x= seq(1,5,1), y = actual_value, group = 1,color = "actual"), alpha = 0.3)+
    theme_bw()+theme(axis.title.x=element_blank(),axis.text.x = element_text(angle = 20, face = "bold"),text = element_text(size = 10),legend.title=element_blank(), legend.position ="none")
  
  return(h)
}

single_gap_missing_data_interpolation_run <- function(data_full, prop_missing, burn_in = 5000,iterations = 15000, dt = 2.05){
  
  d_length = dim(data_full)[1]
  n_gaps = d_length*prop_missing
  yA_to_change = data_full
  
  for(i in 1:(n_gaps)){ 
    k = rdunif(1,(1+(i-1)*floor(d_length/n_gaps)), i*floor(d_length/n_gaps) )
    yA_to_change[k,] =NA    
  }
  
  data_single_pair= data.frame(SisterPairID=double(2*nFrames),Position_1=double(2*nFrames),
                               Frame=integer(2*nFrames),Time=integer(2*nFrames),
                               SisterID=integer(2*nFrames))
  data_single_pair[,1:5]=c(matrix(1,nrow=2*nFrames,ncol=1),cbind(c(yA_to_change[,1],yA_to_change[,2])),rep(c(1:nFrames),times=2),
                           rep(seq(from=0,to=((nFrames-1)*dt),by=dt),times=2),
                           cbind(c(rep(1,times=nFrames),rep(2,times=nFrames))))
 
data_single_pair <-  data_single_pair%>%
    dplyr::group_by(SisterPairID,SisterID) %>%
    dplyr::mutate(Position_1=interpolate_missing_data(Position_1,Time)) %>%
     dplyr::ungroup() 
  
  pairIDs <- unique(data_single_pair$SisterPairID)
  y = prepare_for_stan_format(data_single_pair)
  y_missing = purrr::map(y, is.na) %>% purrr::map(function(x) x[,1] | x[,2])
  for (i in 1:length(y)) {
    y[[i]][y_missing[[i]],] <- 0
  }
  y_missing <- purrr::map(y_missing,as.integer)
  K=max(data_single_pair$Frame)
  
  trajectory_index <- which(pairIDs==1)
  start_end <- purrr::map_df(y_missing,get_start_end_of_nonmissing_data,.id="traj_index") %>%
    mutate(traj_index=as.integer(traj_index)) #use this info to only fit to the existing tracks
  T0 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T0);
  T1 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T1);
  
  y_missing=rep(0,times=nFrames)
  sigma_0 = c(0.0,0.5,0.5,0.0)
  stan_input = list(dt=dt, Frames=K, nStates = 4,
                    y = y[[1]],
                    y_missing = y_missing,
                    sigma0 =sigma_0, #assume each state equally likely initially
                    T0 =T0, #shall i say T0=T0 based on lines 205,206?
                    T1 = T1,cos_phi = cos_phi)
  
  ##Asymmetric metaphase model
  stan_file= here::here('src/stan_reparametrise/metaphase_asymmetric_full_likelihood_reparametrised.stan') #Asymmetry in all variables   (v+,v-,tau)
  m <- stan_model(stan_file) #slightly different way to call a stan model, this compiles the model and the next line does the sampling
  estimate_asym <- sampling(m,
                             data=stan_input,
                             seed = 42,
                             chains = 4,
                             warmup = burn_in,
                             iter = iterations,
                             init_r=0.1,
                             pars = c("f","eta","xi","auxStates","P","aux"),
                             include=FALSE, #avoid saving the params listed above
                             control=list(adapt_delta=0.95,
                                          max_treedepth=12))
  
  bridge1<-bridge_sampler(estimate_asym)
  print(bridge1$logml)
  estimate_single_missing_summary = summary(estimate_asym) 
  post_asym_1=as.matrix( estimate_asym)
  states_pos=which(1*(colnames(post_asym_1)=="sigma_sim[1]")==1)
  nstates= 4
  smo_prob_asymmetric_1=matrix(NA,nrow=nFrames,ncol=nstates)
  for(frames in 1:nFrames){
    for(st in 1:nstates){
      smo_prob_asymmetric_1[frames,st]=sum(1*(post_asym_1[,states_pos+frames-1]==st))/(dim(post_asym_1)[1])
    }
  }
  Mp_states=max.col(as.data.frame(smo_prob_asymmetric_1))
  return(list(estimate_asym, bridge1))
}


double_gap_missing_data_interpolation_run <- function(data_full, prop_missing, burn_in = 5000,iterations = 15000,dt = 2.05){
  
  d_length = dim(data_full)[1]
  n_gaps = floor(d_length*prop_missing/2)
  yA_to_change = data_full
  i=1
  while(i<= n_gaps){ 
    k = rdunif(1,(1+(i-1)*floor(d_length/n_gaps)), i*floor(d_length/n_gaps) )
    print(k)
    if(k<=(d_length-1)){ 
      yA_to_change[k:(k+1),] =NA  
      i=i+1
    } 
  }
  
  data_single_pair= data.frame(SisterPairID=double(2*nFrames),Position_1=double(2*nFrames),
                               Frame=integer(2*nFrames),Time=integer(2*nFrames),
                               SisterID=integer(2*nFrames))
  data_single_pair[,1:5]=c(matrix(1,nrow=2*nFrames,ncol=1),cbind(c(yA_to_change[,1],yA_to_change[,2])),rep(c(1:nFrames),times=2),
                           rep(seq(from=0,to=((nFrames-1)*dt),by=dt),times=2),
                           cbind(c(rep(1,times=nFrames),rep(2,times=nFrames))))
  
  data_single_pair <-  data_single_pair%>%
    dplyr::group_by(SisterPairID,SisterID) %>%
    dplyr::mutate(Position_1=interpolate_missing_data(Position_1,Time)) %>%
    dplyr::ungroup() 
  
  pairIDs <- unique(data_single_pair$SisterPairID)
  y = prepare_for_stan_format(data_single_pair)
  y_missing = purrr::map(y, is.na) %>% purrr::map(function(x) x[,1] | x[,2])
  for (i in 1:length(y)) {
    y[[i]][y_missing[[i]],] <- 0
  }
  y_missing <- purrr::map(y_missing,as.integer)
  K=max(data_single_pair$Frame)
  
  trajectory_index <- which(pairIDs==1)
  start_end <- purrr::map_df(y_missing,get_start_end_of_nonmissing_data,.id="traj_index") %>%
    mutate(traj_index=as.integer(traj_index)) #use this info to only fit to the existing tracks
  T0 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T0);
  T1 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T1);
  
  y_missing=rep(0,times=nFrames)
  sigma_0 = c(0.0,0.5,0.5,0.0)
  stan_input = list(dt=dt, Frames=K, nStates = 4,
                    y = y[[1]],
                    y_missing = y_missing,
                    sigma0 =sigma_0, #assume each state equally likely initially
                    T0 =T0, #shall i say T0=T0 based on lines 205,206?
                    T1 = T1,cos_phi = cos_phi)
  
  ##Asymmetric metaphase model
  stan_file= here::here('src/stan_reparametrise/metaphase_asymmetric_full_likelihood_reparametrised.stan') #Asymmetry in all variables   (v+,v-,tau)
  m <- stan_model(stan_file) #slightly different way to call a stan model, this compiles the model and the next line does the sampling
  estimate_asym <- sampling(m,
                             data=stan_input,
                             seed = 42,
                             chains = 4,
                             warmup = burn_in,
                             iter = iterations,
                             init_r=0.1,
                             pars = c("f","eta","xi","auxStates","P","aux"),
                             include=FALSE, #avoid saving the params listed above
                             control=list(adapt_delta=0.95,
                                          max_treedepth=12))
  
  bridge1<-bridge_sampler(estimate_asym)
  estimate_single_missing_summary = summary(estimate_asym) 
  post_asym_1=as.matrix( estimate_asym)
  states_pos=which(1*(colnames(post_asym_1)=="sigma_sim[1]")==1)
  nstates= 4
  smo_prob_asymmetric_1=matrix(NA,nrow=nFrames,ncol=nstates)
  for(frames in 1:nFrames){
    for(st in 1:nstates){
      smo_prob_asymmetric_1[frames,st]=sum(1*(post_asym_1[,states_pos+frames-1]==st))/(dim(post_asym_1)[1])
    }
  }
  Mp_states=max.col(as.data.frame(smo_prob_asymmetric_1))
  return(list(estimate_asym, bridge1))
  ########################### double Missing point #####################################################################
  ##################################################################################################################
}




triple_gap_missing_data_interpolation_run <- function(data_full, prop_missing, burn_in = 5000,iterations = 15000,dt=2.05){
  
  d_length = dim(data_full)[1]
  n_gaps = floor(d_length*prop_missing/3)
  yA_to_change = data_full
  i=1
  while(i<= n_gaps){ 
    k = rdunif(1,(1+(i-1)*floor(d_length/n_gaps)), i*floor(d_length/n_gaps) )
    if(k<=(d_length-2)){ 
      yA_to_change[k:(k+2),] =NA  
      i=i+1
    } 
  }
  data_single_pair= data.frame(SisterPairID=double(2*nFrames),Position_1=double(2*nFrames),
                               Frame=integer(2*nFrames),Time=integer(2*nFrames),
                               SisterID=integer(2*nFrames))
  data_single_pair[,1:5]=c(matrix(1,nrow=2*nFrames,ncol=1),cbind(c(yA_to_change[,1],yA_to_change[,2])),rep(c(1:nFrames),times=2),
                           rep(seq(from=0,to=((nFrames-1)*dt),by=dt),times=2),
                           cbind(c(rep(1,times=nFrames),rep(2,times=nFrames))))
  
  data_single_pair <-  data_single_pair%>%
    dplyr::group_by(SisterPairID,SisterID) %>%
    dplyr::mutate(Position_1=interpolate_missing_data(Position_1,Time)) %>%
    dplyr::ungroup() 
  
  pairIDs <- unique(data_single_pair$SisterPairID)
  y = prepare_for_stan_format(data_single_pair)
  y_missing = purrr::map(y, is.na) %>% purrr::map(function(x) x[,1] | x[,2])
  for (i in 1:length(y)) {
    y[[i]][y_missing[[i]],] <- 0
  }
  y_missing <- purrr::map(y_missing,as.integer)
  K=max(data_single_pair$Frame)
  
  trajectory_index <- which(pairIDs==1)
  start_end <- purrr::map_df(y_missing,get_start_end_of_nonmissing_data,.id="traj_index") %>%
    mutate(traj_index=as.integer(traj_index)) #use this info to only fit to the existing tracks
  T0 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T0);
  T1 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T1);
  
  y_missing=rep(0,times=nFrames)
  sigma_0 = c(0.0,0.5,0.5,0.0)
  stan_input = list(dt=dt, Frames=K, nStates = 4,
                    y = y[[1]],
                    y_missing = y_missing,
                    sigma0 =sigma_0, #assume each state equally likely initially
                    T0 =T0, #shall i say T0=T0 based on lines 205,206?
                    T1 = T1,cos_phi = cos_phi)
  
  ##Asymmetric metaphase model
  stan_file= here::here('src/stan_reparametrise/metaphase_asymmetric_full_likelihood_reparametrised.stan') #Asymmetry in all variables   (v+,v-,tau)
  m <- stan_model(stan_file) #slightly different way to call a stan model, this compiles the model and the next line does the sampling
  estimate_asym <- sampling(m,
                             data=stan_input,
                             seed = 42,
                             chains = 4,
                             warmup = burn_in,
                             iter = iterations,
                             init_r=0.1,
                             pars = c("f","eta","xi","auxStates","P","aux"),
                             include=FALSE, #avoid saving the params listed above
                             control=list(adapt_delta=0.95,
                                          max_treedepth=12))
  
  bridge1<-bridge_sampler(estimate_asym)
  return(list(estimate_asym, bridge1))
  ############################ triple Missing point #####################################################################
  ##################################################################################################################
}



all_in_a_row_gap_missing_data_interpolation_run <- function(data_full, prop_missing, burn_in = 5000,iterations = 15000,dt = 2.05){
  
  d_length = dim(data_full)[1]
  n_missing = floor(d_length*prop_missing)
  yA_to_change = data_full
  
  k = rdunif(1,1, d_length-2*n_missing)
  yA_to_change[k:(k+n_missing),] =NA    
  
  data_single_pair= data.frame(SisterPairID=double(2*nFrames),Position_1=double(2*nFrames),
                               Frame=integer(2*nFrames),Time=integer(2*nFrames),
                               SisterID=integer(2*nFrames))
  data_single_pair[,1:5]=c(matrix(1,nrow=2*nFrames,ncol=1),cbind(c(yA_to_change[,1],yA_to_change[,2])),rep(c(1:nFrames),times=2),
                           rep(seq(from=0,to=((nFrames-1)*dt),by=dt),times=2),
                           cbind(c(rep(1,times=nFrames),rep(2,times=nFrames))))
  
  data_single_pair <-  data_single_pair%>%
    dplyr::group_by(SisterPairID,SisterID) %>%
    dplyr::mutate(Position_1=interpolate_missing_data(Position_1,Time)) %>%
    dplyr::ungroup() 
  
  pairIDs <- unique(data_single_pair$SisterPairID)
  y = prepare_for_stan_format(data_single_pair)
  y_missing = purrr::map(y, is.na) %>% purrr::map(function(x) x[,1] | x[,2])
  for (i in 1:length(y)) {
    y[[i]][y_missing[[i]],] <- 0
  }
  y_missing <- purrr::map(y_missing,as.integer)
  K=max(data_single_pair$Frame)
  
  trajectory_index <- which(pairIDs==1)
  start_end <- purrr::map_df(y_missing,get_start_end_of_nonmissing_data,.id="traj_index") %>%
    mutate(traj_index=as.integer(traj_index)) #use this info to only fit to the existing tracks
  T0 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T0);
  T1 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T1);
  
  y_missing=rep(0,times=nFrames)
  sigma_0 = c(0.0,0.5,0.5,0.0)
  stan_input = list(dt=dt, Frames=K, nStates = 4,
                    y = y[[1]],
                    y_missing = y_missing,
                    sigma0 =sigma_0, #assume each state equally likely initially
                    T0 =T0, #shall i say T0=T0 based on lines 205,206?
                    T1 = T1,cos_phi = cos_phi)
  
  ##Asymmetric metaphase model
  stan_file= here::here('src/stan_reparametrise/metaphase_asymmetric_full_likelihood_reparametrised.stan') #Asymmetry in all variables   (v+,v-,tau)
  m <- stan_model(stan_file) #slightly different way to call a stan model, this compiles the model and the next line does the sampling
  estimate_asym <- sampling(m,
                             data=stan_input,
                             seed = 42,
                             chains = 4,
                             warmup = burn_in,
                             iter = iterations,
                             init_r=0.1,
                             pars = c("f","eta","xi","auxStates","P","aux"),
                             include=FALSE, #avoid saving the params listed above
                             control=list(adapt_delta=0.95,
                                          max_treedepth=12))
  
  bridge<-bridge_sampler(estimate_asym)

  return(list(estimate_asym, bridge))
  ############################ all in a row Missing point #####################################################################
  ##################################################################################################################
}

