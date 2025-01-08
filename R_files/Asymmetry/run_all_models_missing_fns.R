run_all_models_missing = function(jobset_str, warm_up_iter = 10000, total_iter=20000, warm_up_iter_anaphase = 2000, total_iter_anaphase = 4000, dt=2.05, K=Inf, p_missing_final = 0.2, prop_missing_initial = 0.3){
  #set up
  
  whole_name2 = str_split(jobset_str,"per_frame/")[[1]][2]
  whole_name2 = str_split(whole_name2,".ome.csv")[[1]][1]
  whole_name=str_split(jobset_str,"per_frame/kittracking0")
  split_name=str_split(whole_name[[1]][2],".ome.csv")
  cell_name=paste(split_name[[1]][1], sep="")
  split_name= str_split(cell_name,"-kitjobset")
  cell_name=paste(split_name[[1]][1], sep="")
  
  print(whole_name)
  print(cell_name)
  
  
#  data_single_pair_old = process_jobset_keep_missing(jobset_str,K=Inf, start_from = 0,max_missing = p_missing_final) %>%
#    filter(!is.na(SisterID)) #omit unpaired KTs
#  length(unique(data_single_pair_old$SisterPairID))

  data_single_pair = filter_cell_data(jobset_str, max_anaphase_frames = 45, p_missing_final = p_missing_final, prop_missing = prop_missing_initial)
  length(unique(data_single_pair$SisterPairID))
  
  
  g <- ggplot(data_single_pair, aes(x=Time, y=Position_1,color=factor(SisterID))) +
    geom_line() +
    facet_wrap(.~SisterPairID) +
    theme_bw() +
    theme(legend.position = "None",
          strip.text.x = element_text(size = 8),
          axis.text.x = element_text(angle=90)) +
    labs(x="Time (s)",y="Position (um)")
  print(length(unique(data_single_pair$SisterPairID)))
  #  ggsave(here::here(gsub(""," ",paste("Trajectories_cell_",cell_name,split_name[[1]][2],".pdf"))))
  
  #Sister pairs with no missing data for this particular cell
  Complete_sisters=unique(data_single_pair$SisterPairID)
  
  for(i in Complete_sisters){
    h <- ggplot(data_single_pair %>% filter(SisterPairID==i), aes(x=Time, y=Position_1,color=factor(SisterID))) +
      geom_line() +
      theme_bw() +
      theme(legend.position = "None",
            strip.text.x = element_text(size = 8),
            axis.text.x = element_text(angle=90)) +
      labs(x="Time (s)",y="Position (um)")
    print(h)
  }
  
  g <- ggplot(data_single_pair, aes(x=Time, y=Position_1,color=factor(SisterID))) +
    geom_line() +
    facet_wrap(.~SisterPairID) +
    theme_bw() +
    theme(legend.position = "None",
          strip.text.x = element_text(size = 8),
          axis.text.x = element_text(angle=90)) +
    labs(x="Time (s)",y="Position (um)")
  #   ggsave(here::here(gsub(""," ",paste("Trajectories2_cell_",cell_name,split_name[[1]][2],".pdf"))))
  print("first_point")
  
  p_missing_interpolate = 0.25
  #####################################################################################################
  pairs_to_include = unique(data_single_pair$SisterPairID)
  data_single_pair_interpolate <- process_jobset(jobset_str,K=K,max_missing=p_missing_interpolate) %>%
    filter(!is.na(SisterID)) 
  #omit unpaired KTs 
  data_single_pair_interpolate <- data_single_pair_interpolate %>%filter( SisterPairID %in% (pairs_to_include)) #keep sisters that have pass the criterion of filter_cell_data
  
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
  for(i in pairIDs){
    stan_input = list(dt=dt, T=K,
                      nTracks = 1,
                      y = y_interpolate[[which(pairIDs==i)]],
                      y_missing = y_interpolate_missing[[which(pairIDs==i)]]
    )
    changept_estimate <- stan(file=here::here('~/Documents/GitHub/MetaAnaDynamics/STAN_files/Asymmetry/anaphase_changepoint.stan'),
                              data=stan_input,
                              seed = 42,
                              chains =4,
                              warmup = warm_up_iter_anaphase,
                              iter = total_iter_anaphase,
                              control=list(adapt_delta=0.95,max_treedepth=12))

    #Find the minimum anaphase time for all sisters

    pars_to_plot = c("tau","alpha","beta","a","t_ana")
    color_scheme_set("purple")
    posterior <- as.array(changept_estimate)
    #mcmc_trace(posterior,pars=pars_to_plot)
    ch_es<-as.data.frame(changept_estimate)
    min2tana<-mean(ch_es$t_ana)
    anaphase_times[,j]<-c(min2tana,i)
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
  Mean_anaphase = Mean_anaphase %>% select(Mean_t_ana)%>%mutate(Mean_anaphase, Mean_t_ana = ifelse((Mean_t_ana< median(Mean_anaphase$Mean_t_ana)-100), median(Mean_anaphase$Mean_t_ana)-30,Mean_t_ana ))

 # write_xlsx(Mean_anaphase,here::here(paste0("Mean_anaphase_times_",whole_name2,"_all.xlsx")))
#  write.csv(Mean_anaphase,here::here(paste0("Mean_anaphase_times_",whole_name2,".csv")),row.names=FALSE)


  #Maybe more accurate to go with median instead of minimum to reduce the effect of a possible outlier/false change point.

  min_tana_new =  min(Mean_anaphase$Mean_t_ana)
  data_single_pair <- data_single_pair[(data_single_pair$Time<floor(min_tana_new-30))==TRUE,]
  rm(posterior)
  rm(pars_to_plot)


  df_sisters_more_info = get_df_sisters_more(jobset_str, Mean_anaphase, max_anaphase_frames, p_missing_final, prop_missing_initial)
                                   
    
  write_xlsx(df_sisters_more_info,here::here(paste0("Df_sisters_info_more_",whole_name2,"_all.xlsx")))

  g <- ggplot(data_single_pair, aes(x=Time, y=Position_1,color=factor(SisterID))) +
    geom_line() +
    facet_wrap(.~SisterPairID) +
    theme_bw() +
    theme(legend.position = "None",
          strip.text.x = element_text(size = 8),
          axis.text.x = element_text(angle=90)) +
    labs(x="Time (s)",y="Position (um)")
  g
  ggsave(here::here(gsub(" "," ", paste0("Trajectories_cell_",whole_name2,"_metaphase_only.pdf"))))



  #This makes the autocorellation based on the mean position of two sisters
  df_acf <-df_sisters_more_info %>%
    split(.$SisterPairID) %>%
    map(~acf(.$Mean_X_sisters , main = unique(.$SisterPairID), lag.max = 150, plot = F,na.action=na.pass))%>%
    map_dfr(~data.frame(lag = .$lag, acf = .$acf, ci = qnorm(0.975)/sqrt(.$n.used)), .id = "SisterPairID")


  h = ggplot(df_acf, aes(lag, acf))+  geom_line() +
    #  geom_segment(aes(xend = lag, yend = 0)) +
    geom_hline(aes(yintercept = ci), linetype = "dashed", color = "deepskyblue2") +
    geom_hline(aes(yintercept = -ci), linetype = "dashed", color = "deepskyblue2") +
    facet_wrap(~SisterPairID)+
    theme_bw() +
    theme(legend.position = "None",
          strip.text.x = element_text(size = 8),
          axis.text.x = element_text(angle=90))+
    labs(x="Lags",y="Autocorrelation")
  h
  ggsave(here::here(gsub("","",paste0("Autocorrelations_cell_",whole_name2,"_metaphase_NAN_included.pdf"))))

  #Seconds instead of frames in lags.
  df_acf_seconds <- df_acf %>% select(lag) %>%  mutate(df_acf, seconds = lag*2.04933, .keep = "all")

  h = ggplot(df_acf_seconds, aes(seconds, acf))+  geom_line() +
    #  geom_segment(aes(xend = lag, yend = 0)) +
    geom_hline(aes(yintercept = ci), linetype = "dashed", color = "deepskyblue2") +
    geom_hline(aes(yintercept = -ci), linetype = "dashed", color = "deepskyblue2") +
    facet_wrap(~SisterPairID)+
    theme_bw() +
    theme(legend.position = "None",
          strip.text.x = element_text(size = 8),
          axis.text.x = element_text(angle=90))+
    labs(x="Time(s)",y="Autocorrelation")
  h
  ggsave(here::here(gsub("","",paste0("Autocorrelations_cell_",cell_name,split_name[[1]][2],"_metaphase_NAN_included_seconds.pdf"))))
  write.csv(df_acf_seconds,here::here(paste0("Autocorrelations_seconds_metaphase_nan_included_",whole_name2,".csv")),row.names=FALSE)

  pairIDs <- unique(data_single_pair$SisterPairID)
  y = prepare_for_stan_format(data_single_pair)
  y_missing = purrr::map(y, is.na) %>% purrr::map(function(x) x[,1] | x[,2])
  for (i in 1:length(y)) {
    y[[i]][y_missing[[i]],] <- 10^9
  } #replace missing values high values--- see also in .stan file to treat missing (stan does not like NA values)
  y_missing <- purrr::map(y_missing,as.integer)
  K=max(data_single_pair$Frame)
  ########################Bridge Sampling############################################
  Marginal_likelihoods_bridge <- data.frame(SisterPairID=integer(length(pairIDs)),Full=double(length(pairIDs)), v_minus=double(length(pairIDs)), v_plus=double(length(pairIDs)),v_minus_plus=double(length(pairIDs)),symmetric=double(length(pairIDs)))
  ######################################################################################
  Bayes_Factors_bridge <- data.frame(SisterPairID=integer(5*length(pairIDs)),Full=double(5*length(pairIDs)), v_minus=double(5*length(pairIDs)), v_plus=double(5*length(pairIDs)),v_minus_plus=double(5*length(pairIDs)),symmetric=double(5*length(pairIDs)),post_prob=double(5*length(pairIDs)))
  ###########################################################################################

  j=1
  for(i in pairIDs){

    Marginal_likelihoods_bridge[j,1]=i
    # Elpd_estimate[j,1]=i
    # DF_Loo_cv[2*j-1:2*j,1]=i
    trajectory_index <- which(pairIDs==i)
    start_end <- purrr::map_df(y_missing,get_start_end_of_nonmissing_data,.id="traj_index") %>%
    mutate(traj_index=as.integer(traj_index)) #use this info to only fit to the existing tracks
    cos_theta <- get_cos_theta(data_single_pair,pairIDs[trajectory_index])
    cos_theta[is.na(cos_theta)] <- 0 #replace NAs as stan does not deal with these, should not come into the fit due to T0 and T1
    T0 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T0);
    T1 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T1);
    Frames=K
    sigma0 = c(0.0,0.5,0.5,0.0)
    stan_input = list(dt=dt, Frames=K, nStates = 4,
                      y = y[[trajectory_index]],
                      y_missing = y_missing[[trajectory_index]],
                      sigma0 =sigma0, #initial state probabilities
                      T0 =T0, #shall i say T0=T0 based on lines 205,206?
                      T1 = T1,cos_theta = cos_theta )

    ##Asymmetric metaphase model
    stan_file= here::here("~/Documents/GitHub/MetaAnaDynamics/STAN_files/Asymmetry/metaphase_asymmetric_full_likelihood_reparametrised_missing.stan")
    m <- stan_model(stan_file) #slightly different way to call a stan model, this compiles the model and the next line does the sampling

    estimate_asym1 <- sampling(m,
                                data=stan_input,
                                seed = 42,
                                chains = 4,
                                warmup = warm_up_iter,
                                iter = total_iter,
                                init_r=0.1,
                                pars = c("f","eta","xi","auxStates","P","aux"),
                                include=FALSE, #avoid saving the params listed above
                                control=list(adapt_delta=0.99,
                                             max_treedepth=12))
     name2save2 <- paste("Asymmetric_full_vminus_vplus_tau_",whole_name2,"_pair_",i,".rds",sep="")
     print(name2save2)
     saveRDS(estimate_asym1,name2save2 )
    #
#
   #first check the convergence and ESS etc etc and get rid of the kt with issues
    bridge1<-bridgesampling::bridge_sampler(estimate_asym1)
    print(bridge1$logml)
    #bridge1=bridge_sampler(estimate_asym1,method = "warp3")

    Marginal_likelihoods_bridge[j,2]=bridge1$logml

    stan_file= here::here('/Users/constandinakoki/Documents/GitHub/MetaAnaDynamics/STAN_files/Asymmetry/metaphase_asymmetric_vminus_full_likelihood_reparametrised_missing.stan') #Asymmetry in v- variables

    #note that the file 'src/stan_files/metaphase_3d_v2.stan' has an alternative formulation of the metaphase model using the log_sum_exp operator which should avoid any numerical overflow/underflow but this seems to run several times slower than this version, parameter estimates are similar
    m <- stan_model(stan_file) #slightly different way to call a stan model, this compiles the model and the next line does the sampling
    estimate_asym2 <- sampling(m,data=stan_input,
                               seed = 42,
                               chains = 4,
                               warmup = warm_up_iter,
                               iter = total_iter,
                               init_r=0.1,
                               pars = c("eta","f","xi","auxStates","P","aux"),
                               include=FALSE, #avoid saving the params listed above
                               control=list(adapt_delta=0.99,max_treedepth=12))
    #name2save <- paste("Asymmetric_full_vminus_",cell_name,"pair_",i,"_2.rds",sep="")
    #saveRDS(estimate_asym2,name2save)


    name2save2 <- paste("Asymmetric_full_vminus_",whole_name2,"_pair_",i,".rds",sep="")
    saveRDS(estimate_asym2,name2save2 )
    ############################################################################

    bridge2<-bridgesampling::bridge_sampler(estimate_asym2)
    print(bridge2$logml)
    #error_measures(bridge2)
    #bridge2=bridgesampling::bridge_sampler(estimate_asym2,method = "warp3")
    Marginal_likelihoods_bridge[j,3]=bridge2$logml
    # }

    stan_file= here::here('/Users/constandinakoki/Documents/GitHub/MetaAnaDynamics/STAN_files/Asymmetry/metaphase_asymmetric_vplus_full_likelihood_reparametrised_missing.stan') #Asymmetry in v+ variables         #note that the file 'src/stan_files/metaphase_3d_v2.stan' has an alternative formulation of the metaphase model using the log_sum_exp operator which should avoid any numerical overflow/underflow but this seems to run several times slower than this version, parameter estimates are similar

    m <- stan_model(stan_file) #slightly different way to call a stan model, this compiles the model and the next line does the sampling
    estimate_asym3 <- sampling(m,
                               data=stan_input,
                               seed = 42,
                               chains = 4,
                               warmup = warm_up_iter,
                               iter = total_iter,
                               init_r=0.1,
                               pars = c("eta","f","xi","auxStates","P","aux"),
                               include=FALSE, #avoid saving the params listed above
                               control=list(adapt_delta=0.99,
                                            max_treedepth=12))
    # name2save <- paste("Asymmetric_full_vplus_",cell_name,"pair_",i,"_2.rds",sep="")
    #saveRDS(estimate_asym3,name2save)

    name2save2 <- paste("Asymmetric_full_vplus_",whole_name2,"_pair_",i,".rds",sep="")
    saveRDS(estimate_asym3,name2save2 )
    ############################################################################


    ############################################################################

    bridge3=bridgesampling::bridge_sampler(estimate_asym3)
    print(bridge3$logml)
    #bridge3=bridgesampling::bridge_sampler(estimate_asym3,method = "warp3")
    Marginal_likelihoods_bridge[j,4]=bridge3$logml
    stan_file= here::here('/Users/constandinakoki/Documents/GitHub/MetaAnaDynamics/STAN_files/Asymmetry/metaphase_asymmetric_vminusvplus_full_likelihood_reparametrised_missing.stan') #Asymmetry in v-,v+ variables


    #note that the file 'src/stan_files/metaphase_3d_v2.stan' has an alternative formulation of the metaphase model using the log_sum_exp operator which should avoid any numerical overflow/underflow but this seems to run several times slower than this version, parameter estimates are similar
    m <- stan_model(stan_file) #slightly different way to call a stan model, this compiles the model and the next line does the sampling
    estimate_asym4 <- sampling(m,
                               data=stan_input,
                               seed = 42,
                               chains = 4,
                               warmup = warm_up_iter,
                               iter = total_iter,
                               init_r=0.1,
                               pars = c("eta","f","xi","auxStates","P","aux"),
                               include=FALSE, #avoid saving the params listed above
                               control=list(adapt_delta=0.99,max_treedepth=12))
    # name2save <- paste("Asymmetric_full_vminus_vplus_",cell_name,"pair_",i,"_2.rds",sep="")
    #saveRDS(estimate_asym4,name2save)

    name2save2 <- paste("Asymmetric_full_vminus_vplus_",whole_name2,"_pair_",i,".rds",sep="")
    saveRDS(estimate_asym4,name2save2 )



    # if(max(summary(estimate_asym4)$summary[c(1:17),'Rhat'])<1.1){
    bridge4<-bridgesampling::bridge_sampler(estimate_asym4)
    print(bridge4$logml)
    #bridge4=bridgesampling::bridge_sampler(estimate_asym4,method = "warp3")
    Marginal_likelihoods_bridge[j,5]=bridge4$logml
    #    loo4=loo(estimate_asym4)#}

    ###Symmetric metaphase model
    stan_file= here::here('/Users/constandinakoki/Documents/GitHub/MetaAnaDynamics/STAN_files/Asymmetry/metaphase_symmetric_full_likelihood_missing.stan')#Symmetric with full likelihood for using bridge sampling

    #note that the file 'src/stan_files/metaphase_3d_v2.stan' has an alternative formulation of the metaphase model using the log_sum_exp operator which should avoid any numerical overflow/underflow but this seems to run several times slower than this version, parameter estimates are similar
    m <- stan_model(stan_file) #slightly different way to call a stan model, this compiles the model and the next line does the sampling
    symmetric <- sampling(m,data=stan_input,
                          seed = 42,
                          chains = 4,
                          warmup = warm_up_iter,
                          iter = total_iter,
                          pars = c("eta","f","xi","auxStates","P","aux"),
                          include=FALSE, #avoid saving the params listed above
                          control=list(adapt_delta=0.99,max_treedepth=12))
    #name2save <- paste("Symmetric_",cell_name,"pair_",i,"_2.rds",sep="")
    #ÍsaveRDS(symmetric,name2save)

    name2save2 <- paste("Symmetric_",whole_name2,"_pair_",i,".rds",sep="")
    saveRDS(symmetric,name2save2 )

    # if(max(summary(symmetric)$summary[c(1:17),'Rhat'])<1.1){
    bridge5=bridgesampling::bridge_sampler(symmetric)
    print(bridge5$logml)
    Marginal_likelihoods_bridge[j,6]=bridge5$logml

    Bayes_Factors_bridge[1+(j-1)*5,]=c(i,bridge1$logml,bridgesampling::bf(bridge1,bridge2)[1],bridgesampling::bf(bridge1,bridge3)[1],bridgesampling::bf(bridge1,bridge4)[1],bridgesampling::bf(bridge1,bridge5)[1], post_prob(bridge1,bridge2,bridge3,bridge4,bridge5)[[1]])
    Bayes_Factors_bridge[2+(j-1)*5,]=c(i,bridgesampling::bf(bridge2,bridge1)[1],bridge2$logml,bridgesampling::bf(bridge2,bridge3)[1],bridgesampling::bf(bridge2,bridge4)[1],bridgesampling::bf(bridge2,bridge5)[1], post_prob(bridge1,bridge2,bridge3,bridge4,bridge5)[[2]])
    Bayes_Factors_bridge[3+(j-1)*5,]=c(i,bridgesampling::bf(bridge3,bridge1)[1],bridgesampling::bf(bridge3,bridge2)[1],bridge3$logml,bridgesampling::bf(bridge3,bridge4)[1],bridgesampling::bf(bridge3,bridge5)[1], post_prob(bridge1,bridge2,bridge3,bridge4,bridge5)[[3]])
    Bayes_Factors_bridge[4+(j-1)*5,]=c(i,bridgesampling::bf(bridge4,bridge1)[1],bridgesampling::bf(bridge4,bridge2)[1],bridgesampling::bf(bridge4,bridge3)[1],bridge4$logml,bridgesampling::bf(bridge4,bridge5)[1], post_prob(bridge1,bridge2,bridge3,bridge4,bridge5)[[4]])
    Bayes_Factors_bridge[5+(j-1)*5,]=c(i,bridgesampling::bf(bridge5,bridge1)[1],bridgesampling::bf(bridge5,bridge2)[1],bridgesampling::bf(bridge5,bridge3)[1],bridgesampling::bf(bridge5,bridge4)[1],bridge5$logml, post_prob(bridge1,bridge2,bridge3,bridge4,bridge5)[[5]])

    print(i)
    j=j+1

    write.csv( Marginal_likelihoods_bridge,here::here(paste0("Marginal_likelihoods_bridge_",whole_name2,".csv")), row.names = FALSE)
    write.csv(Bayes_Factors_bridge,here::here(paste0("Bayes_factors_bridge",whole_name2,".csv")), row.names = FALSE)

  }
}
