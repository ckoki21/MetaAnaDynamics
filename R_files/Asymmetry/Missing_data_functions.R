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


# Get mean anaphase times and  to cut the data and then get the missing data statistics.
# If the trajectory has missing data after the anaphase time, then we are not affected by that
#  We also want to know if a trajectory has missing data on both sisters


get_statistics_missing_data  = function(jobset_str, full_cell_name,max_mis_prop, mean_anaphase_path = NA){
  #Mean anaphase_time is the dataframe with the anaphase times as derived by the change-point model
  #we will only use the Average(mean_anaphase_time$Mean_t_ana)
  if(is.na(mean_anaphase_path)){
    K = Inf
    print(K) 
  }else{
    mean_t_ana_df = read.csv(mean_anaphase_path, header =TRUE)
    K = (mean(mean_t_ana_df[,1])-30)/2.05 
    print(K)
    #Equivalent
    #    K = colMeans(mean_t_ana_df)[1]
    #    K = colMeans(mean_t_ana_df)["Mean_t_ana"]
  }
  #CAREFUL: T0 = 0 ALWAYS 
  df_missing = how_much_missing(jobset_str, K  , T0=0,max_mis_prop)
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
      # check 
      #prepeei na skefto pos na to kano. xano mia paratirisi otan vrisko lengths of missing
      #data >1 kai exo mia parapano otan exo mono mia paratirisi. auto symbaine giati exo kai to diff
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
