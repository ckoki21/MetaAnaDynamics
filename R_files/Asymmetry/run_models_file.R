rm(list=ls())
#setwd("/Users/constandinakoki/Library/CloudStorage/OneDrive-UniversityofWarwick/Asymmetric_KTs_folder_mac/CLUSTER/Asymmetric_KTs/")
setwd("~/Asymmetric_models")
library(rstan)
library(dplyr) 
library(purrr)
library(zoo) 
library(here) 
library(ggplot2) 
library(patchwork) 
library(bayesplot)
library(tidybayes)
library(stringr)
library(readr)
library(tidyr)
library(bridgesampling)
library(writexl)
library(truncnorm)
library(mcp)
library(shinystan)

options(mc.cores = parallel::detectCores()) #uses as many cores as you h
rstan::rstan_options(auto_write = TRUE) #tries to avoid recompiling stan code

here::i_am("asymmetric_metaphase_model.Rmd")

source(here::here('R/helper_fns.R'))
source(here::here('R/df_sisters_more.R'))
source(here::here('R/run_all_models_fns.R'))
source(here::here('R/run_all_models_missing_fns.R'))

dt0=2.05
#jobset_str <- here::here("data/OS_LLSM_200818_MC191_Untreated_2.04933s_per_frame/kittracking001-kitjobset_200825_DonaldDuck_auto_v125-OS_LLSM_200818_MC191_Untreated_capture10.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200818_MC191_Untreated_2.04933s_per_frame/kittracking003-kitjobset_200825_DonaldDuck_auto_v125-OS_LLSM_200818_MC191_Untreated_capture12.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200818_MC191_Untreated_2.04933s_per_frame/kittracking005-kitjobset_200825_DonaldDuck_auto_v125-OS_LLSM_200818_MC191_Untreated_capture2.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200818_MC191_Untreated_2.04933s_per_frame/kittracking006-kitjobset_200825_DonaldDuck_auto_v125-OS_LLSM_200818_MC191_Untreated_capture5.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200818_MC191_Untreated_2.04933s_per_frame/kittracking007-kitjobset_200825_DonaldDuck_auto_v125-OS_LLSM_200818_MC191_Untreated_capture6.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200818_MC191_Untreated_2.04933s_per_frame/kittracking008-kitjobset_200825_DonaldDuck_auto_v125-OS_LLSM_200818_MC191_Untreated_capture8.ome.csv")


#jobset_str <- here::here("data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame/kittracking004-kitjobset_200827_DonaldDuck_auto_v125-OS_LLSM_200820_MC191_Untreated_capture13.ome.csv")
##jobset_str <- here::here("data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame/kittracking005-kitjobset_200827_DonaldDuck_auto_v125-OS_LLSM_200820_MC191_Untreated_capture15.ome.csv") #>>>????
##jobset_str <- here::here("data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame/kittracking006-kitjobset_200827_DonaldDuck_auto_v125-OS_LLSM_200820_MC191_Untreated_capture16.ome.csv") #>>>>not working
#jobset_str <- here::here("data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame/kittracking007-kitjobset_200827_DonaldDuck_auto_v125-OS_LLSM_200820_MC191_Untreated_capture17.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame/kittracking009-kitjobset_200827_DonaldDuck_auto_v125-OS_LLSM_200820_MC191_Untreated_capture19.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame/kittracking011-kitjobset_200827_DonaldDuck_auto_v125-OS_LLSM_200820_MC191_Untreated_capture24.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame/kittracking012-kitjobset_200827_DonaldDuck_auto_v125-OS_LLSM_200820_MC191_Untreated_capture25.ome.csv")
##jobset_str <- here::here("data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame/kittracking014-kitjobset_200827_DonaldDuck_auto_v125-OS_LLSM_200820_MC191_Untreated_capture27.ome.csv") #>>>Too many missing data
##jobset_str <- here::here("data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame/kittracking015-kitjobset_200827_DonaldDuck_auto_v125-OS_LLSM_200820_MC191_Untreated_capture28.ome.csv") #>>>Too many missing data
#jobset_str <- here::here("data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame/kittracking017-kitjobset_200827_DonaldDuck_auto_v125-OS_LLSM_200820_MC191_Untreated_capture30.ome.csv")#--------------->missing data
#jobset_str <- here::here("data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame/kittracking019-kitjobset_200827_DonaldDuck_auto_v125-OS_LLSM_200820_MC191_Untreated_capture4.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame/kittracking020-kitjobset_200827_DonaldDuck_auto_v125-OS_LLSM_200820_MC191_Untreated_capture5.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame/kittracking021-kitjobset_200827_DonaldDuck_auto_v125-OS_LLSM_200820_MC191_Untreated_capture6.ome.csv") 
#jobset_str <- here::here("data/OS_LLSM_200820_MC191_Untreated_2.04933s_per_frame/kittracking023-kitjobset_200827_DonaldDuck_auto_v125-OS_LLSM_200820_MC191_Untreated_capture9.ome.csv")




#jobset_str <- here::here("data/OS_LLSM_200929_MC191_Untreated_2.04933s_per_frame/kittracking003-kitjobset_210821_DonaldDuck_auto_v313-OS_LLSM_200929_MC191kittracking015-kitjobset_210821_DonaldDuck_auto_v313-OS_LLSM_200929_MC191_Untreated_capture9_Untreated_capture17.ome.csv")
##?jobset_str <- here::here("data/OS_LLSM_200929_MC191_Untreated_2.04933s_per_frame/kittracking006-kitjobset_210821_DonaldDuck_auto_v313-OS_LLSM_200929_MC191_Untreated_capture21.ome.csv") 
#jobset_str <- here::here("data/OS_LLSM_200929_MC191_Untreated_2.04933s_per_frame/kittracking007-kitjobset_210821_DonaldDuck_auto_v313-OS_LLSM_200929_MC191_Untreated_capture25.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200929_MC191_Untreated_2.04933s_per_frame/kittracking008-kitjobset_210821_DonaldDuck_auto_v313-OS_LLSM_200929_MC191_Untreated_capture26.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200929_MC191_Untreated_2.04933s_per_frame/kittracking010-kitjobset_210821_DonaldDuck_auto_v313-OS_LLSM_200929_MC191_Untreated_capture28.ome.csv")
#?jobset_str <- here::here("data/OS_LLSM_200929_MC191_Untreated_2.04933s_per_frame/kittracking011-kitjobset_210821_DonaldDuck_auto_v313-OS_LLSM_200929_MC191_Untreated_capture3.ome.csv")
#?jobset_str <- here::here("data/OS_LLSM_200929_MC191_Untreated_2.04933s_per_frame/kittracking012-kitjobset_210821_DonaldDuck_auto_v313-OS_LLSM_200929_MC191_Untreated_capture4.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200929_MC191_Untreated_2.04933s_per_frame/kittracking015-kitjobset_210821_DonaldDuck_auto_v313-OS_LLSM_200929_MC191_Untreated_capture9.ome.csv")

#jobset_str <- here::here("data/OS_LLSM_200930_MC191_Untreated_2.04933s_per_frame/kittracking001-kitjobset_210818_DonaldDuck_auto_v312-OS_LLSM_200930_MC191_Untreated_capture11.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200930_MC191_Untreated_2.04933s_per_frame/kittracking002-kitjobset_210818_DonaldDuck_auto_v312-OS_LLSM_200930_MC191_Untreated_capture12.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200930_MC191_Untreated_2.04933s_per_frame/kittracking003-kitjobset_210818_DonaldDuck_auto_v312-OS_LLSM_200930_MC191_Untreated_capture14.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200930_MC191_Untreated_2.04933s_per_frame/kittracking004-kitjobset_210818_DonaldDuck_auto_v312-OS_LLSM_200930_MC191_Untreated_capture15.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200930_MC191_Untreated_2.04933s_per_frame/kittracking006-kitjobset_210818_DonaldDuck_auto_v312-OS_LLSM_200930_MC191_Untreated_capture18.ome.csv")
#jobset_str <- here::here("data/OS_LLSM_200930_MC191_Untreated_2.04933s_per_frame/kittracking007-kitjobset_210818_DonaldDuck_auto_v312-OS_LLSM_200930_MC191_Untreated_capture19.ome.csv") 
##jobset_str <- here::here("data/OS_LLSM_200930_MC191_Untreated_2.04933s_per_frame/kittracking008-kitjobset_210818_DonaldDuck_auto_v312-OS_LLSM_200930_MC191_Untreated_capture3.ome.csv") #only 1 trajectory
#jobset_str <- here::here("data/OS_LLSM_200930_MC191_Untreated_2.04933s_per_frame/kittracking010-kitjobset_210818_DonaldDuck_auto_v312-OS_LLSM_200930_MC191_Untreated_capture5.ome.csv")



#NEW DATA
#####################################Untreated##################################

#NEW DATA
##dt0=2.07 #NOT GOOD TRAJECTORY
##jobset_str <- here::here("data/New_data_per_frame/kittracking001-kitjobset_230404_AI_2023-02-23_Capture12_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv") #~0.92 on average missing data

##dt0=2.14 
##jobset_str <- here::here("data/New_data_per_frame/kittracking001-kitjobset_Al_2023-02-16_Capture10_MC191_488Ndc80EGFP_2-24sec-frames_metaphase_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking002-kitjobset_230404_AI_2023-02-23_Capture14_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")

##dt0=2.07  #BAD
##jobset_str <- here::here("data/New_data_per_frame/kittracking002-kitjobset_Al_2023-02-16_Capture12_MC191_488Ndc80EGFP_2-07sec-frames_metaphase_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking003-kitjobset_230404_AI_2023-02-23_Capture15_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking003-kitjobset_Al_2023-02-16_Capture13_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")

##dt0=2.07 #BAD
##jobset_str <- here::here("data/New_data_per_frame/kittracking004-kitjobset_230404_AI_2023-02-23_Capture7_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")

##dt0=2.07 #BAD
##jobset_str <- here::here("data/New_data_per_frame/kittracking004-kitjobset_Al_2023-02-16_Capture14_MC191_488Ndc80EGFP_2-07sec-frames_metaphase_deconvolved.ome.csv")

##dt0=2.07 #BAD
##jobset_str <- here::here("data/New_data_per_frame/kittracking005-kitjobset_230404_AI_2023-02-23_Capture9_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")


#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking005-kitjobset_Al_2023-02-16_Capture15_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking006-kitjobset_Al_2023-02-16_Capture16_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking007-kitjobset_Al_2023-02-16_Capture17_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.15
#jobset_str <- here::here("data/New_data_per_frame/kittracking008-kitjobset_Al_2023-02-16_Capture4_MC191_488Ndc80EGFP_2-15sec-frames_deconvolved.ome.csv")

#dt0=2.09
#jobset_str <- here::here("data/New_data_per_frame/kittracking009-kitjobset_Al_2023-02-16_Capture6_MC191_488Ndc80EGFP_2-09sec-frames_deconvolved.ome.csv")

##dt0=2.07 #BAD
##jobset_str <- here::here("data/New_data_per_frame/kittracking010-kitjobset_Al_2023-02-16_Capture7_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking011-kitjobset_Al_2023-02-16_Capture8_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking012-kitjobset_Al_2023-02-17_Capture4_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")

#########
#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking013-kitjobset_Al_2023-02-17_Capture5_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking014-kitjobset_Al_2023-02-22_Capture1_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking015-kitjobset_AI_2023-02-22_Capture4_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking016-kitjobset_AI_2023-02-22_Capture5_MC191_488Ndc80EGFP_2-07sec-frames_metaphase_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking017-kitjobset_AI_2023-02-23_Capture1_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking018-kitjobset_AI_2023-02-23_Capture4_MC191_488Ndc80EGFP_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking019-kitjobset_AI_2023-02-23_Capture5_MC191_488Ndc80EGFP_2-07sec-frames_metaphase_deconvolved.ome.csv")


#-------------------------#

#dt0=2.09
#jobset_str <- here::here("data/New_data_per_frame/kittracking003-kitjobset_230406_AI_2023-03-17_Capture4_MC191_488Ndc80EGFP_330nM_Nocodazole_2-09sec-frames_deconvolved.ome.csv")

#dt0 = 1.93
#jobset_str <- here::here("data/New_data_per_frame/kittracking001-kitjobset_230406_AI_2023-03-17_Capture1_MC191_488Ndc80EGFP_330nM_Nocodazole_1-93sec-frames_deconvolved.ome.csv")


##dt0 = 2.07 #BAD
##jobset_str <- here::here("data/New_data_per_frame/kittracking004-kitjobset_230406_AI_2023-03-17_Capture5_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved.ome.csv")

##dt0=2.07 #BAD
##jobset_str <- here::here("data/New_data_per_frame/kittracking008-kitjobset_230406_AI_2023-03-22_Capture12_MC191_488Ndc80EGFP_DMSO_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking007-kitjobset_230406_AI_2023-03-22_Capture11_MC191_488Ndc80EGFP_DMSO_2-07sec-frames_deconvolved.ome.csv")


#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking005-kitjobset_230406_AI_2023-03-17_Capture6_MC191_488Ndc80EGFP_DMSO_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking011-kitjobset_230406_AI_2023-03-24_Capture14_MC191_488Ndc80EGFP_DMSO_2-07sec-frames_deconvolved.ome.csv")


#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking010-kitjobset_230406_AI_2023-03-22_Capture4_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking016-kitjobset_230406_AI_2023-03-24_Capture9_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved.ome.csv")

##dt0=2.07 #BAD
##jobset_str <- here::here("data/New_data_per_frame/kittracking002-kitjobset_230406_AI_2023-03-17_Capture3_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved.ome.csv")

##dt0=2.07 #BAD
##jobset_str <- here::here("data/New_data_per_frame/kittracking006-kitjobset_230406_AI_20 23-03-17_Capture7_MC191_488Ndc80EGFP_DMSO_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.07 #BAD
#jobset_str <- here::here("data/New_data_per_frame/kittracking004-kitjobset_230406_AI_2023-03-17_Capture5_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved.ome.csv")

##dt0=2.07 #BAD
##jobset_str <- here::here("data/New_data_per_frame/kittracking009-kitjobset_230406_AI_2023-03-22_Capture3_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved.ome.csv")


##dt0=2.07 #BAD
##jobset_str <- here::here("data/New_data_per_frame/kittracking012-kitjobset_230406_AI_2023-03-24_Capture2_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved.ome.csv")


#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking015-kitjobset_230406_AI_2023-03-24_Capture8_MC191_488Ndc80EGFP_330nM_Nocodazole_2-09sec-frames_deconvolved.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_data_per_frame/kittracking014-kitjobset_230406_AI_2023-03-24_Capture7_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved.ome.csv")

#dt0=2.07 #BAD
#jobset_str <- here::here("data/New_data_per_frame/kittracking013-kitjobset_230406_AI_2023-03-24_Capture3_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved.ome.csv")


#New nocodazole data
#dt0=2.07
#jobset_str <- here::here("data/New_nocodazole_data_per_frame/kittracking001-kitjobset_230702_2023-05-16_Capture8_AI_2023-05-16_Capture8_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_nocodazole_data_per_frame/kittracking001-kitjobset_230621_2023-04-21_Capture9_AI_2023-04-21_Capture9_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_nocodazole_data_per_frame/kittracking001-kitjobset_230621_2023-04-21-Capture5_AI_2023-04-21_Capture5_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_nocodazole_data_per_frame/kittracking001-kitjobset_230621_2023-05-04_Capture11_AI_2023-05-04_Capture11_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_nocodazole_data_per_frame/kittracking001-kitjobset_230621_2023-05-04_Capture12_AI_2023-05-04_Capture12_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_nocodazole_data_per_frame/kittracking001-kitjobset_230621_2023-05-11_Capture1_AI_2023-05-11_Capture1_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#dt0=1.94
#jobset_str <- here::here("data/New_nocodazole_data_per_frame/kittracking001-kitjobset_230621_2023-05-11_Capture2_AI_2023-05-11_Capture2_MC191_488Ndc80EGFP_330nM_Nocodazole_1-94sec-frames_deconvolved_5days.ome.csv")

#dt0=2.07
#jobset_str <- here::here("data/New_nocodazole_data_per_frame/kittracking001-kitjobset_230621_2023-04-21_Capture8_AI_2023-04-21_Capture8_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")


#############################New data 2#########################################
#dt0=2.07
#jobset_str <- here::here("data/New_data2_per_frame/kittracking001-kitjobset_230615_2023_04_20_Capture10_AI_Capture10_DMSO_2-07sec-frames_5days.ome.csv")


#BAD
#dt0=2.07 #only 1 pair
#jobset_str <- here::here("data/New_data2_per_frame/kittracking002-kitjobset_230722_AI_2023-04-20_Capture3_Nocodazole_2-07sec-frames_t420-t660_5days.ome.csv")


#dt0=1.94 #short ???
#jobset_str <- here::here("data/New_data2_per_frame/kittracking003-kitjobset_230801_AI_2023-04-20_Capture4_Nocodazole_1-94sec-frames_t340-t500_5days.ome.csv")


# dt0=2.07 #is it short?
# jobset_str <- here::here("data/New_data2_per_frame/kittracking004-kitjobset_230722_AI_2023-04-20_Capture6_Nocodazole_2-07sec-frames_t150-t310_5days.ome.csv")


 #dt0=2.07
 #jobset_str <- here::here("data/New_data2_per_frame/kittracking005-kitjobset_230615_AI_2023-04-21_Capture16_DMSO_2-07sec-frames_5days.ome.csv")

# BAD 
# dt0=1.94 #only 3 pairs, short metaphase time
# jobset_str <- here::here("data/New_data2_per_frame/kittracking006-kitjobset_230801_AI_2023-04-21_Capture1_Nocodazole_1-94sec-frames_t450-t600_5days.ome.csv")


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking012-kitjobset_230802_AI_2023-05-04_Capture18_DMSO_2-07sec-frames_5days.ome.csv")


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking013-kitjobset_230802_AI_2023-05-04_Capture1_Nocodazole_2-07sec-frames_t180-t350_5days.ome.csv")


 #dt0=2.07
 #jobset_str <- here::here("data/New_data2_per_frame/kittracking014-kitjobset_230802_AI_2023-05-04_Capture6_Nocodazole_2-07sec-frames_t300-t500_5days.ome.csv")


#BAD
# dt0=2.07 #Not good trajectories
# jobset_str <- here::here("data/New_data2_per_frame/kittracking015-kitjobset_230802_AI_2023-05-04_Capture8_Nocodazole_2-07sec-frames_t150-t350_5days.ome.csv")

#BAD
# dt0=2.07 #Only 6 trajectories
# jobset_str <- here::here("data/New_data2_per_frame/kittracking016-kitjobset_230802_AI_2023-05-10_Capture10_Nocodazole_2-07sec-frames_t240-t420_5days.ome.csv")


 #dt0=2.07
 #jobset_str <- here::here("data/New_data2_per_frame/kittracking017-kitjobset_230803_AI_2023-05-10_Capture11_Nocodazole_2-07sec-frames_t200-t360_5days.ome.csv")


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking018-kitjobset_230802_AI_2023-05-10_Capture13_DMSO_2-07sec-frames_t250-t450_5days.ome.csv")


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking019-kitjobset_230803_AI_2023-05-10_Capture1_Nocodazole_2-07sec-frames_t250-t400_5days.ome.csv")


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking021-kitjobset_230803_AI_2023-05-10_Capture4_Nocodazole_2-07sec-frames_t150-t310_5days.ome.csv")

# BAD
# dt0=2.07 #only 6 trajectories
# jobset_str <- here::here("data/New_data2_per_frame/kittracking022-kitjobset_230803_AI_2023-05-10_Capture4_Nocodazole_2-07sec-frames_t250-t450_5days.ome.csv")


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking023-kitjobset_230803_AI_2023-05-10_Capture5_Nocodazole_2-07sec-frames_t90-t270_5days.ome.csv")


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking024-kitjobset_230728_AI_2023-05-10_Capture6_Nocodazole_2-07sec-frames_t350-t550_5days.ome.csv")

#BAD 
# dt0=2.07 #not one good trajectory
# jobset_str <- here::here("data/New_data2_per_frame/kittracking025-kitjobset_230803_AI_2023-05-10_Capture7_Nocodazole_2-07sec-frames_t250-t450_5days.ome.csv")


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking026-kitjobset_230803_AI_2023-05-11_Capture10_Nocodazole_2-07sec-frames_t1-t250_5days.ome.csv")


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking027-kitjobset_230801_AI_2023-05-11_Capture18_Nocodazole_2-07sec-frames_t100-t300_5days.ome.csv")


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking029-kitjobset_230728_AI_2023-05-11_Capture22_DMSO_2-07sec-frames_flowdec_t90-t290_5days.ome.csv")

#BAD  
# dt0=1.94 #only 8 trajectories
# jobset_str <- here::here("data/New_data2_per_frame/kittracking030-kitjobset_230803_AI_2023-05-11_Capture2_Nocodazole_1-94sec-frames_t250-t400_5days.ome.csv")

#BAD  
# dt0=2.07 #not even one good paired trajectory
# jobset_str <- here::here("data/New_data2_per_frame/kittracking032-kitjobset_230803_AI_2023-05-11_Capture3_Nocodazole_2-07sec-frames_t200-t390_5days.ome.csv")

#BAD  
# dt0=2.07 #only 6 trajectories
# jobset_str <- here::here("data/New_data2_per_frame/kittracking033-kitjobset_230803_AI_2023-05-11_Capture9_Nocodazole_2-07sec-frames_t250-t400_5days.ome.csv")


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking034-kitjobset_230803_AI_2023-05-11_Capture8_Nocodazole_2-07sec-frames_t200-t380_5days.ome.csv")

#BAD?  
# dt0=2.07 #not going to anaphase !!! not good oscillations
# jobset_str <- here::here("data/New_data2_per_frame/kittracking035-kitjobset_230728_AI_2023-05-16_Capture12_Nocodazole_2-07sec-frames_t150-t250_5days.ome.csv")

#BAD 
# dt0=2.07 #Not even one good paired trajectory
# jobset_str <- here::here("data/New_data2_per_frame/kittracking036-kitjobset_230728_AI_2023-05-16_Capture13_DMSO_2-07sec-frames_t400-t540_5days.ome.csv")


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking037-kitjobset_230803_AI_2023-05-16_Capture4_Nocodazole_2-07sec-frames_t200-t330_5days.ome.csv")

#BAD  
# dt0=2.07 #only 4 trajecctories
# jobset_str <- here::here("data/New_data2_per_frame/kittracking038-kitjobset_230803_AI_2023-05-16_Capture5_Nocodazole_2-07sec-frames_t200-t440_5days.ome.csv")


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking040-kitjobset_230801_AI_2023-05-16_Nocodazole_2-07sec-frames_t200-t350_5days.ome.csv")


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking041-kitjobset_230801_AI_2023-05-16_Capture7_Nocodazole_2-07sec-frames_t400-t540_5days.ome.csv")


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking043-kitjobset_230803_AI_2023-05-16_Capture9_Nocodazole_2-07sec-frames_t100-t330_5days.ome.csv")

#BAD  
# dt0=2.07 #only 6 trajectories
# jobset_str <- here::here("data/New_data2_per_frame/kittracking044-kitjobset_230803_AI_2023-05-17_Capture2_Nocodazole_2-07sec-frames_t350-t590_5days.ome.csv")

#BAD
# dt0=2.07 #not even one paired trajectory
# jobset_str <- here::here("data/New_data2_per_frame/kittracking045-kitjobset_230803_AI_2023-05-17_Capture4_Nocodazole_2-07sec-frames_t100-t300_5days.ome.csv")

#BAD  
# dt0=2.07 #only 1 good trajectory
# jobset_str <- here::here("data/New_data2_per_frame/kittracking046-kitjobset_230801_AI_2023-05-17_Capture6_Nocodazole_2-07sec-frames_t300-t500_5days.ome.csv")

###############################New_untrimmed DMSO data####################


# dt0=2.07
# jobset_str <- here::here("data/New_data2_per_frame/kittracking035-kitjobset_240303_All_CFAR_dt2-AI_2023-05-11_Capture22_DMSO_2-07sec-frames_5days.ome.csv")


###############################New_untrimmed Nocodazole data####################


# dt0=2.07 #only one trajectory
# jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking001-kitjobset_240318_All_CFAR_dt2-AI_2023-04-20_Capture3_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

# dt0=1.94 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking002-kitjobset_240318_All_CFAR_dt2-AI_2023-04-20_Capture4_Nocodazole_1-94sec-frames_deconvolved_5days.ome.csv")

#  dt0=2.07 #only 8 trajectories
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking003-kitjobset_240318_All_CFAR_dt2-AI_2023-04-21_Capture1_Nocodazole_1-94sec-frames_deconvolved_5days.ome.csv")


#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking004-kitjobset_240318_All_CFAR_dt2-AI_2023-04-21_Capture5_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")


#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking005-kitjobset_240318_All_CFAR_dt2-AI_2023-04-21_Capture8_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")


# dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking006-kitjobset_240318_All_CFAR_dt2-AI_2023-04-21_Capture9_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#    dt0=2.07 
# jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking007-kitjobset_240318_All_CFAR_dt2-AI_2023-05-04_Capture11_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking008-kitjobset_240318_All_CFAR_dt2-AI_2023-05-04_Capture12_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#    dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking009-kitjobset_240318_All_CFAR_dt2-AI_2023-05-04_Capture1_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#    dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking010-kitjobset_240318_All_CFAR_dt2-AI_2023-05-04_Capture6_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#    dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking011-kitjobset_240318_All_CFAR_dt2-AI_2023-05-04_Capture8_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 #too many missing data -- only 8 trajectories
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking012-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture10_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking013-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture11_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")


#  dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking014-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture1_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking015-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture2_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking016-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture4bottom_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking017-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture4top_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#  dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking018-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture5_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking019-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture10_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 #only 9 trajectories
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking020-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture18_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking021-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture1_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=1.94 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking022-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture2bottom_Nocodazole_1-94sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking023-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture2top_Nocodazole_1-94sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 #too many missing data
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking024-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture3_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
# jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking025-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture8_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 #only 5 trajectories
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking026-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture9_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 #too many missing data
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking027-kitjobset_240318_All_CFAR_dt2-AI_2023-05-16_Capture4_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#  dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking028-kitjobset_240318_All_CFAR_dt2-AI_2023-05-16_Capture5_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#     dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking029-kitjobset_240318_All_CFAR_dt2-AI_2023-05-16_Capture6bottom_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#     dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking030-kitjobset_240318_All_CFAR_dt2-AI_2023-05-16_Capture6top_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking031-kitjobset_240318_All_CFAR_dt2-AI_2023-05-16_Capture8_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#    dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking032-kitjobset_240318_All_CFAR_dt2-AI_2023-05-16_Capture9_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#     dt0=2.07 #only 6 trajectories
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking033-kitjobset_240318_All_CFAR_dt2-AI_2023-05-17_Capture2_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#     dt0=2.07 #not even one good trajectory
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking034-kitjobset_240318_All_CFAR_dt2-AI_2023-05-17_Capture4_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#      dt0=2.07 #not even one good trajectory
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking035-kitjobset_240318_All_CFAR_dt2-AI_2023-05-17_Capture6_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")



##############DMSO UNTRIMMED CFAR DATA####################################

#   dt0=2.58 #short
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking001-kitjobset_240303_All_CFAR_dt2-AI_2022-12-01_Capture11_2-58sec-frames.ome.csv")

#    dt0=2.58 
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking002-kitjobset_240303_All_CFAR_dt2-AI_2022-12-01_Capture15_2-58sec-frames.ome.csv")

#  dt0=1.93 #short
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking003-kitjobset_240303_All_CFAR_dt2-AI_2022-12-01_Capture16_1-93sec-frames.ome.csv")

# dt0=2.15
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking004-kitjobset_240303_All_CFAR_dt2-AI_2022-12-01_Capture2_2-15sec-frames.ome.csv")


#   dt0=2.15 #no anaphase
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking005-kitjobset_240303_All_CFAR_dt2-AI_2022-12-01_Capture7_2-15sec-frames.ome.csv")


#  dt0=2.15
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking006-kitjobset_240303_All_CFAR_dt2-AI_2022-12-01_Capture9_2-15sec-frames.ome.csv")


#  dt0=2.24
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking007-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture10_2-24sec-frames_metaphase_5days.ome.csv")

#   dt0=2.07 #no good trajectories
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking008-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture12_2-07sec-frames_metaphase_5days.ome.csv")


#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking009-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture13_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking010-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture15_2-07sec-frames_5days.ome.csv")

#   dt0=2.07 #no good trajectories
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking011-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture16_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking012-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture17_2-07sec-frames_5days.ome.csv")

#   dt0=2.15
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking013-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture4_2-15sec-frames_5days.ome.csv")


#   dt0=2.09
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking014-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture6_2-09sec-frames_5days.ome.csv")


#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking015-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture8_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking016-kitjobset_240303_All_CFAR_dt2-AI_2023-02-17_Capture4_2-07sec-frames_5days.ome.csv")

 #  dt0=2.07
 #jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking017-kitjobset_240303_All_CFAR_dt2-AI_2023-02-17_Capture5_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking018-kitjobset_240303_All_CFAR_dt2-AI_2023-02-22_Capture1_2-07sec-frames_5days.ome.csv")

#   dt0=2.58
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking019-kitjobset_240303_All_CFAR_dt2-AI_2023-02-22_Capture4_2-58sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking020-kitjobset_240303_All_CFAR_dt2-AI_2023-02-22_Capture5_2-07sec-frames_metaphase_5days.ome.csv")

#   dt0=2.07 #no good trajectories
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking021-kitjobset_240303_All_CFAR_dt2-AI_2023-02-23_Capture12_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking022-kitjobset_240303_All_CFAR_dt2-AI_2023-02-23_Capture14_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking023-kitjobset_240303_All_CFAR_dt2-AI_2023-02-23_Capture15_2-07sec-frames_5days.ome.csv")

#  dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking024-kitjobset_240303_All_CFAR_dt2-AI_2023-02-23_Capture1_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking025-kitjobset_240303_All_CFAR_dt2-AI_2023-02-23_Capture4_2-07sec-frames_5days.ome.csv")

#   dt0=2.07 #only 7 trajectories
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking026-kitjobset_240303_All_CFAR_dt2-AI_2023-02-23_Capture5_2-07sec-frames_metaphase_5days.ome.csv")

 #  dt0=2.07
 #jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking027-kitjobset_240303_All_CFAR_dt2-AI_2023-03-17_Capture6_DMSO_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking028-kitjobset_240303_All_CFAR_dt2-AI_2023-03-22_Capture11_DMSO_2-07sec-frames_5days.ome.csv")

#   dt0=2.07 #no good trajectories
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking029-kitjobset_240303_All_CFAR_dt2-AI_2023-03-22_Capture12_DMSO_2-07sec-frames_5days.ome.csv")


#   dt0=2.07
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking030-kitjobset_240303_All_CFAR_dt2-AI_2023-03-24_Capture14_DMSO_2-07sec-frames_flowdec_5days.ome.csv")

#  dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking031-kitjobset_240303_All_CFAR_dt2-AI_2023-04-20_Capture10_DMSO_2-07sec-frames_5days.ome.csv")

#    dt0=2.07
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking032-kitjobset_240303_All_CFAR_dt2-AI_2023-04-21_Capture16_DMSO_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking033-kitjobset_240303_All_CFAR_dt2-AI_2023-05-04_Capture18_DMSO_2-07sec-frames_5days.ome.csv")

#    dt0=2.07
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking034-kitjobset_240303_All_CFAR_dt2-AI_2023-05-10_Capture13_DMSO_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking035-kitjobset_240303_All_CFAR_dt2-AI_2023-05-11_Capture22_DMSO_2-07sec-frames_flowdec_5days.ome.csv")

#    dt0=2.07 #not good trajectories
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking036-kitjobset_240303_All_CFAR_dt2-AI_2023-05-16_Capture13_DMSO_2-07sec-frames_5days.ome.csv")


run_all_models_missing(jobset_str = here::here(jobset_str), warm_up_iter = 8000, total_iter = 12000, warm_up_iter_anaphase = 3000, total_iter_anaphase = 5000, dt = dt0, K=Inf, p_missing_final = 0.2, prop_missing_initial = 0.3)
