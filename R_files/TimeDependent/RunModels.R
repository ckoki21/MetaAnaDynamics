rm(list=ls())
setwd("~/Temporal_models")
library(rstan)
library(bridgesampling)
library(dplyr) 
library(purrr)
library(stringr)
library(readr)
library(tidyr)
library(zoo)
library(utils)
library(writexl)


rstan::rstan_options(auto_write = TRUE) #tries to avoid recompiling stan code
rstan::rstan_options(javascript=FALSE)
options(mc.cores = parallel::detectCores())
source('R/helper_fns.R')
source('R/df_sisters_more.R')
source('R/RunAllTimeDependent.R')

K=Inf # 

######################### #New_untrimmed Nocodazole data ######################

#dt0=1.94 
# jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking002-kitjobset_240318_All_CFAR_dt2-AI_2023-04-20_Capture4_Nocodazole_1-94sec-frames_deconvolved_5days.ome.csv")
 
 
#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking004-kitjobset_240318_All_CFAR_dt2-AI_2023-04-21_Capture5_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
  
  
#  dt0=2.07 
# jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking005-kitjobset_240318_All_CFAR_dt2-AI_2023-04-21_Capture8_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
 
  
# dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking006-kitjobset_240318_All_CFAR_dt2-AI_2023-04-21_Capture9_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
 
#  dt0=2.07 
#jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking007-kitjobset_240318_All_CFAR_dt2-AI_2023-05-04_Capture11_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
  
#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking008-kitjobset_240318_All_CFAR_dt2-AI_2023-05-04_Capture12_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
  
 #   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking009-kitjobset_240318_All_CFAR_dt2-AI_2023-05-04_Capture1_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
  
#    dt0=2.07 
# jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking010-kitjobset_240318_All_CFAR_dt2-AI_2023-05-04_Capture6_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
  
#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking013-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture11_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
 
#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking014-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture1_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
  
#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking015-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture2_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
  
#   dt0=2.07 
 # jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking016-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture4bottom_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking017-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture4top_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
  
 #  dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking018-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture5_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
  
#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking019-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture10_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
  
#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking021-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture1_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
  
#   dt0=1.94 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking022-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture2bottom_Nocodazole_1-94sec-frames_deconvolved_5days.ome.csv")
  
#   dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking023-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture2top_Nocodazole_1-94sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
# jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking025-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture8_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

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
  

# dt0=1.93 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking036-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-03-17_Capture1_Nocodazole_1-93sec-frames_deconvolved_5days.ome.csv")

# dt0=2.09
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking037-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-03-17_Capture4_Nocodazole_2-09sec-frames_deconvolved_5days.ome.csv")
 
 
# dt0=2.09
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking038-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-03-24_Capture8_Nocodazole_2-09sec-frames_deconvolved_5days.ome.csv")
 
 
# dt0=2.07
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking039-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-03-24_Capture9_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
 

 #  dt0=2.07 #
 # jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking041-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-03-17_Capture3_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#  dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking044-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-03-22_Capture4_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
  
 #  dt0=2.07 
 # jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking045-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-03-24_Capture2_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
  
#  dt0=2.07 
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking046-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-04-20_Capture6_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
  
#   dt0=2.07 # difficult to define anaphase --- only one sister shown or none!
#  jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking047-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-05-10_Capture6_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
 
 #  dt0=2.07 #
 # jobset_str <- here::here("data/New_nocodazole_untrimmed_per_frame/kittracking050-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-03-24_Capture3_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")
  


 ##############DMSO UNTRIMMED CFAR DATA####################################
 
dt0=2.58 
jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking002-kitjobset_240303_All_CFAR_dt2-AI_2022-12-01_Capture15_2-58sec-frames.ome.csv")

# dt0=2.15
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking004-kitjobset_240303_All_CFAR_dt2-AI_2022-12-01_Capture2_2-15sec-frames.ome.csv")
 
 
#  dt0=2.15
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking006-kitjobset_240303_All_CFAR_dt2-AI_2022-12-01_Capture9_2-15sec-frames.ome.csv")
  
  
#  dt0=2.24
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking007-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture10_2-24sec-frames_metaphase_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking009-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture13_2-07sec-frames_5days.ome.csv")
 
#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking010-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture15_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking012-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture17_2-07sec-frames_5days.ome.csv")
 
#   dt0=2.15
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking013-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture4_2-15sec-frames_5days.ome.csv")
 

#   dt0=2.09
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking014-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture6_2-09sec-frames_5days.ome.csv")
 
 
 #  dt0=2.07
#  jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking015-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture8_2-07sec-frames_5days.ome.csv")
 
#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking016-kitjobset_240303_All_CFAR_dt2-AI_2023-02-17_Capture4_2-07sec-frames_5days.ome.csv")
 
#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking017-kitjobset_240303_All_CFAR_dt2-AI_2023-02-17_Capture5_2-07sec-frames_5days.ome.csv")
 
#  dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking018-kitjobset_240303_All_CFAR_dt2-AI_2023-02-22_Capture1_2-07sec-frames_5days.ome.csv")
 
#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking019-kitjobset_240303_All_CFAR_dt2-AI_2023-02-22_Capture4_2-07sec-frames_5days.ome.csv")
 
 #  dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking020-kitjobset_240303_All_CFAR_dt2-AI_2023-02-22_Capture5_2-07sec-frames_metaphase_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking022-kitjobset_240303_All_CFAR_dt2-AI_2023-02-23_Capture14_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking023-kitjobset_240303_All_CFAR_dt2-AI_2023-02-23_Capture15_2-07sec-frames_5days.ome.csv")
 
#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking024-kitjobset_240303_All_CFAR_dt2-AI_2023-02-23_Capture1_2-07sec-frames_5days.ome.csv")
 
#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking025-kitjobset_240303_All_CFAR_dt2-AI_2023-02-23_Capture4_2-07sec-frames_5days.ome.csv")
 
#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking027-kitjobset_240303_All_CFAR_dt2-AI_2023-03-17_Capture6_DMSO_2-07sec-frames_5days.ome.csv")
 
#   dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking028-kitjobset_240303_All_CFAR_dt2-AI_2023-03-22_Capture11_DMSO_2-07sec-frames_5days.ome.csv")
 

 #   dt0=2.07
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking030-kitjobset_240303_All_CFAR_dt2-AI_2023-03-24_Capture14_DMSO_2-07sec-frames_flowdec_5days.ome.csv")
 
#   dt0=2.97
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking031-kitjobset_240303_All_CFAR_dt2-AI_2023-04-20_Capture10_DMSO_2-07sec-frames_5days.ome.csv")
 
#    dt0=2.07
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking032-kitjobset_240303_All_CFAR_dt2-AI_2023-04-21_Capture16_DMSO_2-07sec-frames_5days.ome.csv")
 
 #  dt0=2.07
# jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking033-kitjobset_240303_All_CFAR_dt2-AI_2023-05-04_Capture18_DMSO_2-07sec-frames_5days.ome.csv")
 
#    dt0=2.07
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking034-kitjobset_240303_All_CFAR_dt2-AI_2023-05-10_Capture13_DMSO_2-07sec-frames_5days.ome.csv")
 
#   dt0=2.07
#jobset_str <- here::here("data/DMSO_data_CFAR_per_frame/kittracking035-kitjobset_240303_All_CFAR_dt2-AI_2023-05-11_Capture22_DMSO_2-07sec-frames_flowdec_5days.ome.csv")

print(jobset_str)

TimeDependentModelFunction3(jobset_str = jobset_str,
                           warm_up_iter = 5000, total_iter = 10000, 
                           warm_up_iter_anaphase = 3000, total_iter_anaphase = 5000,
                           dt = dt0, K=Inf,p_missing_final = 0.2, prop_missing_initial =0.3)