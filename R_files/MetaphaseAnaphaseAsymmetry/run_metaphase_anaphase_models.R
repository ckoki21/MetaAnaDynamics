
#  title: "Asymmetric metaphase anaphase model"
rm(list=ls())
setwd("~/Documents/GitHub/MetaAnaDynamics/R_files/MetaphaseAnaphaseAsymmetry/")
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
library(here)
library(bridgesampling)
library(loo)
library(knitr)
library(writexl)
library(truncnorm)
library(shinystan)


rstan::rstan_options(auto_write = TRUE) #tries to avoid recompiling stan code
rstan::rstan_options(javascript=FALSE)
options(mc.cores = parallel::detectCores()) #uses as many cores as you h

here::i_am("run_metaphase_anaphase.R")
source(here::here('helper_fns.R'))
source(here::here('df_sisters_more.R'))


K=Inf # Number of frames

#Load the data
######################### #New_untrimmed Nocodazole data ######################
dt0=1.94 
jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking002-kitjobset_240318_All_CFAR_dt2-AI_2023-04-20_Capture4_Nocodazole_1-94sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking004-kitjobset_240318_All_CFAR_dt2-AI_2023-04-21_Capture5_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

# dt0=2.07 
#jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking005-kitjobset_240318_All_CFAR_dt2-AI_2023-04-21_Capture8_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")


# dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking006-kitjobset_240318_All_CFAR_dt2-AI_2023-04-21_Capture9_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#dt0=2.07 
#jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking007-kitjobset_240318_All_CFAR_dt2-AI_2023-05-04_Capture11_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking008-kitjobset_240318_All_CFAR_dt2-AI_2023-05-04_Capture12_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking009-kitjobset_240318_All_CFAR_dt2-AI_2023-05-04_Capture1_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#    dt0=2.07 
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking010-kitjobset_240318_All_CFAR_dt2-AI_2023-05-04_Capture6_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking013-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture11_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#  dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking014-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture1_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking015-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture2_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking016-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture4bottom_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking017-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture4top_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#  dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking018-kitjobset_240318_All_CFAR_dt2-AI_2023-05-10_Capture5_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking019-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture10_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")


#   dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking021-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture1_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=1.94 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking022-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture2bottom_Nocodazole_1-94sec-frames_deconvolved_5days.ome.csv")

#  dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking023-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture2top_Nocodazole_1-94sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking025-kitjobset_240318_All_CFAR_dt2-AI_2023-05-11_Capture8_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#  dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking028-kitjobset_240318_All_CFAR_dt2-AI_2023-05-16_Capture5_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking029-kitjobset_240318_All_CFAR_dt2-AI_2023-05-16_Capture6bottom_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#    dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking030-kitjobset_240318_All_CFAR_dt2-AI_2023-05-16_Capture6top_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#   dt0=2.07 
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking031-kitjobset_240318_All_CFAR_dt2-AI_2023-05-16_Capture8_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

# dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking032-kitjobset_240318_All_CFAR_dt2-AI_2023-05-16_Capture9_Nocodazole_2-07sec-frames_deconvolved_5days.ome.csv")

#  dt0=1.93 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking036-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-03-17_Capture1_MC191_488Ndc80EGFP_330nM_Nocodazole_1-93sec-frames_deconvolved_5_days_sigma.ome.csv")

# dt0=2.09
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking037-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-03-17_Capture4_MC191_488Ndc80EGFP_330nM_Nocodazole_2-09sec-frames_deconvolved_5_days_sigma.ome.csv")

# dt0=2.09
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking038-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-03-24_Capture8_MC191_488Ndc80EGFP_330nM_Nocodazole_2-09sec-frames_deconvolved_5_days_sigma.ome.csv")

#dt0=2.07
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking039-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-03-24_Capture9_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved_5_days_sigma.ome.csv")

#  dt0=2.07 #
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking041-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-03-17_Capture3_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved_5_days_sigma.ome.csv")

# dt0=2.07 
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking044-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-03-22_Capture4_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved_5_days_sigma.ome.csv")

#   dt0=2.07 
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking045-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-03-24_Capture2_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved_5_days_sigma.ome.csv")

#  dt0=2.07 
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/New_nocodazole_untrimmed_per_frame/kittracking046-kitjobset_240318_Nocodazole_All_CFAR_dt2-AI_2023-04-20_Capture6_MC191_488Ndc80EGFP_330nM_Nocodazole_2-07sec-frames_deconvolved_5_days_sigma.ome.csv")


##############DMSO UNTRIMMED CFAR DATA####################################

#  dt0=2.58 
#jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking002-kitjobset_240303_All_CFAR_dt2-AI_2022-12-01_Capture15_2-58sec-frames.ome.csv")

#  dt0=2.15
#jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking004-kitjobset_240303_All_CFAR_dt2-AI_2022-12-01_Capture2_2-15sec-frames.ome.csv")

#dt0=2.58
#jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking006-kitjobset_240303_All_CFAR_dt2-AI_2022-12-01_Capture9_2-15sec-frames.ome.csv")


#  dt0=2.24
#jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking007-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture10_2-24sec-frames_metaphase_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking009-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture13_2-07sec-frames_5days.ome.csv")

#  dt0=2.07
#jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking010-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture15_2-07sec-frames_5days.ome.csv")


#  dt0=2.07
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking012-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture17_2-07sec-frames_5days.ome.csv")

#  dt0=2.15
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking013-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture4_2-15sec-frames_5days.ome.csv")

#   dt0=2.09
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking014-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture6_2-09sec-frames_5days.ome.csv")

# dt0=2.07
#jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking015-kitjobset_240303_All_CFAR_dt2-AI_2023-02-16_Capture8_2-07sec-frames_5days.ome.csv")

#  dt0=2.07
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking016-kitjobset_240303_All_CFAR_dt2-AI_2023-02-17_Capture4_2-07sec-frames_5days.ome.csv")

# dt0=2.07
#  jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking017-kitjobset_240303_All_CFAR_dt2-AI_2023-02-17_Capture5_2-07sec-frames_5days.ome.csv")

#  dt0=2.07
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking018-kitjobset_240303_All_CFAR_dt2-AI_2023-02-22_Capture1_2-07sec-frames_5days.ome.csv")

#   dt0=2.58
#jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking019-kitjobset_240303_All_CFAR_dt2-AI_2023-02-22_Capture4_2-58sec-frames_5days.ome.csv")

#   dt0=2.07
#jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking020-kitjobset_240303_All_CFAR_dt2-AI_2023-02-22_Capture5_2-07sec-frames_metaphase_5days.ome.csv")

#  dt0=2.07
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking022-kitjobset_240303_All_CFAR_dt2-AI_2023-02-23_Capture14_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking023-kitjobset_240303_All_CFAR_dt2-AI_2023-02-23_Capture15_2-07sec-frames_5days.ome.csv")

#  dt0=2.07
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking024-kitjobset_240303_All_CFAR_dt2-AI_2023-02-23_Capture1_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking025-kitjobset_240303_All_CFAR_dt2-AI_2023-02-23_Capture4_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking027-kitjobset_240303_All_CFAR_dt2-AI_2023-03-17_Capture6_DMSO_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking028-kitjobset_240303_All_CFAR_dt2-AI_2023-03-22_Capture11_DMSO_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
#jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking030-kitjobset_240303_All_CFAR_dt2-AI_2023-03-24_Capture14_DMSO_2-07sec-frames_flowdec_5days.ome.csv")

#  dt0=2.07
# jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking031-kitjobset_240303_All_CFAR_dt2-AI_2023-04-20_Capture10_DMSO_2-07sec-frames_5days.ome.csv")

#dt0=2.07
#jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking032-kitjobset_240303_All_CFAR_dt2-AI_2023-04-21_Capture16_DMSO_2-07sec-frames_5days.ome.csv")

# dt0=2.07
#jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking033-kitjobset_240303_All_CFAR_dt2-AI_2023-05-04_Capture18_DMSO_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
#jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking034-kitjobset_240303_All_CFAR_dt2-AI_2023-05-10_Capture13_DMSO_2-07sec-frames_5days.ome.csv")

#   dt0=2.07
#jobset_str <- here::here("~/Documents/GitHub/MetaAnaDynamics/data/DMSO_data_CFAR_per_frame/kittracking035-kitjobset_240303_All_CFAR_dt2-AI_2023-05-11_Capture22_DMSO_2-07sec-frames_flowdec_5days.ome.csv")


#Read the cell-trajectories. keep the ones with no more than 25% missing data
#Create various plots
#Bring the trajectory data into format for stan runs
proportional_missing = 1
whole_name=str_split(jobset_str,"per_frame/kittracking0")
split_name=str_split(whole_name[[1]][2],".ome.csv")
cell_name=paste(split_name[[1]][1], sep="")
split_name= str_split(cell_name,"-kitjobset")
cell_name=paste(split_name[[1]][1], sep="")

Missing_positions <- how_much_missing(jobset_str, max_mis_prop = proportional_missing)
Missing_positions2 = Missing_positions %>% select( SisterPairID,SisterID,proportionNaN, which_missing) %>% mutate(new = as.character(which_missing), .keep = "unused")
#write_xlsx(Missing_positions2,here::here(paste0("Missingpositions_",wl,".xlsx")))

whole_name2 = str_split(jobset_str,"per_frame/")[[1]][2]
whole_name2 = str_split(whole_name2,".ome.csv")[[1]][1]

all_data = process_jobset_keep_missing(jobset_str,K=K,max_missing=1) %>%
  filter(!is.na(SisterID))
length(unique(all_data$SisterPairID))

################################################################################

data_single_pair <- process_jobset_keep_missing(jobset_str,K=Inf, start_from = 0,max_missing = proportional_missing) %>%
  filter(!is.na(SisterID)) #omit unpaired KTs
length(unique(data_single_pair$SisterPairID))

yz_plot <- ggplot(data_single_pair, aes(x=Position_2, y=Position_3,color=factor(SisterPairID))) +
  geom_line(alpha = 0.4)+ 
  geom_point(alpha = 0.7, size = 0.3)+ 
  theme_bw() +
  theme(legend.position = "None",
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(angle=90)) +theme_bw()+
  labs(x="Y-axis",y="Z-axis")+scale_color_discrete(name = "Sister Pair ID")
# ggsave(here::here(gsub("","",paste0("Plots_new_data/xy_including_missing_",cell_name,split_name[[1]][2],".pdf"))),widt0h = 22, height = 20)

g <- ggplot(data_single_pair, aes(x=Time, y=Position_1,color=factor(SisterID))) +
  geom_line() +
  facet_wrap(.~SisterPairID) +
  theme_bw() +
  theme(legend.position = "None",
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(angle=90)) +
  labs(x="Time (s)",y="Position (um)")
g


g <- ggplot(data_single_pair, aes(x=Time, y=Position_1,color=factor(SisterID))) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None",
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(angle=90)) +
  labs(x="Time (s)",y="Position (um)")
g


p_missing_final = 0.2
prop_missing_initial =0.3
p_missing_interpolate = 0.25
data_single_pair_old = process_jobset_keep_missing(jobset_str,K=Inf, start_from = 0,max_missing = p_missing_interpolate) %>%
  filter(!is.na(SisterID)) #omit unpaired KTs
length(unique(data_single_pair_old$SisterPairID))


data_single_pair = filter_cell_data(jobset_str, max_anaphase_frames = 45, p_missing_final = p_missing_final, prop_missing = prop_missing_initial)
length(unique(data_single_pair$SisterPairID))


#  data_single_pair %>% group_by(Time,Frame,SisterPairID) %>% mutate(Ymean = mean(Position_2), Zmean = mean(Position_3))%>%filter(SisterPairID ==36) %>% ggplot(aes(x= Ymean, y=Zmean,color=Frame)) +
# #  geom_line(alpha = 0.4)+ 
#   geom_point(alpha = 0.7)+ 
#   theme_bw() +
#   theme(legend.position = "None",
#         strip.text.x = element_text(size = 8),
#         axis.text.x = element_text(angle=90)) +theme_bw()+
#   labs(x="Y-axis",y="Z-axis")+scale_color_discrete(name = "Sister Pair ID")
# 


# y_time_series <- ggplot(data_single_pair, aes(x= Time, y=Position_2,color=factor(SisterID))) + facet_wrap(~SisterPairID)+
#   geom_line()+ 
#   theme_bw() +
#   theme(legend.position = "None",
#         strip.text.x = element_text(size = 8),
#         axis.text.x = element_text(angle=90)) +theme_bw()+
#   labs(x="Time (s)",y="Y_axis")+scale_color_discrete(name = "Sister Pair ID")
# 
#  y_time_series
#  
#   yz_plot <- data_single_pair %>% group_by(Time,Frame,SisterPairID) %>% mutate(Ymean = mean(Position_2), Zmean = mean(Position_3)) %>% ggplot(aes(x= Ymean, y=Zmean,color=factor(SisterPairID))) +
#   geom_line(alpha = 0.4)+ 
#   geom_point(alpha = 0.7, size= 0.3)+ 
#   theme_bw() +
#   theme(legend.position = "None",
#         strip.text.x = element_text(size = 8),
#         axis.text.x = element_text(angle=90)) +theme_bw()+
#   labs(x="Y-axis",y="Z-axis")+scale_color_discrete(name = "Sister Pair ID")
# # ggsave(here::here(gsub("","",paste0("Plots_new_data/xy_0.2missing_",cell_name,split_name[[1]][2],".pdf"))),widt0h = 22, height = 20)
#  


g <- ggplot(data_single_pair, aes(x=Time, y=Position_1,color=factor(SisterID))) +
  geom_line() +
  facet_wrap(.~SisterPairID) +
  theme_bw() +
  theme(legend.position = "None",
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(angle=90)) +
  labs(x="Time (s)",y="Position (um)")
g

# Center each series around the mean of its SisterPairID
data_single_pair1 <- data_single_pair %>% filter(SisterID == 1)%>%
  group_by(SisterPairID) %>%
  mutate(Position_1_centered = Position_1 - mean(Position_1, na.rm = TRUE))

data_single_pair2 <- data_single_pair %>% filter(SisterID == 2)%>%
  group_by(SisterPairID) %>%
  mutate(Position_1_centered = Position_1 - mean(Position_1, na.rm = TRUE))

data_1 = data_single_pair %>% filter(SisterID == 1) 
data_2 = data_single_pair %>% filter(SisterID == 2) 


#Sister pairs with no missing data for this particular cell
Complete_sisters=unique(data_single_pair$SisterPairID)
df_test = cbind(data_1$SisterPairID, data_1$Position_1,data_2$Position_1, data_1$Frame, data_1$Time)
colnames(df_test)=c("SisterPairID","X1","X2","Frame","Time")
df_test = as.data.frame(df_test)
df_test = df_test %>% rowwise() %>% mutate(Mean_position = mean(c(X1, X2)))

#This makes the autocorellation based on the mean position of two sisters
df_acf <-df_test %>%
  split(.$SisterPairID) %>% 
  map(~acf(.$Mean_position , main = unique(.$SisterPairID), lag.max = 150, plot = F,na.action=na.pass))%>%
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
#ggsave(here::here(gsub("","",paste0("Autocorrelations_cell_",cell_name,split_name[[1]][2],"NAN_included.pdf"))))


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
g


#we consider a changepoint model of anaphase onset that ignores the dynamics in metaphase, as described in Armond et al 2019 bioRxiv, https://doi.org/10.1101/582379. This model looks at the distance between a kinetochore pair. In metaphase, it assumes this is constant, and in anaphase it assumes this distance grows linearly. 

#We will fit the changepoint model to the experimental data. We want to to find the minimum anaphase time point for all kt-s in the cell 
pairs_to_include = unique(data_single_pair$SisterPairID)
data_single_pair_interpolate <- process_jobset(jobset_str,K=K,max_missing=1) %>%
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
  stan_input = list(dt=dt0, T=K,
                    nTracks = 1,
                    y = y_interpolate[[which(pairIDs==i)]],
                    y_missing = y_interpolate_missing[[which(pairIDs==i)]]
  )
  changept_estimate <- stan(file=here::here('~/Documents/GitHub/MetaAnaDynamics/STAN_files/MetaphaseAnaphaseAsymmetry/anaphase_changepoint.stan'),
                            data=stan_input,
                            seed = 42,
                            chains =4,
                            warmup = 3000,
                            iter = 6000,
                            control=list(adapt_delta=0.95,max_treedepth=12))
  #Find the minimum anaphase time for all sisters
  
  #pars_to_plot = c("tau","alpha","beta","a","t_ana")
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
  print(min2tana)
  j<-j+1
}

min_tana = min_tana/j   #This is the mean anaphase time from all sisters (or we could take the median. We have noticed that there are some sisters
#with outliers anaphase times and hence minimum is not the right choice.)

#Another way to go is to replace the outlier with the median time and then find the minimum anaphase time. 
Mean_anaphase = data.frame("Mean_t_ana"= colMeans(df_t_ana),"SisterPairID" = as.numeric(names(df_t_ana))) #Dataframe with mean anaphase times for all sister pairs.
#Remove outliers. I define them as the values which are < than (median(anaphase_times)/2)
Mean_anaphase = Mean_anaphase %>% select(Mean_t_ana)%>%mutate(Mean_anaphase, Mean_t_ana = ifelse((Mean_t_ana< median(Mean_anaphase$Mean_t_ana)-100), median(Mean_anaphase$Mean_t_ana)-30,Mean_t_ana ))
#Maybe more accurate to go with median instead of minimum to reduce the effect of a possible outlier/false change point.
t_ana_input0 =  median(Mean_anaphase$Mean_t_ana)


# Extracting the "2-07sec" part
extracted <- str_extract(jobset_str, "\\d+-\\d+sec")
# Removing 'sec' and replacing '-' with '.'
dt0 <- as.numeric(str_replace(str_remove(extracted, "sec"), "-", "."))
print(dt0)
t_ana_input0 = dt0*rough_anaphase_frame(jobset_str) #Another way to find the rough anaphase time, instead of the change point model above. These two methods give similar results.
print(t_ana_input0)


ggplot(data_single_pair, aes(x=Time-t_ana_input0, y=Position_1,color=factor(SisterID))) +
  geom_line() +
  facet_wrap(.~SisterPairID) +
  geom_vline(aes(xintercept = 0+30), linetype = "dashed", color = "gray40")+
  geom_vline(aes(xintercept = 0+60), linetype = "dashed", color = "gray40")+
  geom_vline(aes(xintercept = 0+90), linetype = "dashed", color = "gray40")+
  geom_vline(aes(xintercept = 0+120), linetype = "dashed", color = "gray40")+
  geom_vline(aes(xintercept = 0+150), linetype = "dashed", color = "gray40")+
  geom_vline(aes(xintercept = 0+180), linetype = "dashed", color = "gray40")+
  geom_vline(aes(xintercept = 0+210), linetype = "dashed", color = "gray40")+
  geom_vline(aes(xintercept = 0), color = "gray20")+
  theme_bw() +
  theme(legend.position = "None",
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(angle=90)) +
  labs(x="Time (s)",y="Position (um)")

ggsave(here::here(gsub("","",paste0("Trajectories_t_ana_facet_wrap_",cell_name,split_name[[1]][2],".pdf"))))


#Run METAPHASE ANAPHASE MODEL


pairIDs <- unique(data_single_pair$SisterPairID)
y = prepare_for_stan_format(data_single_pair)
y_missing = purrr::map(y, is.na) %>% purrr::map(function(x) x[,1] | x[,2])
for (i in 1:length(y)) {
  y[[i]][y_missing[[i]],] <- 10^9
} #replace missing values high values--- see also in .stan file to treat missing (stan does not like NA values)ß
y_missing <- purrr::map(y_missing,as.integer)
K=max(data_single_pair$Frame)

divergent_chains_sisters_warm = matrix(NA,nrow=length(pairIDs) ,ncol=6)     
divergent_chains_sisters_no_warm = matrix(NA,nrow=length(pairIDs) ,ncol=6)
j=1
model=1
warm_up_iter = 7000
total_iter = 14000
for(i in pairIDs){
  trajectory_index <- which(pairIDs==i)
  start_end <- purrr::map_df(y_missing,get_start_end_of_nonmissing_data,.id="traj_index") %>%
    mutate(traj_index=as.integer(traj_index)) #use this info to only fit to the existing tracks
  cos_theta <- get_cos_theta(data_single_pair,pairIDs[trajectory_index])
  cos_theta[is.na(cos_theta)] <- 0 #replace NAs as stan does not deal with these, should not come into the fit due to T0 and T1
  T0 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T0);
  T1 = start_end %>% filter(traj_index==trajectory_index) %>% pull(T1);
  Frames=K
  sigma0 = c(0.0,0.5,0.5,0.0,0.0)
  stan_input = list(dt=dt0, Frames=K, nStates = 5,
                    y = y[[trajectory_index]],
                    y_missing = y_missing[[trajectory_index]],
                    t_ana_input=t_ana_input0,
                    sigma0 =sigma0, #initial state probabilities
                    T0 =T0, 
                    T1 = T1,cos_theta = cos_theta)
  stan_file= here::here('~/Documents/GitHub/MetaAnaDynamics/STAN_files/MetaphaseAnaphaseAsymmetry/asymmetric_metaphase_anaphase_reparametrised_missing_data.stan')
  m <- stan_model(stan_file) 
  asymmetric_metaphase_anaphase <- sampling(m,data=stan_input,
                                            seed = 25,
                                            chains = 4,
                                            warmup = warm_up_iter,
                                            iter = total_iter,
                                            #                          init_r = 0.001,
                                            pars = c("eta","f","xi"),
                                            include=FALSE, #avoid saving the params listed above
                                            control=list(adapt_delta=0.95,max_treedepth=12))
  name2save <- paste("Asymmetric_metaphase_anaphase_",whole_name2,"_pair_",i,".rds",sep="")
  saveRDS(asymmetric_metaphase_anaphase,name2save)
  print(i)
  print(j)
  j=j+1
}

