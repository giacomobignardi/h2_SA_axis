#Author: Giacomo Bignardi
#Description: load average HCP matrices and compute group-level functional connectomes (FC) averages (templates)
#Program: load_FC ------------------------------------------------------------------------------------------------------------------------------
# load packages
library(tidyverse)
library(readr)
library(psych)

# clear workspace
rm(list = ls())

# set Open Access working directories
wd_oa = getwd()
wd_oa_scripts = "10_github"
# wdNOA_ImageOutput = "05_Figures"

# set not Open Access working directories
wd_noa_output = paste0(substr(
  getwd(),
  0,
  nchar(getwd())-nchar("10_github")-1
),  "/11_noa_output")
wd_noa_data = paste0(substr(
  getwd(),
  0,
  nchar(getwd())-nchar("10_github")-1
),"/03_rawData/noa")
wd_noa_rsfmri = paste0(substr(
  getwd(),
  0,
  nchar(getwd())-nchar("10_github")-1
),"/03_rawData/private/rsfMRIavg")

# load
# inclusion index (based on subject inclusion in Valk et al., 2022 [Nat. Comms])
Valk_sub = read_csv(sprintf("%s/subj_inclusion.csv",wd_noa_data), col_names = F) %>% rename(Subject = X1)

# Load
FC_name = list.files(sprintf("%s/day_avg", wd_noa_rsfmri))
FC_name_d1 = list.files(sprintf("%s/day1", wd_noa_rsfmri))
FC_name_d2 = list.files(sprintf("%s/day2", wd_noa_rsfmri))

# Compute mean FC for template
FC_m = matrix(0, 419, 419)
isna = 0
FC_name_index = FC_name[(substr(FC_name,1,6) %in% Valk_sub$Subject)]

# FC:AVERAGE####
# create a long dataframe with the FC matrix per subject
for (i in 1:(length(FC_name_index))){
  FC_i = read_csv(sprintf("%s/day_avg/%s", wd_noa_rsfmri,FC_name_index[i]), col_names = F)
  if(is.na(FC_i$X1[1])){
    isna = isna +1
  }else{
    FC_i[FC_i < 0] = 0
    FC_m =  FC_m + FC_i
  }
}

# compute the template mean
FC_m_f = data.frame(FC_m/(length(FC_name_index)-isna))

# FC:AVERAGE DAY1####
# Compute mean FC for template at day 1
FC_m_d1 = matrix(0, 419, 419)
isna = 0
FC_name_index_d1 = FC_name_d1[(substr(FC_name_d1,1,6) %in% Valk_sub$Subject)]

# create a long df with the FC matrix per subject
for (i in 1:(length(FC_name_index_d1))){
  FC_i = read_csv(sprintf("%s/day1/%s", wd_noa_rsfmri,FC_name_index_d1[i]), col_names = F)
  if(is.na(FC_i$X1[1])){
    isna = isna +1
  }else{
    FC_i[FC_i < 0] = 0
    FC_m_d1 =  FC_m_d1 + FC_i
  }
}

# compute the template mean
FC_m_d1_f = data.frame(FC_m_d1/(length(FC_name_index_d1)-isna))

# FC:AVERAGE DAY2####
# Compute mean FC for template at day 2
FC_m_d2 = matrix(0, 419, 419)
isna = 0
FC_name_index_d2 = FC_name_d2[(substr(FC_name_d2,1,6) %in% Valk_sub$Subject)]

# create a long df with the FC matrix per subject
for (i in 1:(length(FC_name_index_d2))){
  FC_i = read_csv(sprintf("%s/day2/%s", wd_noa_rsfmri,FC_name_index_d2[i]), col_names = F)
  if(is.na(FC_i$X1[1])){
    isna = isna +1
  }else{
    FC_i[FC_i < 0] = 0
    FC_m_d2 =  FC_m_d2 + FC_i
  }
}

# compute the template mean
FC_m_d2_f = data.frame(FC_m_d2/(length(FC_name_index_d2)-isna))
wdOA_output = "03_outputs/processedData"

# SAVE####
# save FC averages
write_csv(FC_m_f,sprintf("%s/%s/02_FC_m.csv",wdNOA,wdNOA_output))
write_csv(FC_m_d1_f,sprintf("%s/%s/02_FC_m_t1.csv",wdNOA,wdNOA_output))
write_csv(FC_m_d2_f,sprintf("%s/%s/02_FC_m_t2.csv",wdNOA,wdNOA_output))