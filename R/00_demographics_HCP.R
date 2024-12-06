#Author: Giacomo Bignardi
#Description: sample descriptive of the HCP
#------------------------------------------------------------------------------------------------------------------------------
# load packages
library(tidyverse)
library(readr)

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

# load
# inclusion index (based on subject inclusion in Valk et al., 2022 [Nat. Comms])
Valk_sub = read_csv(sprintf("%s/subj_inclusion.csv",wd_noa_data), col_names = F) %>% rename(Subject = X1)

# load data unrestricted and restricted HCP data
HCP  = read_csv(sprintf("%s/unrestricted_giaco_6_25_2021_3_50_21.csv",wd_noa_data))
HCP_restricited  = read_csv(sprintf("%s/RESTRICTED_giaco_8_13_2021_11_47_39.csv",wd_noa_data))
HCP_all = merge(HCP,HCP_restricited, by = "Subject", all = T)

# FULL SAMPLE####
# full sample (inclusion criteria Valk et al. - see above)
HCP_exc = HCP_all %>% filter(Subject %in% Valk_sub$Subject)

# total N of subjects
nrow(HCP_exc)

# gender
table(HCP_exc$Gender)

# age
HCP_exc %>% reframe(mean = mean(Age_in_Yrs), 
                    sd = sd(Age_in_Yrs),
                    min =  min(Age_in_Yrs),
                    max =  max(Age_in_Yrs))

# SAMPLE EXCLUDING TWINS####
# select only twins with genotyped ID
HCP_not_twin = HCP_exc %>% filter(ZygositySR == "NotTwin")

# total n of subjects
nrow(HCP_not_twin)

# gender
table(HCP_not_twin$Gender)

# age
HCP_not_twin %>% reframe(mean = mean(Age_in_Yrs), 
                    sd = sd(Age_in_Yrs),
                    min =  min(Age_in_Yrs),
                    max =  max(Age_in_Yrs))

## SAVE sample excluding twins####
write_csv(HCP_not_twin%>% select(Subject), sprintf("%s/00_nottwin_ids.csv",wd_noa_output))

## twin sample
# IDs of twins with mislabeled or unknown zygosity
# mislabeled DZ
HCP_notDZ_index = HCP_exc %>% 
  filter(ZygositySR=='NotMZ' & ZygosityGT=='MZ') %>% 
  pull(Subject)

# mislabeled DZ
HCP_uknownSR_index = HCP_exc %>% 
  filter(is.na(ZygositySR) & !is.na(ZygosityGT)) %>% 
  pull(Subject)


HCP_twin = HCP_exc %>% 
  filter(!is.na(ZygosityGT)) %>% 
  filter(!(Subject %in% c(HCP_notDZ_index,HCP_uknownSR_index)))
table(HCP_twin$Gender,HCP_twin$ZygosityGT)
HCP_twin %>% reframe(mean = mean(Age_in_Yrs),
                     sd = sd(Age_in_Yrs),
                     min =  min(Age_in_Yrs),
                     max =  max(Age_in_Yrs))

## SAVE twin sample####
write_csv(HCP_twin%>% select(Subject), sprintf("%s/00_twin_ids.csv",wd_noa_output))
