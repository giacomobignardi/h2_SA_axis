#Author: Giacomo Bignardi
#Program: Median Absolute Deviation (MAD) analysis------------------------------------------------------------------------------------------------------------------------------
#for simplicity we copy previous scripts and minimal adapt them to run MAD-based analyses
#clean working environment 
rm(list = ls())
library(tidyverse)
library(reshape2)
library(patchwork)
library(lavaan)

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

# load measurement error model
source(sprintf("%s/R/functions/meermo/SMEM.R", wd_oa))

# Microstructural intensity
MPmi_i  =  read_csv(sprintf("%s/HCP_S1200_MPC_400.csv",wd_noa_data)) %>% rename(Sub = "Var1")

# Geodesic distances
GD_i  =   read_csv(sprintf("%s/01_GD.csv",wd_noa_output))

# Functional gradient
FC_g1_i_d1 =  read_csv(sprintf("%s/03_GFC_i_d1.csv",wd_noa_output))
FC_g1_i_d2 =  read_csv(sprintf("%s/03_GFC_i_d2.csv",wd_noa_output))

# inclusion index 
Twin_sub =  read_csv(sprintf("%s/00_twin_ids.csv",wd_noa_output))

# load dataFrames
HCP  = read_csv(sprintf("%s/unrestricted_giaco_6_25_2021_3_50_21.csv",wd_noa_data))
HCP_restricited  = read_csv(sprintf("%s/RESTRICTED_giaco_8_13_2021_11_47_39.csv",wd_noa_data))
HCP_all = merge(HCP,HCP_restricited, by = "Subject", all = T)

# FULL SAMPLE####
# full sample (inclusion criteria Valk et al., 2022 Nat Comm and not twins)
HCP_twin = HCP_all %>% 
  select(Family_ID,Subject, Gender, Age_in_Yrs,ZygosityGT,FS_IntraCranial_Vol) %>% 
  filter(Subject %in% Twin_sub$Subject)%>% 
  mutate(Sib_ID = ifelse(duplicated(Family_ID),2,1), #create a twinship order
         ZygosityGT = fct_recode(ZygosityGT, "1" = "MZ", "2" = "DZ"))
# INDIVIDUALS####
# get individual values for mp
MPmi_i = MPmi_i%>%
  filter(Sub %in% Twin_sub$Subject) %>%
  pivot_longer(names_to = "Parcel", values_to = "t1w/t2w_mi", c(2:ncol(MPmi_i)))%>%mutate(Parcel = substr(Parcel,6,nchar(Parcel)))%>%
  rename(value = `t1w/t2w_mi`) %>%
  mutate(Parcel = as.numeric(Parcel)) %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  rename(MP = value)

# get individual values for gd
GD_i = GD_i %>%
  filter(Sub %in% Twin_sub$Subject) %>%
  group_by(Parcel) %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  rename(GD = value)

# bind dataframe to be able to fit a measurement error model later
FC_gradients = rbind(FC_g1_i_d1, FC_g1_i_d2)
FC_gradients = FC_gradients%>%
  filter(Sub %in% Twin_sub$Subject) %>% #filter for twins
  select(Sub,Parcel,`0`, Session)%>%
  rename(value = `0`) %>% 
  pivot_wider(names_from = Session, values_from = value) %>% 
  rename(G1_d1 = `1.5`, G1_d2 = `2.5`)

# put all of them togheter
SFs_i_list = list(MPmi_i, GD_i,FC_gradients)

# merge all data frames in list
SFs_i = SFs_i_list %>% reduce(full_join, by=c("Sub",'Parcel'))
SFs = merge(SFs_i, HCP_all %>% rename(Sub = Subject), by = ("Sub"), all = T)

# merge
SFB_i_list = list(SFs_i,HCP_twin %>% rename(Sub = Subject))
SF_twin_i_final = SFB_i_list %>% reduce(full_join, by=c('Sub'))

# compute variances for within level analysis
SF_twin_i_final_within = SF_twin_i_final %>% 
  group_by(Sub,Family_ID,Gender,Age_in_Yrs,ZygosityGT,Sib_ID,FS_IntraCranial_Vol) %>% 
  reframe(MP = mad(MP),
          GD = mad(GD),
          G1_d1 = mad(G1_d1),
          G1_d2 = mad(G1_d2)) 

# create a list of parcel to run models on later
parcel_list = unique(FC_gradients%>%select(Parcel))

SF_twin_i_final_within_res = umx::umx_residualize(c("G1_d1","G1_d2","MP","GD"), cov = c("Age_in_Yrs","Gender"), data = SF_twin_i_final_within)


# rename for sem function
SFB_i_final = SF_twin_i_final_within_res %>% 
  select(Family_ID,Sub,Sib_ID, ZygosityGT, Gender, Age_in_Yrs,G1_d1,G1_d2,MP,GD,FS_IntraCranial_Vol) %>% 
  select(-Sub) %>% 
  mutate(ICV = (FS_IntraCranial_Vol)) %>% 
  ungroup() %>% 
  #prepare for twin modeling
  pivot_wider(names_from = Sib_ID, values_from = Gender:ICV) %>% 
  rename(P1_r1_T1=G1_d1_1,
         P1_r1_T2=G1_d1_2,
         P1_r2_T1=G1_d2_1,   
         P1_r2_T2=G1_d2_2,   
         P2_T1=MP_1,
         P2_T2=MP_2,
         P3_T1=GD_1,
         P3_T2=GD_2,
  ) %>% 
  as.data.frame() %>% 
  mutate(age = Age_in_Yrs_1,#note that the age of the tiwn is the same
         sex = Gender_1,#note that we only have same-sex twins
         zyg = ZygosityGT) %>% 
  select(-c(Family_ID, Gender_1,Gender_2,Age_in_Yrs_1,Age_in_Yrs_2,ZygosityGT))

# MAD ####
# source CTD model specification
source(sprintf("%s/R/functions/lavaantwda/CFM_CPMEM_AE_2g3p.R", wd_oa))

# estimate h2 
# correlated factor solution (multivariate via Direct Symmetric Approach)
cfm_cpmem_AE_fit = sem(cfm_cpmem_AE_model,
                           data = SFB_i_final, #we do this as it is exogenous
                           group = "zyg",
                           missing = "ML",
                           group.label= c(1,2),
                       std.ov = T)
  
# extract summary statistics
cfm_cpmem_AE_sumy = summary(cfm_cpmem_AE_fit, ci = T, fit.measures = T)
fitmeasures = fitmeasures(cfm_cpmem_AE_fit)[c("cfi","rmsea")]

# calculate standardized coefficients
stand_AE = standardizedSolution(cfm_cpmem_AE_fit)

# genetic correlations (sanity check these shold be the same)
cfm_cpmem_AE_sumy$pe%>% filter(grepl("rA_P",label)) %>% distinct(label,.keep_all = T)
stand_AE %>% filter(grepl("covA",label)) %>% distinct(label,.keep_all = T)

# environemntal correlations
cfm_cpmem_AE_sumy$pe%>% filter(grepl("rE_P",label)) %>% distinct(label,.keep_all = T)
stand_AE %>% filter(grepl("covE",label)) %>% distinct(label,.keep_all = T)

# h^2 (note that standardized and not should differ, as one include age and sex in the denominator, the other don't)
cfm_cpmem_AE_sumy$pe%>% filter(grepl("h2",label)) %>% distinct(label,.keep_all = T)
stand_AE %>% filter(grepl("=~",op) & grepl("A",lhs) ) %>% distinct(lhs,.keep_all = T) %>% mutate(h2 = est.std^2)

# e^2
stand_AE %>% filter(grepl("=~",op) & grepl("E",lhs) ) %>% distinct(lhs,.keep_all = T) %>% mutate(h2 = est.std^2)

# mem
stand_AE %>% filter(grepl("f1",label))
stand_AE %>% filter(grepl("var_intra", label))

# MAD with ICV as additonal covariate ####
# two step residualization (ICV is correlated between twin within a pair, would need extra care to be included as a manifest covariate)
SF_twin_i_final_within_res = umx::umx_residualize(c("G1_d1","G1_d2","MP","GD"), cov = c("FS_IntraCranial_Vol"), data = SF_twin_i_final_within_res)

# rename for sem function
SFB_i_final_res = SF_twin_i_final_within_res %>% 
  select(Family_ID,Sub,Sib_ID, ZygosityGT, Gender, Age_in_Yrs,G1_d1,G1_d2,MP,GD,FS_IntraCranial_Vol) %>% 
  select(-Sub) %>% 
  mutate(ICV = (FS_IntraCranial_Vol)) %>% 
  ungroup() %>% 
  #prepare for twin modeling
  pivot_wider(names_from = Sib_ID, values_from = Gender:ICV) %>% 
  rename(P1_r1_T1=G1_d1_1,
         P1_r1_T2=G1_d1_2,
         P1_r2_T1=G1_d2_1,   
         P1_r2_T2=G1_d2_2,   
         P2_T1=MP_1,
         P2_T2=MP_2,
         P3_T1=GD_1,
         P3_T2=GD_2,
  ) %>% 
  as.data.frame() %>% 
  mutate(zyg = ZygosityGT) %>% 
  select(-c(Family_ID, Gender_1,Gender_2,Age_in_Yrs_1,Age_in_Yrs_2,ZygosityGT))

# estimate h2 
# correlated factor solution (multivariate via Direct Symmetric Approach)
cfm_cpmem_AE_fit_res = sem(cfm_cpmem_AE_model,
                       data = SFB_i_final_res,
                       group = "zyg",
                       missing = "ML",
                       group.label= c(1,2),
                       std.ov = T)

# extract summary statistics
cfm_cpmem_AE_sumy_res = summary(cfm_cpmem_AE_fit_res, ci = T, fit.measures = T)
fitmeasures_res = fitmeasures(cfm_cpmem_AE_fit_res)[c("cfi","rmsea")]

# calculate standardized coefficients
stand_AE_res= standardizedSolution(cfm_cpmem_AE_fit_res)

# genetic correlations (sanity check these shold be the same)
cfm_cpmem_AE_sumy_res$pe%>% filter(grepl("rA_P",label)) %>% distinct(label,.keep_all = T)
stand_AE_res%>% filter(grepl("covA",label)) %>% distinct(label,.keep_all = T)

# environemntal correlations
cfm_cpmem_AE_sumy_res$pe%>% filter(grepl("rE_P",label)) %>% distinct(label,.keep_all = T)
stand_AE_res %>% filter(grepl("covE",label)) %>% distinct(label,.keep_all = T)

# h^2 (note that standardized and not should differ, as one include age and sex in the denominator, the other don't)
cfm_cpmem_AE_sumy_res$pe%>% filter(grepl("h2",label)) %>% distinct(label,.keep_all = T)
stand_AE_res %>% filter(grepl("=~",op) & grepl("A",lhs) ) %>% distinct(lhs,.keep_all = T) %>% mutate(h2 = est.std^2)

# e^2
stand_AE_res %>% filter(grepl("=~",op) & grepl("E",lhs) ) %>% distinct(lhs,.keep_all = T) %>% mutate(h2 = est.std^2)

# mem
stand_AE_res %>% filter(grepl("f1",label))
stand_AE_res %>% filter(grepl("var_intra", label))
