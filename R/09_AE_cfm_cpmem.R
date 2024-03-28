#Author: Giacomo Bignardi
#Description: estimate heritability of and genetic association between structure and function
#Program: Multivariate h2 twin measurement model------------------------------------------------------------------------------------------------------------------------------
#clean working environment 
rm(list = ls())
library(tidyverse)
library(reshape2)
library(patchwork)
library(ggseg)
library(ggseg3d)
library(ggsegSchaefer)
library(ggridges)
library(lavaan)

# clear workspace
rm(list = ls())

# set Open Access working directories
wdOA = getwd()
wdOA_scripts = "02_scripts"
wdNOA_ImageOutput = "05_Figures"


# set not Open Access working directories
wdNOA = getwd()
wdNOA_Data = "/01_input"
wdNOA_output = "/03_outputs/processedData"

wdNOA_rawData = paste0(substr(
  getwd(),
  0,
  nchar(getwd())-nchar("04_analysis_OA")-1
),"/03_rawData/private")

# Microstructural intensity
MPmi_i  =  read_csv(sprintf("%s/t1t2w/HCP_S1200_MPC_400.csv",wdNOA_rawData)) %>% rename(Sub = "Var1")

# Geodesic distances
GD_i  =   read_csv(sprintf("%s/%s/01_GD.csv",wdNOA,wdNOA_output))

# Functional gradient
FC_g1_i_d1 =  read_csv(sprintf("%s/%s/03_GFC_i_d1.csv", wdNOA,wdNOA_output))
FC_g1_i_d2 =  read_csv(sprintf("%s/%s/03_GFC_i_d2.csv", wdNOA,wdNOA_output))

# inclusion index 
Twin_sub =  read_csv(sprintf("%s/%s/00_twin_ids.csv", wdNOA,wdNOA_output))

# cortical types and Yeo functional network annotations
annotations = read_csv(sprintf("%s/%s/merged_YeoKongMesulamTypes.csv", wdNOA,wdNOA_Data))

# tidy yeo 7 annotations 
network_7_yeo = c()
for(i in 1:400){
  net_i = c()
  net_i = strsplit(annotations$label_Yeo_7, "_")[[i]][3]
  network_7_yeo = rbind(network_7_yeo,net_i)
}

# network_7_yeo = c(network_7_yeo, rep("Other", 19))
annotations$label_Yeo7_short = network_7_yeo[,1]
annotations = annotations %>% rename(Parcel = parcel_Yeo)

# inclusion index 
Twin_sub =  read_csv(sprintf("%s/%s/00_twin_ids.csv", wdNOA,wdNOA_output))

# load dataFrames
HCP  = read_csv(sprintf("%s/HCP/unrestricted_giaco_6_25_2021_3_50_21.csv",wdNOA_rawData))
HCP_restricited  = read_csv(sprintf("%s/HCP/RESTRICTED_giaco_8_13_2021_11_47_39.csv",wdNOA_rawData))
HCP_all = merge(HCP,HCP_restricited, by = "Subject", all = T)

# FULL SAMPLE####
# full sample (inclusion criteria Valk et al., 2022 Nat Comm and not twins)
HCP_twin = HCP_all %>% 
  select(Family_ID,Subject, Gender, Age_in_Yrs,ZygosityGT) %>% 
  filter(Subject %in% Twin_sub$Subject)%>% 
  mutate(Sib_ID = ifelse(duplicated(Family_ID),2,1), #create a twinship order
         ZygosityGT = fct_recode(ZygosityGT, "1" = "MZ", "2" = "DZ"))
#INDIVIDUALS####
#get individual values for mp
MPmi_i = MPmi_i%>%
  filter(Sub %in% Twin_sub$Subject) %>%
  pivot_longer(names_to = "Parcel", values_to = "t1w/t2w_mi", c(2:ncol(MPmi_i)))%>%mutate(Parcel = substr(Parcel,6,nchar(Parcel)))%>%
  rename(value = `t1w/t2w_mi`) %>%
  mutate(Parcel = as.numeric(Parcel)) %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  rename(MP = value)

#get individual values for gd
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
SFs_i_final_list = list(SFs_i, annotations)
SFs_i_final = SFs_i_final_list %>% reduce(full_join, by=c('Parcel'))%>% 
  mutate(label_Yeo7_short = factor(label_Yeo7_short, levels = rev(c("Default", 
                                                                    "Limbic", 
                                                                    "Cont",
                                                                    "SalVentAttn",
                                                                    "DorsAttn",
                                                                    "SomMot",
                                                                    "Vis"))))

SFs = merge(SFs_i_final, HCP_all %>% rename(Sub = Subject), by = ("Sub"), all = T)

# create a list of parcel to run models on later
parcel_list = unique(FC_gradients%>%select(Parcel))

# merge
SFB_i_list = list(SFs_i_final,HCP_twin %>% rename(Sub = Subject))
SF_twin_i_final = SFB_i_list %>% reduce(full_join, by=c('Sub'))

# create a list of parcel to run models on later
parcel_list = unique(FC_gradients%>%select(Parcel))

# rename for sem function
SFB_i_final = SF_twin_i_final %>% 
  select(Family_ID,Sub,Parcel,Sib_ID, ZygosityGT, Gender, Age_in_Yrs,G1_d1,G1_d2,MP, GD) %>% 
  select(-Sub) %>% 
  group_by(Parcel) %>% 
  ungroup() %>% 
  #prepare for twin modeling
  pivot_wider(names_from = Sib_ID, values_from = Gender:GD) %>% 
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
  mutate(age_T1 = Age_in_Yrs_1,#note that the age of the twin is the same
         sex_T1 = Gender_1,#note that we only have same-sex twins
         age_T2 = Age_in_Yrs_1,
         sex_T2 = Gender_1,
         zyg = ZygosityGT) %>% 
  select(-c(Family_ID, Gender_1,Gender_2,Age_in_Yrs_1,Age_in_Yrs_2,ZygosityGT))

# Twin multivariate model with measurement error ####
# source CTD model specification
source(sprintf("%s/%s/functions/lavaantwda/CFM_cpmem_AE_2g3p.R", wdOA,wdOA_scripts))

# estimate h2 across parcels
cfm_cpmem_AE_sumy = c()
for(i in parcel_list$Parcel){
  # reinitialize at the beginning of each loop
  cfm_cpmem_AE_fit_i = c()
  cfm_cpmem_AE_sumy_i = c()
  fitmeasures_i= c()
  SFB_i_final_i = c()
  
  # Including manifest covariates in multigroup SEM can make things complicated (e.g., unequal variance across groups, total phenotypic variance estimates)
  # so for twin modelin we apply a classic residualization procedures
  #  two step residualization for age and sex 
  SFB_i_final_i = umx::umx_residualize(c("P1_r1","P1_r2", "P2", "P3"), 
                                       cov = c("age","sex"), 
                                       suffixes=c("_T1","_T2"),
                                       data = SFB_i_final %>% filter(Parcel== i))
  
  # correlated factor solution (multivariate via Direct Symmetric Approach)
  cfm_cpmem_AE_fit_i = sem(cfm_cpmem_AE_model, # common pathway measurement error twin multivariate model
                          data = SFB_i_final_i, # fit data to each parcel i
                          group = "zyg",
                          missing = "ML",
                          std.ov = T, # standardize all observed variables first
                          group.label= c(1,2))
  
  # extract summary statistics
  cfm_cpmem_AE_sumy_i = summary(cfm_cpmem_AE_fit_i, ci = T)
  fitmeasures_i = fitmeasures(cfm_cpmem_AE_fit_i,c("cfi","rmsea"))

  # bind to previous estimates
  cfm_cpmem_AE_sumy = rbind(cfm_cpmem_AE_sumy,
                               cbind(cfm_cpmem_AE_sumy_i$pe %>% 
                                       mutate(Parcel = i,
                                              cfi = fitmeasures_i[1],
                                              rmsea = fitmeasures_i[2])
                               )
                            )
}

## Supplementary File 4 ####
# SF4A: parameters in strucutral models
SF4A = cfm_cpmem_AE_sumy %>%                                 
  filter((grepl("covA",label) | grepl("covE",label)| grepl("var",label)|grepl("f1",label)|grepl("mean",label))) %>% 
  group_by(Parcel) %>% 
           distinct(label, .keep_all = T) %>% 
  mutate(label = recode(label,
                        mean_P1_r1 = "mean_FCG11",
                        mean_P1_r2 = "mean_FCG12",
                        mean_P2 = "mean_T1w/T2wmi",
                        mean_P3 = "GD",
                        f1 = "λ",
                        var_intra_r1 = "σ1intra^2",
                        var_intra_r2 = "σ2intra^2",
                        varA_P1 = "σAinter^2",
                        varA_P2 = "σAmi^2",
                        varA_P3 = "σAGD^2",
                        varE_P1 = "σEinter^2",
                        varE_P2 = "σEmi^2",
                        varE_P3 = "σEGD^2",
                        covA_P12 = "σA_mig1",
                        covA_P13 = "σA_gdg1",
                        covA_P23 = "σA_migd",
                        covE_P12 = "σE_mig1",
                        covE_P13 = "σE_gdg1",
                        covE_P23 = "σE_migd",
                        )) %>% 
  filter(!label %in% c("varP1","varP2","varP3"))

SF4A_annot = merge(SF4A, annotations %>% select(Parcel,label_Yeo_7,label_Yeo7_short), by = "Parcel")%>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  select(-c(lhs,op,rhs,block,group))

# SF4B: genetic correlations between microstrucutre and functional gradient loadings
SF4B = cfm_cpmem_AE_sumy %>%                                 
  filter(grepl("rA_P12",label)) %>% 
  group_by(Parcel) %>% 
  distinct(label, .keep_all = T) %>% 
  mutate(label = recode(label,rA_P12 = "rA_mig1" ))

SF4B_annot = merge(SF4B, annotations %>% select(Parcel,label_Yeo_7,label_Yeo7_short), by = "Parcel")%>% 
  mutate(p_adjust = p.adjust(pvalue, method ="bonferroni"),
         across(where(is.numeric), round, 3)) %>%
  select(-c(lhs,op,rhs,block,group))

# SF4C: genetic correlations between geodesic distances and functional gradient loadings
SF4C = cfm_cpmem_AE_sumy %>%                                 
  filter(grepl("rA_P13",label)) %>% 
  group_by(Parcel) %>% 
  distinct(label, .keep_all = T) %>% 
  mutate(label = recode(label,rA_P12 = "rA_mig1" ))

SF4C_annot = merge(SF4C, annotations %>% select(Parcel,label_Yeo_7,label_Yeo7_short), by = "Parcel")%>% 
  mutate(p_adjust = p.adjust(pvalue, method ="bonferroni"),
         across(where(is.numeric), round, 3)) %>%
  select(-c(lhs,op,rhs,block,group))

## save####
write_csv(SF4A_annot, sprintf("%s/%s/09_SF4A.csv",wdNOA,wdNOA_output))
write_csv(SF4B_annot, sprintf("%s/%s/09_SF4B.csv",wdNOA,wdNOA_output))
write_csv(SF4C_annot, sprintf("%s/%s/09_SF4C.csv",wdNOA,wdNOA_output))

cfm_cpmem_AE_sumy_short =  cfm_cpmem_AE_sumy %>%                                 
  filter((grepl("rE",label) | grepl("rA",label) | grepl("h",label)) & !(grepl("var",label)|grepl("bh",label)))

#each parcel has 9 (6 cor and 3 h2) parameter estimates
nrow(cfm_cpmem_AE_sumy_short)/9

# retain only good fitting model (following Teeuw et al., 2021 Neuroimage)
cfs_fltr_fit = cfm_cpmem_AE_sumy_short %>% 
  filter(rmsea < 0.08, cfi > .90)

# how many with good fit
nrow(cfs_fltr_fit)/9

# exclude heywood cases
cfs_fltr = cfs_fltr_fit %>%   
  filter(est < 1 | est > -1)

# how many left
n_pacerl_model = nrow(cfs_fltr)/9

# merge with annotation to stratify across cortical networks
AE_fltr_annot = merge(cfs_fltr, annotations, by = "Parcel") %>%   mutate(label_Yeo7_short = factor(label_Yeo7_short, levels = rev(c("Default", 
                                                                                                                                    "Limbic", 
                                                                                                                                    "Cont",
                                                                                                                                    "SalVentAttn",
                                                                                                                                    "DorsAttn",
                                                                                                                                    "SomMot",
                                                                                                                                    "Vis"))))
# correct for multiple comparisions
# NOTE P1P2 = inter FCg1 <-> t1t2w
# NOTE P1P3 = inter FCg1 <-> GD
est_AE_annot = AE_fltr_annot %>% 
  group_by(label) %>% 
  mutate(p_adjust = p.adjust(pvalue, method ="bonferroni"))

## save####
write_csv(est_AE_annot, sprintf("%s/%s/09_CPMEM_output.csv",wdNOA,wdNOA_output))

# Heritability estimates
library(ggtext)

# average heritability
est_AE_annot %>% 
  filter(grepl("h2",label)) %>% 
  group_by(lhs) %>% 
  reframe(mean(est), sd(est))

# display percentage of significant genetic and environemntal correlations 
# amongst all the possible FCg1 <-> t1t2w (that is twice the number of parcels)
(est_AE_annot %>% 
  filter(p_adjust<.05 & label == "rA_P12" & est <0) %>% nrow() +
est_AE_annot %>% 
  filter(p_adjust<.05 & label == "rA_P12" & est >0) %>% nrow() +
est_AE_annot %>% 
  filter(p_adjust<.05 & label == "rE_P12" & est <0) %>% nrow() +
est_AE_annot %>% 
  filter(p_adjust<.05 & label == "rE_P12" & est >0) %>% nrow()) /
  (n_pacerl_model*2)

# focus on genetic correlations
(est_AE_annot %>% 
    filter(p_adjust<.05 & label == "rA_P12") %>% nrow())/
  (n_pacerl_model)


# display percentage of significant genetic and environemntal correlations 
# amongst all the possible FCg1 <-> GD (that is twice the number of parcels)
(est_AE_annot %>% 
    filter(p_adjust<.05 & label == "rA_P13" & est <0) %>% nrow() +
    est_AE_annot %>% 
    filter(p_adjust<.05 & label == "rA_P13" & est >0) %>% nrow()) /
  (n_pacerl_model)

(est_AE_annot %>% 
    filter(p_adjust<.05 & label == "rE_P13" & est <0) %>% nrow() +
    est_AE_annot %>% 
    filter(p_adjust<.05 & label == "rE_P13" & est >0) %>% nrow()) /
  (n_pacerl_model)


# focus on genetic correlations
(est_AE_annot %>% 
    filter(p_adjust<.05 & label == "rA_P13" & est <0) %>% nrow())/
  (n_pacerl_model)

(est_AE_annot %>% 
    filter(p_adjust<.05 & label == "rA_P13" & est >0) %>% nrow()) /
  (n_pacerl_model)


# average magnitude of genetic correlations
est_AE_annot %>% 
    filter(p_adjust<.05 & label == "rA_P13" & est <0) %>% 
  summarise(mean= psych::fisherz2r(mean(psych::fisherz(est))),
            sd = psych::fisherz2r(sd(psych::fisherz(est))))

est_AE_annot %>% 
  filter(p_adjust<.05 & label == "rA_P13" & est >0) %>% 
  summarise(mean= psych::fisherz2r(mean(psych::fisherz(est))),
            sd = psych::fisherz2r(sd(psych::fisherz(est))))


# focus on enviromental correlations
(est_AE_annot %>% 
    filter(p_adjust<.05 & label == "rE_P13" & est <0) %>% 
  nrow() +
  est_AE_annot %>% 
    filter(p_adjust<.05 & label == "rE_P13" & est >0) %>% 
  nrow()) /
  (n_pacerl_model)


# plot associations on the cortex
#Figure 3B-C####
# template surface 
associations2cortex = est_AE_annot %>% 
  filter(label == "rA_P13" & p_adjust<.05) %>% 
  ungroup() %>% 
  select(Parcel, label_Yeo_7, est)

min(associations2cortex$est, na.rm = T)
max(associations2cortex$est, na.rm = T)
# adapt to plot
schaefer7_association_GD = schaefer7_400_3d %>% #
  filter(hemi == "left")%>%
  mutate(ggseg_3d = map(ggseg_3d, ~ .x %>% 
                          full_join(associations2cortex%>%
                                      filter(Parcel < 201) %>% #left hemisphere
                                      rename(region = label_Yeo_7) %>% 
                                      select(region,est),
                                    by = "region")))
# plot parcelations using ggseg3d
association_GD_left_medial = 
  ggseg3d(
    atlas = schaefer7_association_GD,
    colour = "est", text = "est",
    palette =  c( "#4d48a0"= min(associations2cortex$est, na.rm = T), "white"  = 0, "#a65a53" =max(associations2cortex$est, na.rm = T)),
    surface = "inflated",
    hemisphere = "left")%>% 
  remove_axes()%>%
  pan_camera("left medial")

association_GD_left_lateral = 
  ggseg3d(
    atlas = schaefer7_association_GD,
    colour = "est", text = "est",
    palette =  c( "#4d48a0"= min(associations2cortex$est, na.rm = T), "white"  = 0, "#a65a53" =max(associations2cortex$est, na.rm = T)),
    surface = "inflated",
    hemisphere = "left")%>% 
  remove_axes()%>%
  pan_camera("left lateral")

#plot association with microstructure
associations2cortex_mi = est_AE_annot %>% 
filter(label == "rA_P12" & p_adjust<.05) %>% 
  ungroup() %>% 
  select(Parcel, label_Yeo_7, est)

est_AE_annot %>% 
  filter(label == "rA_P12") %>% 
  pull(est)


min(est_AE_annot %>% 
      filter(label == "rA_P12") %>% 
      pull(est)
    , na.rm = T)
max(est_AE_annot %>% 
      filter(label == "rA_P12") %>% 
      pull(est), na.rm = T)

# adapt to plot
schaefer7_association_mi = schaefer7_400_3d %>% #
  filter(hemi == "left")%>%
  mutate(ggseg_3d = map(ggseg_3d, ~ .x %>% 
                          full_join(associations2cortex_mi%>%
                                      filter(Parcel < 201) %>% #left hemisphere
                                      rename(region = label_Yeo_7) %>% 
                                      select(region,est),
                                    by = "region")))
# plot parcelations using ggseg3d
association_mi_left_medial = 
  ggseg3d(
    atlas = schaefer7_association_mi,
    colour = "est", text = "est",
    palette =  c( "#4d48a0"= min(associations2cortex$est, na.rm = T), "white"  = 0, "#a65a53" =max(associations2cortex$est, na.rm = T)),
    surface = "inflated",
    hemisphere = "left")%>% 
  remove_axes()%>%
  pan_camera("left medial")

association_mi_left_lateral = 
  ggseg3d(
    atlas = schaefer7_association_mi,
    colour = "est", text = "est",
    palette =  c( "#4d48a0"= min(associations2cortex$est, na.rm = T), "white"  = 0, "#a65a53" =max(associations2cortex$est, na.rm = T)),
    surface = "inflated",
    hemisphere = "left")%>% 
  remove_axes()%>%
  pan_camera("left lateral")


