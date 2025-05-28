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

# Microstructural intensity
MPmi_i  =  read_csv(sprintf("%s/HCP_S1200_MPC_400.csv",wd_noa_data)) %>% rename(Sub = "Var1")

# Geodesic distances
GD_i  =   read_csv(sprintf("%s/01_GD.csv",wd_noa_output))

# Functional gradient
FC_g1_i_d1 =  read_csv(sprintf("%s/03_GFC_i_d1.csv",wd_noa_output))
FC_g1_i_d2 =  read_csv(sprintf("%s/03_GFC_i_d2.csv",wd_noa_output))

# inclusion index 
Twin_sub =  read_csv(sprintf("%s/00_twin_ids.csv",wd_noa_output))

# cortical types and Yeo functional network annotations
annotations = read_csv(sprintf("%s/cortical_types.csv", wd_noa_data))

# tidy yeo 7 annotations 
network_7_yeo = c()
for(i in 1:400){
  net_i = c()
  net_i = strsplit(annotations$label_Yeo_7, "_")[[i]][3]
  network_7_yeo = rbind(network_7_yeo,net_i)
}

# network_7_yeo = c(network_7_yeo, rep("Other", 19))
annotations$label_Yeo7_short = network_7_yeo[,1]
annotations = annotations %>% rename(Parcel = parcel_Yeo) %>%   
  mutate(label_Yeo7_short = factor(label_Yeo7_short, levels = rev(c("Default",    
                                                                    "Cont",
                                                                    "Limbic", 
                                                                    "SalVentAttn",
                                                                    "DorsAttn",
                                                                    "SomMot",
                                                                    "Vis")))) %>%   
  mutate(label_Yeo7_short = fct_recode(label_Yeo7_short, 
                                       "Default" =  "Default",
                                       "Frontoparietal" =  "Cont",
                                       "Limbic" =  "Limbic",
                                       "Ventral attention" =  "SalVentAttn",
                                       "Dorsal attention" =  "DorsAttn",
                                       "Somatomotor" =  "SomMot",
                                       "Visual" =  "Vis"))

# inclusion index 
Twin_sub =  read_csv(sprintf("%s/00_twin_ids.csv",wd_noa_output))

# load dataFrames
HCP  = read_csv(sprintf("%s/unrestricted_giaco_6_25_2021_3_50_21.csv",wd_noa_data))
HCP_restricited  = read_csv(sprintf("%s/RESTRICTED_giaco_8_13_2021_11_47_39.csv",wd_noa_data))
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
SFs_i_final = SFs_i_final_list %>% reduce(full_join, by=c('Parcel'))

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
  select(-c(Gender_1,Gender_2,Age_in_Yrs_1,Age_in_Yrs_2,ZygosityGT))

#n of complete pairs for FGC1 
SFB_i_final %>% 
  filter(!is.na(P1_r1_T1) & !is.na(P1_r1_T2) & !is.na(P1_r2_T1) & !is.na(P1_r2_T2)) %>% 
  group_by(zyg) %>% 
  distinct(Family_ID) %>% 
  reframe(table(zyg))

#n of complete pairs for mp
SFB_i_final %>% 
  filter(!is.na(P2_T1) & !is.na(P2_T2)) %>% 
  group_by(zyg) %>% 
  distinct(Family_ID) %>% 
  reframe(table(zyg))

#n of complete pairs for gd
SFB_i_final %>% 
  filter(!is.na(P3_T1) & !is.na(P3_T2)) %>% 
  group_by(zyg) %>% 
  distinct(Family_ID) %>% 
  reframe(table(zyg))

# Twin multivariate model with measurement error ####
# source CTD model specification
source(sprintf("%s/R/functions/lavaantwda/CFM_cpmem_AE_2g3p.R", wd_oa))

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
write_csv(SF4A_annot, sprintf("%s/SI/SFILE3.csv",wd_oa))
write_csv(SF4B_annot, sprintf("%s/SI/SFILE4.csv",wd_oa))
write_csv(SF4C_annot, sprintf("%s/SI/SFILE5.csv",wd_oa))

#plot h2 - microstructure
mi_h2 = cfm_cpmem_AE_sumy %>% 
  filter(label == "h2_P2") %>% 
  merge(annotations %>% select(Parcel,label_Yeo7_short), by = "Parcel") %>% 
  ggplot(aes(label_Yeo7_short,  est, fill = label_Yeo7_short))+
  geom_boxplot() +
  #geom_hline(yintercept = mean(cpmem_AE_fltr_annot$est), linetype = "dashed")+
  ylim(0,1) +
  labs(y = expression(paste(italic(h)[twin]^2)),
       x = "Yeo-Krienen 7 networks",
       subtitle = "")+
  theme_classic(base_size = 10)+
  scale_fill_manual(values = c('#9F53AA',
                               '#7A9ABD',
                               '#3d8043',
                               '#b584cf',
                               '#F4FEC8',
                               '#e8a633',
                               '#D8707A'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),legend.position = "none")

#plot h2 - geodesic distance
gd_h2 = cfm_cpmem_AE_sumy %>% 
  filter(label == "h2_P3") %>% 
  merge(annotations %>% select(Parcel,label_Yeo7_short), by = "Parcel") %>% 
  ggplot(aes(label_Yeo7_short,  est, fill = label_Yeo7_short))+
  geom_boxplot() +
  #geom_hline(yintercept = mean(cpmem_AE_fltr_annot$est), linetype = "dashed")+
  ylim(0,1) +
  labs(y = expression(paste(italic(h)[twin]^2)),
       x = "Yeo-Krienen 7 networks",
       subtitle = "")+
  theme_classic(base_size = 10)+
  scale_fill_manual(values = c('#9F53AA',
                               '#7A9ABD',
                               '#3d8043',
                               '#b584cf',
                               '#F4FEC8',
                               '#e8a633',
                               '#D8707A'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),legend.position = "none")

# correlations
cfm_cpmem_AE_sumy_short =  cfm_cpmem_AE_sumy %>%                                 
  filter((grepl("rE",label) | grepl("rA",label) | grepl("h",label)) & !(grepl("var",label)|grepl("bh",label)))



#each parcel has 9 (6 cor and 3 h2) parameter estimates
nrow(cfm_cpmem_AE_sumy_short)/9

# retain only good fitting model (following Teeuw et al., 2021 Neuroimage)
cfs_fltr_fit = cfm_cpmem_AE_sumy_short %>% 
  filter(rmsea < 0.08, cfi > .90)

# how many with good fit
n_parcel_model = nrow(cfs_fltr_fit)/9

# exclude Heywood cases
cfs_fltr_fit %>% filter(est > 1 | est < -1)
# the parcel which was already identify (114, limbic network) produce out of bounds h2  
parcel_exc <- cfs_fltr_fit %>% filter(est > 1 & grepl("h2",label)) %>% pull(Parcel)

# merge with annotation to stratify across cortical networks
AE_fltr_annot = merge(cfs_fltr_fit, annotations, by = "Parcel") 

# correct for multiple comparisions
# NOTE P1P2 = inter FCg1 <-> t1t2w
# NOTE P1P3 = inter FCg1 <-> GD
est_AE_annot = AE_fltr_annot %>% 
  group_by(label) %>% 
  mutate(p_adjust = p.adjust(pvalue, method ="bonferroni"))

## save####
write_csv(est_AE_annot, sprintf("%s/09_CPMEM_output.csv",wd_noa_output))

# average heritability
est_AE_annot %>% 
  filter(grepl("h2",label)) %>% 
  filter(Parcel != parcel_exc) %>% #we remove parcel 114 as it produced out of bounds estimates for the measurement error model
  group_by(lhs) %>% 
  reframe(mean(est), sd(est))

# display percentage of significant genetic and environemntal correlations 
# amongst all the possible FCg1 <-> t1t2w (that is twice the number of parcels)
(est_AE_annot %>% 
  filter(Parcel != parcel_exc) %>%
  filter(p_adjust<.05 & label == "rA_P12" & est <0) %>% nrow() +
est_AE_annot %>% 
  filter(Parcel != parcel_exc) %>%
  filter(p_adjust<.05 & label == "rA_P12" & est >0) %>% nrow() +
est_AE_annot %>% 
  filter(Parcel != parcel_exc) %>%
  filter(p_adjust<.05 & label == "rE_P12" & est <0) %>% nrow() +
est_AE_annot %>% 
  filter(Parcel != parcel_exc) %>%
  filter(p_adjust<.05 & label == "rE_P12" & est >0) %>% nrow()) /
  ((n_parcel_model-1)*2)

# focus on genetic correlations
(est_AE_annot %>% 
    filter(Parcel != parcel_exc) %>%
    filter(p_adjust<.05 & label == "rA_P12") %>% nrow())/
  (n_parcel_model-1)

# display percentage of significant genetic and environemntal correlations 
# amongst all the possible FCg1 <-> GD (that is twice the number of parcels)
(est_AE_annot %>% 
    filter(Parcel != parcel_exc) %>%
    filter(p_adjust<.05 & label == "rA_P13" & est <0) %>% nrow() +
    est_AE_annot %>% 
    filter(Parcel != parcel_exc) %>%
    filter(p_adjust<.05 & label == "rA_P13" & est >0) %>% nrow()) /
  (n_parcel_model-1)

(est_AE_annot %>% 
    filter(Parcel != parcel_exc) %>%
    filter(p_adjust<.05 & label == "rE_P13" & est <0) %>% nrow() +
    est_AE_annot %>% 
    filter(Parcel != parcel_exc) %>%
    filter(p_adjust<.05 & label == "rE_P13" & est >0) %>% nrow()) /
  (n_parcel_model-1)

# focus on genetic correlations
(est_AE_annot %>% 
    filter(Parcel != parcel_exc) %>%
    filter(p_adjust<.05 & label == "rA_P13" & est <0) %>% nrow())/
  (n_parcel_model-1)

(est_AE_annot %>% 
    filter(Parcel != parcel_exc) %>%
    filter(p_adjust<.05 & label == "rA_P13" & est >0) %>% nrow()) /
  (n_parcel_model-1)

# average magnitude of genetic correlations
est_AE_annot %>% 
  filter(Parcel != parcel_exc) %>%
    filter(p_adjust<.05 & label == "rA_P13" & est <0) %>% 
  summarise(mean= psych::fisherz2r(mean(psych::fisherz(est))),
            sd = psych::fisherz2r(sd(psych::fisherz(est))))

est_AE_annot %>% 
  filter(Parcel != parcel_exc) %>%
  filter(p_adjust<.05 & label == "rA_P13" & est >0) %>% 
  summarise(mean= psych::fisherz2r(mean(psych::fisherz(est))),
            sd = psych::fisherz2r(sd(psych::fisherz(est))))


# focus on enviromental correlations
(est_AE_annot %>% 
    filter(Parcel != parcel_exc) %>%
    filter(p_adjust<.05 & label == "rE_P13" & est <0) %>% 
  nrow() +
  est_AE_annot %>% 
    filter(Parcel != parcel_exc) %>%
    filter(p_adjust<.05 & label == "rE_P13" & est >0) %>% 
  nrow()) /
  (n_parcel_model-1)


# plot associations on the cortex
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

# plot in 2D
cort_rA = est_AE_annot %>%
  filter((grepl("rA",label))) %>% 
  filter(Parcel!=parcel_exc) %>% 
  #first we re annotate the correlations names
  mutate(label = factor(label, levels = c("rA_P12","rA_P13","rA_P23")),
         label = fct_recode( label, 
                             "Microstructural intensity" = "rA_P12",
                             "Geodesic distances" = "rA_P13",
                             "remove"  = "rA_P23")) %>% 
  filter(label !="remove") %>%
  rename(region = label_Yeo_7,
         cor = label) %>% 
  mutate(est = ifelse(p_adjust<.05,est,NA)) %>% 
  ggplot() +
  facet_grid(col = vars(cor)) +
  geom_brain(atlas = schaefer7_400, 
             position = position_brain(hemi ~ side), #to stack them
             aes(fill = est),
             color = "#414a4c",
             size = .5) +
  scale_fill_distiller(palette = "Reds", direction  = 1)+
  labs(fill = expression(italic(r)[A]))+
  theme_classic(base_size = 10)+
  theme(axis.text  = element_blank(),
        axis.line  = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank())

cort_rE = est_AE_annot %>%
  filter((grepl("rE",label))) %>% 
  filter(Parcel!=parcel_exc) %>% 
  #first we re annotate the correlations names
  mutate(label = factor(label, levels = c("rE_P12","rE_P13","rE_P23")),
         label = fct_recode( label, 
                             "Microstructural intensity" = "rE_P12",
                             "Geodesic distances" = "rE_P13",
                             "remove"  = "rE_P23")) %>% 
  filter(label !="remove") %>%
  rename(region = label_Yeo_7,
         cor = label) %>% 
  mutate(est = ifelse(p_adjust<.05,est,NA)) %>% 
  ggplot() +
  facet_grid(col = vars(cor)) +
  geom_brain(atlas = schaefer7_400, 
             position = position_brain(hemi ~ side), #to stack them
             aes(fill = est),
             color = "#414a4c",
             size = .5) +
  scale_fill_distiller(palette = "Blues", direction  = 1)+
  labs(fill = expression(italic(r)[E]))+
  theme_classic(base_size = 10) +
  theme(axis.text  = element_blank(),
        axis.line  = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank())
  

# Figure 5####
pdf(sprintf("%s/Figures/09_Fig_5.pdf",wd_oa),
    width=9, 
    height=5)
(((mi_h2/gd_h2) + plot_layout(axes = "collect"))|(cort_rA/cort_rE) + plot_layout(axes = "collect")) +plot_layout(widths = c(1,2.5), heights = c(1,1)) + plot_annotation(tag_levels = "A")
dev.off()