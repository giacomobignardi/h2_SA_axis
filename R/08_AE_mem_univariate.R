#Author: Giacomo Bignardi
#Description: estimate association between strucutre and function
#Program: Univariate h2 twin with measurement error model------------------------------------------------------------------------------------------------------------------------------
# load packages
library(tidyverse)
library(readr)
library(lavaan)

# clear workspace
rm(list = ls())
set.seed(42)

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
# Functional gradient
FC_g1_i_d1 =  read_csv(sprintf("%s/03_GFC_i_d1.csv",wd_noa_output))
FC_g1_i_d2 =  read_csv(sprintf("%s/03_GFC_i_d2.csv",wd_noa_output))
FC_g1_i =  read_csv(sprintf("%s/03_GFC_i.csv",wd_noa_output))

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
annotations = annotations %>% rename(Parcel = parcel_Yeo)%>%   
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


# load dataFrames
HCP  = read_csv(sprintf("%s/unrestricted_giaco_6_25_2021_3_50_21.csv",wd_noa_data))
HCP_restricited  = read_csv(sprintf("%s/RESTRICTED_giaco_8_13_2021_11_47_39.csv",wd_noa_data))
HCP_all = merge(HCP,HCP_restricited, by = "Subject", all = T)

# FULL SAMPLE####
# full sample (inclusion criteria Valk et al., 2022 Nat Comm. excluding twins)
HCP_twin = HCP_all %>% 
  select(Family_ID,Subject, Gender, Age_in_Yrs,ZygosityGT) %>% 
  filter(Subject %in% Twin_sub$Subject)%>% 
  mutate(Sib_ID = ifelse(duplicated(Family_ID),2,1), #create a twinship order
         ZygosityGT = fct_recode(ZygosityGT, "1" = "MZ", "2" = "DZ"))

# INDIVIDUALS####
# bind dataframe to be able to fit a measurement error model later
FC_gradients = rbind(FC_g1_i_d1, FC_g1_i_d2, FC_g1_i)
FC_gradients = FC_gradients%>%
  filter(Sub %in% Twin_sub$Subject) %>% 
  select(Sub,Parcel,`0`, Session)%>%
  rename(value = `0`) %>% 
  pivot_wider(names_from = Session, values_from = value) %>% 
  rename(G1_d1 = `1.5`, G1_d2 = `2.5`, G1 = `0`)

# merge all data frames in list
SFs_i_final_list = list(FC_gradients, annotations)
SFs_i_final = SFs_i_final_list %>% reduce(full_join, by=c('Parcel'))

# merge with aesthetic sensitivity data
SFB_i_list = list(SFs_i_final,HCP_twin %>% rename(Sub = Subject))
SF_twin_i_final = SFB_i_list %>% reduce(full_join, by=c('Sub'))

# create a list of parcel to run models on later
parcel_list = unique(FC_gradients%>%select(Parcel))

# rename for sem function
SFB_i_final = SF_twin_i_final %>% 
  select(Family_ID,Sub,Parcel,Sib_ID, ZygosityGT, Gender, Age_in_Yrs,G1_d1,G1_d2,G1) %>% 
  select(-Sub) %>% 
  group_by(Parcel) %>% 
  ungroup() %>% 
# prepare for twin modeling
  pivot_wider(names_from = Sib_ID, values_from = Gender:G1) %>% 
  rename(P1_r1_T1=G1_d1_1,
         P1_r1_T2=G1_d1_2,
         P1_r2_T1=G1_d2_1,   
         P1_r2_T2=G1_d2_2,   
         P_T1=G1_1,
         P_T2=G1_2) %>% 
  as.data.frame() %>% 
  mutate(age_T1 = Age_in_Yrs_1,#note that the age of the twin is the same
         age_T2 = Age_in_Yrs_1,
         sex_T1 = Gender_1,#note that we only have same-sex twins so sex (not gender) is the same
         sex_T2 = Gender_1,
         zyg = ZygosityGT) %>% 
  select(-c(Family_ID, Gender_1,Gender_2,Age_in_Yrs_1,Age_in_Yrs_2,ZygosityGT))

# Classical Twin Design ####
# source CTD model specification
source(sprintf("%s/R/functions/lavaantwda/ADE_2g1p.R", wd_oa))

# estimate h2 across parcels
AE_sumy = c()
for(i in parcel_list$Parcel){
  # reinitialize at the beginning of each loop
  SFB_i_final_i = c()
  AE_fit_i = c()
  AE_sumy_i = c()
  AE_fit_i= c()
  
  # Including manifest covariates in multigroup SEM can make things complicated (e.g., unequal variance across groups, total phenotypic variance estimates)
  # so for twin modelin we apply a classic residualization procedures
  #  two step residualization for age and sex 
  SFB_i_final_i = umx::umx_residualize(c("P"), 
                                       cov = c("age","sex"), 
                                       suffixes=c("_T1","_T2"),
                                       data = SFB_i_final %>% filter(Parcel== i))
  
  # correlated factor solution (multivariate via Direct Symmetric Approach)
  AE_fit_i = sem(AE_model,
                 data = SFB_i_final_i,
                 group = "zyg",
                 missing = "ML",
                 std.ov = T, # standardize all observed variables first
                 group.label= c(1,2))
                         
  # extract summary statistics
  AE_sumy_i = summary(AE_fit_i, ci = T)
  AE_fit_i = fitmeasures(AE_fit_i,c("cfi","rmsea"))
  
  # bind to previous estimates
  AE_sumy = rbind(AE_sumy,
                        cbind(AE_sumy_i$pe %>%
                                filter(label %in% c("h2_P")) %>%
                                mutate(Parcel = i,
                                       cfi = AE_fit_i[1],
                                       rmsea = AE_fit_i[2])
                        )
  )
}

# retain only good fitting model (following Teeuw et al., 2021 Neuroimage)
AE_fltr = AE_sumy %>% filter(rmsea < 0.08, cfi > .90)
# note that 398 parcel meet goodness of fit criteria
nrow(AE_fltr)
# add information for which parcel the model did not fit well
AE_sumy %>% filter(rmsea >= 0.08, cfi <= .90) %>% pull(Parcel)
#[1]  61 251

# merge with annotation to stratify across cortical networks
AE_fltr_annot = merge(AE_fltr, annotations, by = "Parcel")

# vizualize heritability stratified by cortical networks
box_plot_AE = AE_fltr_annot %>% 
  ggplot(aes(label_Yeo7_short,  est, fill = label_Yeo7_short))+
  geom_boxplot() +
  geom_hline(yintercept = mean(AE_fltr_annot$est), linetype = "dashed")+
  ylim(0,1) +
  labs(y = expression(italic(h)[twin]^2),
       x = "Yeo-Krienen 7 networks")+
  theme_classic()+
  scale_fill_manual(values = c('#9F53AA',
                                        '#7A9ABD',
                                        '#3d8043',
                                        '#b584cf',
                                        '#F4FEC8',
                                        '#e8a633',
                                        '#D8707A'))+
                                          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),legend.position = "none")


# Common Pathway Measurement Error Model ####
# source MEM model specification
# Note that the model specification follows the model specification for a common pathway model
# but with only two observed variables
source(sprintf("%s/R/functions/lavaantwda/CPMEM_2g1p.R", wd_oa))

# estimate h2 across parcels
cpmem_AE_sumy = c()
for(i in parcel_list$Parcel){
  # reinitialize at the beginning of each loop
  cpmem_AE_fit_i = c()
  cpmem_AE_sumy_i = c()
  fitmeasures_i= c()
  SFB_i_final_i = c()
  
  # Including manifest covariates in multigroup SEM can make things complicated (e.g., unequal variance across groups, total phenotypic variance estimates)
  # so for twin modelin we apply a classic residualization procedures
  #  two step residualization for age and sex 
  SFB_i_final_i = umx::umx_residualize(c("P1_r1","P1_r2"), 
                                       cov = c("age","sex"), 
                                       suffixes=c("_T1","_T2"),
                                       data = SFB_i_final %>% filter(Parcel== i))
  
  
  # correlated factor solution (multivariate via Direct Symmetric Approach)
  cpmem_AE_fit_i = sem(cpmem_AE_model,
                               data = SFB_i_final_i,
                               group = "zyg",
                               missing = "ML",
                               std.ov = T, # standardize all observed variables first
                               group.label= c(1,2))
  
  # extract summary statistics
  cpmem_AE_sumy_i = summary(cpmem_AE_fit_i, ci = T)
  fitmeasures_i = fitmeasures(cpmem_AE_fit_i,c("cfi","rmsea"))

  # bind to previous estimates
  cpmem_AE_sumy = rbind(cpmem_AE_sumy,
                                cbind(cpmem_AE_sumy_i$pe %>%
                                        filter(label %in% c("h2_P1")) %>%
                                        mutate(Parcel = i,
                                               cfi = fitmeasures_i[1],
                                               rmsea = fitmeasures_i[2])
                                )
  )
}

# retain only good fitting model (following Teeuw et al., 2021 Neuroimage)
cpmem_AE_fltr = cpmem_AE_sumy %>% filter(rmsea < 0.08, cfi > .90)
# 395 good model fitting parcels
nrow(cpmem_AE_fltr)
# add information for which parcel the model did not fit well
cpmem_AE_sumy %>% filter(rmsea >= 0.08 | cfi <= .90) %>% pull(Parcel)
# 100 152 163 214 235

# exclude Heywood cases
cpmem_AE_fltr = cpmem_AE_fltr %>% filter(est < 1, est > -1)
# 1 Heywood cases
cpmem_AE_sumy %>% filter(est >= 1 | est <= -1) %>% pull(Parcel)
#[1] 114

# merge with annotation to stratify across cortical networks
cpmem_AE_fltr_annot = merge(cpmem_AE_fltr, annotations, by = "Parcel")


# box plot of univariate heritability (with mem)
box_plot_cpmem_AE = cpmem_AE_fltr_annot %>% 
  ggplot(aes(label_Yeo7_short,  est, fill = label_Yeo7_short))+
  geom_boxplot() +
  geom_hline(yintercept = mean(cpmem_AE_fltr_annot$est), linetype = "dashed")+
  ylim(0,1) +
  labs(y = expression(paste("deattenuated ", italic(h)[twin]^2)),
       x = "Yeo-Krienen 7 networks")+
  theme_classic()+
  scale_fill_manual(values = c('#9F53AA',
                                        '#7A9ABD',
                                        '#3d8043',
                                        '#b584cf',
                                        '#F4FEC8',
                                        '#e8a633',
                                        '#D8707A'))+
                                          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),legend.position = "none")

# comparison between the two models####
# 1 goodness of fit
nrow(cpmem_AE_fltr)-nrow(AE_fltr)

# average increase in heritability
((mean(cpmem_AE_fltr_annot$est) - mean(AE_fltr_annot$est))/(mean(AE_fltr_annot$est)))*100

sd(cpmem_AE_fltr_annot$est)
sd(AE_fltr_annot$est)
# visually
library(patchwork)

# point wise increase in heritability
# merge parcel that display sufficient model fit in both models
AE_comparison = merge(cpmem_AE_fltr_annot, AE_fltr_annot, by = c("Parcel"))
scatter_plot = AE_comparison %>% 
  ggplot(aes(est.y,est.x)) + 
  geom_hline(yintercept = mean(cpmem_AE_fltr_annot$est), linetype = "dashed")+
  geom_hline(yintercept = mean(AE_fltr_annot$est), linetype = "dashed", color = "lightgray")+
  geom_abline(a = 0, b = 1, linetype = "dashed")+
  geom_point(aes(fill = label_Yeo7_short.x), pch = 21, color = "black", size = 2) + 
  scale_fill_manual(values = c('#9F53AA',
                                         '#7A9ABD',
                                         '#3d8043',
                                         '#b584cf',
                                         '#F4FEC8',
                                         '#e8a633',
                                         '#D8707A'))+
                                           
  xlim(0,1)+
  ylim(0,1)+
  labs(x = expression(italic(h)[twin]^2),
       y = expression(paste("deattenuated ", italic(h)[twin]^2)),
       fill = "Yeo-Krienen 7")+
  theme_classic()

## save FIG. 4####
pdf(sprintf("%s/Figures/08_Fig_4.pdf",wd_oa), 
    width=9, 
    height=3)
(box_plot_AE|box_plot_cpmem_AE|scatter_plot) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
dev.off()
