#Author: Giacomo Bignardi
#Description: estimate association between structure and function between-individuals correcting for measurement error
#Program: Strucutral Measurement Error Equation Modeling------------------------------------------------------------------------------------------------------------------------------
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

# load measurement error model
source(sprintf("%s/R/functions/meermo/SMEM.R", wd_oa))

# load
# Microstructural intensity
MPmi_i = read_csv(sprintf("%s/HCP_S1200_MPC_400.csv",wd_noa_data)) %>% rename(Sub = "Var1")

# Geodesic distances
GD_i = read_csv(sprintf("%s/01_GD.csv",wd_noa_output))

# Functional gradient
FC_g1_i_d1 = read_csv(sprintf("%s/03_GFC_i_d1.csv",wd_noa_output))
FC_g1_i_d2 = read_csv(sprintf("%s/03_GFC_i_d2.csv",wd_noa_output))

# inclusion index 
notTwin_sub = read_csv(sprintf("%s/00_nottwin_ids.csv",wd_noa_output))

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
annotations = annotations %>% rename(Parcel = parcel_Yeo)

# load demographics
HCP  = read_csv(sprintf("%s/unrestricted_giaco_6_25_2021_3_50_21.csv",wd_noa_data))
HCP_restricited  = read_csv(sprintf("%s/RESTRICTED_giaco_8_13_2021_11_47_39.csv",wd_noa_data))
HCP_all = merge(HCP,HCP_restricited, by = "Subject", all = T) %>% 
  select(Subject,Gender,Age_in_Yrs) %>% 
  filter(Subject %in% notTwin_sub$Subject)

# INDIVIDUALS####
# get individual values for mp
MPmi_i = MPmi_i%>%
  filter(Sub %in% notTwin_sub$Subject) %>%
  pivot_longer(names_to = "Parcel", values_to = "t1w/t2w_mi", c(2:ncol(MPmi_i)))%>%mutate(Parcel = substr(Parcel,6,nchar(Parcel)))%>%
  rename(value = `t1w/t2w_mi`) %>%
  mutate(Parcel = as.numeric(Parcel)) %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  rename(MP = value)

# get individual values for gd
GD_i = GD_i %>%
  filter(Sub %in% notTwin_sub$Subject) %>%
  group_by(Parcel) %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  rename(GD = value)

# bind dataframe to be able to fit a measurement error model later
FC_gradients = rbind(FC_g1_i_d1, FC_g1_i_d2)
FC_gradients = FC_gradients%>%
  filter(Sub %in% notTwin_sub$Subject) %>% #filter for twins
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

# SEM####
# rename for sem function
SFs = SFs %>% 
  rename(P1_r1 = G1_d1, P1_r2 =G1_d2,
         P2 = MP,
         P3 = GD,
         sex = Gender,
         age = Age_in_Yrs)
set.seed(42)

# store estimates for each parcel
est_smem = c()

# loop mem function over parcels
for( i in parcel_list$Parcel){
  # reinitialize
  est_smem_i = c()
  SFs_i = c()
  
  # residualise for sex and age
  SFs_i = umx::umx_residualize(c("P1_r1","P1_r2","P2","P3"), 
                       cov = c("age","sex"), 
                       data = SFs %>% filter(Parcel== i))
  
  # fit the measurement error model to parcel i
  smem_fit = sem(smem, data = SFs_i, std.ov = T, missing = "ML")
  
  # extract standardized estimated for the parameters
  smem_parm_i = standardizedSolution(smem_fit)
 
  # extract correlations
  est_smem_i = smem_parm_i %>% 
    mutate(Parcel = i)
  
  # bind across parcels
  est_smem = rbind(est_smem,est_smem_i)
}

## Supplementary File 3 ####
# SF3a: parameters in mem
SF3A = est_smem %>%                                 
  filter(grepl("var_intra",label) | grepl("f1",label) | grepl("rP",label)) %>% 
  group_by(Parcel) %>% 
  distinct(label, .keep_all = T) %>% 
  mutate(label = recode(label,
                        var_intra_r1 = "σ1intra^2",
                        var_intra_r2 = "σ2intra^2",
                        f1 = "λ",
                        rP1P2 = "σ_mig1",
                        rP1P3 = "σ_gdg1",
                        rP2P3 = "σ_migd"
  ))

SF3A_annot = merge(SF3A, annotations %>% select(Parcel,label_Yeo_7,label_Yeo7_short), by = "Parcel")%>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  select(-c(lhs,op,rhs,))


# SF2b: parameters in mem
write_csv(SF3A_annot, sprintf("%s/SI/SFILE2.csv",wd_oa))


# merge with annotation to stratify across cortical networks
est_smem_annot = merge(est_smem %>%  filter(grepl("rP",label)) %>%  group_by(Parcel) %>% 
                         distinct(label, .keep_all = T) %>% ungroup() %>%  select(label, est.std ,pvalue, ci.lower, ci.upper, Parcel) , annotations, by = "Parcel") %>% mutate(label_Yeo7_short = factor(label_Yeo7_short, levels = rev(c("Default", 
                                                                                                                                           "Limbic", 
                                                                                                                                           "Cont",
                                                                                                                                           "SalVentAttn",
                                                                                                                                           "DorsAttn",
                                                                                                                                           "SomMot",
                                                                                                                                           "Vis"))))



# correct for multiple comparisions
# NOTE P1P2 = inter FCg1 <-> t1t2w
# NOTE P1P3 = inter FCg1 <-> GD
est_smem_annot = est_smem_annot %>% 
  group_by(label) %>% 
  mutate(p_adjust = p.adjust(pvalue, method ="bonferroni"))

# display total significant parcels
n_negative_sig_t1wt2w = est_smem_annot %>% 
  filter(p_adjust<.05 & label == "rP1P2" & est.std <0) %>% nrow()
est_smem_annot %>% 
  filter(p_adjust<.05 & label == "rP1P2" & est.std >0) %>% nrow()
n_positive_sig = est_smem_annot %>% 
  filter(p_adjust<.05 & label == "rP1P3" & est.std >0) %>% nrow()
n_negative_sig = est_smem_annot %>% 
  filter(p_adjust<.05 & label == "rP1P3" & est.std <0) %>% nrow()
est_smem_annot %>% 
  filter(p_adjust<.05 & label == "rP2P3" & est.std >0) %>% nrow()
est_smem_annot %>% 
  filter(p_adjust<.05 & label == "rP1P2" & est.std <0) %>% 
  select(Parcel,est.std,class_types,label_Yeo7_short)
names(est_smem_annot)

# percentage of significant parcels
((n_negative_sig_t1wt2w) / 400)*100 #1.75
((n_positive_sig + n_negative_sig) / 400)*100 #56.75

# clustering of results based on functional network
est_smem_annot %>% 
  filter(label == "rP1P3" & p_adjust<.05) %>% 
  group_by(label_Yeo7_short) %>% 
  summarise(mean_r = psych::fisherz2r(mean(psych::fisherz(est.std))))

# one sample t
t_vis = t.test(psych::fisherz(est_smem_annot %>% filter(label_Yeo7_short == "Vis" & label == "rP1P3" & p_adjust<.05) %>% pull(est.std)))
t_som = t.test(psych::fisherz(est_smem_annot %>% filter(label_Yeo7_short == "SomMot" & label == "rP1P3" & p_adjust<.05) %>% pull(est.std)))
t_dor = t.test(psych::fisherz(est_smem_annot %>% filter(label_Yeo7_short == "DorsAttn" & label == "rP1P3" & p_adjust<.05) %>% pull(est.std)))
t_sal = t.test(psych::fisherz(est_smem_annot %>% filter(label_Yeo7_short == "SalVentAttn" & label == "rP1P3" & p_adjust<.05) %>% pull(est.std)))
t_lim = t.test(psych::fisherz(est_smem_annot %>% filter(label_Yeo7_short == "Limbic" & label == "rP1P3" & p_adjust<.05) %>% pull(est.std)))
t_con = t.test(psych::fisherz(est_smem_annot %>% filter(label_Yeo7_short == "Cont" & label == "rP1P3" & p_adjust<.05) %>% pull(est.std)))
t_dmn = t.test(psych::fisherz(est_smem_annot %>% filter(label_Yeo7_short == "Default" & label == "rP1P3" & p_adjust<.05) %>% pull(est.std)))

# conservative approach (adjust for Bonferroni)
p_t = p.adjust(c(t_vis$p.value, t_som$p.value, t_dor$p.value, t_sal$p.value, t_lim$p.value, t_con$p.value, t_dmn$p.value), method = "bonferroni")
p_t <.05
t_vis; t_som; t_sal; t_dmn

# plot associations on the cortex
## Figure 2F###
# template surface 
associations2cortex = est_smem_annot %>% 
  filter(label == "rP1P3" & p_adjust<.05) %>% 
  ungroup() %>% 
  select(Parcel, label_Yeo_7, est.std)

# range of parcel-wise estimates
min(associations2cortex$est.std, na.rm = T)
max(associations2cortex$est.std, na.rm = T)

# adapt to plot
schaefer7_association_GD = schaefer7_400_3d %>% #
  filter(hemi == "left")%>%
  mutate(ggseg_3d = map(ggseg_3d, ~ .x %>% 
                          full_join(associations2cortex%>%
                                      filter(Parcel < 201) %>% #left hemisphere
                                      rename(region = label_Yeo_7) %>% 
                                      select(region,est.std),
                                    by = "region")))

# plot parcelations using ggseg3d
association_GD_left_medial = 
  ggseg3d(
    atlas = schaefer7_association_GD,
    colour = "est.std", text = "est.std",
    palette =  c( "#4d48a0"= min(associations2cortex$est.std, na.rm = T), "white"  = 0, "#a65a53" =max(associations2cortex$est.std, na.rm = T)),
    surface = "inflated",
    hemisphere = "left")%>% 
  remove_axes()%>%
  pan_camera("left medial")

association_GD_left_lateral = 
  ggseg3d(
    atlas = schaefer7_association_GD,
    colour = "est.std", text = "est.std",
    palette =  c( "#4d48a0"= min(associations2cortex$est.std, na.rm = T), "white"  = 0, "#a65a53" =max(associations2cortex$est.std, na.rm = T)),
    surface = "inflated",
    hemisphere = "left")%>% 
  remove_axes()%>%
  pan_camera("left lateral")


# associaiton with microstructure
associations2cortex_mi = est_smem_annot %>% 
  filter(label == "rP1P2" & p_adjust<.05) %>% 
  ungroup() %>% 
  select(Parcel, label_Yeo_7, est.std)

# range of parcel-wise estimates
min(associations2cortex_mi$est.std, na.rm = T)
max(associations2cortex_mi$est.std, na.rm = T)

# adapt to plot
schaefer7_association_mi = schaefer7_400_3d %>% #
  filter(hemi == "left")%>%
  mutate(ggseg_3d = map(ggseg_3d, ~ .x %>% 
                          full_join(associations2cortex_mi%>%
                                      filter(Parcel < 201) %>% #left hemisphere
                                      rename(region = label_Yeo_7) %>% 
                                      select(region,est.std),
                                    by = "region")))

# plot parcelations using ggseg3d
association_mi_left_medial = 
  ggseg3d(
    atlas = schaefer7_association_mi,
    colour = "est.std", text = "est.std",
    palette =  c( "#4d48a0"= min(associations2cortex_mi$est.std, na.rm = T), "white"  = 0),
    surface = "inflated",
    hemisphere = "left")%>% 
  remove_axes()%>%
  pan_camera("left medial")

association_mi_left_lateral = 
  ggseg3d(
    atlas = schaefer7_association_mi,
    colour = "est.std", text = "est.std",
    palette =  c( "#4d48a0"= min(associations2cortex_mi$est.std, na.rm = T), "white"  = 0),
    surface = "inflated",
    hemisphere = "left")%>% 
  remove_axes()%>%
  pan_camera("left lateral")



# FIGURE 3B (2D)####
#here we make a simpler 2d figure
cort_r = est_smem_annot %>%
  #first we re annotate the correlations names
  mutate(label = fct_recode( label, 
                            "Microstructural intensity" = "rP1P2",
                            "Geodesic distances" = "rP1P3",
                            "remove"  = "rP2P3")) %>% 
  filter(label !="remove") %>%
  rename(region = label_Yeo_7,
         cor = label) %>% 
  mutate(est.std = ifelse(p_adjust<.05,est.std,NA)) %>% 
  ggplot() +
  facet_grid(~cor) +
  geom_brain(atlas = schaefer7_400, 
             position = position_brain(hemi ~ side), #to stack them
             aes(fill = est.std),
             color = "#414a4c",
             size = .5) +
  scale_fill_distiller(palette = "RdBu")+
  labs(fill = expression(italic(r)[p]))+
  theme_void(base_size = 20)+
  theme(legend.position = "bottom",
        legend.title.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

## save FIG. 3B####
pdf(sprintf("%s/Figures/07_Fig_3B.pdf",wd_oa),
    width = 10,
    height = 5 )
cort_r
dev.off()

# SUPPLEMENTARY ####
# estimate correlations without partitioning intra-individual variability out to check for overall improvement
# Load average Functional gradient
FC_g1_i =  read_csv(sprintf("%s/03_GFC_i.csv",wd_noa_output))

#get individual values for the functional gradient
FC_g1_i = FC_g1_i%>%
  filter(Sub %in% notTwin_sub$Subject) %>%
  select(Sub,Parcel,`0`)%>%
  rename(G1 = `0`)

# put all of them togheter
SFs_averaged_i_list = list(MPmi_i, GD_i,FC_g1_i)

# merge all data frames in list
SFs_averaged_i = SFs_averaged_i_list %>% reduce(full_join, by=c("Sub",'Parcel'))
SFs_averaged_i_final_list = list(SFs_averaged_i, annotations)
SFs_averaged_i_final = SFs_averaged_i_final_list %>% reduce(full_join, by=c('Parcel'))%>% 
  mutate(label_Yeo7_short = factor(label_Yeo7_short, levels = rev(c("Default", 
                                                                    "Limbic", 
                                                                    "Cont",
                                                                    "SalVentAttn",
                                                                    "DorsAttn",
                                                                    "SomMot",
                                                                    "Vis"))))

# CORRELATIONS####
# rename
SFs_averaged_i_final = SFs_averaged_i_final %>% 
  rename(P1 = G1,
         P2 = MP,
         P3 = GD)
set.seed(42)

# store estimates for each parcel
est_r = c()

# loop simple correlations over parcels
for(i in parcel_list$Parcel){
  # reinitialize
  est_r_i = c()
  SFs_averaged_i_final_i = c()
  
  # scale
  SFs_averaged_i_final_i = SFs_averaged_i_final %>% 
    filter(Parcel==i) %>% 
    mutate(P1 = ((P1-mean(P1,na.rm = T))/sd(P1,na.rm = T)),
           P2 = ((P2-mean(P2,na.rm = T))/sd(P2,na.rm = T)),
           P3 = ((P3-mean(P3,na.rm = T))/sd(P3,na.rm = T)))
  
  # fit the measurement error model to parcel i
  rP1P2_fit = cor.test(SFs_averaged_i_final_i$P1,SFs_averaged_i_final_i$P2)
  rP1P3_fit = cor.test(SFs_averaged_i_final_i$P1,SFs_averaged_i_final_i$P3)
  rP2P3_fit = cor.test(SFs_averaged_i_final_i$P2,SFs_averaged_i_final_i$P3)
  
  rP1P2_est = data.frame(Parcel = i,
                                label = "rP1P2",
                                est = rP1P2_fit$estimate, 
                                pvalue = rP1P2_fit$p.value)
  rP1P3_est = data.frame(Parcel = i,
                                label = "rP1P3",
                                est = rP1P3_fit$estimate, 
                                pvalue = rP1P3_fit$p.value)
  rP2P3_est = data.frame(Parcel = i,
                                label = "rP2P3",
                                est = rP2P3_fit$estimate, 
                                pvalue = rP2P3_fit$p.value)
  

  est_r_i = rbind(rP1P2_est,rP1P3_est, rP2P3_est)

  
  # bind across parcels
  est_r = rbind(est_r,est_r_i)
}

# merge with annotation to stratify across cortical networks
est_r_annot = merge(est_r, annotations, by = "Parcel") %>%   mutate(label_Yeo7_short = factor(label_Yeo7_short, levels = rev(c("Default", 
                                                                                                                                     "Limbic", 
                                                                                                                                     "Cont",
                                                                                                                                     "SalVentAttn",
                                                                                                                                     "DorsAttn",
                                                                                                                                     "SomMot",
                                                                                                                                     "Vis"))))

# correct for multiple comparisions
# NOTE P1P2 = inter FCg1 <-> t1t2w
# NOTE P1P3 = inter FCg1 <-> GD
est_r_annot = est_r_annot %>% 
  group_by(label) %>% 
  mutate(p_adjust = p.adjust(pvalue, method ="bonferroni"))

## save####
write_csv(est_smem, sprintf("%s/07_SMEM_output.csv",wd_noa_output))

# plot overall improvement
est = merge(est_smem %>% rename(est_smem = est.std), est_r %>% rename(est_r = est), by = c("Parcel","label"))
est_annot = merge(est, annotations, by = "Parcel") %>%   mutate(label_Yeo7_short = factor(label_Yeo7_short, levels = rev(c("Default", 
                                                                                                                               "Limbic", 
                                                                                                                               "Cont",
                                                                                                                               "SalVentAttn",
                                                                                                                               "DorsAttn",
                                                                                                                               "SomMot",
                                                                                                                               "Vis"))))

## SFIG. S1####
scatter_plot = est_annot %>% 
  filter(label == "rP1P3") %>% 
  ggplot(aes(est_r,est_smem)) + 
  geom_abline(a = 0, b = 1, linetype = "dashed")+
  geom_point(aes(fill = label_Yeo7_short), pch = 21, color = "black", size = 2.0) + 
  geom_smooth(method = "lm", se = F, color = "lightGray")+
  scale_fill_manual(values = c('#9F53AA',
                                         '#7A9ABD',
                                         '#3d8043',
                                         '#b584cf',
                                         '#e8a633',
                                         '#F4FEC8',
                                         '#D8707A'))+
  xlim(-1,1)+
  ylim(-1,1)+
  labs(x = "estimates pearson r")+
  labs(y = "estimates smem")+
  labs(fill = "Yeo-Krienen\n 7 networks")+
  theme_classic()

## save FIG. S1####
pdf(sprintf("%s/Figures/07_Fig_S2.pdf",wd_oa), 
    width=4.5, 
    height=3)
scatter_plot 
dev.off()

# range of improvement
range(est_annot$est_smem-est_annot$est_r)
