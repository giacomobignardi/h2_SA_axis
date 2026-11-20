#Author: Giacomo Bignardi
#Description: estimate association between structure and function between-individuals correcting for measurement error
#Program: Strucutral Measurement Error Equation Modeling------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(reshape2)
library(patchwork)
library(ggseg)
library(ggsegSchaefer)
library(ggridges)
library(lavaan)

# clear workspace
rm(list = ls())
g1_tag = "_signmatched"   # "" = Procrustes-aligned (main); "_signmatched" / "_ordermatched" = supplement (R2.2)

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
FC_g1_i_d1 = read_csv(sprintf("%s/03_GFC_i_d1%s.csv",wd_noa_output, g1_tag))
FC_g1_i_d2 = read_csv(sprintf("%s/03_GFC_i_d2%s.csv",wd_noa_output, g1_tag))

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
SFs_i_final = SFs_i_final_list %>% reduce(full_join, by=c('Parcel'))

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
  mutate(across(where(is.numeric), round, 5)) %>% 
  select(-c(lhs, op, rhs))

# SF2: parameters in the structural measurement error model
if (g1_tag == "") write_csv(SF3A_annot, sprintf("%s/SI/SFILE2.csv",wd_oa))

# merge with annotation to stratify across cortical networks
est_smem_annot = merge(est_smem %>%  filter(grepl("rP",label)) %>%  group_by(Parcel) %>% 
                         distinct(label, .keep_all = T) %>% ungroup() %>%  select(label, est.std ,pvalue, ci.lower, ci.upper, Parcel) , annotations, by = "Parcel")

# correct for multiple comparisons (Bonferroni: main text; BH-FDR, more lenient: supplement)
# NOTE P1P2 = inter FCg1 <-> t1t2w
# NOTE P1P3 = inter FCg1 <-> GD
est_smem_annot = est_smem_annot %>% 
  group_by(label) %>% 
  mutate(p_adjust = p.adjust(pvalue, method = "bonferroni"),
         p_fdr    = p.adjust(pvalue, method = "BH"))

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
((n_negative_sig_t1wt2w) / 400)*100 #2.25
((n_positive_sig + n_negative_sig) / 400)*100 #57.5

# counts per correction (Bonferroni: main text; FDR and uncorrected: supplement)
est_smem_annot %>% group_by(label) %>%
  summarise(n_bonf = sum(p_adjust < .05), n_fdr = sum(p_fdr < .05), n_uncorr = sum(pvalue < .05),
            n_bonf_pos = sum(p_adjust < .05 & est.std > 0), n_bonf_neg = sum(p_adjust < .05 & est.std < 0),
            median_r = median(est.std), median_abs_r = median(abs(est.std)))

# network clustering of GD-G1 associations: test ALL parcels within each network
net_t = est_smem_annot %>% ungroup() %>% filter(label == "rP1P3") %>%
  group_by(label_Yeo7_short) %>%
  summarise(n      = n(),
            mean_r = psych::fisherz2r(mean(psych::fisherz(est.std))),
            t      = unname(t.test(psych::fisherz(est.std))$statistic),
            df     = unname(t.test(psych::fisherz(est.std))$parameter),
            p      = t.test(psych::fisherz(est.std))$p.value, .groups = "drop") %>%
  mutate(p_bonf = p.adjust(p, "bonferroni"))
net_t

# plot associations on the cortex
# FIGURE 3B (2D)####
cort_r_data = est_smem_annot %>% ungroup() %>%
  mutate(label = fct_recode(label,
                            "Microstructural intensity" = "rP1P2",
                            "Geodesic distances" = "rP1P3",
                            "remove" = "rP2P3")) %>%
  filter(label != "remove") %>%
  rename(region = label_Yeo_7, cor = label) %>%
  mutate(est.std = ifelse(p_adjust < .05, est.std, NA)) %>%
  select(region, cor, est.std) %>%
  group_by(cor)

cort_r = ggplot() +
  geom_brain(atlas = schaefer7_400(),
             data = cort_r_data,
             mapping = aes(fill = est.std),
             position = position_brain(hemi ~ view), #to stack them
             color = "#414a4c",
             linewidth = .5) +
  facet_grid(~cor) +
  scale_fill_distiller(palette = "RdBu") +
  labs(fill = expression(italic(r)[p])) +
  theme_void(base_size = 14) +
  theme(legend.title.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

# SUPPLEMENTARY ####
# estimate correlations without partitioning intra-individual variability out to check for overall improvement
# Load average Functional gradient
FC_g1_i =  read_csv(sprintf("%s/03_GFC_i%s.csv",wd_noa_output, g1_tag))

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
SFs_averaged_i_final = SFs_averaged_i_final_list %>% reduce(full_join, by=c('Parcel'))
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
est_r_annot = merge(est_r, annotations, by = "Parcel")

# correct for multiple comparisons
# NOTE P1P2 = inter FCg1 <-> t1t2w
# NOTE P1P3 = inter FCg1 <-> GD
est_r_annot = est_r_annot %>% 
  group_by(label) %>% 
  mutate(p_adjust = p.adjust(pvalue, method = "bonferroni"))

# symmetric analysis without measurement-error correction (R3.3)
est_r_annot %>% group_by(label) %>% summarise(n_bonf = sum(p_adjust < .05), median_r = median(est))

## save####
write_csv(est_smem, sprintf("%s/07_SMEM_output%s.csv",wd_noa_output, g1_tag))

# plot overall improvement
est = merge(est_smem %>% rename(est_smem = est.std), est_r %>% rename(est_r = est), by = c("Parcel","label"))
est_annot = merge(est, annotations, by = "Parcel")

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
                                         '#F4FEC8',
                                         '#e8a633',
                                         '#D8707A'))+
  xlim(-1,1)+
  ylim(-1,1)+
  labs(x = expression(italic(r)[gd-g1]))+
  labs(y =  expression(paste("deattenuated ",italic(r)[gd-g1])))+
  labs(fill = "Yeo-Krienen 7")+
  theme_classic(base_size = 14)


## save FIG. 3####
pdf(sprintf("%s/07_Fig_3%s.pdf", if (g1_tag == "") paste0(wd_oa, "/Figures") else wd_noa_output, g1_tag),
    width=9.5, 
    height=2.75)
(cort_r|scatter_plot) +plot_layout(widths = c(2.5,1)) + plot_annotation(tag_levels = "A")
dev.off()

# range of improvement
range(est_annot$est_smem-est_annot$est_r)


