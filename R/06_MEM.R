#Author: Giacomo Bignardi
#Description: measurement error model to account for test-retest data in FCG1
#Program: measurement error model------------------------------------------------------------------------------------------------------------------------------
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

# load measurement error model (meermo directory)
source(sprintf("%s/R/functions/meermo/mem.R", wd_oa))

# load
# Functional gradient
FC_g1_i_d1 =  read_csv(sprintf("%s/03_GFC_i_d1.csv",wd_noa_output))
FC_g1_i_d2 =  read_csv(sprintf("%s/03_GFC_i_d2.csv",wd_noa_output))

# inclusion index 
notTwin_sub =  read_csv(sprintf("%s/00_nottwin_ids.csv",wd_noa_output))

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

                                                                                                                              
# FUNCTIONAL GRADIENTS####
# bind dataframe to be able to fit a measurement error model later
FC_gradients = rbind(FC_g1_i_d1, FC_g1_i_d2)
FC_gradients = FC_gradients%>%
  filter(Sub %in% notTwin_sub$Subject) %>% #filter for twins
  select(Sub,Parcel,`0`, Session)%>%
  rename(value = `0`) %>% 
  pivot_wider(names_from = Session, values_from = value) %>% 
  rename(G1_d1 = `1.5`, G1_d2 = `2.5`)

# create a list of parcel to run models on later
parcel_list = unique(FC_gradients%>%select(Parcel))

# MEASURMENT ERROR MODEL####
# rename for mem function
FC_gradients = FC_gradients %>% 
  rename(P1_r1 = G1_d1, P1_r2 =G1_d2)
set.seed(42)

# store estimates for each parcel
est_mem = c()
mem_parm = c()
# loop mem function over parcels
for( i in parcel_list$Parcel){
  
  # reinitialize
  est_mem_i = c()
  FC_gradients_i = c()
  mem_fitmeasures_i = c()

  # fit the measurement error model to parcel i
  mem_fit_i = cfa(mem, 
                  data = FC_gradients %>% filter(Parcel==i), 
                  std.ov = T)
  
  # extract parameters
  mem_parm_i = parameterestimates(mem_fit_i)
  mem_fitmeasures_i = fitmeasures(mem_fit_i,c("cfi", "rmsea"))
  
  # extract ICCs and inter-individual variance
  est_mem_i = mem_parm_i %>% 
    filter(lhs == "ICC2" | lhs == "ICC2k" | label == "var_inter") %>% 
    select(label, est) %>% 
    pivot_wider(names_from = "label", values_from = est) %>% 
    mutate(Parcel = i,
           cfi = mem_fitmeasures_i[1],
           rmsea = mem_fitmeasures_i[2])
  
  # bind across parcels
  est_mem = rbind(est_mem,est_mem_i)
  mem_parm = rbind(mem_parm,mem_parm_i %>% mutate(Parcel = i))
}

## Supplementary File 2 ####
# SF2c: parameters in mem
SF2C = mem_parm %>%                                 
  filter(grepl("var_intra",label) | grepl("var_inter",label)) %>% 
  mutate(label = recode(label,
                        var_intra_r1 = "σ1intra^2",
                        var_intra_r2 = "σ2intra^2",
                        var_inter = "σinter^2",
                        mean_P3 = "GD"
  ))

SF2C_annot = merge(SF2C, annotations %>% select(Parcel,label_Yeo_7,label_Yeo7_short), by = "Parcel")%>% 
  mutate(across(where(is.numeric), round, 5)) %>% 
  select(-c(lhs,op,rhs, pvalue,z))

# SF2c: parameters in mem
write_csv(SF2C_annot, sprintf("%s/SI/SFILE1.csv",wd_oa))

# merge with annotation to stratify across cortical networks
est_mem_annot = merge(est_mem, annotations, by = "Parcel")

# FIGURE 2D ####
box_plot_intra = est_mem_annot %>% 
  ggplot(aes(label_Yeo7_short, 1-ICC2k, fill = label_Yeo7_short))+
  geom_boxplot() +
  geom_hline(yintercept = mean(1-est_mem$ICC2k), linetype = "dashed")+
  ylim(0,1) +
  labs(y = "proportion of \nintra-individual variance",
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

box_plot_inter = est_mem_annot %>% 
  ggplot(aes(label_Yeo7_short, ICC2k, fill = label_Yeo7_short))+
  geom_boxplot() +
  geom_hline(yintercept = mean(est_mem$ICC2k), linetype = "dashed")+
  ylim(0,1) +
  labs(y = "proportion of \ninter-individual variance",
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

## save FIG. 2D####
pdf(sprintf("%s/Figures/06_Fig_2D.pdf",wd_oa), 
     width=6, 
     height=3)
box_plot_intra + box_plot_inter
dev.off()

# clustering of results based on functional network
est_mem_annot %>% 
  group_by(label_Yeo7_short) %>% 
  summarise(mean_ICC2k = psych::fisherz2r(mean(psych::fisherz(ICC2k))))

# statistical differences stratfied across functional network
aov_mod = aov(zICC2k ~ as.factor(label_Yeo7_short),est_mem_annot %>% mutate(zICC2k =psych::fisherz(ICC2k)))
summary(aov_mod)
aov_hsd = TukeyHSD(aov_mod)

#report anova results
report::report(aov_mod)

# range of intra-individual variances
est_mem_annot %>% 
 reframe(range(ICC2k))
est_mem_annot %>% slice(which.min(ICC2k))
est_mem_annot %>% slice(which.max(ICC2k))

# remember that correlation coefficient is a biased estimate for the true correlation
# bias is a function of the reliability of the measure (or inter-individual variance) -> rox,oy = rtx,ty*sqrt(rxxryy)
# so the lower bound for the bias in the correlation with the functional gradient, assuming perfect reliability of strucure
# range of lower bound for bias(r-observed, r-true) 
sqrt(1*mean(est_mem$ICC2k))

# subsample two individuals to make the conceptual inter- intra-individual
# distinction present in Figure 2C
set.seed(42)
sub1 = sample(FC_gradients$Sub,1) 
sub2 = sample(FC_gradients$Sub,1)

#FIGURE 2C####
FC_gradients_ann <- merge(FC_gradients,annotations, "Parcel")
G1_s1_d1 <- FC_gradients_ann %>%
  filter(Sub == sub1) %>%
  rename(region = label_Yeo_7) %>% 
  ggplot() +
  geom_brain(atlas = schaefer7_400, 
             side= "lateral",
             hemi = "left",
             aes(fill = P1_r1),
             color = "#414a4c",
             size = .5,
             show.legend = FALSE) +
  scale_fill_viridis_c()+
  theme_void()

G1_s1_d2 <- FC_gradients_ann %>%
  filter(Sub == sub1) %>%
  rename(region = label_Yeo_7) %>% 
  ggplot() +
  geom_brain(atlas = schaefer7_400, 
             side= "lateral",
             hemi = "left",
             aes(fill = P1_r2),
             color = "#414a4c",
             size = .5,
             show.legend = FALSE) +
  scale_fill_viridis_c()+
  theme_void()+
  theme(legend.title.position = "none")

G1_s2_d1 <- FC_gradients_ann %>%
  filter(Sub == sub2) %>%
  rename(region = label_Yeo_7) %>% 
  ggplot() +
  geom_brain(atlas = schaefer7_400, 
             side= "lateral",
             hemi = "left",
             aes(fill = P1_r1),
             color = "#414a4c",
             size = .5,
             show.legend = FALSE) +
  scale_fill_viridis_c()+
  theme_void()+
  theme(legend.title.position = "none")

G1_s2_d2 <- FC_gradients_ann %>%
  filter(Sub == sub2) %>%
  rename(region = label_Yeo_7) %>% 
  ggplot() +
  geom_brain(atlas = schaefer7_400, 
             side= "lateral",
             hemi = "left",
             aes(fill = P1_r2),
             color = "#414a4c",
             size = .5,
             show.legend = FALSE) +
  scale_fill_viridis_c()+
  theme_void()+
  theme(legend.title.position = "none")

## save FIGURE 2C####
pdf(sprintf("%s/Figures/06_Fig_2C.pdf",wd_oa),
    width = 3,
    height = 3 )
(G1_s1_d1 + G1_s2_d1) / (G1_s1_d2|G1_s2_d2)
dev.off()