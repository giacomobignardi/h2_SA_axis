#Author: Giacomo Bignardi
#Description: move to the individual level of analysis
#Program: group level analyses------------------------------------------------------------------------------------------------------------------------------
#clean working environment 
rm(list = ls())
library(tidyverse)
library(reshape2)
library(patchwork)
library(ggseg)
library(ggseg3d)
library(ggsegSchaefer)
library(ggridges)

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


#load
# Microstructural intensity
MPmi_i  =  read_csv(sprintf("%s/HCP_S1200_MPC_400.csv",wd_noa_data)) %>% rename(Sub = "Var1")

# Geodesic distances
GD_i  =   read_csv(sprintf("%s/01_GD.csv",wd_noa_output))

# Functional gradient
FC_g1_i =  read_csv(sprintf("%s/03_GFC_i.csv",wd_noa_output))

# Functional gradient(group-level)
FC_g1 =  read_csv(sprintf("%s/03_FC_m_G.csv",wd_noa_output))

#inclusion index 
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
annotations$label_Yeo7_short = network_7_yeo[,1]
annotations = annotations %>% rename(Parcel = parcel_Yeo)

#INDIVIDUALS####
#get individual values for mp
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

# get individual values for g1
FC_g1_i = FC_g1_i%>%
  filter(Sub %in% notTwin_sub$Subject) %>%
  select(Sub,Parcel,`0`)%>%
  rename(G1 = `0`)

# tidy FCG1 for ease of interpretation
FC_g1 = FC_g1 %>%
  rename(FC_G1 = `0`) %>% 
  mutate(Parcel = c(1:400)) %>% 
  select(Parcel, FC_G1)

# put all of them togheter
SFs_i_list = list(MPmi_i, GD_i,FC_g1_i)

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

# put  FCG1 group level and annotation all together
SFs_list = list(FC_g1, annotations)

# merge all data frames in list
SFs = SFs_list %>% 
  reduce(full_join, by='Parcel') %>% 
  mutate(label_Yeo7_short = factor(label_Yeo7_short, levels = rev(c("Default", 
                                                                    "Limbic", 
                                                                    "Cont",
                                                                    "SalVentAttn",
                                                                    "DorsAttn",
                                                                    "SomMot",
                                                                    "Vis"))))
#FIG. 2A####
var_FCG1_p = SFs_i_final %>%
  ggplot(aes(Parcel, scale(G1))) +
  scale_fill_manual(values = c('#9F53AA',
                                        '#7A9ABD',
                                        '#3d8043',
                                        '#b584cf',
                                        '#e8a633',
                                        '#F4FEC8',
                                        '#D8707A'))+
                                          geom_point(aes(fill = label_Yeo7_short),color = "white",pch = 21) +
  geom_point(data= SFs %>%  rename(G1 = FC_G1), aes(Parcel, scale(G1),fill = label_Yeo7_short),color = "black",pch = 21) +
  labs(fill = "Yeo-Krienen\n7 networks",
       x = "Schaefer parcel",
       y = bquote(zFC[G1]),
       subtitle = "")+
  theme_classic() +
  theme(legend.position = "none")


## save FIG. 2A####
tiff(sprintf("%s/Figures/05_Fig_2A.tiff",wd_oa), 
     units="in", 
     width=7.75, 
     height=2, 
     res=300)
var_FCG1_p
dev.off()

