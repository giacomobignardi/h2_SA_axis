#Author: Giacomo Bignardi
#Description: explore average trends of S-A axis of organization across different neurobiological properties 
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

# load
# Microstructural intensity
MPmi_i  =  read_csv(sprintf("%s/t1t2w/HCP_S1200_MPC_400.csv",wdNOA_rawData)) %>% rename(Sub = "Var1")

# Geodesic distances
GD_i  =   read_csv(sprintf("%s/%s/01_GD.csv",wdNOA,wdNOA_output))

# Functional gradient
FC_g1 =  read_csv(sprintf("%s/%s/03_FC_m_G.csv", wdNOA,wdNOA_output))

#inclusion index 
notTwin_sub =  read_csv(sprintf("%s/%s/00_nottwin_ids.csv", wdNOA,wdNOA_output))

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

# GROUP-LEVEL####
## t1w/t2w####
# obtain average t1w/t2w values across participants
MPmi =  MPmi_i%>%
  filter(Sub %in% notTwin_sub$Subject) %>% 
  pivot_longer(names_to = "Parcel", values_to = "value", c(2:ncol(MPmi_i)))%>%
  mutate(Parcel = substr(Parcel,6,nchar(Parcel)))%>%
  mutate(Parcel = as.numeric(Parcel)) %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  group_by(Parcel) %>%
  summarise(`t1w/t2w_mi` = mean(value, na.rm = T))

## GD####
# obtain average GD values across participants
GD = GD_i %>% 
  filter(Sub %in% notTwin_sub$Subject) %>% 
  group_by(Parcel) %>% 
  mutate(value = ifelse(value == 0, NA, value)) %>%
  summarise(GD = mean(value, na.rm = T)) 

## FCG1####
# tidy FCG1 for ease of interpretation
FC_g1 = FC_g1 %>%
  rename(FC_G1 = `0`) %>% 
  mutate(Parcel = c(1:400)) %>% 
  select(Parcel, FC_G1)

# put MP,GD and FCG1 all together
SFs_list = list(MPmi, GD, FC_g1, annotations %>% rename(Parcel = parcel_Yeo)) #ADD classes and yeo_7###

# merge all data frames in list
SFs = SFs_list %>% 
  reduce(full_join, by='Parcel') %>% 
  mutate(label_Yeo7_short = factor(label_Yeo7_short, levels = rev(c("Default", 
                                                                 "Limbic", 
                                                                 "Cont",
                                                                 "SalVentAttn",
                                                                 "DorsAttn",
                                                                 "SomMot",
                                                                 "Vis"))),
         #if needed during revision 
         classes_6_types = recode(class_types, 
                                  '7' = "Agranular",
                                  '6' = "Disgranular",
                                  '5' = "Eulaminate-I",
                                  '4' = "Eulaminate-II",
                                  '3' = "Eulaminate-III",
                                  '2' = "Koniocortical",
                                  '1' = "mask",
                                  '8' = "Other")) 
  

# FIG 1 A-B####
## A: t1wt2w####
# template surface 
t1t2w_template_graph = SFs %>%
  mutate(`t1w/t2w_mi` = scale(`t1w/t2w_mi`)) %>% #z score for consistency
  mutate(`t1w/t2w_mi`  = ifelse(classes_6_types == "mask", NA, `t1w/t2w_mi` )) %>%  #remove mask from graph
  select(Parcel, label_Yeo_7, `t1w/t2w_mi`)

# range of parcel-wise values
min(t1t2w_template_graph$`t1w/t2w_mi`, na.rm = T)
max(t1t2w_template_graph$`t1w/t2w_mi`, na.rm = T)

# adapt to plot
schaefer7_400_t1wt2w = schaefer7_400_3d %>% #
  filter(hemi == "left")%>%
  mutate(ggseg_3d = map(ggseg_3d, ~ .x %>% 
                          full_join(t1t2w_template_graph%>%
                                      filter(Parcel < 201) %>% #left hemisphere
                                      rename(region = label_Yeo_7) %>% 
                                      select(region,`t1w/t2w_mi`),
                                    by = "region")))
# plot parcelations using ggseg3d
t1wt2w_left_medial = 
  ggseg3d(
    atlas = schaefer7_400_t1wt2w,
    colour = "t1w/t2w_mi", text = "t1w/t2w_mi",
    palette =  c( "#0D0887FF"= min(t1t2w_template_graph$`t1w/t2w_mi`, na.rm = T), "#CC4678FF"  = 0, "#F0F921FF" =max(t1t2w_template_graph$`t1w/t2w_mi`, na.rm = T)),
    surface = "inflated",
    hemisphere = "left")%>% 
  remove_axes()%>%
  pan_camera("left medial")
t1wt2w_left_lateral = 
  ggseg3d(
    atlas = schaefer7_400_t1wt2w,
    colour = "t1w/t2w_mi", text = "t1w/t2w_mi",
    palette =  c( "#0D0887FF"= min(t1t2w_template_graph$`t1w/t2w_mi`, na.rm = T), "#CC4678FF"  = 0, "#F0F921FF" =max(t1t2w_template_graph$`t1w/t2w_mi`, na.rm = T)),
    surface = "inflated",
    hemisphere = "left")%>% 
  remove_axes()%>%
  pan_camera("left lateral")

## B: GD####
# template surface 
GD_template_graph = SFs %>%
  mutate(GD = scale(GD)) %>% #z score for consistenc 
  select(Parcel, label_Yeo_7, GD)

# min and max for legend
min(GD_template_graph$GD)
max(GD_template_graph$GD)

# adapt to plot
schaefer7_400_GD = schaefer7_400_3d %>% #
  filter(hemi == "left")%>%
  mutate(ggseg_3d = map(ggseg_3d, ~ .x %>% 
                          full_join(GD_template_graph%>%
                                      filter(Parcel < 201) %>% 
                                      rename(region = label_Yeo_7) %>% 
                                      select(region,GD),
                                    by = "region")))
# plot with ggseg3d
GD_left_medial =
  ggseg3d(
    atlas = schaefer7_400_GD,
    colour = "GD", text = "GD",
    palette =  c( "#f0c27b"= min(GD_template_graph$GD, na.rm = T), "#9D6B63"  = 0, "#4b1248" =max(GD_template_graph$GD, na.rm = T)),
    surface = "inflated",
    hemisphere = "left")%>%
  remove_axes()%>%
  pan_camera("left medial")
GD_left_lateral =
  ggseg3d(
    atlas = schaefer7_400_GD,
    colour = "GD", text = "GD",
    palette =  c( "#f0c27b"= min(GD_template_graph$GD, na.rm = T), "#9D6B63"  = 0, "#4b1248" =max(GD_template_graph$GD, na.rm = T)),
    surface = "inflated",
    hemisphere = "left")%>%
  remove_axes()%>%
  pan_camera("left lateral")

# FIGURE 1D####
# correlation between Microstructural intensity and Functional gradient
## left: t1w/t2w####
corplot_left = SFs %>% 
  mutate(zFC_G1 = scale(FC_G1), `zt1w/t2w_mi` = scale(`t1w/t2w_mi`)) %>%
  ggplot(aes(`zt1w/t2w_mi`, zFC_G1)) +
  scale_fill_manual(values = c('#9F53AA',
                                        '#7A9ABD',
                                        '#3d8043',
                                        '#b584cf',
                                        '#e8a633',
                                        '#F4FEC8',
                                        '#D8707A'))+
                                          geom_point(aes(fill = label_Yeo7_short), pch = 21, size = 2, color = "lightGray") +
  geom_smooth(method = "lm", color = "gray")+
  theme_classic() +
  labs(fill = "Yeo-Krienen\n7 networks",
       x = bquote(zT1w/T2w[mi]),
       y = bquote(zFC[G1]),
       subtitle = "")+
  ylim(c(-1.7, 1.9))

## right : GD####            
corplot_right =SFs %>% mutate(zGD = scale(GD), zFC_G1 = scale(FC_G1)) %>% 
  ggplot(aes(zGD, zFC_G1)) +
  scale_fill_manual(values = c('#9F53AA',
                                        '#7A9ABD',
                                        '#3d8043',
                                        '#b584cf',
                                        '#e8a633',
                                        '#F4FEC8',
                                        '#D8707A'))+
                                          geom_point(aes(fill = label_Yeo7_short), pch = 21, size = 2, color ="lightGray") +
  geom_smooth(method = "lm", color = "gray")+
  theme_classic() +
  labs(fill = "Yeo-Krienen\n7 networks",
       y = bquote(zFC[G1]),
       subtitle = "")+
  ylim(c(-1.7, 1.9))

## save FIG. 1D####
pdf(sprintf("%s/%s/04_fig1D.pdf",wdOA,wdNOA_ImageOutput),
    width = 5,
    height =2 )
(corplot_left|corplot_right) + plot_layout(guides = "collect")  & theme(legend.position = "none")
dev.off()

#correlation between different S-F axis
cormat  = round(cor(SFs %>% select(`t1w/t2w_mi`,GD,FC_G1), method = "spearman"),2)
cor_mp = cor.test(SFs$`t1w/t2w_mi`,SFs$FC_G1, method = "spearman")
cor_gd = cor.test(SFs$`GD`,SFs$FC_G1, method = "spearman")

# ##Plot correlation matrix####
# get_lower_tri = function(cormat){
#   cormat[upper.tri(cormat, diag = T)] <- NA
#   return(cormat)
# }
# # melted_cormat =  melt(get_lower_tri(cormat), na.rm = TRUE)
# # melted_cormat = melted_cormat %>% mutate(
# #   Var1 = factor(Var1, levels = rev(c("t1w/t2w_mi", "MPC_g1", "GD", "FC_G1"))),
# #   Var2 = factor(Var2, levels = c("t1w/t2w_mi", "MPC_g1", "GD", "FC_G1")))
# 
# melted_cormat =  melt(get_lower_tri(cormat), na.rm = TRUE)
# melted_cormat = melted_cormat %>% mutate(
#   Var1 = factor(Var1, levels = rev(c("t1w/t2w_mi", "GD", "FC_G1"))),
#   Var2 = factor(Var2, levels = c("t1w/t2w_mi", "GD", "FC_G1")))
# 
# 
# p_corplot = ggplot(data = melted_cormat, aes(Var1,Var2, fill = value))+
#   geom_tile(color = "white")+
#   scale_fill_viridis_c(option = "inferno", limits = c(-1, 1))+
#   theme_classic()+ 
#   coord_fixed()+
#   labs(fill = "r",
#        x = "",
#        y = "")+
#   geom_text(aes(Var1, Var2, label = value), color = "white", size = 3.5)
# 

# ####NEED TO RE-WRITE FOR CONSISTENCY####
# ACE_Lavaan_mp_cov = read_csv(sprintf("%s/%s/02.2_ACE_Lavaan_mp_cov.csv", wdOA,wdOA_output)) %>% 
#   select(parcel, classes_6_types, label_Yeo_7, label_Yeo7_short) %>% 
#   rename(Parcel = parcel) %>% 
#   mutate(classes_6_types  = factor(classes_6_types, levels =  rev(c("Agranular",
#                                                                     "Disgranular",
#                                                                     "Eulaminate-I",
#                                                                     "Eulaminate-II",
#                                                                     "Eulaminate-III",
#                                                                     "Koniocortical",
#                                                                     "mask",
#                                                                     "Other")))) %>% 
#   mutate(label_Yeo7_short  = factor(label_Yeo7_short, levels =  c("Vis",
#                                                                     "SomMot",
#                                                                     "DorsAttn",
#                                                                     "SalVentAttn",
#                                                                     "Cont",
#                                                                     "Limbic",
#                                                                     "Default")))
# 
# 
# brain_cortical_types =  ggplot(ACE_Lavaan_mp_cov %>% 
#                                  filter(classes_6_types != "mask") %>% 
#                                  rename(region = label_Yeo_7)%>%
#                                  select(region,  classes_6_types) %>%
#                                  as.data.frame()) + geom_brain(atlas = ggsegSchaefer::schaefer7_400,
#                                                                aes(fill = classes_6_types), size=.5,  position = position_brain())+
#   scale_fill_manual(values = 
#                       c(
#                     "#7D8BAE",
#                     "#83A7B0",
#                     "#A6BA9F",
#                     "#DDC592",
#                     "#FFD0C5",
#                     "#FFE2FD",
#                     "#808080"))+
#                       labs(fill = "cortical types")+
#   theme_void() +
#   theme(legend.position = "bottom")
# 
# 
# brain_functional_networks =  ggplot(ACE_Lavaan_mp_cov %>% 
#                                  rename(region = label_Yeo_7)%>%
#                                  select(region,  label_Yeo7_short) %>%
#                                  as.data.frame()) + geom_brain(atlas = ggsegSchaefer::schaefer7_400,
#                                                                aes(fill = label_Yeo7_short), size=.5,  position = position_brain())+
#   scale_fill_manual(values = c('#9F53AA',
#                                         '#7A9ABD',
#                                         '#3d8043',
#                                         '#b584cf',
#                                         '#e8a633',
#                                         '#F4FEC8',
#                                         '#D8707A'))+
#                                           labs(fill = "Yeo-Krienen\n7 Networks")+
#   theme_void() +
#   theme(legend.position = "bottom")
# 
# 
# 
# (brain_functional_networks/brain_cortical_types) & theme(legend.position = "bottom")
# 
# 
# 
# pdf(sprintf("%s/05_images/image/17_S1_yeo_corticaltypes.pdf",wdNOA),
#     width = 7,
#     height =7 )
# (brain_functional_networks/brain_cortical_types) & theme(legend.position = "bottom")
# dev.off()
# 
# ggseg3d(atlas = schaefer7_400_3d, surface = "inflated", hemisphere = "left") %>% 
#   remove_axes()%>%
#   pan_camera("left medial")

#FIGURE 1 A-C####
##t1wt2w####
#template surface 
t1t2w_template_graph = SFs %>%
  mutate(`t1w/t2w_mi` = scale(`t1w/t2w_mi`)) %>% #z score for consistency
  mutate(`t1w/t2w_mi`  = ifelse(classes_6_types == "mask", NA, `t1w/t2w_mi` )) %>%  #remove mask from graph %>% 
  select(Parcel, label_Yeo_7, `t1w/t2w_mi`)
#adapt to plot
schaefer7_400_t1wt2w = schaefer7_400_3d %>% #
  filter(hemi == "left")%>%
  mutate(ggseg_3d = map(ggseg_3d, ~ .x %>% 
                          full_join(t1t2w_template_graph%>%
                                      filter(Parcel < 201) %>% 
                                      rename(region = label_Yeo_7) %>% 
                                      select(region,`t1w/t2w_mi`),
                                    by = "region")))
#plot, scheafer parcelations ggseg3d
t1wt2w_left_medial = 
  ggseg3d(
    atlas = schaefer7_400_t1wt2w,
    colour = "t1w/t2w_mi", text = "t1w/t2w_mi",
    palette =  c( "#0D0887FF"= min(t1t2w_template_graph$`t1w/t2w_mi`, na.rm = T), "#CC4678FF"  = 0, "#F0F921FF" =max(t1t2w_template_graph$`t1w/t2w_mi`, na.rm = T)),
    surface = "inflated",
    hemisphere = "left")%>% 
  remove_axes()%>%
  pan_camera("left medial")
t1wt2w_left_lateral = 
  ggseg3d(
    atlas = schaefer7_400_t1wt2w,
    colour = "t1w/t2w_mi", text = "t1w/t2w_mi",
    palette =  c( "#0D0887FF"= min(t1t2w_template_graph$`t1w/t2w_mi`, na.rm = T), "#CC4678FF"  = 0, "#F0F921FF" =max(t1t2w_template_graph$`t1w/t2w_mi`, na.rm = T)),
    surface = "inflated",
    hemisphere = "left")%>% 
  remove_axes()%>%
  pan_camera("left lateral")

##GD####
#template surface 
GD_template_graph = SFs %>%
  mutate(GD = scale(GD)) %>% #z score for consistenc 
  select(Parcel, label_Yeo_7, GD)
#adapt to plot
schaefer7_400_GD = schaefer7_400_3d %>% #
  filter(hemi == "left")%>%
  mutate(ggseg_3d = map(ggseg_3d, ~ .x %>% 
                          full_join(GD_template_graph%>%
                                      filter(Parcel < 201) %>% 
                                      rename(region = label_Yeo_7) %>% 
                                      select(region,GD),
                                    by = "region")))
#plot, shcafer parcelations ggseg3d
GD_left_medial =
  ggseg3d(
    atlas = schaefer7_400_GD,
    colour = "GD", text = "GD",
    palette =  c( "#f0c27b"= min(GD_template_graph$GD, na.rm = T), "#9D6B63"  = 0, "#4b1248" =max(GD_template_graph$GD, na.rm = T)),
    surface = "inflated",
    hemisphere = "left")%>%
  remove_axes()%>%
  pan_camera("left medial")
GD_left_lateral =
  ggseg3d(
    atlas = schaefer7_400_GD,
    colour = "GD", text = "GD",
    palette =  c( "#f0c27b"= min(GD_template_graph$GD, na.rm = T), "#9D6B63"  = 0, "#4b1248" =max(GD_template_graph$GD, na.rm = T)),
    surface = "inflated",
    hemisphere = "left")%>%
  remove_axes()%>%
  pan_camera("left lateral")