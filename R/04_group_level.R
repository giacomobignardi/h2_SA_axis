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
# Microstructural intensity
MPmi_i  =  read_csv(sprintf("%s/HCP_S1200_MPC_400.csv",wd_noa_data)) %>% rename(Sub = "Var1")

# Geodesic distances
GD_i  =   read_csv(sprintf("%s/01_GD.csv",wd_noa_output))

# Functional gradient
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

## A (2D): t1wt2w####
# This automatise the process and makes making figures easier
t1wt2w_left = t1t2w_template_graph %>%
  rename(region = label_Yeo_7) %>% 
  ggplot() +
  geom_brain(atlas = schaefer7_400, 
             hemi = "left",
             position = position_brain(side  + hemi ~ .), #to stack them
             aes(fill = `t1w/t2w_mi`),
             color = "#414a4c",
             size = .5,
             hemisphere = "left") +
  scale_fill_viridis_c(option = "plasma")+
  labs(subtitle = expression(zT1w/T2w[mi]),
       fill = "Association-Sensorimotor")+
  theme_void(base_size = 14)+
  theme(legend.position = "bottom",
        legend.title.position = "top",
        element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

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

## B (2D): t1wt2w####
# This automatise the process and makes making figures easier

GD_left = GD_template_graph %>%
  rename(region = label_Yeo_7) %>% 
  ggplot() +
  geom_brain(atlas = schaefer7_400, 
             hemi = "left",
             position = position_brain(side  + hemi ~ .), #to stack them
             aes(fill = GD),
             color = "#414a4c",
             size = .5,
             hemisphere = "left") +
  scico::scale_fill_scico(palette = "bilbao", direction = -1)+
  labs(#title = "Geodesic distances",
       subtitle = "zGD",
       fill = "Sensorimotor-Association")+
  theme_void(base_size = 14)+
  theme(legend.position = "bottom",
        legend.title.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))


## C (2D): FCG1####
# template surface 
G1_template_graph = SFs %>%
  mutate(G1 = scale(FC_G1)) %>% #z score for consistenc 
  select(Parcel, label_Yeo_7, G1)


G1_left = G1_template_graph %>%
  rename(region = label_Yeo_7) %>% 
  ggplot() +
  geom_brain(atlas = schaefer7_400, 
             hemi = "left",
             position = position_brain(side  + hemi ~ .), #to stack them
             aes(fill = G1),
             color = "#414a4c",
             size = .5,
             hemisphere = "left") +
  scale_fill_viridis_c()+
  labs(#title = "Functional gradient",
       subtitle = expression(zFC[G1]),
       fill = "Sensorimotor-Association")+
  theme_void(base_size = 14)+
  theme(legend.position = "bottom",
        legend.title.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))

## save FIG. 1A####
pdf(sprintf("%s/Figures/04_Fig_1AC.pdf",wd_oa),
    width = 12,
    height = 4 )
(t1wt2w_left + plot_spacer() + GD_left + plot_spacer() + G1_left + plot_layout(widths = c(4, 1 ,4, 1, 4))) 
dev.off()

## save legend for FIG. 1####
pdf(sprintf("%s/Figures/04_legend1.pdf",wd_oa),
    width = 4,
    height = 4 )
plot(schaefer7_400, color = "gray",hemi = "left",position = position_brain(side  + hemi ~ .), size = .5,) + 
  theme_void() +
  theme(legend.position = "none", plot.title = element_blank())
dev.off()


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
pdf(sprintf("%s/Figures/04_Fig_1D.pdf",wd_oa),
    width = 5,
    height =2 )
(corplot_left|corplot_right) + plot_layout(guides = "collect")  & theme(legend.position = "none")
dev.off()

#correlation between different S-F axis
cormat  = round(cor(SFs %>% select(`t1w/t2w_mi`,GD,FC_G1), method = "spearman"),2)
cor_mp = cor.test(SFs$`t1w/t2w_mi`,SFs$FC_G1, method = "spearman")
cor_gd = cor.test(SFs$`GD`,SFs$FC_G1, method = "spearman")
