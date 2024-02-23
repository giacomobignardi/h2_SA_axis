#Author: Giacomo Bignardi
#Description: estimate overlap between two subsamples of HCP
#Program: Multivariate h2 twin measurement model------------------------------------------------------------------------------------------------------------------------------
#clean working environment 
rm(list = ls())
library(tidyverse)

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

# Unrelated sample
SEM_output  =  read_csv(sprintf("%s/%s/07_SMEM_output.csv",wdNOA,wdNOA_output)) %>% 
  filter(label == "rP1P3") %>% 
  mutate(p_adj = p.adjust(pvalue, method = "bonferroni")) %>% 
  filter(p_adj <.05)
  
# Twin sample
CPMEM_output_rA  =  read_csv(sprintf("%s/%s/08_CPMEM_output.csv",wdNOA,wdNOA_output)) %>% 
  filter((label == "rA_P13" )) %>% 
  distinct(Parcel, .keep_all = T)%>% 
  group_by(label) %>% 
  mutate(p_adj = p.adjust(pvalue, method = "bonferroni"))%>% 
  filter(p_adj <.05)

CPMEM_output_rE  =  read_csv(sprintf("%s/%s/08_CPMEM_output.csv",wdNOA,wdNOA_output)) %>% 
  filter((label == "rE_P13")) %>% 
  distinct(Parcel, .keep_all = T)%>% 
  mutate(p_adj = p.adjust(pvalue, method = "bonferroni"))%>% 
  filter(p_adj <.05)

sig_parcel_rAE = unique(c(CPMEM_output_rA %>% pull(Parcel), CPMEM_output_rE %>% pull(Parcel)))

# overlapp 
overlapp = sum((SEM_output %>% pull(Parcel)) %in% (sig_parcel_rAE))/length(sig_parcel_rAE)
length(sig_parcel_rAE)
length(SEM_output %>% pull(Parcel))
