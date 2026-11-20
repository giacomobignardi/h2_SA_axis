#Author: Giacomo Bignardi
#Description: estimate overlap between two subsamples of HCP
#Program: Multivariate h2 twin measurement model------------------------------------------------------------------------------------------------------------------------------
#clean working environment 
library(tidyverse)

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

# Unrelated sample
SEM_output  =  read_csv(sprintf("%s/07_SMEM_output.csv",wd_noa_output)) %>% 
  filter(label == "rP1P3") %>% 
  mutate(p_adj = p.adjust(pvalue, method = "bonferroni")) %>% 
  filter(p_adj <.05)
  
# Twin sample
CPMEM_output_rA  =  read_csv(sprintf("%s/09_CPMEM_output.csv",wd_noa_output)) %>% 
  filter((label == "rA_P13" )) %>% 
  distinct(Parcel, .keep_all = T)%>% 
  group_by(label) %>% 
  mutate(p_adj = p.adjust(pvalue, method = "bonferroni"))%>% 
  filter(p_adj <.05, admissible, Parcel != 114)   # same exclusions as in 09

CPMEM_output_rE  =  read_csv(sprintf("%s/09_CPMEM_output.csv",wd_noa_output)) %>% 
  filter((label == "rE_P13")) %>% 
  distinct(Parcel, .keep_all = T)%>% 
  mutate(p_adj = p.adjust(pvalue, method = "bonferroni"))%>% 
  filter(p_adj <.05, admissible, Parcel != 114)   # same exclusions as in 09

sig_parcel_rAE = unique(c(CPMEM_output_rA %>% pull(Parcel), CPMEM_output_rE %>% pull(Parcel)))

# proportion of twin-significant parcels that are also significant in the non-twin sample
s = sum(sig_parcel_rAE %in% (SEM_output %>% pull(Parcel))) / length(sig_parcel_rAE)

# overlap summary
c(n_nontwin_sig = nrow(SEM_output),
  n_twin_sig    = length(sig_parcel_rAE),
  n_overlap     = sum(sig_parcel_rAE %in% (SEM_output %>% pull(Parcel))),
  pct_overlap   = 100 * s)
