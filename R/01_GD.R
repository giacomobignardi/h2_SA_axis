#Author: Giacomo Bignardi
#Description: Load and select Geodesic Distances
#------------------------------------------------------------------------------------------------------------------------------
# clean working environment 
rm(list = ls())
library(tidyverse)

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
# Geodesic Distances from Valk S.####
GD = R.matlab::readMat(sprintf("%s/hcp_gd.mat",wd_noa_data))
GD_ids = R.matlab::readMat(sprintf("%s/HCPids.mat",wd_noa_data))

# inclusion index (based on subject inclusion in Valk et al., 2022 [Nat. Comms])
Valk_sub = read_csv(sprintf("%s/subj_inclusion.csv",wd_noa_data), col_names = F) %>% rename(Subject = X1)

# Prepare data for later analysis
GD_id= which(GD_ids$iDs %in% Valk_sub$Subject)
GD_df = cbind(Sub = Valk_sub$Subject,GD$gd.ind[GD_id, 1:400]) %>% as.data.frame()
GD_df_long = GD_df%>%pivot_longer(names_to = "Parcel", 
                                  values_to = "gd", 
                                  c(2:ncol(GD_df))) %>%
  mutate(Parcel = substr(Parcel,2,nchar(Parcel)))%>%
  rename(value = gd) %>% 
  mutate(Parcel = as.numeric(Parcel)) %>% 
  mutate(Parcel = Parcel - 1)

# SAVE geodesic distances####
write_csv(GD_df_long, sprintf("%s/01_GD.csv",wd_noa_output))
