#Author: Giacomo Bignardi
#Description: Load and select Geodesic Distances
#------------------------------------------------------------------------------------------------------------------------------
# clean working environment 
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

# load
# Geodesic Distances from Valk S.####
GD = R.matlab::readMat(sprintf("%s/GD/hcp_gd.mat",wdNOA_rawData))
GD_ids = R.matlab::readMat(sprintf("%s/GD/HCPids.mat",wdNOA_rawData))

# inclusion index (based on subject inclusion in Valk et al., 2022 [Nat. Comms])
Valk_sub = read_csv(sprintf("%s%s/subjects_for_giaco.csv",wdNOA,wdNOA_Data), col_names = F) %>% rename(Subject = X1)

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
write_csv(GD_df_long, sprintf("%s/%s/01_GD.csv",wdNOA,wdNOA_output))
