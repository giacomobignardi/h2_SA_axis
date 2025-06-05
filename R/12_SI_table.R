#Author: Giacomo Bignardi
# adapted from https://github.com/giacomobignardi/model-etfd-iamdate/blob/main/SI/06_table5_est.csv
# prepare supplementary file
library(openxlsx)
library(tidyverse)
rm(list = ls())

# set open access working directories
wd_oa = getwd()

# load sup tables
st1 <- read.csv(sprintf("%s/SI/SFILE1.csv", wd_oa))
st2 <- read.csv(sprintf("%s/SI/SFILE2.csv", wd_oa))
st3 <- read.csv(sprintf("%s/SI/SFILE3.csv", wd_oa))
st4 <- read.csv(sprintf("%s/SI/SFILE4.csv", wd_oa))
st5 <- read.csv(sprintf("%s/SI/SFILE5.csv", wd_oa))

# rename labels where needed
st5$label <- "rA_gdg1"

# st1 decomposition of g1 variance in inter- and intra-individual variance
# st2 multivariate measurement error sem
# st3 multivariate measurement error twin-informed sem
# st4 genetic correlations
# st5 environemntal correlations

# create a reference table
st0 <- data.frame(sheet = c(paste0("S", seq(0,5,1))),
                  content = c(
                    "Overview of sheet contents: provides information about the contents of all other sheets",
                    "G1 (Variance Decomposition): Splits G1 variance into inter- and intra-individual components – see Model (Fig. 2C), Estimates (Fig. 2D)",
                    "Multivariate SEM (Measurement Error): Structural equation model accounting for measurement error – see Model (Fig. 7A), Estimates (Fig. 3A)",
                    "Multivariate SEM (Twin-Informed, Measurement Error): Twin-informed SEM incorporating measurement error – see Model (Fig. 7B–C), Estimates (Fig. 4B, Fig. 5A–B)",
                    "Genetic Correlations: Microstructure & Functional Gradient: Correlations modeled between microstructure and functional gradient – see Model (Fig. 7C), Estimates (Fig. 5C)",
                    "Genetic Correlations: Geodesic Distance & Functional Gradient: Correlations modeled between geodesic distances and functional gradient – see Model (Fig. 7C), Estimates (Fig. 5C)"
                  ))
st_list = list(st0, st1, st2, st3, st4, st5)

# create a workbook
wb <- createWorkbook()

# create a list for sheet names
sn <- c("0.overview", 
        "1.measure_model", 
        "2.phenotypic_model", 
        "3.twin_model", 
        "4.rA_mig1", 
        "5.rA_gdg1")

# add each correlation matrix as a sheet
for (i in 1:length(st_list)) {
  sheet_name <- sn[i]
  addWorksheet(wb, sheet_name)  # Add a new sheet
  writeData(wb, sheet = sheet_name, x = st_list[[i]]) 
}

# add note to sup 5
note_sup5 <- "Notes. Lables match the parameter names in Figures"
writeData(wb, "1.measure_model", note_sup5, startCol = 1, startRow =  nrow(st1)+5)
writeData(wb, "2.phenotypic_model", note_sup5, startCol = 1, startRow =  nrow(st2)+5)
writeData(wb, "3.twin_model", note_sup5, startCol = 1, startRow =  nrow(st3)+5)

# save Supplementary file
saveWorkbook(wb,sprintf("%s/SI/Supplementary_Table.xlsx", wd_oa), overwrite = T)
