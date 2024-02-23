# structural measurement error model to estimate the relatinshio
# between a phenotype P1 with two repeated measure r1, r2
# and a phenotype P2
smem <- "
  # means
    P1_r1~c(mean_P1_r1)*1 
    P1_r2~c(mean_P1_r2)*1 
    P2~c(mean_P2)*1 
    P3~c(mean_P3)*1 
  # measurement error model
    P1_inter =~ NA*P1_r1 + f1*P1_r1 + f1*P1_r2
  # intra-individual variance 
    P1_r1~~var_intra_r1*P1_r1
    P1_r2~~var_intra_r2*P1_r2
  # inter-individual variance 
    P1_inter~~1*P1_inter
    P1_r1~~0*P1_r2
  # variances
    P2~~var_P2*P2
    P3~~var_P3*P3
  # correlations
    P1_inter~~rP1P2*P2
    P1_inter~~rP1P3*P3
    P2~~rP2P3*P3
  "