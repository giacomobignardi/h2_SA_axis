# measurement error model for a phenotype P with two repeated measure r1, r2
mem <- "
  # measurement error model
    P1_inter =~ 1*P1_r1 + 1*P1_r2
  # intra-individual variance 
    P1_r1~~var_intra_r1*P1_r1
    P1_r2~~var_intra_r2*P1_r2
  # inter-individual variance 
    P1_inter~~var_inter*P1_inter
    P1_r1~~0*P1_r2
  # Intra-Class-Correlation (2,1)
    ICC2 := var_inter / (var_inter + (var_intra_r1+ var_intra_r2)/2)
  # Intra-Class-Correlation (2,K) sperman-brown correction
    k := 2
    ICC2k := (k*ICC2)/(1 +(k-1)*ICC2)
  "
# Spearman-Brown prediction (2*ICC_2_full)/(1 +(2-1)*ICC_2_full)