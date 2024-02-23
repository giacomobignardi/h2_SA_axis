#AE UNIVARIATE COMMON PATHWAY MEASUREMENT ERROR MODEL####
# Version: V.0.1
# Author: Giacomo Bignardi
# structural equation model specification for 2 groups
# 1: same-sex monozygotic 
# 2: same-sex dizygotic
#AE MODEL####
cpmem_AE_model <- "
    # means
      P1_r1_T1~c(mean_P1_r1)*1
      P1_r2_T1~c(mean_P1_r2)*1
       
      P1_r1_T2~c(mean_P1_r1)*1
      P1_r2_T2~c(mean_P1_r2)*1
      
    # measurement error model
      P1_inter_T1 =~ NA*P1_r1_T1 + f1*P1_r1_T1 + f1*P1_r2_T1
      P1_inter_T2 =~ NA*P1_r1_T2 + f1*P1_r1_T2 + f1*P1_r2_T2
      
    # intra-individual variance 
      P1_r1_T1~~var_intra_r1*P1_r1_T1
      P1_r2_T1~~var_intra_r2*P1_r2_T1
      P1_r1_T2~~var_intra_r1*P1_r1_T2
      P1_r2_T2~~var_intra_r2*P1_r2_T2
      
    # set residual correlation to 0
      P1_r1_T1~~0*P1_r2_T1
      P1_r1_T2~~0*P1_r2_T2
   
    # AE components
      A_P1_T1=~ P1_inter_T1
      A_P1_T2=~ P1_inter_T2

      E_P1_T1=~ P1_inter_T1
      E_P1_T2=~ P1_inter_T2

    # variances
      A_P1_T1 ~~ varA_P1*A_P1_T1
      A_P1_T2 ~~ varA_P1*A_P1_T2
    
      E_P1_T1 ~~ varE_P1*E_P1_T1
      E_P1_T2 ~~ varE_P1*E_P1_T2

    # scale factor variance to be equal to 1
      varA_P1 + varE_P1 == 1
      
    # fix variances to 0
      P1_inter_T1~~0*P1_inter_T1
      P1_inter_T2~~0*P1_inter_T2
    
    # constraints (twin pair covariances, p = phenotype)
      cov_MZ_P1 == varA_P1
      cov_DZ_P1 == .5*varA_P1

    # covariances cross-twin wihtin-trait
      A_P1_T1 ~~ c((cov_MZ_P1),(cov_DZ_P1))*A_P1_T2
      E_P1_T1 ~~ 0*E_P1_T2

    # set every other covariance to zero
    # within-trait A*G  
      A_P1_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2
      A_P1_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2

    # definitions:estimates
    # phenotypic variance
      varP1 := 1
    
    # compute standard deviations
      sdA_P1 := sqrt(varA_P1)
      sdE_P1 := sqrt(varE_P1)

    # twin-h2
      h2_P1 := varA_P1"
