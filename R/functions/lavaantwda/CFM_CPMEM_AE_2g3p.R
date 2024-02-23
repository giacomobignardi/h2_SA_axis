#AE MULTIVARIATE COMMON PATHWAY MEASUREMENT ERROR MODEL####
# Version: V.0.1
# Author: Giacomo Bignardi
# structural equation model specification for 2 groups
# 1: same-sex monozygotic 
# 2: same-sex dizygotic
#AE MODEL####
# cfm: common factor model
cfm_cpmem_AE_model <-"
    # means
      P1_r1_T1~c(mean_P1_r1)*1 
      P1_r2_T1~c(mean_P1_r2)*1 
      P1_r1_T2~c(mean_P1_r1)*1 
      P1_r2_T2~c(mean_P1_r2)*1 
  
      P2_T1~c(mean_P2)*1 
      P2_T2~c(mean_P2)*1 
      P3_T1~c(mean_P3)*1 
      P3_T2~c(mean_P3)*1 
    
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
      
    # scale factor variance to be equal to 1
      varA_P1 + varE_P1 == 1

    #AE components
      A_P1_T1=~ P1_inter_T1
      A_P1_T2=~ P1_inter_T2
      A_P2_T1=~ P2_T1
      A_P2_T2=~ P2_T2
      A_P3_T1=~ P3_T1
      A_P3_T2=~ P3_T2   
  
      E_P1_T1=~ P1_inter_T1
      E_P1_T2=~ P1_inter_T2
      E_P2_T1=~ P2_T1
      E_P2_T2=~ P2_T2
      E_P3_T1=~ P3_T1
      E_P3_T2=~ P3_T2   

    #variances
      A_P1_T1 ~~ varA_P1*A_P1_T1
      A_P1_T2 ~~ varA_P1*A_P1_T2
      A_P2_T1 ~~ varA_P2*A_P2_T1
      A_P2_T2 ~~ varA_P2*A_P2_T2
      A_P3_T1 ~~ varA_P3*A_P3_T1
      A_P3_T2 ~~ varA_P3*A_P3_T2
  
      E_P1_T1 ~~ varE_P1*E_P1_T1
      E_P1_T2 ~~ varE_P1*E_P1_T2
      E_P2_T1 ~~ varE_P2*E_P2_T1
      E_P2_T2 ~~ varE_P2*E_P2_T2
      E_P3_T1 ~~ varE_P3*E_P3_T1
      E_P3_T2 ~~ varE_P3*E_P3_T2

    #fix remaining variance to 0
      P1_inter_T1~~0*P1_inter_T1
      P1_inter_T2~~0*P1_inter_T2
      P2_T1~~0*P2_T1
      P2_T2~~0*P2_T2
      P3_T1~~0*P3_T1
      P3_T2~~0*P3_T2

    #constraints (twin pair covariances, p = phenotype)
      cov_MZ_P1 == varA_P1
      cov_MZ_P2 == varA_P2
      cov_MZ_P3 == varA_P3
  
      cov_DZ_P1 == .5*varA_P1
      cov_DZ_P2 == .5*varA_P2
      cov_DZ_P3 == .5*varA_P3
    
    #covariances cross-twin wihtin-trait
      A_P1_T1 ~~ c((cov_MZ_P1),(cov_DZ_P1))*A_P1_T2
      A_P2_T1 ~~ c((cov_MZ_P2),(cov_DZ_P2))*A_P2_T2
      A_P3_T1 ~~ c((cov_MZ_P3),(cov_DZ_P3))*A_P3_T2 
  
      E_P1_T1 ~~ 0*E_P1_T2
      E_P2_T1 ~~ 0*E_P2_T2
      E_P3_T1 ~~ 0*E_P3_T2

    #covariances within-twin cross-trait
      A_P1_T1 ~~ (covA_P12)*A_P2_T1 + (covA_P13)*A_P3_T1 
      A_P2_T1 ~~ (covA_P23)*A_P3_T1 
      
      E_P1_T1 ~~ (covE_P12)*E_P2_T1 + (covE_P13)*E_P3_T1 
      E_P2_T1 ~~ (covE_P23)*E_P3_T1 
      
      A_P1_T2 ~~ (covA_P12)*A_P2_T2 + (covA_P13)*A_P3_T2 
      A_P2_T2 ~~ (covA_P23)*A_P3_T2 
      
      E_P1_T2 ~~ (covE_P12)*E_P2_T2 + (covE_P13)*E_P3_T2
      E_P2_T2 ~~ (covE_P23)*E_P3_T2 
   
    #constrains
      cov_MZ_P12 == covA_P12
      cov_MZ_P13 == covA_P13
      cov_MZ_P23 == covA_P23
      
      cov_DZ_P12 == .5*covA_P12
      cov_DZ_P13 == .5*covA_P13
      cov_DZ_P23 == .5*covA_P23
    
    #covariances cross-twin cross-trait
      A_P1_T1 ~~ c((cov_MZ_P12),(cov_DZ_P12))*A_P2_T2 + c((cov_MZ_P13),(cov_DZ_P13))*A_P3_T2 
      A_P2_T1 ~~ c((cov_MZ_P12),(cov_DZ_P12))*A_P1_T2 + c((cov_MZ_P23),(cov_DZ_P23))*A_P3_T2 
      A_P3_T1 ~~ c((cov_MZ_P13),(cov_DZ_P13))*A_P1_T2 + c((cov_MZ_P23),(cov_DZ_P23))*A_P2_T2 
  
      E_P1_T1 ~~ 0*E_P2_T2 + 0*E_P3_T2 
      E_P2_T1 ~~ 0*E_P1_T2 + 0*E_P3_T2 
      E_P3_T1 ~~ 0*E_P1_T2 + 0*E_P2_T2 

    #set every other covariance to zero
    #within-trait A*G  
      A_P1_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2
      A_P1_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2
      
      A_P2_T1 ~~ 0*E_P2_T1 + 0*E_P2_T2
      A_P2_T2 ~~ 0*E_P2_T1 + 0*E_P2_T2
          
      A_P3_T1 ~~ 0*E_P3_T1 + 0*E_P3_T2
      A_P3_T2 ~~ 0*E_P3_T1 + 0*E_P3_T2

    
    #cross-trait A*G
      A_P1_T1 ~~ 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2  
      A_P1_T2 ~~ 0*E_P2_T1 + 0*E_P2_T2 + 0*E_P3_T1 + 0*E_P3_T2 
  
      A_P2_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2 
      A_P2_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P3_T1 + 0*E_P3_T2 
      
      A_P3_T1 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2 
      A_P3_T2 ~~ 0*E_P1_T1 + 0*E_P1_T2 + 0*E_P2_T1 + 0*E_P2_T2
    
    # definitions:estimates
    # phenotypic variance
      varP1 := (varA_P1+varE_P1)
      varP2 := (varA_P2+varE_P2)
      varP3 := (varA_P3+varE_P3)

    # compute genetic covariance (bivariate heritabilities - optional)
      # bh2_P12 := covA_P12/ (covA_P12 + covE_P12)
      # bh2_P13 := covA_P13/ (covA_P13 + covE_P13)
      # bh2_P23 := covA_P23/ (covA_P23 + covE_P23)

    # compute standard deviations (A)
      sdA_P1 := sqrt(varA_P1)
      sdA_P2 := sqrt(varA_P2)
      sdA_P3 := sqrt(varA_P3)

    # compute genetic correlations
      rA_P12 := covA_P12/(sdA_P1*sdA_P2)
      rA_P13 := covA_P13/(sdA_P1*sdA_P3)
      rA_P23 := covA_P23/(sdA_P2*sdA_P3)

    # compute standard deviations (E)
      sdE_P1 := sqrt(varE_P1)
      sdE_P2 := sqrt(varE_P2)
      sdE_P3 := sqrt(varE_P3)

    # compute environmental correlations
      rE_P12 := covE_P12/(sdE_P1*sdE_P2)
      rE_P13 := covE_P13/(sdE_P1*sdE_P3)
      rE_P23 := covE_P23/(sdE_P2*sdE_P3)

    # h2-twin
      h2_P1 := varA_P1/varP1
      h2_P2 := varA_P2/varP2
      h2_P3 := varA_P3/varP3
  "
