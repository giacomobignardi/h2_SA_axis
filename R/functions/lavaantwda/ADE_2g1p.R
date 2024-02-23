#ADE MODEL####
# Version: V.0.1
# Author: Giacomo Bignardi
# structural equation model specification for 2 groups
# 1: same-sex monozygotic (MZ)
# 2: same-sex dizygotic (DZ)
#ADE MODEL####
#full ADE model
ADE_model <-"
  # estimate intercepts (means)
    P_T1 ~ c(mean)*1 
    P_T2 ~ c(mean)*1 
  
  # ADE components
    A_T1=~ P_T1
    A_T2=~ P_T2
    D_T1 =~ P_T1
    D_T2 =~ P_T2
    E_T1 =~ P_T1
    E_T2 =~ P_T2
    
  # ADE variances
    A_T1 ~~ varA*A_T1
    A_T2 ~~ varA*A_T2
    D_T1 ~~ varD*D_T1
    D_T2 ~~ varD*D_T2
    E_T1 ~~ varE*E_T1
    E_T2 ~~ varE*E_T2
  
  # fix remaining P variance to 0
    P_T1~~0*P_T1
    P_T2~~0*P_T2
  
  # CTD twin informed constraints
    covMZA == varA 
    covDZA == .5*varA 
    covMZD == varD
    covDZD == .25*varD
  
  # covariances
  # covariances cross-twin wihtin-trait
    A_T1 ~~ c((covMZA),(covDZA))*A_T2
    D_T1 ~~ c((covMZD),(covDZD))*D_T2
    E_T1 ~~ 0*E_T2
  
  # set every other covariance to zero
    A_T1 ~~ 0*E_T1 + 0*E_T2 + 0*D_T1 + 0*D_T2
    A_T2 ~~ 0*E_T1 + 0*E_T2 + 0*D_T1 + 0*D_T2
    D_T1 ~~ 0*E_T1 + 0*E_T2
    D_T2 ~~ 0*E_T1 + 0*E_T2
  
  # calculate summary output
    varP := (varA+varD+varE)
    sdA := sqrt(varA)
    sdD := sqrt(varD)
    sdE := sqrt(varE)
  
  # narrow-sense twin-h2
    h2_P := varA/varP
    d2 := varD/varP
    e2 := varE/varP
"
#AE MODEL####
#constrained AE model
AE_model <-"
  # estimate intercepts (means)
    P_T1 ~ c(mean)*1 
    P_T2 ~ c(mean)*1 
  
  # ADE components
    A_T1 =~ P_T1
    A_T2 =~ P_T2
    E_T1 =~ P_T1
    E_T2 =~ P_T2
  
  # ADE variances
    A_T1 ~~ varA*A_T1
    A_T2 ~~ varA*A_T2
    E_T1 ~~ varE*E_T1
    E_T2 ~~ varE*E_T2
  
  # fix remaining P variance to 0
    P_T1~~0*P_T1
    P_T2~~0*P_T2
  
  # CTD twin informed constraints
    covMZA == varA 
    covDZA == .5*varA 
  
  # covariances
  # covariances cross-twin wihtin-trait
    A_T1 ~~ c((covMZA),(covDZA))*A_T2
    E_T1 ~~ 0*E_T2
  
  #set every other covariance to zero
    A_T1 ~~ 0*E_T1 + 0*E_T2 
    A_T2 ~~ 0*E_T1 + 0*E_T2
  
  #calculate summary output
    varP := (varA+varE)
    sdA := sqrt(varA)
    sdE := sqrt(varE)
  
  #narrow-sense twin-h2
    h2_P := varA/varP
    e2 := varE/varP
"
#E MODEL####
#constrained E model
E_model <-"
  # estimate intercepts (means) 
    P_T1 ~ c(mean)*1 
    P_T2 ~ c(mean)*1 
  
  # ADE components
    E_T1 =~ P_T1
    E_T2 =~ P_T2
  
  # ADE variances
    E_T1 ~~ varE*E_T1
    E_T2 ~~ varE*E_T2
  
  # fix remaining P variance to 0
    P_T1~~0*P_T1
    P_T2~~0*P_T2
  
  # covariances
  # covariances cross-twin wihtin-trait
    E_T1 ~~ 0*E_T2
  
  # calculate summary output
    varP := varE
    sdE := sqrt(varE)
    E := varE/varP
"