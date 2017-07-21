

coelicolor_ODE<- function(time, state, theta){
  ### Parameters
  
  ###AdpA
  beta_AdpA <- theta["beta_AdpA"]
  k_cat_AdpA <- theta["k_cat_AdpA"]
  gamma_AdpA <- theta["gamma_AdpA"]
  k1_AdpA <- theta["k1_AdpA"]
  k2_AdpA <- theta["k2_AdpA"]
  sigma_AdpA <- theta["sigma_AdpA"]
  n1 <- theta["n1"]
  n2 <- theta["n2"]
  
  
  ###BldA
  gamma_BldA <- theta["gamma_BldA"]
  k1_BldA <- theta["k1_BldA"]
  sigma_BldA <- theta["sigma_BldA"]
  p <- theta["p"]
  
  ###BldN 
  gamma_BldN <- theta["BldN"]
  k1_BldN <- theta["BldN"]
  sigma_BldN <- theta["sigma_BldN"]
  h <- theta["h"]
  
  ###States
  AdpA <- state["AdpA"]
  BldA <- state["BldA"]
  BldN <- state["BldN"]
  
  ### ODE's 
  dAdpA <- (beta_AdpA*AdpA^n1)/(k1_AdpA^n1 + AdpA^n1) + (gamma_AdpA*AdpA)^n2/(k2_AdpA^n2 + AdpA^n2)- sigma_AdpA* AdpA 
  dBldA <- (gamma_BldA*BldA^p)/(k1_BldA^p + BldA^p) - sigma_BldA* BldA
  dBldN <- (gamma_BldN*BldN^h)/(k1_BldN^h + BldN^h) - sigma_BldN* BldN
  
  return(list(c(dAdpA, dBldA)))
  
  
  
  
}