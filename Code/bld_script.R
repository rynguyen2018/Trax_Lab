library(deSolve)

coelicolor_ODE<- function(time, state, theta){
  ### Parameters
  
  ###AdpA
  beta_AdpA <- theta["beta_AdpA"]
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

  
  ###States
  AdpA <- state["AdpA"]
  BldA <- state["BldA"]
  
  ### ODE's 
  dAdpA <- (beta_AdpA*AdpA^n1)/(k1_AdpA^n1 + AdpA^n1) + (gamma_AdpA*AdpA)^n2/(k2_AdpA^n2 + AdpA^n2)- sigma_AdpA* AdpA 
  dBldA <- (gamma_BldA*BldA^p)/(k1_BldA^p + BldA^p) - sigma_BldA* BldA
  
  return(list(c(dAdpA, dBldA)))
  
  
}

theta <- c(beta_AdpA= 0.001, k1_AdpA= 0.005, k2_AdpA= 0.003, gamma_adpA= 0.1, sigma_AdpA=1/120, n1=1, n2=1, 
            gamma_bldA= 0.2, k1_BldA=0.004, sigma_bldA= 1/140  )
initState <- c(AdpA= 0.005,BldA= 0.003 )
trajModel <- data.frame(ode(y=initState, times=time, func=coelicolor_ODE, 
                            parms=theta, method = "ode45"))