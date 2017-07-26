library(deSolve)
library(ggplot2)
setwd("C:/Users/Ryan/Desktop/gitcode/Trax_Lab/Code")
#setwd("C:/Users/Death Star/Desktop/Trax_Lab/Code")
points<- read.csv("adpA-experimental.csv", header=TRUE)


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
  dAdpA <- (beta_AdpA*AdpA^n1)/(k1_AdpA^n1 + AdpA^n1) + (gamma_AdpA*AdpA^n2)/(k2_AdpA^n2 + AdpA^n2)- sigma_AdpA* AdpA 
  dBldA <- (gamma_BldA*BldA^p)/(k1_BldA^1 + BldA^p) - sigma_BldA* BldA
  
  return(list(c(dAdpA, dBldA)))
  
  
}

#time<- points$time 
#ODEtime<- points$time
ODEtime<- c(1:15)
theta <- c(beta_AdpA= 0.1, k1_AdpA= 0.005, k2_AdpA= 0.003, gamma_AdpA= 0.1, sigma_AdpA=41.7, n1=1, n2=1, 
            gamma_BldA= 0.2, k1_BldA=0.004, sigma_BldA= 40,p=1  )
initState <- c(AdpA= 0.001,BldA= 0.005 )

trajModel <- data.frame(ode(y=initState, times=ODEtime, func=coelicolor_ODE, 
                            parms=theta, method = "ode45"))

trajAdpA<- data.frame(trajModel$AdpA)
ODEtime_data <- data.frame(trajModel$time)
ggplot(trajAdpA, aes(x=ODEtime_data, y=trajAdpA)) +
  geom_line(aes(y = trajAdpA, col = "AdpA"), size = 1.2) 


logPrior <- function(theta) {
  logPriorbeta_AdpA <- dunif(theta[["beta_AdpA"]], min = 1*10^-5, max = 1*10^-2, log = TRUE)
  logPriork1_AdpA <- dunif(theta[["k1_AdpA"]], min = 1*10^-5, max = 0.006, log = TRUE)
  logPriork2_AdpA <- dunif(theta[["k2_AdpA"]], min = 1*10^-5, max = 0.006, log = TRUE)
  logPriorgamma_AdpA <- dunif(theta[["gamma_AdpA"]], min = 1*10^-5, max = 1*10^-2, log = TRUE)
  logPriorsigma_AdpA <- dunif(theta[["sigma_AdpA"]], min = 1*10^-5, max = 1*10^-2, log = TRUE)
  logPriorn1 <-  dunif(theta[["n1"]], min = 1*10^-5, max = 1*10^-2, log = TRUE)
  logPriorn2 <- dunif(theta[["n2"]], min = 1*10^-5, max = 1*10^-2, log = TRUE)
  
  logPriorgamma_BldA <- dunif(theta[["gamma_BldA"]], min = 1*10^-5, max = 1*10^-2, log = TRUE)
  logPriork1_BldA <- dunif(theta[["k1_BldA"]], min = 1*10^-5, max = 0.006, log = TRUE)
  logPriorsigma_BldA <- dunif(theta[["sigma_BldA"]], min = 1*10^-5, max = 1*10^-2, log = TRUE)
  logPriorp <- dunif(theta[["p"]], min = 1*10^-5, max = 1*10^-2, log = TRUE)
  
  
  return(logPriorbeta_AdpA+ logPriork1_AdpA +logPriork2_AdpA + logPriorgamma_AdpA +logPriorgamma_AdpA 
         +logPriorsigma_AdpA+logPriorn1 +logPriorgamma_BldA+logPriork1_BldA +  logPriorsigma_BldA + logPriorp)
}