library(deSolve)
library(ggplot2)
#setwd("C:/Users/Ryan/Desktop/gitcode/Trax_Lab/Code")
setwd("C:/Users/Death Star/Desktop/Trax_Lab/Code")
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

time<- points$time 
adpAExpression<- points$IMV.intensity.Mean.Value.

###Prior Function
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

###Likelihood function for a single data point 
pointLogLike <- function(i, expressionData, expressionModel){
 #Fluorescence is observed through a Poisson process
  poissonLike <-dpois(x= expressionData[i], expressionModel[i], log=TRUE) 
  return (poissonLike)
  
}

## Likelihood function for all data points:
trajLogLike <- function(time, expressionData, theta, initState) {
  trajModel <- data.frame(ode(y=initState, times=ODEtime, func=coelicolor_ODE, 
                              parms=theta, method = "ode45"))
  expressionModel <- trajModel$AdpA
  logLike <- 0
  for (i in 1:length(ODEtime)) {
    logLike <- logLike + pointLogLike(i, expressionData , expressionModel )
  }
  return(logLike)
}
## Posterior function:
logPosterior <- function(time, BurdenData, theta, initState) {
  ## Calculate the log prior (logPrior) for the vector of model
  ## parameters (theta).
  logPrior <- logPrior(theta)
  
  ## Calculate the log likelihood (logLike) of the data given theta, the
  ## incidence data (IncData), and the initial values of the state
  ## variables (initState).
  logLike <- trajLogLike(time, expressionData, theta, initState)
  
  logPosterior <- logPrior + logLike
  return(logPosterior)
}

## Posterior function for Metropolis_Hastings algorithm (only depends on theta): 
logPosteriorMH <- function(MHparams) {
  return(logPosterior(ODEtime, adpAExpression, 
                      theta= c(MHparams), 
                      initState=c(AdpA= 0.001,BldA= 0.005 ))) 

}

mcmcMH <- function(posterior, initTheta, proposalSD, numIterations) {
  
  # Evaluate the function "posterior" at "initTheta", and assign to a
  # variable called posteriorThetaCurrent.
  posteriorThetaCurrent <- logPosteriorMH(initTheta)
  
  # Initialise variables to store the current value of theta, the
  # vector of sample values, and the number of accepted proposals.
  thetaCurrent <- initTheta 
  samples <- initTheta
  accepted <- 0
  
  # Run the MCMC algorithm for numIterations interations.
  for (i in 1:numIterations) {
    
    # Draw a new theta from a Gaussian proposal distribution and
    # assign this to a variable called thetaProposed.
    thetaProposed <- rnorm(n= length(thetaCurrent), mean= thetaCurrent, sd=proposalSD)
    
    # Assign names to the thetaProposed vector.
    names(thetaProposed) <- names(thetaCurrent)
    
    # Evaluate the log) posterior function at the proposed theta
    # value and assign to a variable called 
    # posteriorThetaProposed.
    posteriorThetaProposed <- logPosteriorMH(thetaProposed) 
    # Compute the Metropolis-Hastings (log) acceptance
    # probability and assign to a variable called
    # logAcceptance.
    logAcceptance <- posteriorThetaProposed-posteriorThetaCurrent
    
    # Draw a random number uniformly-distributed between 0 and 1
    # using "runif" and assign to a variable called randNum.
    randNum <- runif(1,min=0,max=1)
    
    # Use the random number and the acceptance probability to 
    # determine if thetaProposed will be accepted.
    if (randNum< exp(logAcceptance)) {
      
      # If accepted, change the current value of theta to the
      # proposed value of theta.
      thetaCurrent <- thetaProposed
      
      # And update the current value of the posterior 
      # function.
      posteriorThetaCurrent <- posteriorThetaProposed
      
      # And update number of accepted proposals.
      accepted <- accepted+1
      
    }
    
    # Add the current theta to the vector of samples.
    samples <- c(samples, thetaCurrent)
    
    cat("iteration:", i, "chain:", thetaCurrent,
        "acceptance rate:", accepted / i, "\n")
  }
  return(samples)
}


# Running the MCMC algorithm to vary the parameters R0 and D:
mcmcTrace <- mcmcMH(posterior = logPosteriorMH, # posterior distribution
                    initTheta =c(beta_AdpA= 0.1, k1_AdpA= 0.005, k2_AdpA= 0.003, gamma_AdpA= 0.1, sigma_AdpA=41.7, n1=1, n2=1, gamma_BldA= 0.2, k1_BldA=0.004, sigma_BldA= 40,p=1  ), # intial parameter guess
                    proposalSD = c(2.6*10^-5, 3.2*10^-4 ), # standard deviations of 
                    # parameters for Gaussian proposal distribution
                    numIterations = 5000) # number of iterations

trace <- matrix(mcmcTrace, ncol = 2, byrow = T)

library(coda)
trace <- mcmc(trace)
plot(trace)
summary(trace)

#Burn in 
traceBurn <- trace[-(1:1000),]
traceBurn <- mcmc(traceBurn)
plot(traceBurn)
summary(traceBurn)







#ODEtime<- points$time
#ODEtime<- c(1:15)
#theta <- c(beta_AdpA= 0.1, k1_AdpA= 0.005, k2_AdpA= 0.003, gamma_AdpA= 0.1, sigma_AdpA=41.7, n1=1, n2=1, 
#            gamma_BldA= 0.2, k1_BldA=0.004, sigma_BldA= 40,p=1  )
#initState <- c(AdpA= 0.001,BldA= 0.005 )

#trajModel <- data.frame(ode(y=initState, times=ODEtime, func=coelicolor_ODE, 
#                            parms=theta, method = "ode45"))

#trajAdpA<- data.frame(trajModel$AdpA)
#ODEtime_data <- data.frame(trajModel$time)
#ggplot(trajAdpA, aes(x=ODEtime_data, y=trajAdpA)) +
#  geom_line(aes(y = trajAdpA, col = "AdpA"), size = 1.2) 

