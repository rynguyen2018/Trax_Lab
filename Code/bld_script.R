library(deSolve)
library(ggplot2)
#setwd("C:/Users/Ryan/Desktop/gitcode/Trax_Lab/Code")
setwd("C:/Users/Death Star/Desktop/Trax_Lab/Code")
points<- read.csv("adpA-experimental.csv", header=TRUE)

adpAExpression <- points$Rounded.IMV.s
ODEtime<- seq(from =1, to=points$time[length(points$time)], by=0.01)

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

  dAdpA<- (beta_AdpA*AdpA^n1)/(k1_AdpA^n1 + AdpA^n1) + (gamma_AdpA*BldA^n2)/(k2_AdpA^n2 + BldA^n2)- sigma_AdpA* AdpA 
  #print(k1_AdpA^n1)
  dBldA <- (gamma_BldA*AdpA^p)/(k1_BldA^p + AdpA^p) - sigma_BldA* BldA
  return(list(c(dAdpA, dBldA)))
  
  
}

theta<- c(beta_AdpA= 0.989, gamma_AdpA= 0.500579864, k1_AdpA= 0.04955542, k2_AdpA= 0.02996290,  sigma_AdpA=1/120, n1=0.781525209, n2=0.88210864, gamma_BldA= 0.200629732, k1_BldA=0.03989260, sigma_BldA= 1/120,p= 5, shape_parameter= 1.052736575  )
        
logPrior <- function(theta) {
  logPriorbeta_AdpA <- dunif(theta[["beta_AdpA"]], min = 1*10^-7, max = 10, log = TRUE)
  
  logPriorgamma_AdpA <- dunif(theta[["gamma_AdpA"]], min = 1*10^-7, max = 10, log = TRUE)
  logPriork1_AdpA <- dunif(theta[["k1_AdpA"]], min = 1*10^-5, max = 10-1, log = TRUE)
  logPriork2_AdpA <- dunif(theta[["k2_AdpA"]], min = 1*10^-5, max = 10^-1, log = TRUE)
  
  logPriorsigma_AdpA <- dunif(theta[["sigma_AdpA"]], min = 10^-7, max = 1, log = TRUE)
  logPriorn1 <-  dunif(theta[["n1"]], min = -20, max = 20, log = TRUE)
  logPriorn2 <- dunif(theta[["n2"]], min = -20, max = 20, log = TRUE)
  
  logPriorgamma_BldA <- dunif(theta[["gamma_BldA"]], min = 1*10^-7, max = 1*10, log = TRUE)
  logPriork1_BldA <- dunif(theta[["k1_BldA"]], min = 1*10^-5, max = 10^-1, log = TRUE)
  logPriorsigma_BldA <- dunif(theta[["sigma_BldA"]],  min = 10^-7, max = 1, log = TRUE)
  logPriorp <- dunif(theta[["p"]], min = -20, max = 20, log = TRUE)
  logpriorshape_parameter<- dunif(theta[["shape_parameter"]], min = 0, max = 30, log = TRUE)
  
  return(logPriorbeta_AdpA+ logPriork1_AdpA +logPriork2_AdpA + logPriorgamma_AdpA +logPriorgamma_AdpA 
  +logPriorsigma_AdpA+logPriorn1 +logPriorgamma_BldA+logPriork1_BldA +  logPriorsigma_BldA + logPriorp+ logpriorshape_parameter)
}

###Likelihood function for a single data point 
pointLogLike <- function(i, expressionData, expressionModel, theta){
  #Fluorescence is observed through a  process
  nbLike <-dnbinom(x= expressionData[i], size=theta[["shape_parameter"]],mu= expressionModel[i], log=TRUE)
  return (nbLike)
  
}

## Likelihood function for all data points:
trajLogLike <- function(time, expressionData, theta, initState) {
  
  trajModel <- data.frame(ode(y=initState, times=time, func=coelicolor_ODE, 
                              parms=theta, method = "ode45"))
  expressionModel <- trajModel$AdpA
  logLike <- 0
  for (i in 1:length(time)) {
    if(time[i]%%1 ==0){
      logLike <- logLike + pointLogLike(time[i], expressionData , expressionModel, theta )
    }
  }
  return(logLike)
}

## Posterior function:
logPosterior <- function(time, expressionData, theta, initState) {
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
  return(logPosterior(time= ODEtime, expressionData= adpAExpression, 
                      theta= c(MHparams), 
                      initState=c(AdpA= 1000,BldA= 1000 ))) 
  
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


initState<- c(AdpA= 1000, BldA=1000)

theta<- c(beta_AdpA= 0.989, gamma_AdpA= 0.500579864, k1_AdpA= 0.04955542, k2_AdpA= 0.02996290,  sigma_AdpA=1/120, n1=0.781525209, n2=0.88210864, gamma_BldA= 0.200629732, k1_BldA=0.03989260, sigma_BldA= 1/120,p= 5, shape_parameter= 1.052736575  )


# Running the MCMC algorithm to vary the parameters R0 and D:
mcmcTrace <- mcmcMH(posterior = logPosteriorMH, # posterior distribution
                    initTheta = theta, # intial parameter guess
                    proposalSD = c(10^-2,            10^-5 ,                         10^-5 ,             10^-3 ,          10^-4 ,                   0.25 ,    0.25 ,                          10^-4 ,             10^-4 ,          10^-4,                 0.5,      0.3 ), # standard deviations of # parameters for Gaussian proposal distribution
                    numIterations = 20000) # number of iterations


trace <- matrix(mcmcTrace, ncol = length(theta), byrow = T)
trace <- mcmc(trace)
plot(trace)
summary(trace)

## Remove the first 2000 iterations to allow for burn-in:
traceBurn <- trace[-(1:2000),]
traceBurn <- mcmc(traceBurn)
plot(traceBurn)
summary(traceBurn)































