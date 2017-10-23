library(deSolve)
library(ggplot2)
library(parallel)
setwd("/Users/Echo_Base/Desktop/Trax_code/Code")
#setwd("C:/Users/Death Star/Desktop/Trax_Lab/Code")
#setwd("/Users/ryannguyen/Desktop/Trax_code/Code")
points<- read.csv("adpA-experimental_rev2.csv", header=TRUE)
#dyn.load("gene_circuit.dll")
dyn.load("gene_circuit.so")
adpAExpression <- points$Concentration_micromolar[0:35]
ODEtime<- points$Time[0:35]#seq(from =1, to=points$Time[length(points$Time)], by=0.0005)

#I love computational things because it makes me sad when I do stupid things. Swag
logPrior <- function(theta) {
  reject_value<- -99999999

  logPriorbeta_AdpA <- max(reject_value, dunif(theta[["beta_AdpA"]], min = 1*10^-7, max = 10^5, log = TRUE))
  
  logPriorgamma_AdpA <- dunif(theta[["gamma_AdpA"]], min = 1*10^-7, max = 10^5, log = TRUE)
  if(theta[["k1_AdpA"]]<0.00005){
    logPriork1_AdpA<- reject_value
  }else{
    logPriork1_AdpA <- dnorm(log(theta[["k1_AdpA"]]), mean= 4.605, sd = 2.7, log = TRUE)
  }
  if(theta[["k2_AdpA"]]<0.000005){
    logPriork2_AdpA<- reject_value
  }else{
    logPriork2_AdpA <- dnorm(log(theta[["k1_AdpA"]]), mean= 4.605, sd = 2.7, log = TRUE)
  }
  logPriorn1 <-  dunif(theta[["n1"]], min = -20, max = 20, log = TRUE)
  logPriorn2 <- dunif(theta[["n2"]], min = -20, max = 20, log = TRUE)
  logPriorgamma_BldA <- dunif(theta[["gamma_BldA"]], min = 1*10^-7, max = 1*10^5, log = TRUE)
  if(theta[["k1_BldA"]]<0.000005){
    logPriork1_BldA<- reject_value
  }else{
    logPriork1_BldA <- dnorm(log(theta[["k1_BldA"]]), mean= 4.605, sd = 2.7, log = TRUE)
  }
  logPriorp <- dunif(theta[["p"]], min = -20, max = 20, log = TRUE)
  logpriorshape_parameter<- dunif(theta[["shape_parameter"]], min = 0.1, max = 1000, log = TRUE)
  return(logPriorbeta_AdpA+ logPriorgamma_AdpA + logPriork1_AdpA +logPriork2_AdpA
         +logPriorn1+logPriorn2 +logPriorgamma_BldA+logPriork1_BldA +   logPriorp +logpriorshape_parameter)
}

###Likelihood function for a single data point
pointLogLike <- function(i, expressionData, expressionModel, theta){
  #Fluorescence is observed through a negative binomial process with non integer counts. Thus, you first sample from a gamma distribution and then use the result as the mean for a poisson process
  if(i==0){
    i<-i+1
  }
  nbLike <-dgamma(x= expressionData[i], shape=theta[["shape_parameter"]], scale= expressionModel[i]/theta[["shape_parameter"]] ,log=TRUE)
  if(is.na(nbLike)){
    return(0)
  }
  return(nbLike)
}

## Likelihood function for all data points:
trajLogLike <- function(time, expressionData, theta, initState) {
  trajModel <- data.frame(ode(y=initState, times=time, func="derivs",parms=theta[0:9], dllname= "gene_circuit", initfunc = "initmod", nout = 2))
  
  expressionModel <- trajModel$AdpA
  logLike <- 0
  for (i in time) {
      logLike <- logLike + pointLogLike(i/1800, expressionData , expressionModel, theta )
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
                      initState=c(AdpA= 1, BldA=1.5)))
  
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
    #print(thetaProposed)
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
        "acceptance rate:", accepted /i, "\n")
  }
  #print(accepted/numIterations)
  return(samples)
}


initState<- c(AdpA= 1, BldA=1.5) #initial concentration in micromolar

theta_vec <- list()

#theta1<- c(beta_AdpA= 90, gamma_AdpA=320, k1_AdpA= 178, k2_AdpA= 90,  sigma_AdpA=3.5*10^-3, n1=1.2, n2=5, gamma_BldA= 203, k1_BldA=200, sigma_BldA= 3.5*10^-3,p= 5, sigma_adpAchange= 2*10^-2, sigma_bldAchange= 4*10^-2, shape_parameter= 32)
#theta2<- c(beta_AdpA= 100, gamma_AdpA=200, k1_AdpA= 100, k2_AdpA= 90,  sigma_AdpA=5*10^-3, n1=1.2, n2=1.68, gamma_BldA= 200, k1_BldA=200, sigma_BldA= 5*10^-2,p= 10, sigma_adpAchange= 4*10^-2, sigma_bldAchange= 0.3, shape_parameter= 3)


no_cores<- detectCores()-2
for(i in 1:no_cores){
  theta<- c(beta_AdpA= runif(1, min= 40, max= 100), gamma_AdpA=runif(1, min= 40, max= 100), k1_AdpA= runif(1, min= 9*10^-2, max= 10), k2_AdpA= runif(1, min= 9*10^-2, max= 10), n1=runif(1, min=-3, max= 4), n2=runif(1, min= -4, max= 5), gamma_BldA= runif(1, min= 20, max= 100), k1_BldA=runif(1, min= 9*10^-2, max= 10),p= runif(1, min= -2, max= 7), shape_parameter= 32)
  theta_vec<- c(theta_vec, list(theta))
}

print(theta_vec)
cl<-makeCluster(no_cores)
parallel::clusterSetRNGStream(cl=cl,iseed=NULL)

mcmcTrace<- mclapply(X= theta_vec, FUN=function(theta){mcmcMH(posterior = logPosteriorMH, # posterior distribution
                                                             initTheta = theta, # intial parameter guess
                                                             proposalSD = c(4*10^-2, 4*10^-2, 5*10^-2, 3*10^-2, 6*10^-2, 3*10^-2, 4*10^-2, 5*10^-2, 5*10^-2,2*10^-1), # standard deviations of # parameters for Gaussian proposal distribution
                                                             numIterations = 500000)}) # number of iterations 


stopCluster(cl)
rm(cl)

library(coda)

saveRDS(object=mcmcTrace,file="mcmc_out_10_23_17.rds")


#traceBurn <- trace[-(1:1000),]
#traceBurn <- mcmc(traceBurn)
#plot(traceBurn)
#summary(traceBurn)
#autocorr.plot(traceBurn)
#library(coda)
#trace <- mcmc(trace)
#plot(trace)
#summary(trace)
