library(deSolve)
library(ggplot2)
#setwd("/Users/Echo_Base/Desktop/Trax_code/Code")
setwd("C:/Users/Death Star/Desktop/Trax_Lab/Code")
#setwd("/Users/Ryan/Desktop/gitcode/Trax_Lab/Code")
points<- read.csv("adpA-experimental_rev2.csv", header=TRUE)
dyn.load("gene_circuit.dll")
#dyn.load("gene_circuit.so")
adpAExpression <- points$Concentration
ODEtime<- points$Time#seq(from =1, to=points$Time[length(points$Time)], by=0.0005)


##
## post-10-mclapply.hack.R
##
## Nathan VanHoudnos
## nathanvan AT northwestern FULL STOP edu
## July 14, 2014
## Last Edit:  August 26, 2014
##
## A script to implement a hackish version of
## parallel:mclapply() on Windows machines.
## On Linux or Mac, the script has no effect
## beyond loading the parallel library.

require(parallel)

## Define the hack
mclapply.hack <- function(..., mc.cores=NULL) {
  ## Create a cluster
  if( is.null(mc.cores) ) {
    size.of.list <- length(list(...)[[1]])
    mc.cores <- min(size.of.list, detectCores())
  }
  ## N.B. setting outfile to blank redirects output to
  ##      the master console, as is the default with
  ##      mclapply() on Linux / Mac
  cl <- makeCluster( mc.cores, outfile="" )
  
  ## Find out the names of the loaded packages
  loaded.package.names <- c(
    ## Base packages
    sessionInfo()$basePkgs,
    ## Additional packages
    names( sessionInfo()$otherPkgs ))
  
  tryCatch( {
    
    ## Copy over all of the objects within scope to
    ## all clusters.
    this.env <- environment()
    while( identical( this.env, globalenv() ) == FALSE ) {
      clusterExport(cl,
                    ls(all.names=TRUE, env=this.env),
                    envir=this.env)
      this.env <- parent.env(environment())
    }
    clusterExport(cl,
                  ls(all.names=TRUE, env=globalenv()),
                  envir=globalenv())
    
    ## Load the libraries on all the clusters
    ## N.B. length(cl) returns the number of clusters
    parLapply( cl, 1:length(cl), function(xx){
      lapply(loaded.package.names, function(yy) {
        require(yy , character.only=TRUE)})
    })
    
    ## Run the lapply in parallel
    return( parLapply( cl, ...) )
  }, finally = {
    ## Stop the cluster
    stopCluster(cl)
  })
}

## If the OS is Windows, set mclapply to the
## the hackish version. Otherwise, leave the
## definition alone.
mclapply <- switch( Sys.info()[['sysname']],
                    Windows = {mclapply.hack},
                    Linux   = {mclapply},
                    Darwin  = {mclapply})

## end post-10-mclapply.hack.R

#I love computational things because it makes me sad when I do stupid things. Swag
logPrior <- function(theta) {
    logPriorbeta_AdpA <- dunif(theta[["beta_AdpA"]], min = 1*10^-7, max = 10^5, log = TRUE)
    
    logPriorgamma_AdpA <- dunif(theta[["gamma_AdpA"]], min = 1*10^-7, max = 10^5, log = TRUE)
    logPriork1_AdpA <- dnorm(log(abs(theta[["k1_AdpA"]])), mean= 4.605, sd = 2.7, log = TRUE)
    logPriork2_AdpA <- dnorm(log(abs(theta[["k2_AdpA"]])), mean = 4.605, sd = 2.7, log = TRUE)
    logPriorsigma_AdpA <- dunif(theta[["sigma_AdpA"]], min = 10^-7, max = 1, log = TRUE)
    logPriorn1 <-  dunif(theta[["n1"]], min = -20, max = 20, log = TRUE)
    logPriorn2 <- dunif(theta[["n2"]], min = -20, max = 20, log = TRUE)
    
    logPriorgamma_BldA <- dunif(theta[["gamma_BldA"]], min = 1*10^-7, max = 1*10^5, log = TRUE)
    logPriork1_BldA <- dnorm(log(abs(theta[["k1_BldA"]])), mean = 4.605, sd = 2.7, log = TRUE)
    logPriorsigma_BldA <- dunif(theta[["sigma_BldA"]],  min = 10^-7, max = 1, log = TRUE)
    logPriorp <- dunif(theta[["p"]], min = -20, max = 20, log = TRUE)
    logPriorsigma_AdpA_change <- dunif(theta[["sigma_adpAchange"]], min= 10^-7, max=1, log= TRUE)
    logPrior_bldA_change <- dunif(theta[["sigma_bldAchange"]], min= 10^-7, max=1, log= TRUE)
    logpriorshape_parameter<- dunif(theta[["shape_parameter"]], min = 0.1, max = 1000, log = TRUE)
    return(logPriorbeta_AdpA+ logPriorgamma_AdpA + logPriork1_AdpA +logPriork2_AdpA
    +logPriorsigma_AdpA+logPriorn1+logPriorn2 +logPriorgamma_BldA+logPriork1_BldA +  logPriorsigma_BldA + logPriorp+logPriorsigma_AdpA_change+logPrior_bldA_change +logpriorshape_parameter)
}

###Likelihood function for a single data point
pointLogLike <- function(i, expressionData, expressionModel, theta){
    #Fluorescence is observed through a negative binomial process with non integer counts. Thus, you first sample from a gamma distribution and then use the result as the mean for a poisson process
    nbLike <-dgamma(x= expressionData[i], shape=theta[["shape_parameter"]], scale= expressionModel[i]/theta[["shape_parameter"]] ,log=TRUE)
    if(is.na(nbLike)){
        return(0)
    }
    return(nbLike)
}

## Likelihood function for all data points:
trajLogLike <- function(time, expressionData, theta, initState) {
    trajModel <- data.frame(ode(y=initState, times=time, func="derivs",parms=theta[0:13], dllname= "gene_circuit", initfunc = "initmod", nout = 2,events = list(func="event",time=28) ))
    expressionModel <- trajModel$AdpA
    logLike <- 0
    for (i in 1:length(time)) {
        logLike <- logLike + pointLogLike(time[i], expressionData , expressionModel, theta )
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
    initState=c(AdpA= 0.00001, BldA=0.000015)))
    
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
        thetaProposed <- round(rnorm(n= length(thetaCurrent), mean= thetaCurrent, sd=proposalSD),5)
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
        "acceptance rate:", accepted / i, "\n")
    }
    return(samples)
}


initState<- c(AdpA= 0.00001, BldA=0.000015)

theta<- c(beta_AdpA= 90, gamma_AdpA=320, k1_AdpA= 178, k2_AdpA= 90,  sigma_AdpA=3.5*10^-3, n1=1.2, n2=5, gamma_BldA= 203, k1_BldA=200, sigma_BldA= 3.5*10^-3,p= 5, sigma_adpAchange= 2*10^-2, sigma_bldAchange= 4*10^-2, shape_parameter= 32)


#cl <- parallel::makeCluster(spec=detectCores()-10,type="PSOCK")
#parallel::clusterSetRNGStream(cl=cl,iseed=NULL)


#clusterExport(cl, list("adpAExpression", "ODEtime", "initState", "logPrior", "theta", "logPosteriorMH", "logPosterior", "pointLogLike", "trajLogLike", "mcmcMH"))

mcmcTrace <-mcmcMH(posterior = logPosteriorMH, # posterior distribution
initTheta = theta, # intial parameter guess
proposalSD = c(4*10^-1, 4*10^-1, 5*10^-1, 5*10^-1, 5*10^-4, 6*10^-3, 2*10^-2, 4*10^-1, 5*10^-1, 5*10^-4,5*10^-1, 7*10^-4,5*10^-3,2*10^-3, 2*10^-1), # standard deviations of # parameters for Gaussian proposal distribution
numIterations = 350000) # number of iterations
