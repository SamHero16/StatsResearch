  ### Simulation 1: SLR
  
  library(spectralGP)
  library(mgcv)
  library(broom)
  library(ggplot2)
  library(dplyr)
  
  ##Settings##
  
  #Fully adjustable
  
  L = 5 ## Number of samples taken
  
  N = 16384 ## Sample Size
  
  M = 50 ## Number of TPRS 
  
  ## Desired Outcomes 
  b0 = 0
  b1 = 1
  
  D = 512 ## Dimension of Population
  
  
  ## Gaussian Process, c(.1,2) set arbitrarily
  gp=gp(c(D,D),matern.specdens,c(.1,2))
  simulate(gp)
  
  ##Confounding surface
  confounder = gp(c(D,D),matern.specdens,c(.15,2))
  simulate(confounder)
  
  
  ##Setup TPRS: Create the full set, then increment later using loops
  ##TPRS are size of population, will take sample later
  o1 = gp$omega$omega1
  o2 = gp$omega$omega2
  o1o2 = data.frame(o1,o2)
  smoother = s(o1,o2,k = M)
  Smoothframe = smooth.construct(smoother, data = o1o2, knots = NULL)
  TPRS = data.frame(Smoothframe$X)
  
  ##rearrange columns for easy looping
  TPRS <- TPRS %>%
    select(M-2:M, everything())
  
  ##Create population
  WD = D*D
  popx = gp$process + .5*confounder$process
  popy = popx + .5*confounder$process + rnorm(WD,0,1)
  
  #############
  
  #Matrix for storage
  mat <- matrix(ncol = M -2)
  
  
  
  #Simulation Loop: Sample data, try 3:10 TPRS's, record MSE.
  
  for(i in 1:L){
    ##Create Sample set: important because data is related by location
    sample = sample(1:WD,16384)
    ##Create sample
    x = popx[sample]
    
    ##Create sample
    
    y = popy[sample]
    
    
    ##Create results vector for this Sample
    vecResults = c()
    
   
    ##Loop through different number of TPRS
    for(i in 3:M){
      workingSplines = data.matrix(TPRS[sample,1:i]) #Switched back to matrix because data frame wasnt working
      
      #Fit Linear Regression Model and collect results
      model = tidy(lm(y~x + workingSplines) )
      vecResults = append(vecResults, model$estimate[2]) 
    }
    
  ##Add data to matrix 
  mat = rbind(mat,vecResults)
   
   
   
  }
  #Take off first useless row
  mat = mat[-1,]
  
  
  ##MSE CALC
  MSEbynumspline = c()
  for(i in 1:ncol(mat)){
    mse = mean((mat[,i] - 1)^2)
    MSEbynumspline[i] = mse
  }
  
  #Resulting Information and loose summary of findings
  lm(y~x)
  MSEbynumspline
  plot(colMeans(mat))
  plot(MSEbynumspline)
  min(colMeans(mat))
  min(MSEbynumspline)
  
  
#propose to Study actual time effectiveness of creating the splines as side project. 
  