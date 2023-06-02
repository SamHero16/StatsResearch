### Simulation 1

library(spectralGP)
library(mgcv)
library(broom)
library(ggplot2)
library(dplyr)

##Settings##

L = 5 ## Number of samples taken

N = 16384 ## Sample Size 

M = 10 ## Number of TPRS 

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


##Setup TPRS: Create the full set of 10, then increment down in loops

#Make TPRS the size of the sample
ocreate = gp(c(128,128),matern.specdens,c(1,1))
o1 = ocreate$omega$omega1
o2 = ocreate$omega$omega2
smoother = s(o1,o2,k = M)
o1o2 = data.frame(o1,o2)
Smoothframe = smooth.construct(smoother, data = o1o2, knots = NULL)
TPRS = data.frame(Smoothframe$X)

##rearrange columns for easy looping
##NEED TO ABSTRACT
TPRS <- TPRS %>%
  select(X8, X9,X10, everything())

##Create population
WD = D*D
popx = gp$process + confounder$process
popy = popx + confounder$process + rnorm(WD,0,1)

#############




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
    workingSplines = data.matrix(TPRS[,1:i]) #Switched back to matrix because data frame wasnt working
    
    #Fit Linear Regression Model
    model = tidy(lm(y~x + workingSplines)) # not working: does not help regression AT ALL
    append(vecResults, model$estimate[2]). # not working
    
  }
  
  
 
  
  
  
}



