##Poisson Regression 
##(via rounding the dependent variable)

library(spectralGP)
library(mgcv)
library(broom)
library(ggplot2)
library(dplyr)

##Settings##

L = 5 ## Number of samples taken

N = 1000 ## Sample Size

M = 25## Number of TPRS 

## Desired Outcomes 
b0 = 0
b1 = 1

D = 512 ## Dimension of Population




## Gaussian Process, c(.1,2) set arbitrarily
gp=gp(c(D,D),matern.specdens,c(.1,2))
simulate(gp)

##Confounding surface// change first entry of vector for adjustment
confounder = gp(c(D,D),matern.specdens,c(.25,2))
simulate(confounder)


##Setup TPRS: Create the full set, then increment later using loops
##TPRS are size of population, will take sample later
o1 = gp$omega$omega1
o2 = gp$omega$omega2
o1o2 = data.frame(o1,o2)
smoother = s(o1,o2,k = M+1)
Smoothframe = smooth.construct(smoother, data = o1o2, knots = NULL)
TPRS = data.frame(Smoothframe$X)

##rearrange columns for easy looping
TPRS <- TPRS[,c(M:(M+1),1:(M-2))]



##Create population

WD = D*D

popx = gp$process + .5*confounder$process
popy = popx + .5*confounder$process

#Matrix's for storage

estimates <- matrix(nrow = L, ncol = M - 2)
coverage <- matrix(nrow = L, ncol = M - 2)
standardError <- matrix(nrow = L, ncol = M - 2)
b0estimates <- matrix(nrow = L, ncol = M - 2)

#Betas
#coverage,yes or no
#Standard error


#Simulation Loop: Sample data, try 3:M TPRS's, record MSE.

for(i in 1:L){
  ##Create Sample set: important because data is related by location
  sample = sample(1:WD,N)
  
  ##Create sample
  x = popx[sample]
  y = rpois(WD,lambda = exp(-2 + popy))[sample]
  
  

  ##Create results vector for this sample
  tempEstimates = numeric(M-2)
  tempCoverage = numeric(M-2)
  tempStandardError = numeric(M-2)
  tempb0 = numeric(M-2)
  
  
  
  
  ##Loop through different number of TPRS
  for(j in 3:M){
    workingSplines = data.matrix(TPRS[sample,1:j])
    
    #Fit Linear Regression Model and collect results
    model = glm(y~x + workingSplines,family = "poisson")
    
    
    tempEstimates[j-2] = (tidy(model)$estimate[2])
    
    tempb0[j-2] = (tidy(model)$estimate[1])
    
    
    
    tempCI = tidy(confint(model))
    tempCoverage[j-2] = tempCI$x[,"2.5 %"][2] < b1 && tempCI$x[,"97.5 %"][2] > b1
    
    tempStandardError[j-2] = (tidy(model)$std.error[2])
  }
  
  ##Add data to matrix 
  estimates[i,] = tempEstimates
  coverage[i,] = tempCoverage
  standardError[i,] = tempStandardError
  b0estimates[i,] = tempb0
  
  
}

#Rename columns for clarity
colnames(estimates) = c(3:M)
colnames(coverage) = c(3:M)
colnames(standardError) = c(3:M)



#Calculations

#MSE CALC
RMSE = apply(estimates,2,function(w) sqrt(mean((w-b1)^2)))
#Means of estimates
columnMeans = colMeans(estimates)
#Standard Deviation of the distribution of Estimates
sdOfEstimates = apply(estimates, 2, sd)
#bias
bias = apply(estimates, 2,function(w) mean(w-b1))
#percent Coverage
percentCoveage = apply(coverage,2,mean)



ggplot() + geom_point(aes(x = 3:M,y =  columnMeans),size = .5) + geom_path(aes(3:M, columnMeans))+
  labs( x = "TPRS degrees of freedom", y = "Point estimate of b1", title = "Point Estimate's of b1",subtitle =  "by TPRS degrees of freedom") + 
  geom_hline(yintercept = b1)

ggplot() + geom_point(aes(x = 3:M,y =  colMeans(b0estimates)),size = .5) + geom_path(aes(3:M, colMeans(b0estimates)))+
  labs( x = "TPRS degrees of freedom", y = "Point estimate of b0", title = "Point Estimate's of b0",subtitle =  "by TPRS degrees of freedom") 



summary(model)
#LOOK AT INTERCEPT



