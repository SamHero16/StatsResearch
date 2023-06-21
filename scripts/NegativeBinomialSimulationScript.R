##Negative Binomial Regression with variable specdens. 


library(spectralGP)
library(mgcv)
library(broom)
library(ggplot2)
library(dplyr)
library(readr)
library(MASS)



##Settings

#Number of samples taken
L = 10

#Sample Size
N = 1000 

#Number of TPRS 
M = 15

#Desired Outcomes 
b0 = 0
b1 = 1

#Dimension of Population
D = 512 
WD = D*D


#Spectral Density for Gaussian Procces and Confounder ###VARIABLE
g = .25
c = .25

###########



##Gaussian Process, c(.1,2) set arbitrarily
gp=gp(c(D,D),matern.specdens,c(g,2))
simulate(gp)

##Confounding surface// change first entry of vector for adjustment
confounder = gp(c(D,D),matern.specdens,c(c,2))
simulate(confounder)


##Setup TPRS
o1 = gp$omega$omega1
o2 = gp$omega$omega2
o1o2 = data.frame(o1,o2)
smoother = s(o1,o2,k = M+1)
Smoothframe = smooth.construct(smoother, data = o1o2, knots = NULL)
TPRS = data.frame(Smoothframe$X)
TPRS <- TPRS[,c(M:(M+1),1:(M-2))]



##Create population
popx = gp$process + .5*confounder$process
popy = popx + .5*confounder$process


##Matrix's for storage
estimates <- matrix(nrow = L, ncol = M - 2)
coverage <- matrix(nrow = L, ncol = M - 2)
standardError <- matrix(nrow = L, ncol = M - 2)
b0estimates <- matrix(nrow = L, ncol = M - 2)




##Simulation Loop: Sample data, try 3:M TPRS's, record MSE.
for(i in 1:L){
  sample = sample(1:WD,N)
  x = popx[sample]
  y = rpois(WD,lambda = exp(-2 + popy))[sample]
  
  
  #Temporary storage variables
  tempEstimates = numeric(M-2)
  tempCoverage = numeric(M-2)
  tempStandardError = numeric(M-2)
  tempb0 = numeric(M-2)
  
  
  
  
  #Loop through different number of TPRS
  for(j in 3:M){
    Splines = data.matrix(TPRS[sample,1:j])
    
    model = glm.nb(y~x + Splines)
    
    tempEstimates[j-2] = (tidy(model)$estimate[2])
    tempb0[j-2] = (tidy(model)$estimate[1])
    
    tempCI = confint(model)
    tempCoverage[j-2] = tempCI["x","2.5 %"] < b1 && tempCI["x","97.5 %"] > b1
    
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
colnames(b0estimates) = c(3:M)



write_csv(data.frame(coverage),file=paste0("RawData/",format(Sys.time(), "%d%B%Y-%H:%M-G:"),g,"-C:",c, "-CoverageMatrixNegBinom.csv"))
write_csv(data.frame(estimates),file=paste0("RawData/",format(Sys.time(), "%d%B%Y-%H:%M-G:"),g,"-C:",c, "-EstimatesMatrixNegBinom.csv"))
write_csv(data.frame(standardError),file=paste0("RawData/",format(Sys.time(), "%d%B%Y-%H:%M-G:"),g,"-C:",c, "-StandardErrorMatrixNegBinom.csv"))
write_csv(data.frame(b0estimates),file=paste0("RawData/",format(Sys.time(), "%d%B%Y-%H:%M-G:"),g,"-C:",c, "-b0estimatesMatrixNegBinom.csv"))






