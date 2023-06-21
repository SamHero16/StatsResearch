##Poisson Regression with variable specdens. 


library(spectralGP)
library(mgcv)
library(broom)
library(ggplot2)
library(dplyr)
library(readr)



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
g = .1
c = .1

###########



##Gaussian Process, c(.1,2) set arbitrarily
gp1=gp(c(D,D),matern.specdens,c(g,2))
simulate(gp1)

gp2 = gp(c(D,D),matern.specdens,c(g,2))
simulate(gp2)

##Confounding surface// change first entry of vector for adjustment
confounder = gp(c(D,D),matern.specdens,c(c,2))
simulate(confounder)


##Setup TPRS
o1 = gp1$omega$omega1
o2 = gp1$omega$omega2
o1o2 = data.frame(o1,o2)
smoother = s(o1,o2,k = M+1)
Smoothframe = smooth.construct(smoother, data = o1o2, knots = NULL)
TPRS = data.frame(Smoothframe$X)
TPRS <- TPRS[,c(M:(M+1),1:(M-2))]



##Create population
popx = gp1$process + .5*confounder$process
popx2 = + gp2$process + .5*confounder$process
popy = popx + popx2 + .5*confounder$process


##Matrix's for storage
estimatesx <- matrix(nrow = L, ncol = M - 2)
coveragex <- matrix(nrow = L, ncol = M - 2)
standardErrorx <- matrix(nrow = L, ncol = M - 2)
estimatesx2 <- matrix(nrow = L, ncol = M - 2)
coveragex2 <- matrix(nrow = L, ncol = M - 2)
standardErrorx2 <- matrix(nrow = L, ncol = M - 2)






##Simulation Loop: Sample data, try 3:M TPRS's, record MSE.
for(i in 1:L){
  sample = sample(1:WD,N)
  x = popx[sample]
  x2 = popx2[sample]
  y = rpois(WD,lambda = exp(-2 + popy))[sample]
  
  
  
  
  #Temporary storage variables
  tempEstimatesx = numeric(M-2)
  tempEstimatesx2 = numeric(M-2)
  tempCoveragex = numeric(M-2)
  tempCoveragex2 = numeric(M-2)
  tempStandardErrorx = numeric(M-2)
  tempStandardErrorx2 = numeric(M-2)
  
  
  
  
  #Loop through different number of TPRS
  for(j in 3:M){
    Splines = data.matrix(TPRS[sample,1:j])
    
    model = glm(y~x2 + x + Splines,family = "poisson")
    
    tempEstimatesx[j-2] = (tidy(model)$estimate[2])
    tempEstimatesx2[j-2] = (tidy(model)$estimate[3])
    
    
    tempCI = confint(model)
    tempCoveragex[j-2] = tempCI["x","2.5 %"] < b1 && tempCI["x","97.5 %"] > b1
    tempCoveragex[j-2] = tempCI["x2","2.5 %"] < b1 && tempCI["x2","97.5 %"] > b1
    tempStandardErrorx[j-2] = (tidy(model)$std.error[2])
    tempStandardErrorx2[j-2] = (tidy(model)$std.error[3])
  }
  
  ##Add data to matrix 
  estimatesx[i,] = tempEstimatesx
  coveragex[i,] = tempCoveragex
  standardErrorx[i,] = tempStandardErrorx
  estimatesx2[i,] = tempEstimatesx2
  coveragex2[i,] = tempCoveragex2
  standardErrorx2[i,] = tempStandardErrorx2
 
  
  
}

#Rename columns for clarity
colnames(estimatesx) = c(3:M)
colnames(coveragex) = c(3:M)
colnames(standardErrorx) = c(3:M)
colnames(estimatesx2) = c(3:M)
colnames(coveragex2) = c(3:M)
colnames(standardErrorx2) = c(3:M)




write_csv(data.frame(coveragex),file=paste0("RawData/",format(Sys.time(), "%d%B%Y-%H:%M-G:"),g,"-C:",c, "-XCoverageMatrixTwoPredictorsPoisson.csv"))
write_csv(data.frame(estimatesx),file=paste0("RawData/",format(Sys.time(), "%d%B%Y-%H:%M-G:"),g,"-C:",c, "-XEstimatesMatrixTwoPredictorsPoisson.csv"))
write_csv(data.frame(standardErrorx),file=paste0("RawData/",format(Sys.time(), "%d%B%Y-%H:%M-G:"),g,"-C:",c, "-XStandardErrorMatrixTwoPredictorsPoisson.csv"))

write_csv(data.frame(coveragex2),file=paste0("RawData/",format(Sys.time(), "%d%B%Y-%H:%M-G:"),g,"-C:",c, "-X2CoverageMatrixTwoPredictorsPoisson.csv"))
write_csv(data.frame(estimatesx2),file=paste0("RawData/",format(Sys.time(), "%d%B%Y-%H:%M-G:"),g,"-C:",c, "-X2EstimatesMatrixTwoPredictorsPoisson.csv"))
write_csv(data.frame(standardErrorx2),file=paste0("RawData/",format(Sys.time(), "%d%B%Y-%H:%M-G:"),g,"-C:",c, "-X2StandardErrorMatrixTwoPredictorsPoisson.csv"))






