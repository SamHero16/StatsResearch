##Poisson Regression 
##(via rounding the dependent variable)

library(spectralGP)
library(mgcv)
library(broom)
library(ggplot2)
library(dplyr)

##Settings##

L = 10 ## Number of samples taken

N = 1000 ## Sample Size

M = 50## Number of TPRS 

## Desired Outcomes tbd

D = 512 ## Dimension of Population


## Gaussian Process, c(.1,2) set arbitrarily
gp=gp(c(D,D),matern.specdens,c(.1,2))
simulate(gp)

##Confounding surface// change first entry of vector for adjustment
confounder = gp(c(D,D),matern.specdens,c(.5,2))
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



##Create population and desired coeffecients.

#Create desired Poisson model and log coeffiecents. 
WD = D*D

desiredx = gp$process
desiredy = ceiling(desiredx) + ceiling(rnorm(N,0,1))
min = -1*min(desiredy)
desiredy = desiredx + min
desiredmodel = tidy(glm(desiredy~desiredx,family = "poisson"))
b0 = desiredmodel$estimate[1]
b1 = desiredmodel$estimate[2]



#Actual counfounded population creation
popx = gp$process + ceiling(.5*confounder$process) 
popy = ceiling(popx) + ceiling(.5*confounder$process) 
min2 = -1*min(popy)
popy = popy + min2






##Storage
mat <- matrix(nrow = L, ncol = M - 2)
CI <- matrix(0,nrow = 2, ncol = M-2)
colnames(CI) = c(3:M)



##Main loop

for(i in 1:L){
  ##Create Sample set: important because data is related by location
  sample = sample(1:WD,N)
  ##Create sample
  x = popx[sample]
  
  ##Create sample
  ##abs is for rare case thas popy is 0 and rnorm is -1, happens approx less than once in a sample of 1000
  y = abs(popy[sample] + ceiling(rnorm(N,0,1)))
  
  
  
  
  
  
  ##Create results vector for this sample
  observation = numeric(M-2)
  
  
  
 z = 1
  
  ##Loop through different number of TPRS
  for(j in 3:M){
    workingSplines = data.matrix(TPRS[sample,1:j]) #Switched back to matrix because data frame wasnt working
    #Fit Linear Regression Model and collect results
    model = glm(y~x + workingSplines, family = "poisson")
    observation[j-2] = (tidy(model)$estimate[2])
    tempCI = tidy(confint(model, level=0.95)) #THIS TAKES VERY LONG
    
     ##First row is lower bound and second row is upper bound
    CI[z] = CI[z] + tempCI$x[,"2.5 %"][2]
    CI[z+1]= CI[z+1] + tempCI$x[,"97.5 %"][2]
    z = z + 2
    
  }
  
  ##Add data to matrix 
  mat[i,] = observation
  
  

}

colnames(mat) = c(3:M)
CI = CI/L


ggplot() + geom_point(aes(x = 3:M,y = colMeans(mat)),size = .5) + geom_path(aes(3:M,colMeans(mat)))+
  geom_point(aes(x = 3:M, y = CI[1,],colour = "red"),size = .5) + geom_path(aes(3:M,CI[1,],colour = "red"))+
  geom_point(aes(x = 3:M, y = CI[2,],colour = "red"),size = .5 ) + geom_path(aes(3:M,CI[2,],colour = "red"))+
  labs( x = "TPRS degrees of freedom", y = "Point estimate") + 
  geom_hline(yintercept = b1)




#notes: 
#creating CI's for general linearized models takes insanely long
#The deisred outcome is b1 with no confounding right? What else would it be?









