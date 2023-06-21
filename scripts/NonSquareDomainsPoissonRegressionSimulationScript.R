##Poisson Regression with non-square domains. 


library(spectralGP)
library(mgcv)
library(broom)
library(ggplot2)
library(dplyr)
library(readr)


##ADJUSTMENT TO SHAPE OF THE SURFACE... Check out the ggplot after I make my GP

##Make the data a circle
circlify = function(q1){
  points = expand.gpgrid(q) - .5
  filledCircle = matrix(ncol = 2,nrow = nrow(points))
  
  for(a in 1:nrow(points)){
    if(.5^2 >= points[a,1]^2 + points[a,2]^2){ #Equation of circle only gives corrdinates in a .5 radius circle
      filledCircle[a,1] = points[a,1]
      filledCircle[a,2] = points[a,2]
    }
  }
  
  filledCircle = data.frame(na.omit(filledCircle) +.5)
  circleValues = predict(q,newdata = filledCircle)
  cbind(filledCircle,circleValues)
}

#Make data a rectangle
randomRect = function(q2) {
  
  points = expand.gpgrid(q2)
  filledRect = matrix(ncol = 2,nrow = nrow(points))
  xL = runif(1,0,.5)
  xR = runif(1,.5,1)
  yL = runif(1,0,.5)
  yR = runif(1,.5,1)
  
  for(a in 1:nrow(points)){
  if(xR>points[a,1] && points[a,1]>xL && yR>points[a,2] && points[a,2]>yL){
    filledRect[a,1] = points[a,1]
    filledRect[a,2] = points[a,2]
    
  }
  
  }

  filledRect = data.frame(na.omit(filledRect))
  RectValues = predict(q2,newdata = filledRect)
  cbind(filledRect,RectValues)
}

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
D = 1024
WD = D*D


#Spectral Density for Gaussian Procces and Confounder
g = .1
c = .25

###########



##Gaussian Process, c(.1,2) set arbitrarily
gp=gp(c(D,D),matern.specdens,c(g,2))
simulate(gp)


    #Make the surface not a square 
modGP = circlify(gp)
ggplot() + geom_point(aes(x=modGP[,1],y=modGP[,2],col = modGP[,3])) #just to show it works




##Confounding surface// change first entry of vector for adjustment
confounder = gp(c(D,D),matern.specdens,c(c,2))
simulate(confounder)

    #Make confounder not a square
modConfounder = circlify(confounder)

  

##Setup TPRS (#Still a square and original size??)
o1 = gp$omega$omega1
o2 = gp$omega$omega2
o1o2 = data.frame(o1,o2)
smoother = s(o1,o2,k = M+1)
Smoothframe = smooth.construct(smoother, data = o1o2, knots = NULL)
TPRS = data.frame(Smoothframe$X)
TPRS <- TPRS[,c(M:(M+1),1:(M-2))]



##Create population
popx = modGP[,3] + modConfounder[,3]
popy = popx + modConfounder[,3]



##Matrix's for storage
estimates <- matrix(nrow = L, ncol = M - 2)
coverage <- matrix(nrow = L, ncol = M - 2)
standardError <- matrix(nrow = L, ncol = M - 2)
b0estimates <- matrix(nrow = L, ncol = M - 2)




##Simulation Loop: Sample data, try 3:M TPRS's, record MSE.
for(i in 1:L){
  sample = sample(1:length(popx),N)
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
    
    model = glm(y~x + Splines,family = "poisson")
    
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








