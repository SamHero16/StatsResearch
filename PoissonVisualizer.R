#Poisson Visualizer


L = 10 ## Number of samples taken

N = 1000 ## Sample Size

M = 50 ## Number of TPRS 

## Desired Outcomes 
b0 = 0
b1 = 1

D = 512 ## Dimension of Population


##Create population and desired coeffecients.
WD = D*D
desiredx = gp$process
desiredy = ceiling(desiredx) + ceiling(rnorm(N,0,1))
min = -1*min(desiredy)
desiredy = desiredy + min
desiredmodel = glm(desiredy~desiredx,family = "poisson")
b0 = tidy(desiredmodel)$estimate[1]
b1 = tidy(desiredmodel)$estimate[2]
newy = predict(desiredmodel,list(wt = desiredx),type = "response")


popx = gp$process + ceiling(.5*confounder$process) 
popy = ceiling(popx) + ceiling(.5*confounder$process) + ceiling(rnorm(N,0,1))
min2 = -1*min(popy)
popy = popy + min2
badmodel = glm(popy~popx,family = "poisson")
badb0 = tidy(badmodel)$estimate[1]
badb1 = tidy(badmodel)$estimate[2]
newbady = predict(badmodel,list(wt = popx),type = "response")

plot(popx,popy)
lines(popx,newbady,col = "blue")


plot(desiredx,desiredy)
lines(desiredx,newy,col = "red")
lines(popx,newbady,col = "blue")
