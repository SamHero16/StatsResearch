#Poisson Visualizer



D = 256


##Create population
WD = D*D
gp=gp(c(D,D),matern.specdens,c(.1,2))
simulate(gp)

##Confounding surface
confounder = gp(c(D,D),matern.specdens,c(.5,2))
simulate(confounder)
conf = ceiling(.5*confounder$process) 


#Normal Poisson Variable creation and poisson model fitting
normalx = gp$process
normaly = ceiling(normalx) + ceiling(rnorm(WD,0,1))
min = -1*min(normaly)
normaly = normaly + min

normalmodel = glm(normaly~normalx,family = "poisson")
b0 = tidy(normalmodel)$estimate[1]
b1 = tidy(normalmodel)$estimate[2]
poissonlinenormal = predict(normalmodel,list(wt = normalx),type = "response")




#Confounding Poisson Variable creation and poisson model fitting

cx = gp$process + conf 
cy = ceiling(cx) + conf + ceiling(rnorm(WD,0,1))
min2 = -1*min(cy)
cy = cy + min2


confoundedmodel = glm(cy~cx,family = "poisson")
cb0 = tidy(confoundedmodel)$estimate[1]
cb1 = tidy(confoundedmodel)$estimate[2]
poissonLineConfoundingOnConfounding = predict(confoundedmodel,list(wt = cx),type = "response")




#Adjusted Poisson Variable creation and poisson model fitting
range = seq(range(cx)[1],range(cx)[2],.1)
adjustedmodel = glm(cy~cx + conf,family = "poisson")
ab0 = tidy(adjustedmodel)$estimate[1]
ab1 = tidy(adjustedmodel)$estimate[2]
adjustedLineOnConfounding = predict(adjustedmodel,type = "response",newdata = list(cx = cx))




#Plots

  #Plot no confounding
plot(normalx,normaly)
lines(normalx,poissonlinenormal,col = "red")

  #Plot confounding
plot(cx,cy)
lines(cx,poissonLineConfoundingOnConfounding, col = "blue")

  #Plot adjusted
summary(adjustedmodel)
plot(cx,cy)
lines(normalx,adjustedLineOnConfounding, col = "green")##This does not make sense in my brain


  #collected plot shown on confounding variables. 
plot(cx,cy)
lines(normalx,poissonlinenormal,col = "red")
lines(cx,poissonLineConfoundingOnConfounding, col = "blue")
##lines(normalx,adjustedLineOnConfounding, col = "green")
summary(adjustedmodel)
