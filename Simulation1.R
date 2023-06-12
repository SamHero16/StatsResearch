    ### Simulation 1: SLR
    
    library(spectralGP)
    library(mgcv)
    library(broom)
    library(ggplot2)
    library(dplyr)
    
    ##Settings##
    
    
    
    
    
    L = 10 ## Number of samples taken
    
    N = 1000 ## Sample Size
    
    M = 50 ## Number of TPRS 
    
    ## Desired Outcomes 
    b0 = 0
    b1 = 1
    
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
    
    ##Create population
    WD = D*D
    popx = gp$process + .5*confounder$process
    popy = popx + .5*confounder$process 
    
    #############
    
    #Matrix's for storage
    mat <- matrix(nrow = L, ncol = M - 2)
    CI <- matrix(0,nrow = 2, ncol = M-2)
    colnames(CI) = c(3:M)
    
    
    #Simulation Loop: Sample data, try 3:10 TPRS's, record MSE.
    
    for(i in 1:L){
      ##Create Sample set: important because data is related by location
      sample = sample(1:WD,N)
      ##Create sample
      x = popx[sample]
      
      ##Create sample
      
      y = popy[sample] + rnorm(N,0,1)
      
      
      ##Create results vector for this sample
      observation = numeric(M-2)
      
      
      
      z = 1
      ##Loop through different number of TPRS
      for(j in 3:M){
        workingSplines = data.matrix(TPRS[sample,1:j]) #Switched back to matrix because data frame wasnt working
        #Fit Linear Regression Model and collect results
        model = lm(y~x + workingSplines)
        observation[j-2] = (tidy(model)$estimate[2])
        tempCI = tidy(confint(model, level=0.95))
        
        ## First row is lower bound and second row is upper bound
        CI[z] = CI[z] + tempCI$x[,"2.5 %"][2]
        CI[z+1]= CI[z+1] + tempCI$x[,"97.5 %"][2]
        z = z + 2
        
      }
      
    ##Add data to matrix 
    mat[i,] = observation
    
     
     
     
    }
    #Rename columns for clarity

    colnames(mat) = c(3:M)
    
    #Find average of the CI's
    CI = CI/L
    
    
    
    ##MSE CALC
    MSEbynumspline = apply(mat,2,function(w) sqrt(mean((w-b1)^2)))
   
    
    #Resulting Information and loose summary of findings
    
      #Finding CI on MSE, dont know if this is correct
    
    rmse_interval <- function(rmse, deg_free, p_lower = 0.025, p_upper = 0.975){
      tibble(.pred_lower = sqrt(deg_free / qchisq(p_upper, df = deg_free)) * rmse,
             .pred_upper = sqrt(deg_free / qchisq(p_lower, df = deg_free)) * rmse)
    }
    
    rmseint = rmse_interval(MSEbynumspline,M)
    
    
    
    ggplot() + 
      geom_point(aes(x = 3:M,y = MSEbynumspline),size = .75) +
      geom_point(aes(x = 3:M, y = rmseint$.pred_lower,colour = "red"),size = .1) +
      geom_point(aes(x = 3:M, rmseint$.pred_upper,colour = "red"),size = .1) + 
      labs( x = "TPRS degrees of freedom", y = "Mean Squared Error")
    
    ##ggsave("/Users/samherold/Desktop/StatsResearch/MSEplot.csv")
    
    
    ggplot() + geom_point(aes(x = 3:M,y = colMeans(mat)),size = .5) + geom_path(aes(3:M,colMeans(mat)))+
      geom_point(aes(x = 3:M, y = CI[1,],colour = "red"),size = .5) + geom_path(aes(3:M,CI[1,],colour = "red"))+
      geom_point(aes(x = 3:M, y = CI[2,],colour = "red"),size = .5 ) + geom_path(aes(3:M,CI[2,],colour = "red"))+
      labs( x = "TPRS degrees of freedom", y = "Point estimate") + 
      geom_hline(yintercept = 1)
    
    ##ggsave("/Users/samherold/Desktop/StatsResearch/PointEstimatePlot.csv")
  
    
    ##Saving results 
    
      #####How do I add date to this? paste0 was acting wierd..
      ##why is my filepath wrong
    ##write.csv(mat,file="/Users/samherold/Desktop/StatsResearch/Estimates.csv")
    
    
    
  
    