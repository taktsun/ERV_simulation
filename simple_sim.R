library(betapart)
library(plyr)
library(tidyverse)
library(ppcor)

#======================
# to do
#======================

# 1. Make number of ER strategies as a hyperparameter

#======================
# simulation setup
#======================

n <- 10 #number of observation (time points per "participant")
ERn <- 3 #number of ER strategies. (Single) Bray-curtis cannot handle ERn = 1 for now
simn <- 2 #number of simulations ("participant")
scalemin <- 0
scalemax <- 100
rounding <- TRUE # whether or not round scores to integers
set.seed(1999)

# curve to sample the starting value of ER strategy A[1]
resample <- function(n){
  # different curves can be set here other than the normal curve
  return (rnorm(n, mean = 50, sd = 15))
} 

# varying simulation parameters 
hyperparameters<- data.frame(
  # higher auto correlation, lower ERV expected
  autocorrA = c(0.2,0.5,0.8), 
  #higher ratioA (% directly explained by A linearly), lower ERV expected
  ratioA = c(0.2,0.5,0.8),
  #higher independence, higher ERV expected
  independence = c(0.2,0.5,0.8)
)


# ==================================
# define functions
# ==================================
funsim <- function(autocorrA,ratioA,independence){
    simOutput <- data.frame(n = n, autocorrA = autocorrA, ratioA = ratioA,
                          independence=independence)
  
  #---- simulating ER strategies
  
  dfSim <- data.frame(a = 1:n)
  dfSim$a[1] <- resample(1)
  for (i in 2:n){
    dfSim$a[i] <- dfSim$a[i-1]*autocorrA + (1-autocorrA)*resample(1)
  }
  

  #Create other strategies
  if(ERn>1){
    for (i in 2:ERn){
      dfSim[letters[i]] <- dfSim$a*ratioA + independence*(1-ratioA)*resample(n) 
    }
  }
  
  #limit the timeseries to min/max
  dfSim[] <- apply(dfSim, 2, function(x) ifelse(x < scalemin , scalemin, x))
  dfSim[] <- apply(dfSim, 2, function(x) ifelse(x > scalemax , scalemax, x))
  
  # rounding
  if(rounding){
    dfSim[] <- round(dfSim,0)
  }
                                 
  #---- ER Variability candidate indices
  
  dfNew <- dfSim
  
  # Momentary SD
  dfNew$sd <- base::apply(dfSim,1,FUN = sd)
  
  # Momentary Euclidean distance
  
  dfNew$edist<-0
  for (i in 1:n){
    dfNew$edist[i] <- (dist(rbind(dfSim[i,],base::lapply(dfSim[,],FUN = mean))))
  }
  
  # Momentary Bray-Curtis dissimilarity
  if (ERn > 1){
      # there appears to be some limitation on the bray.part on handling 1 ER stgy only
      # will throw 0 and warning message if there is only 1 ER stgy
      for (i in 1:n){
      tempres <- bray.part(rbind(dfSim[i,],base::lapply(dfSim[,],FUN = mean)))
      dfNew$bray[i] <- tempres$bray[1]
      dfNew$bray.bal[i] <- tempres$bray.bal[1]
      dfNew$bray.gra[i] <- tempres$bray.gra[1]
    }
  }else{
    dfNew$bray<-0
    dfNew$bray.bal<-0
    dfNew$bray.gra<-0
  }

  # Multi-site Bray-Curtis dissimilarity
  resBray<- beta.multi.abund(dfSim)
  
  
  #---- simOutput
  
  simOutput$em_autocorrA<-acf(dfSim$a, plot = FALSE)$acf[2] #empirical auto correlation
  simOutput$em_cor<- ifelse(ERn > 1,mean(cor(dfSim)[1,2:length(cor(dfSim))^0.5]),NA) #empirical correlation
  simOutput$mean_sd <- mean(dfNew$sd)
  simOutput$mean_edist <- mean(dfNew$edist)
  simOutput$mean_bray <- mean(dfNew$bray)
  simOutput$mean_bray.bal <- mean(dfNew$bray.bal)
  simOutput$mean_bray.gra <- mean(dfNew$bray.gra)
  simOutput$multibray.bal <- resBray$beta.BRAY.BAL 
  simOutput$multibray.gra <- resBray$beta.BRAY.GRA 
  simOutput$multibray <- resBray$beta.BRAY
  
  return(simOutput)

}

# function to return partial correlations ofdissimlarity measures to the 3 parameters 
returnfit <- function(measure,dfreturn){
  dftemp <- data.frame(dfreturn[,measure],dfreturn$autocorrA, dfreturn$ratioA,dfreturn$independence)
  return (pcor(dftemp)$estimate[1,])
}


# ====================================
# actually running the simulation
# ====================================

siminput <- expand.grid(hyperparameters, stringsAsFactors = FALSE)

resOutput <- data.frame()
for (i in 1:nrow(siminput)){
  resSim<-rdply(simn,funsim(siminput$autocorrA[i],siminput$ratioA[i],siminput$independence[i]))
  #resSim <- simOutput  
  if(length(resOutput)==0){
      resOutput <- resSim
    }else{
      resOutput <- rbind.data.frame(resOutput,resSim)
    }
}


# ============
# calculate the fit
# ============

#write.csv(resOutput,"sim.csv")
#resExample<- read.csv("sim.csv")
#returnfit("multibray",resExample)

returnfit("multibray",resOutput)
returnfit("mean_sd",resOutput)
returnfit("mean_edist",resOutput)
returnfit("mean_bray",resOutput)
