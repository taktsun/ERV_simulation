library(betapart)
library(plyr)
library(tidyverse)
library(ppcor)
library(vegan)

#======================
# to do
#======================

# 1. Make number of ER strategies as a hyperparameter
# 2. Make the starting value (in "resample" function) as a hyperparameter

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
  
  
  # standardize
  dfnorm_chord <- decostand(dfSim,"norm")
  dfnorm_logchord <- decostand(log1p(dfSim),"norm")
  dfnorm_chi <- decostand(dfSim,"chi.square")
  dfnorm_hel <- decostand(dfSim,"hel")
  
  # Momentary Euclidean distance, manhattan, chi, hellinger, jaccard
  for (i in 1:n){
    
    x <- dfSim[i,]
    y <- colMeans(dfSim)
    #dfNew$edist[i] <- (dist(rbind(x,y)))
    dfNew$edist[i] <- dis_from_moment(dist(dfSim),i)
    #dfNew$manhattan[i] <- sum(abs(x -y))
    dfNew$manhattan[i] <- dis_from_moment(vegdist(dfSim,method = "manhattan"),i)
    xy <-x + y 
    y. <-y / sum(y) 
    x. <-x / sum(x) 
    
    dfNew$chordnormed[i] <- dis_from_moment(vegdist(dfnorm_chord, method = "euc"),i)
    #dfNew$chordnormed[i] <- vegdist(rbind(dfnorm_chord[i,],colMeans(dfnorm_chord)),"euc")
    dfNew$logchordnormed[i] <- dis_from_moment(vegdist(dfnorm_logchord, method = "euc"),i)
    #dfNew$logchordnormed[i] <- vegdist(rbind(dfnorm_logchord[i,],colMeans(dfnorm_logchord)),"euc")
    
    #!!! migrating to from_moment calculations...
    dfNew$chinormed[i] <- dis_from_moment(vegdist(dfnorm_chi, method = "euc"), i)
    #dfNew$chinormed[i] <- vegdist(rbind(dfnorm_chi[i,],colMeans(dfnorm_chi)),"euc")
    dfNew$hellinger[i] <- sqrt(1/2 * sum(sqrt(x) -sqrt(y))^2)
    
    #!!! the below is actually not correct. should get a matrix and average instead
    dfNew$hellingernormed[i] <- dis_from_moment(vegdist(dfnorm_hel, method = "euc"),i) 
    #dfNew$hellingernormed[i] <- vegdist(rbind(dfnorm_hel[i,],colMeans(dfnorm_hel)),"euc") 
    
    dfNew$chord[i] <- dis_from_moment(vegdist(dfSim,method="chord"),i)
    #dfNew$chord[i] <- vegdist(rbind(x,y),method="chord")[1]
    dfNew$jaccard[i] <- dis_from_moment(vegdist(dfSim,method="jaccard"),i)
    #dfNew$jaccard[i] <- vegdist(rbind(x,y),method="jaccard")[1]
    dfNew$chisq[i] <- dis_from_moment(vegdist(dfSim,method="chisq"),i)
    #dfNew$chisq[i] <- vegdist(rbind(x,y),method="chisq")[1]
    dfNew$kulczynski[i] <- dis_from_moment(vegdist(dfSim,method="kulczynski"),i)
    #dfNew$kulczynski[i] <- vegdist(rbind(x,y),method="kulczynski")[1]
    dfNew$brayveg[i] <- dis_from_moment(vegdist(dfSim,method="bray"),i)
    #dfNew$brayveg[i] <- vegdist(rbind(x,y),method="bray")[1]
    
  }
  
  # Momentary Bray-Curtis dissimilarity by beta.part
  if (ERn > 1){
      # there appears to be some limitation on the bray.part on handling 1 ER stgy only
      # will throw 0 and warning message if there is only 1 ER stgy
      for (i in 1:n){
      tempres <- bray.part(dfSim)
      dfNew$bray[i] <- dis_from_moment(tempres$bray,i)
      dfNew$bray.bal[i] <- dis_from_moment(tempres$bray.bal,i)
      dfNew$bray.gra[i] <- dis_from_moment(tempres$bray.gra,i)
      #old codes that compares moment with mean values
      #tempres <- bray.part(rbind(dfSim[i,],base::lapply(dfSim[,],FUN = mean)))
      #dfNew$bray[i] <- tempres$bray[1]
      #dfNew$bray.bal[i] <- tempres$bray.bal[1]
      #dfNew$bray.gra[i] <- tempres$bray.gra[1]
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
  simOutput$mean_manhattan <- mean(dfNew$manhattan)
  simOutput$mean_chinormed <- mean(dfNew$chinormed)
  simOutput$mean_hellinger <- mean(dfNew$hellinger)
  simOutput$mean_hellingernormed <- mean(dfNew$hellingernormed)
  simOutput$mean_jaccard <- mean(dfNew$jaccard)
  simOutput$mean_chisq <- mean(dfNew$chisq)
  simOutput$mean_logchordnormed <- mean(dfNew$logchordnormed)
  simOutput$mean_chord <- mean(dfNew$chord)
  simOutput$mean_chordnormed <- mean(dfNew$chordnormed)
  simOutput$mean_kulczynski <- mean(dfNew$kulczynski)
  simOutput$mean_brayveg <- mean(dfNew$brayveg)
  simOutput$mean_bray <- mean(dfNew$bray)
  simOutput$mean_bray.bal <- mean(dfNew$bray.bal)
  simOutput$mean_bray.gra <- mean(dfNew$bray.gra)
  simOutput$multibray.bal <- resBray$beta.BRAY.BAL 
  simOutput$multibray.gra <- resBray$beta.BRAY.GRA 
  simOutput$multibray <- resBray$beta.BRAY
  
  return(simOutput)

}

#function to calculate the mean dissimilarity from one obs to all other obs, input dist
dis_from_moment <- function(distobj,obs){
  nobs <- n
  return (mean(as.matrix(distobj)[obs,])*nobs/(nobs-1))
}

# function to return partial correlations ofdissimlarity measures to the 3 parameters 
returnfit <- function(measure,dfreturn){
  dftemp <- data.frame(dfreturn[,measure],dfreturn$autocorrA, dfreturn$ratioA,dfreturn$independence)
  if (is.na(dfreturn[1,measure])){
    return (data.frame(NA,NA,NA,NA))
    }else{
    return (pcor(dftemp)$estimate[1,])
  }
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
returninput <- c("multibray","mean_sd","mean_edist", 
                 "mean_manhattan","mean_hellinger",
                 "mean_hellingernormed",
                 "mean_jaccard",
                 "mean_chinormed", "mean_chisq",
                 "mean_logchordnormed","mean_chord","mean_chordnormed",
                 "mean_kulczynski",
                 "mean_bray","mean_brayveg")
resSummary <- data.frame()
for (i in 1:length(returninput)){
  if(length(resSummary)==0){
    resSummary <- data.frame(as.list(returnfit(returninput[i],resOutput)))
  }else{
    tempdf <- data.frame(as.list(returnfit(returninput[i],resOutput)))
    names(tempdf) <- names(resSummary)
    resSummary <- rbind.data.frame(resSummary,tempdf)
  }
  resSummary[i,1] <- returninput[i]
}

resSummary
