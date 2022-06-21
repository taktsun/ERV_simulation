library(betapart)
library(plyr)
library(tidyverse)
library(ppcor)
library(vegan)
library(entropy)

#======================
# notes/to-do
#======================

# 1. ERn = 1 case can be handled by adding a blank strategy (set ERblankn>=1)
# 2. Explore if there's a way out other than choosing between an unwanted
    #-ve association with resample_mean or +ve association with scalemax.
    #E.g., Extra adjustment, say, add ERn and scalemax to denominator of Euclidean distance?
# 3. Produce data visualizations to see if relationships are linear or not
# 4. Add RQA

#======================
# simulation setup
#======================

n <- 10 #number of observation (time points per participant)
ERblankn <- 0 #number of "unused" (always 0) ER strategies
simn <- 2 #number of simulations (participant)
scalemin <- 0
rounding <- TRUE # whether or not round scores to integers
firstrowzero <- FALSE # testing purpose: set TRUE to force first row to be all 0
testnonoverlap <- TRUE # testing purpose: set TRUE to make [1,1] = 0 and [2,2] = 0
zerotransform <- TRUE # whether or not to replace 0 with 0.0001. TRUE recommended because 0 are not true 0
set.seed(1999)


# varying simulation parameters 
siminput <- expand_grid(
  #higher independence, higher ERV expected
  independence = c(0.2,0.8),
  #higher ratioA (% directly explained by A linearly), lower ERV expected
  ratioA = c(0.2,0.5),
  # higher auto correlation, lower ERV expected
  autocorrA = c(0.25,0.75), 
  #expects no relationship with resample mean (starting value of ER)
  resample_mean = c(0.5),
  #higher resample SD, higher ERV expected
  resample_sd = c(0.15), #c(0.15,0.35),
  #Number of ER strategies: expects no relationship
  ERn = c(1),
  #max of scale: expects no relationship
  scalemax = c(100)
)


# ==================================
# define functions
# ==================================
funsim <- function(autocorrA,ratioA,independence, resample_mean, resample_sd, ERn, scalemax){
    simOutput <- data.frame(n = n, autocorrA = autocorrA, ratioA = ratioA,
                          independence=independence,
                          resample_mean = resample_mean,
                          resample_sd = resample_sd,
                          ERn = ERn,
                          scalemax = scalemax)
  
  # curve to sample the starting value of ER strategy A[1]
  resample <- function(n){
    # different curves can be set here other than the normal curve
    return (rnorm(n, mean = resample_mean*scalemax, sd = resample_sd*scalemax))
  } 
    
    
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

  #Create "blank" strategies
  if(ERblankn>0){
    for (i in (ERn+1):(ERn+ERblankn)){
      dfSim[letters[i]] <- 0
    }
  }
  
  #limit the timeseries to min/max
  dfSim[] <- apply(dfSim, 2, function(x) ifelse(x < scalemin , scalemin, x))
  dfSim[] <- apply(dfSim, 2, function(x) ifelse(x > scalemax , scalemax, x))
  
  # testing purpose: first row to zero
  if(firstrowzero){
    dfSim[1,] <- 0
  }
  
  # testing purpose: non overlapping ER strategy use
  if(testnonoverlap & ((ERn+ERblankn)>1)){
    dfSim[1,1] <- 0
    dfSim[2,2] <- 0
  }
  
  # rounding
  if(rounding){
    dfSim[] <- round(dfSim,0)
  }
    
  # replace 0 with 0.0001 because 0 are not true zero (interval scale not ratio scale)
  if (zerotransform){
  dfSim[dfSim == 0] <- 0.0001
  }
  
  #---- ER Variability candidate indices
  
  dfNew <- dfSim
  
  # Momentary SD
  dfNew$sd <- base::apply(dfSim,1,FUN = sd)
  
  
  # standardize
  dfnorm_chord <- decostand(dfSim,"norm")
  dfnorm_logchord <- decostand(log1p(dfSim),"norm")
  dfnorm_hel <- decostand(dfSim,"hel")
  dfnorm_prob <- decostand(dfSim,"total")
  if (ERblankn == 0  | zerotransform == TRUE){
    dfnorm_chi <- decostand(dfSim,"chi.square")
  }
  
  # KL divergence 
  
  tempdist <- c()
  for (k in 1:(nrow(dfnorm_prob)-1)){
    for (j in (k+1):nrow(dfnorm_prob)){
      tempdist <- append(tempdist,(KL.plugin(dfnorm_prob[k,],dfnorm_prob[j,])))
    }
  }
  tempdist <- as.matrix(tempdist)
  mat.KLdiv <- matrix(0, nrow = nrow(dfnorm_prob), ncol = nrow(dfnorm_prob))
  mat.KLdiv[lower.tri(mat.KLdiv, diag = FALSE)] <- tempdist
  mat.KLdiv[upper.tri(mat.KLdiv)] <- t(mat.KLdiv)[upper.tri(mat.KLdiv)]
  
  # Other matrices
  
  mat.euc <- dist(dfSim)
  mat.manhattan <- vegdist(dfSim,method = "manhattan")
  mat.chord <- vegdist(dfnorm_chord, method = "euc") # or vegdist(dfSim,method="chord")
  mat.logchord <- vegdist(dfnorm_logchord, method = "euc")
  mat.chisq <- vegdist(dfnorm_chi, method = "euc") # or vegdist(dfSim,method="chisq")
  mat.hel <- vegdist(dfnorm_hel, method = "euc")
  mat.jaccard <- vegdist(dfSim,method="jaccard")
  mat.kulczynski <- vegdist(dfSim,method="kulczynski")
  mat.bray <- vegdist(dfSim,method="bray")
  
  
  # Momentary ERV: different measures 
  for (i in 1:n){
    
    # x <- dfSim[i,]
    # y <- colMeans(dfSim)
    # xy <-x + y 
    # y. <-y / sum(y) 
    # x. <-x / sum(x) 
    
    dfNew$edist[i] <- dis_from_moment(mat.euc,i)
    dfNew$manhattan[i] <- dis_from_moment(mat.manhattan,i)

    
    

    # **zero row handling**
    # chord, logchord, chisq, hellinger appear to give okay results 
    # when 0 are replaced by 0.0001
    
    dfNew$chordnormed[i] <- dis_from_moment(mat.chord,i)
    dfNew$logchord[i] <- dis_from_moment(mat.logchord,i)
    
    # chi sq cannot handle blank strategies
    if (ERblankn > 0 & zerotransform==FALSE){
      #dfNew$chinormed[i] <- NA
      dfNew$chisq[i] <- NA
    }else{
      #dfNew$chinormed[i] <- dis_from_moment(mat.chisq, i)
      dfNew$chisq[i] <- dis_from_moment(mat.chisq,i)
    }
    

    dfNew$hellinger[i] <- dis_from_moment(mat.hel,i) 
    dfNew$chord[i] <- dis_from_moment(mat.chord,i)

    # **zero row handling**
    #jaccard & bray approaches 1 when all elements in a row approaches 0
    #kulczynski approaches 0.5 when all elements in a row approaches 0
    
    dfNew$jaccard[i] <- dis_from_moment(mat.jaccard,i)
    dfNew$kulczynski[i] <- dis_from_moment(mat.kulczynski,i)
    dfNew$brayveg[i] <- dis_from_moment(mat.bray,i)
    
    # KL Divergence
    
    dfNew$KLdiv[i] <- dis_from_moment(mat.KLdiv,i)
    
    # the below are old lines: distance between moment and mean 
    #   is different from mean distance of the moment to all other moments
    #   might be needed to handle the complimentary zero problem (0,1,0) vs (1,0,1)

    #dfNew$edist[i] <- (dist(rbind(x,y)))
    #dfNew$manhattan[i] <- sum(abs(x -y))
    #dfNew$chordnormed[i] <- vegdist(rbind(dfnorm_chord[i,],colMeans(dfnorm_chord)),"euc")
    #dfNew$logchordnormed[i] <- vegdist(rbind(dfnorm_logchord[i,],colMeans(dfnorm_logchord)),"euc")
    #dfNew$chinormed[i] <- vegdist(rbind(dfnorm_chi[i,],colMeans(dfnorm_chi)),"euc")
    #dfNew$hellinger[i] <- sqrt(1/2 * sum(sqrt(x) -sqrt(y))^2)
    #dfNew$hellingernormed[i] <- vegdist(rbind(dfnorm_hel[i,],colMeans(dfnorm_hel)),"euc") 
    #dfNew$chord[i] <- vegdist(rbind(x,y),method="chord")[1]
    #dfNew$jaccard[i] <- vegdist(rbind(x,y),method="jaccard")[1]
    #dfNew$chisq[i] <- vegdist(rbind(x,y),method="chisq")[1]
    #dfNew$kulczynski[i] <- vegdist(rbind(x,y),method="kulczynski")[1]
    #dfNew$brayveg[i] <- vegdist(rbind(x,y),method="bray")[1]
    
  }
  
  # Momentary Bray-Curtis dissimilarity by beta.part
  
  resbraypart <- tempres <- bray.part(dfSim)
  
  if ((ERn+ERblankn) > 0){
      # there appears to be some limitation on the bray.part on handling 1 ER stgy only
      # will throw 0 and warning message if there is only 1 ER stgy
      for (i in 1:n){
      dfNew$bray[i] <- dis_from_moment(resbraypart$bray,i)
      dfNew$bray.bal[i] <- dis_from_moment(resbraypart$bray.bal,i)
      dfNew$bray.gra[i] <- dis_from_moment(resbraypart$bray.gra,i)
      #old codes that compares moment with mean values
      #resbraypart <- bray.part(rbind(dfSim[i,],base::lapply(dfSim[,],FUN = mean)))
      #dfNew$bray[i] <- resbraypart$bray[1]
      #dfNew$bray.bal[i] <- resbraypart$bray.bal[1]
      #dfNew$bray.gra[i] <- resbraypart$bray.gra[1]
      }
  }else{
    dfNew$bray<-0
    dfNew$bray.bal<-0
    dfNew$bray.gra<-0
  }
  
  
  
  # Multi-site Bray-Curtis dissimilarity
  resBray<- beta.multi.abund(dfSim)
  
  
  #---- simOutput
  if ((ERn+ERblankn)>1){
  tempmodel <- lm (a ~ ., data = dfSim)
  }
  simOutput$em_autocorrA<-acf(dfSim$a, plot = FALSE)$acf[2] #empirical auto correlation
  simOutput$em_cor<- ifelse(ERn > 1,mean(cor(dfSim)[1,2:length(cor(dfSim))^0.5]),NA) #empirical correlation
  simOutput$em_r <- ifelse((ERn+ERblankn) > 1,sqrt(summary(tempmodel)$r.squared),NA)
  simOutput$mean_sd <- mean(dfNew$sd) 
  simOutput$mean_edist <- mean(dfNew$edist) 
  simOutput$mean_manhattan <- mean(dfNew$manhattan)
  #simOutput$mean_chinormed <- mean(dfNew$chinormed)
  simOutput$mean_hellinger <- mean(dfNew$hellinger)
  simOutput$mean_jaccard <- mean(dfNew$jaccard)
  simOutput$mean_chisq <- mean(dfNew$chisq) 
  simOutput$mean_logchord <- mean(dfNew$logchord) 
  simOutput$mean_chord <- mean(dfNew$chord)
  #simOutput$mean_chordnormed <- mean(dfNew$chordnormed)
  simOutput$mean_kulczynski <- mean(dfNew$kulczynski)
  simOutput$mean_KLdiv <- mean(dfNew$KLdiv)
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
  dftemp <- data.frame(measure = dfreturn[,measure],
                       autocorrA = dfreturn$autocorrA,
                       ratioA = dfreturn$ratioA,
                       # emautocorr = dfreturn$em_autocorrA,
                       # emcor = dfreturn$em_r,
                       independence = dfreturn$independence,
                       resample_mean = dfreturn$resample_mean,
                       resample_sd = dfreturn$resample_sd,
                       ERn = dfreturn$ERn,
                       scalemax = dfreturn$scalemax)
  if (is.na(dfreturn[1,measure])){
    tempfitres <-  (data.frame(NA,NA,NA,NA,NA,NA,NA,NA))
    }else{
    tempfitres <- (pcor(dftemp)$estimate[1,])
    }
  #renaming is needed because if sim parameters doesn't change the pcor does not return names
  names(tempfitres) <- names(dftemp)
  return(tempfitres)
}

# ====================================
# actually running the simulation
# ====================================


resOutput <- data.frame()
for (i in 1:nrow(siminput)){
  resSim<-rdply(simn,funsim(siminput$autocorrA[i],siminput$ratioA[i],
                            siminput$independence[i],
                            siminput$resample_mean[i],
                            siminput$resample_sd[i],
                            siminput$ERn[i],
                            siminput$scalemax[i]))
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


returninput <- c("multibray","mean_sd","mean_edist", 
                 "mean_manhattan","mean_hellinger",
                 "mean_jaccard",
                 "mean_chisq", 
                 "mean_logchord","mean_chord",
                 "mean_kulczynski",
                 "mean_KLdiv",
                 "mean_bray")
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
