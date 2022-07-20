library(betapart) # partitioned bray-curtis dissimilarity
#library(dplyr)
library(ppcor) # partial correlation
library(vegan) # dissimilarity measures
library(entropy) # KL divergence
library(tsDyn) # package for VAR.sim

options(scipen=999)
#======================
# notes/to-do
#======================

# 1. ERn = 1 case can be handled by adding a blank strategy (set ERblankn>=1)
# 2. Explore if there's a way out other than choosing between an unwanted
    #-ve association with ER_mean or +ve association with scalemax.
    #E.g., Extra adjustment, say, add ERn and scalemax to denominator of Euclidean distance?
# 3. Produce data visualizations to see if relationships are linear or not
# 4. Add RQA
# 5. Handle missing rows (in empirical dataset)

#======================
# simulation setup
#======================

# study design
n <- 2 #number of observation (time points per participant)
simn <- 3 #number of simulations (participant)
scalemin <- 0

# other settings and testing conditions
ERblankn <- 1 # number of "unused" (always 0) ER strategies. Always set >=1 (between-SD needs at least 2 to calculate)
zerotransform <- TRUE # whether or not to replace 0 with 0.0001. TRUE recommended because 0 are not true 0
rounding <- TRUE # whether or not round scores to integers
successivecomp <- TRUE # is successive comparisons made? (one-to-many is in by default)
firstrowzero <- FALSE # testing purpose: set TRUE to force first row to be all 0
testnonoverlap <- FALSE # testing purpose: set TRUE to make [1,1] = 0 and [2,2] = 0
set.seed(1999)

# varying simulation parameters
siminput <- expand.grid(
  # higher meanshift, less rank changes occur, lower ERV expected
  meanshift = c(0.0,2.0),
  # higher auto correlation, lower ERV expected
  autocorr = c(0.25,0.75),
  correlation = c(0),
  #  mean ER endorsement: expects no relationship
  ER_mean = c(0.2,0.3,0.4),
  # higher within-strategy SD, higher ERV expected
  ER_withinSD = c(0.14,0.20,0.26), #c(0.15,0.35),
  # Number of ER strategies: expects no relationship
  ERn = c(1,2,5),
  # max of scale: expects no relationship
  scalemax = c(10,100)
)


# ==================================
# define functions
# ==================================
funsim <- function(i.siminput){
    autocorr = siminput$autocorr[i.siminput]
    meanshift = siminput$meanshift[i.siminput]
    correlation = siminput$correlation[i.siminput]
    ER_mean = siminput$ER_mean[i.siminput]
    ER_withinSD = siminput$ER_withinSD[i.siminput]
    ERn = siminput$ERn[i.siminput]
    scalemax = siminput$scalemax[i.siminput]


  #---- simulating ER strategies
    # create an alternating -1,1,... vector for mean-shifting other strategies
    signvector <- sign(rnorm(1))*as.vector(rbind(rep(1,ERn), rep(-1,ERn)))
    signvector[1] <- 0
    signvector <- signvector[1:ERn]
    # adjust it so that the overall mean (of all strategies) remain as first specified
    signvector <- signvector-sum(signvector)/ERn
    # repeat signvector n times to fit the # of observations
    signvector <- rep(signvector, each = n)

    # # ----------------------------------------------------
    # # Using the VAR.sim method...
    # # VAR.sim must create multivariate time series.
    # if (ERn >1){
    #   tmpERn <- ERn
    # }else{
    #   tmpERn <- 2
    # }
    # CV <- diag(1, tmpERn,tmpERn)
    # CV[CV == 0] <- correlation
    # B1 <- diag(autocorr, tmpERn,tmpERn)
    # # B1[B1 == 0] <- 0 # can replace by cross-lagged paremeter if needed
    # if (ERn >1){
    #   dfSim <- (VAR.sim(B=B1, n=n, include="none",varcov=CV))
    # }else{
    #   # VAR.sim must create multivariate time series;
    #   # cut back to 1 col
    #   dfSim <- (VAR.sim(B=B1, n=n, include="none",varcov=CV)[,1])
    # }
    # # ----------------------------------------------------


  # ----------------------------------------------------
  # Using mvrnorm and meanshift to create strategies
  # Generate main strategy: by mvrnorm.
  # (because by arima.sim the SD is higher than specified and by the old way the SD is lower than specified.)
  tmp.r <- matrix(autocorr, n, n)
  tmp.r <- tmp.r^abs(row(tmp.r)-col(tmp.r))
  tmp.dist <- mvrnorm(1, rep(0,n), tmp.r)
  dfSim <- tmp.dist

  #Create other strategies
  if(ERn>1){
    for (i in 2:ERn){
      tmp.ER <-  (tmp.dist + rnorm(n,0,1))*0.5*sqrt(2)
      dfSim <- cbind(dfSim, tmp.ER)
     }
  }
  # ----------------------------------------------------


  # apply meanshift adjustment
    dfSim <- dfSim + signvector*meanshift
  # Scale up the simulated data to match the mean, SD and scalemax parameters
   dfSim <- dfSim*ER_withinSD*scalemax
   dfSim <- dfSim+ER_mean*scalemax

  #Create "blank" strategies
  if(ERblankn>0){
    dfSim<- cbind(dfSim,matrix(0,n,ERblankn,dimnames = list(1:n,letters[(ERn+1):(ERn+ERblankn)])))
    # for (i in (ERn+1):(ERn+ERblankn)){
    #   dfSim[letters[i]] <- 0
    # }
  }
  colnames(dfSim) <- letters[1:(ERn+ERblankn)]


  #limit the timeseries to min/max
  # CJ: Probably remove this, because it relates to measurement
  # CJ: But in general, boolean indexing is faster than ifelse()
  dfSim[dfSim < scalemin] <- scalemin
  dfSim[dfSim > scalemax] <- scalemax


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
  # CJ: This relates to measurement, so remove for now
  if(rounding){
    dfSim <- round(dfSim,0)
  }

  # replace 0 with 0.0001 because 0 are not true zero (interval scale not ratio scale)
  # CJ: This problem will probably go away if you're no longer rounding
  # CJ: Instead of adding a small constant here, you can either:
  # CJ: Create a wrapper for your metric functions that adds a small constant if necessary, OR create a wrapper that returns NA if the calculation fails
  if (zerotransform){
  dfSim[dfSim == 0] <- 0.0001
  }

  return(dfSim)
}

funcal <- function(i.siminput,dfSim){

  #---- ER Variability candidate indices

  dfNew <- dfSim

  # Momentary SD
  dfNew <- cbind(dfNew,sd = apply(dfSim,1,sd))

  # standardize
  dfnorm_chord <- decostand(dfSim,"norm")
  dfnorm_logchord <- decostand(log1p(dfSim),"norm")
  dfnorm_hellinger <- decostand(dfSim,"hellinger")
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

  resbraypart <- bray.part(dfSim)

  mat.euclidean <- as.matrix(dist(dfSim))
  mat.manhattan <- as.matrix(vegdist(dfSim,method = "manhattan"))
  mat.chord <- as.matrix(vegdist(dfnorm_chord, method = "euclidean")) # or vegdist(dfSim,method="chord")
  mat.logchord <- as.matrix(vegdist(dfnorm_logchord, method = "euclidean"))
  mat.chisq <- as.matrix(vegdist(dfnorm_chi, method = "euclidean")) # or vegdist(dfSim,method="chisq")
  mat.hellinger <- as.matrix(vegdist(dfnorm_hellinger, method = "euclidean"))
  mat.jaccard <- as.matrix(vegdist(dfSim,method="jaccard"))
  mat.kulczynski <- as.matrix(vegdist(dfSim,method="kulczynski"))
  mat.brayveg <- as.matrix(vegdist(dfSim,method="bray"))
  mat.braypart.all <- as.matrix(resbraypart$bray)
  mat.braypart.bal <- as.matrix(resbraypart$bray.bal)
  mat.braypart.gra <- as.matrix(resbraypart$bray.gra)

  # if successive comparison, only observations 2:n are used because the 1st doesn't have a previous comparison
  if (successivecomp){
    istart <- 2
  }else{
    istart <- 1
  }


  # Momentary ERV: different measures

  mom.euclidean <- apply(mat.euclidean,1,mean)*n/(n-1)
  mom.manhattan  <- apply(mat.manhattan,1,mean)*n/(n-1)
  mom.chord <- apply(mat.chord,1,mean)*n/(n-1)
  mom.chisq  <- apply(mat.chisq,1,mean)*n/(n-1)
  mom.logchord <- apply(mat.logchord,1,mean)*n/(n-1)
  mom.hellinger <- apply(mat.hellinger,1,mean)*n/(n-1)
  mom.jaccard<- apply(mat.jaccard,1,mean)*n/(n-1)
  mom.kulczynski <- apply(mat.kulczynski,1,mean)*n/(n-1)
  mom.brayveg <- apply(mat.brayveg,1,mean)*n/(n-1)
  mom.braypart.all  <- apply(mat.braypart.all,1,mean)*n/(n-1)
  mom.braypart.bal <- apply(mat.braypart.bal,1,mean)*n/(n-1)
  mom.braypart.gra <- apply(mat.braypart.gra,1,mean)*n/(n-1)
  mom.KLdiv <- apply(mat.KLdiv,1,mean)*n/(n-1)

  if (successivecomp){

    suc.euclidean <- dis_suc_vector(mat.euclidean)
    suc.manhattan <- dis_suc_vector(mat.manhattan)
    suc.chord <- dis_suc_vector(mat.chord)
    suc.chisq <- dis_suc_vector(mat.chisq)
    suc.logchord <- dis_suc_vector(mat.logchord)
    suc.hellinger <- dis_suc_vector(mat.hellinger)
    suc.jaccard <- dis_suc_vector(mat.jaccard)
    suc.kulczynski <- dis_suc_vector(mat.kulczynski)
    suc.brayveg <- dis_suc_vector(mat.brayveg)
    suc.braypart.all <- dis_suc_vector(mat.braypart.all)
    suc.braypart.bal <- dis_suc_vector(mat.braypart.bal)
    suc.braypart.gra <- dis_suc_vector(mat.braypart.gra)
    suc.KLdiv <- dis_suc_vector(mat.KLdiv)

    mom.euclidean <- (mom.euclidean+suc.euclidean)/2
    mom.manhattan  <-(mom.manhattan+suc.manhattan)/2
    mom.chord <- (mom.chord + suc.chord)/2
    mom.chisq  <- (mom.chisq + suc.chisq)/2
    mom.logchord <- (mom.logchord + suc.logchord)/2
    mom.hellinger <- (mom.hellinger + suc.hellinger)/2
    mom.jaccard<- (mom.jaccard + suc.jaccard)/2
    mom.kulczynski <- (mom.kulczynski + suc.kulczynski)/2
    mom.brayveg <- (mom.brayveg + suc.brayveg)/2
    mom.braypart.all  <- (mom.braypart.all+suc.braypart.all)/2
    mom.braypart.bal <- (mom.braypart.bal+suc.braypart.bal)/2
    mom.braypart.gra <- (mom.braypart.gra+suc.braypart.gra)/2
    mom.KLdiv <- (mom.KLdiv+suc.KLdiv)/2
  }


  dfNew <- cbind(dfNew,euclidean = mom.euclidean)
  dfNew <- cbind(dfNew,manhattan = mom.manhattan)
  dfNew <- cbind(dfNew,chord = mom.chord)
  dfNew <- cbind(dfNew,chisq = mom.chisq)
  dfNew <- cbind(dfNew,logchord = mom.logchord)
  dfNew <- cbind(dfNew,hellinger = mom.hellinger)
  dfNew <- cbind(dfNew,jaccard = mom.jaccard)
  dfNew <- cbind(dfNew,kulczynski = mom.kulczynski)
  dfNew <- cbind(dfNew,brayveg = mom.brayveg)
  dfNew <- cbind(dfNew,braypart.all = mom.braypart.all)
  dfNew <- cbind(dfNew,braypart.bal = mom.braypart.bal)
  dfNew <- cbind(dfNew,braypart.gra = mom.braypart.gra)
  dfNew <- cbind(dfNew,KLdiv = mom.KLdiv)

  # Multi-site Bray-Curtis dissimilarity
  resBray<- beta.multi.abund(dfSim)

  #---- simOutput

  simOutput <- c(em_autocorr <- acf(dfSim[,"a"], plot = FALSE)$acf[2],
                 em_cor<- tryCatch(mean(cor(dfSim[,1:(ncol(dfSim)-1)])[1,2:length(cor(dfSim[,1:(ncol(dfSim)-1)]))^0.5]), error=function(err) NA),
                 em_r <- tryCatch(sqrt(summary(lm(a ~ ., data = as.data.frame(dfSim)))$r.squared), error=function(err) NA),
                 mean_sd <- mean(dfNew[istart:n,"sd"]) ,
                 mean_euclidean <- mean(dfNew[istart:n,"euclidean"]) ,
                 mean_manhattan <-mean(dfNew[istart:n,"manhattan"]) ,
                 mean_hellinger <- mean(dfNew[istart:n,"hellinger"]) ,
                 mean_jaccard <- mean(dfNew[istart:n,"jaccard"]) ,
                 mean_chisq <-mean(dfNew[istart:n,"chisq"]) ,
                 mean_logchord <-mean(dfNew[istart:n,"logchord"]) ,
                 mean_chord <-mean(dfNew[istart:n,"chord"]) ,
                 mean_kulczynski <-mean(dfNew[istart:n,"kulczynski"]) ,
                 mean_KLdiv <-mean(dfNew[istart:n,"KLdiv"]) ,
                 mean_brayveg <- mean(dfNew[istart:n,"brayveg"]) ,
                 mean_bray.all <-mean(dfNew[istart:n,"braypart.all"]) ,
                 mean_bray.bal <-mean(dfNew[istart:n,"braypart.bal"]) ,
                 mean_bray.gra <-mean(dfNew[istart:n,"braypart.gra"]) ,
                 multibray.bal <- resBray$beta.BRAY.BAL,
                 multibray.gra <- resBray$beta.BRAY.GRA,
                 multibray.all <- resBray$beta.BRAY
                 )

    return(simOutput)

}

simulatecalculate <- function(i.siminput){
  return (funcal(i.siminput,funsim(i.siminput)))
}

#function to calculate the mean dissimilarity from one obs to all other obs, input dist
dis_from_moment <- function(distobj,obs){
  return (mean(as.matrix(distobj)[obs,])*n/(n-1))
}

dis_suc_vector <- function(distmat){
  return (c(0,(distmat[row(distmat) == col(distmat) + 1])))
}

# function to return partial correlations ofdissimlarity measures to the 3 parameters
returnfit <- function(measure,dfreturn,testtype){
  dftemp <- data.frame(measure = dfreturn[,measure],
                       autocorr = dfreturn$autocorr,
                       meanshift = dfreturn$meanshift,
                       correlation = dfreturn$correlation,
                       # emautocorr = dfreturn$em_autocorr,
                       # emcor = dfreturn$em_r,
                       ER_mean = dfreturn$ER_mean,
                       ER_withinSD = dfreturn$ER_withinSD,
                       ERn = dfreturn$ERn,
                       scalemax = dfreturn$scalemax)
    dftemp <- Filter(function(x)(length(unique(x))>1), dftemp)
    if (is.na(dfreturn[1,measure])){
    tempfitres <-  (data.frame(rep(NA, length(dftemp))))
    }else{
      if (testtype=="pcor"){
      tempfitres <- (pcor(dftemp)$estimate[1,])
      }
      if (testtype=="cor"){
      tempfitres <- (cor(dftemp)[1,])
      }
    }
  #renaming is needed because if sim parameters doesn't change the pcor does not return names
  names(tempfitres) <- names(dftemp)
  return(tempfitres)
}

# ====================================
# actually running the simulation
# ====================================


resOutput <- c()
for (i in 1:nrow(siminput)){
  resSim<-t(replicate(simn,simulatecalculate(i)))
  if(length(resOutput)==0){
      resOutput <- resSim
    }else{
      resOutput <- rbind(resOutput,resSim)
  }
}
colnames(resOutput) <- c("em_autocorr",
                         "em_cor",
                         "em_r",
                         "mean_sd",
                         "mean_euclidean",
                         "mean_manhattan",
                         "mean_hellinger",
                         "mean_jaccard",
                         "mean_chisq",
                         "mean_logchord",
                         "mean_chord",
                         "mean_kulczynski",
                         "mean_KLdiv",
                         "mean_brayveg",
                         "mean_braypart.all",
                         "mean_braypart.bal",
                         "mean_braypart.gra",
                         "multibray.bal",
                         "multibray.gra",
                         "multibray.all"
)

#bind the simulation parameters to result
resInfo <- siminput[rep(seq_len(nrow(siminput)), each = simn), ]
resOutput <- cbind(resInfo,resOutput)

# ============================================================
# calculate the fit with partial-correlation and correlation
# ============================================================


returninput <- c("mean_sd",
                 "mean_euclidean",
                 "mean_brayveg",
                 "multibray.all",
                 "mean_hellinger",
                 "mean_jaccard",
                 "mean_manhattan",
                 "mean_chisq",
                 "mean_chord",
                 "mean_logchord",
                 "mean_kulczynski",
                 "mean_KLdiv"
                 )
resSummary_partialcor <- data.frame()
for (i in 1:length(returninput)){
  if(length(resSummary_partialcor)==0){
    resSummary_partialcor <- data.frame(as.list(returnfit(returninput[i],resOutput,"pcor")))
  }else{
    tempdf <- data.frame(as.list(returnfit(returninput[i],resOutput,"pcor")))
    names(tempdf) <- names(resSummary_partialcor)
    resSummary_partialcor <- rbind.data.frame(resSummary_partialcor,tempdf)
  }
  resSummary_partialcor[i,1] <- returninput[i]
}

resSummary_cor <- data.frame()
for (i in 1:length(returninput)){
  if(length(resSummary_cor)==0){
    resSummary_cor <- data.frame(as.list(returnfit(returninput[i],resOutput, "cor")))
  }else{
    tempdf <- data.frame(as.list(returnfit(returninput[i],resOutput, "cor")))
    names(tempdf) <- names(resSummary_cor)
    resSummary_cor <- rbind.data.frame(resSummary_cor,tempdf)
  }
  resSummary_cor[i,1] <- returninput[i]
}


resSummary_partialcor
resSummary_cor
