library(betapart)
library(plyr)
library(ppcor)
library(vegan)
library(entropy)
# library(faux) #doesn't need it anymore if not using mvrnorm
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
simn <- 2 #number of simulations (participant)
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
  meanshift = c(0.0,2),
  # higher auto correlation, lower ERV expected
  autocorr = c(0.25,0.75),
  # correlation: expects no relationship
  correlation = c(0.0,0.5),
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
    # Using the VAR.sim method...
    # VAR.sim must create multivariate time series.
    if (ERn >1){
      tmpERn <- ERn
    }else{
      tmpERn <- 2
    }
    CV <- diag(1, tmpERn,tmpERn)
    CV[CV == 0] <- correlation
    B1 <- diag(autocorr, tmpERn,tmpERn)
    # B1[B1 == 0] <- 0 # can replace by cross-lagged paremeter if needed

    # create an alternating -1,1,... vector for mean-shifting other strategies
    signvector <- sign(rnorm(1))*as.vector(rbind(rep(1,ERn), rep(-1,ERn)))
    signvector[1] <- 0
    signvector <- signvector[1:ERn]
    # adjust it so that the overall mean (of all strategies) remain as first specified
    signvector <- signvector-sum(signvector)/ERn
    # repeat it n times to fit the # of observations
    signvector <- rep(signvector, each = n)

    if (ERn >1){
      dfSim <- (VAR.sim(B=B1, n=n, include="none",varcov=CV))
    }else{
      # VAR.sim must create multivariate time series;
      # cut back to 1 col
      dfSim <- (VAR.sim(B=B1, n=n, include="none",varcov=CV)[,1])
    }

    # apply meanshift adjustment
    dfSim <- dfSim + signvector*meanshift
    # Scale up to within-strategy SD and ER mean endorsement
    dfSim <- dfSim*ER_withinSD*scalemax
    dfSim <- dfSim+ER_mean*scalemax


  # Using mvrnorm and meanshift to create strategies

  # dfSim <- data.frame(a = 1:n)
  # # Generate main strategy: by mvrnorm.
  # # (because by arima.sim the SD is higher than specified and by the old way the SD is lower than specified.)
  # tmp.r <- matrix(autocorr, n, n)
  # tmp.r <- tmp.r^abs(row(tmp.r)-col(tmp.r))
  # tmp.dist <- mvrnorm(1, rep(0,n), tmp.r)
  # dfSim$a <- tmp.dist * ER_withinSD*scalemax + ER_mean*scalemax
  #
  # #Create other strategies
  # if(ERn>1){
  #   for (i in 2:ERn){
  #     dfSim[letters[i]] <- scalemax*(   ER_mean +
  #                                         # varies relative to ratings of 1st ER
  #                                         (tmp.dist + rnorm(n,0,1))*0.5*sqrt(2)*ER_withinSD +
  #                                        # mean shift, random direction per person
  #                                        (sign(rnorm(1))*meanshift*ER_withinSD)
  #                                    )
  #   }
  # }

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
  dfSim <- apply(dfSim, 2, function(x) ifelse(x < scalemin , scalemin, x))
  dfSim <- apply(dfSim, 2, function(x) ifelse(x > scalemax , scalemax, x))

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
  # if(rounding){
  #   dfSim[] <- round(dfSim,0)
  # }

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
  resbraypart <- bray.part(dfSim)


  # Momentary ERV: different measures
  for (i in 1:n){

    if (successivecomp & i>1){
      suc.edist <- as.matrix(mat.euc)[i,i-1]
      suc.manhattan <- as.matrix(mat.manhattan)[i,i-1]
      suc.chord <- as.matrix(mat.chord)[i,i-1]
      suc.chisq <- as.matrix(mat.chisq)[i,i-1]
      suc.logchord <- as.matrix(mat.logchord)[i,i-1]
      suc.hel <- as.matrix(mat.hel)[i,i-1]
      suc.jaccard <- as.matrix(mat.jaccard)[i,i-1]
      suc.kulczynski <- as.matrix(mat.kulczynski)[i,i-1]
      suc.brayveg <- as.matrix(mat.bray)[i,i-1]
      suc.bray <- as.matrix(resbraypart$bray)[i,i-1]
      suc.bray.bal <- as.matrix(resbraypart$bray.bal)[i,i-1]
      suc.bray.gra <- as.matrix(resbraypart$bray.gra)[i,i-1]
    }else{
      suc.edist <- 0
      suc.manhattan <- 0
      suc.chord <- 0
      suc.chisq <- 0
      suc.logchord <- 0
      suc.hel <- 0
      suc.jaccard <- 0
      suc.kulczynski <- 0
      suc.brayveg <- 0
      suc.bray <- 0
      suc.bray.bal <- 0
      suc.bray.gra <- 0

          }

    dfNew$edist[i] <- dis_from_moment(mat.euc,i) + suc.edist
    dfNew$manhattan[i] <- dis_from_moment(mat.manhattan,i) + suc.manhattan

    # **zero row handling**
    # chord, logchord, chisq, hellinger appear to give okay results
    # when 0 are replaced by 0.0001

    dfNew$logchord[i] <- dis_from_moment(mat.logchord,i) + suc.logchord

    # chi sq cannot handle blank strategies
    if (ERblankn > 0 & zerotransform==FALSE){
      dfNew$chisq[i] <- NA
    }else{
      dfNew$chisq[i] <- dis_from_moment(mat.chisq,i) + suc.chisq
    }


    dfNew$hellinger[i] <- dis_from_moment(mat.hel,i) + suc.hel
    dfNew$chord[i] <- dis_from_moment(mat.chord,i) + suc.chord

    # **zero row handling**
    #jaccard & bray approaches 1 when all elements in a row approaches 0
    #kulczynski approaches 0.5 when all elements in a row approaches 0

    dfNew$jaccard[i] <- dis_from_moment(mat.jaccard,i) + suc.jaccard
    dfNew$kulczynski[i] <- dis_from_moment(mat.kulczynski,i) + suc.kulczynski
    dfNew$brayveg[i] <- dis_from_moment(mat.bray,i) + suc.brayveg

    # KL Divergence

    dfNew$KLdiv[i] <- dis_from_moment(mat.KLdiv,i)

    # Momentary Bray-Curtis dissimilarity by beta.part
    # there appears to be some limitation on the bray.part on handling 1 ER stgy only
    # will throw 0 and warning message if there is only 1 ER stgy
    dfNew$bray[i] <- dis_from_moment(resbraypart$bray,i) + suc.bray
    dfNew$bray.bal[i] <- dis_from_moment(resbraypart$bray.bal,i) + suc.bray.bal
    dfNew$bray.gra[i] <- dis_from_moment(resbraypart$bray.gra,i) + suc.bray.gra


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
    #resbraypart <- bray.part(rbind(dfSim[i,],base::lapply(dfSim[,],FUN = mean)))
    #dfNew$bray[i] <- resbraypart$bray[1]
    #dfNew$bray.bal[i] <- resbraypart$bray.bal[1]
    #dfNew$bray.gra[i] <- resbraypart$bray.gra[1]

  }

  # Multi-site Bray-Curtis dissimilarity
  resBray<- beta.multi.abund(dfSim)

  if (successivecomp){
    istart <- 2
  }else{
    istart <- 1
  }

  #---- simOutput
  #empirical auto correlation
  simOutput <- data.frame(em_autocorr = NA)
  simOutput$em_autocorr<-acf(dfSim[,"a"], plot = FALSE)$acf[2]
  #has to exclude the blank column otherwise cor throws error
  simOutput$em_cor<- tryCatch(mean(cor(dfSim[1:(ncol(dfSim)-1)])[1,2:length(cor(dfSim[1:(ncol(dfSim)-1)]))^0.5]), error=function(err) NA)
  simOutput$em_r <- tryCatch(sqrt(summary(lm(a ~ ., data = dfSim))$r.squared), error=function(err) NA)
  simOutput$mean_sd <- mean(dfNew$sd[istart:n])
  simOutput$mean_edist <- mean(dfNew$edist[istart:n])
  simOutput$mean_manhattan <- mean(dfNew$manhattan[istart:n])
  #simOutput$mean_chinormed <- mean(dfNew$chinormed)
  simOutput$mean_hellinger <- mean(dfNew$hellinger[istart:n])
  simOutput$mean_jaccard <- mean(dfNew$jaccard[istart:n])
  simOutput$mean_chisq <- mean(dfNew$chisq[istart:n])
  simOutput$mean_logchord <- mean(dfNew$logchord[istart:n])
  simOutput$mean_chord <- mean(dfNew$chord[istart:n])
  simOutput$mean_kulczynski <- mean(dfNew$kulczynski[istart:n])
  simOutput$mean_KLdiv <- mean(dfNew$KLdiv[istart:n])
  simOutput$mean_brayveg <- mean(dfNew$brayveg[istart:n])
  simOutput$mean_bray <- mean(dfNew$bray[istart:n])
  simOutput$mean_bray.bal <- mean(dfNew$bray.bal[istart:n])
  simOutput$mean_bray.gra <- mean(dfNew$bray.gra[istart:n])
  simOutput$multibray.bal <- resBray$beta.BRAY.BAL
  simOutput$multibray.gra <- resBray$beta.BRAY.GRA
  simOutput$multibray <- resBray$beta.BRAY

  return(simOutput)

}

simulatecalculate <- function(i.siminput){
  return (funcal(i.siminput,funsim(i.siminput)))
}

#function to calculate the mean dissimilarity from one obs to all other obs, input dist
dis_from_moment <- function(distobj,obs){
  nobs <- n
  return (mean(as.matrix(distobj)[obs,])*nobs/(nobs-1))
}

# function to return partial correlations ofdissimlarity measures to the 3 parameters
returnfit <- function(measure,dfreturn){
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
  resSim<-rdply(simn,simulatecalculate(i))
  if(length(resOutput)==0){
      resOutput <- resSim
    }else{
      resOutput <- rbind.data.frame(resOutput,resSim)
  }
}

#bind the simulation parameters to result
resInfo <- siminput[rep(seq_len(nrow(siminput)), each = simn), ]
resOutput <- cbind(resInfo,resOutput)
# ============
# calculate the fit
# ============


returninput <- c("mean_sd",
                 "mean_edist",
                 "mean_bray",
                 "multibray",
                 "mean_hellinger",
                 "mean_jaccard",
                 "mean_manhattan",
                 "mean_chisq",
                 "mean_chord",
                 "mean_logchord",
                 "mean_kulczynski",
                 "mean_KLdiv"
                 )
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
