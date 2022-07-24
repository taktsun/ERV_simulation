library(betapart)
# library(plyr)
# library(tidyverse)
library(ppcor)
library(vegan)
library(entropy)
library(faux)

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
n <- 50 #number of observation (time points per participant)
simn <- 1 #number of simulations (participant)
scalemin <- 0

# other settings andtesting conditions
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
  autoregressive = c(0.25,0.75),
  #  mean ER endorsement: expects no relationship
  ER_mean = c(0.2,0.3,0.4),
  # higher within-strategy SD, higher ERV expected
  ER_withinSD = c(0.14,0.20,0.26), #c(0.15,0.35),
  # Number of ER strategies: expects no relationship
  # CJ: Minimum ER to 2 right now, some debugging necessary for 1
  ERn = c(2,3,5),
  # max of scale: expects no relationship
  scalemax = c(10,100)
)


# ==================================
# define functions
# ==================================
simulate_data <- function(n = 50, ERn = 2, autoregressive = 1, cross = 0, ER_withinSD = 1, ER_mean = 0, ...){
  # Assign autoregressive regression parameter to diagonal
  Bmat <- diag(autoregressive, ERn)
  # Assign cross-strategy regression coefficient to off-diagonal
  Bmat[lower.tri(Bmat)|upper.tri(Bmat)] <- cross
  # Generate data
  out <- VAR.sim(B = Bmat, n = n, lag = 1, include = "none", varcov = diag(ER_withinSD, nrow(Bmat)))
  # Center to zero
  out <- out - matrix(colMeans(out), ncol = ncol(out), nrow = nrow(out), byrow = TRUE)
  # Center to desired mean
  out + matrix(ER_mean, ncol = ncol(out), nrow = nrow(out))
}

funsim <- function(autoregressive,meanshift, ER_mean, ER_withinSD, ERn, scalemax){
  #---- simulating ER strategies
  cl <- match.call()
  cl[[1]] <- quote(simulate_data)
  dfSim <- eval.parent(cl)

  #---- ER Variability candidate indices
  # CJ: Avoid assignment and making copies of things
  # dfNew <- dfSim

  # Momentary SD
  mom_sd <- base::apply(dfSim,1,FUN = sd)


  # standardize
  # CJ: Avoid unnecessary dependencies; this is just basic algebra, so we don't need decostand. It's much slower
  #dfnorm_chord1 <- decostand(dfSim,"normalize")
  dfnorm_chord1 <- dfSim / sqrt(rowSums(dfSim^2))

  # dfnorm_logchord <- decostand(log1p(dfSim),"norm")
  dfSim_log1p <- log1p(dfSim)
  dfnorm_logchord <- dfSim_log1p / sqrt(rowSums(dfSim_log1p^2))

  # dfnorm_prob <- decostand(dfSim,"total")
  dfSim_rowsum <- rowSums(dfSim)
  dfSim_min <- min(dfSim)
  dfnorm_prob <- dfSim / pmax(dfSim_min, dfSim_rowsum)

  # dfnorm_hel <- decostand(dfSim,"hel")
  dfnorm_hel <- sqrt(dfnorm_prob)

  # CJ: Not sure what happens here tbh
  if (ERblankn == 0  | zerotransform == TRUE){
    # dfnorm_chi <- decostand(dfSim,"chi.square")
    dfnorm_chi <-   sqrt(sum(dfSim)) * dfSim/outer(pmax(dfSim_min, dfSim_rowsum), sqrt(colSums(dfSim)))
  }

  # KL divergence
  # CJ: Not sure what happens here and if it's correct. Either way, growing a vector in nested for-loops is very slow, so try to avoid that. The operation here can probably be executed as a single matrix algebra statement
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
  # CJ: I'm not sure how we should summarize the Euclidean distance,
  #     so here I'm taking the median of the row-to-row distance
  mat.euc <- median(diag(as.matrix(dist(dfSim))[-1, -nrow(dfSim)]))

  #mat.manhattan <- vegdist(dfSim,method = "manhattan")
  # CJ: Base R dist() has manhattan too!
  mat.euc <- median(diag(as.matrix(dist(dfSim, method = "manhattan"))[-1, -nrow(dfSim)]))
  # CJ: Same principle applies, we probably don't need vegdist
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
  # CJ: Try to simplify this and just output a vector with all of the metrics
  if ((ERn+ERblankn)>1){
  tempmodel <- lm (a ~ ., data = dfSim)
  }
  simOutput$em_autoregressive<-acf(dfSim$a, plot = FALSE)$acf[2] #empirical auto correlation
  simOutput$em_cor<- ifelse(ERn > 1,mean(cor(dfSim)[1,2:length(cor(dfSim))^0.5]),NA) #empirical correlation
  simOutput$em_r <- ifelse((ERn+ERblankn) > 1,sqrt(summary(tempmodel)$r.squared),NA)
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

#function to calculate the mean dissimilarity from one obs to all other obs, input dist
dis_from_moment <- function(distobj,obs){
  nobs <- n
  return (mean(as.matrix(distobj)[obs,])*nobs/(nobs-1))
}

# function to return partial correlations ofdissimlarity measures to the 3 parameters
returnfit <- function(measure,dfreturn){
  dftemp <- data.frame(measure = dfreturn[,measure],
                       autoregressive = dfreturn$autoregressive,
                       meanshift = dfreturn$meanshift,
                       # emautoregressive = dfreturn$em_autoregressive,
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
  resSim<-rdply(simn,funsim(siminput$autoregressive[i],siminput$meanshift[i],
                            siminput$ER_mean[i],
                            siminput$ER_withinSD[i],
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
