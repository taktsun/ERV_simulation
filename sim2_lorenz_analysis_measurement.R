library(nonlinearTseries)
library(betapart)
library(vegan)
library(data.table)
setDTthreads(threads = 1) # so it won't take up the parallel processing for running simulation
library(ppcor) # partial correlation
library(doSNOW)
library(parallel)

# How was this script amended to include the sensitivity analyses?
# Most changes were made in calcdis: there, the rounding and missingness conditions
# were added so that the same simulated dataset would go through those conditions
# and have ER variability being calculated.
# Wrappers MCARbyrow and calcdis_inside also added to make these changes more concise.

#======================
# simulation setup
#======================
# set seed to reproduce exact results
  seed = 1999
  set.seed(seed)

# simulation paramaters
  siminput <- expand.grid(
    # reps
    rep = 1000,
    # nobs of time series
    n = c(30,70,100),
    # Number of ER strategies
    ERn = c(3,6,9)
    )
# assign the a seed (dependent on the earlier seed set) to each simulation repetition
  siminput$seed <- sample(1:.Machine$integer.max, nrow(siminput))


# ==================================
# define functions
# ==================================
# self-defined functions
  source("func_indices.R")
  specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))



# a list of names that correspond to the indices output specified in calcdis
  list_metrics <- paste0(c(rep("orig_",6),
                           rep("m1",6),
                           rep("m2",6),
                           rep("m3",6),
                           rep("m4",6),
                           rep("m5",6),
                           rep("r4",6),
                           rep("r3",6),
                           rep("r2",6),
                           rep("r1",6),
                           rep("r0",6)
  ),
  c(rep(c("_repl.suc",
          "_nest.suc",
          "_full.suc",
          "_repl.amm",
          "_nest.amm",
          "_full.amm"),11)))

calcdis_inside <- function(matx, outname){
  tempbray.suc <- calc.bray.suc(matx)
  tempbray.amm <- calc.bray.amm(matx)

  out <- c(tempbray.suc[1],        # Bray-Curtis dissimilarity: replacement
           tempbray.suc[2],        # Bray-Curtis dissimilarity: nestedness
           tempbray.suc[3],        # Bray-Curtis dissimilarity: Full index
           tempbray.amm[1],        # Bray-Curtis dissimilarity: replacement
           tempbray.amm[2],        # Bray-Curtis dissimilarity: nestedness
           tempbray.amm[3]        # Bray-Curtis dissimilarity: Full index
  )
  names(out) <- paste0(outname,c(rep(".suc.",3),rep(".amm.",3)),names(out))
  out
}
MCARbyrow <- function(matx, prob){
  matx <- delete_MCAR(matx, prob, 1)
  matx[!!rowSums(is.na(matx)),] <- NA
  matx
}
# function to calculate dissimilarity indices from simulated datasets
  calcdis <- function(matx, returnparam = FALSE){
    out <- c(calcdis_inside(matx,"m0rI"),
             calcdis_inside(MCARbyrow(matx,0.1),"m1rI"),
             calcdis_inside(MCARbyrow(matx,0.2),"m2rI"),
             calcdis_inside(MCARbyrow(matx,0.3),"m3rI"),
             calcdis_inside(MCARbyrow(matx,0.4),"m4rI"),
             calcdis_inside(MCARbyrow(matx,0.5),"m5rI"),
             calcdis_inside(round(matx,4),"m0r4"),
             calcdis_inside(round(matx,3),"m0r3"),
             calcdis_inside(round(matx,2),"m0r2"),
             calcdis_inside(round(matx,1),"m0r1"),
             calcdis_inside(round(matx,0),"m0r0")
    )




    # # missingness 0.1
    # tempbray.suc <- calc.bray.suc(matx)
    # tempbray.amm <- calc.bray.amm(matx)
    # out <- c(out,
    #          sucm1rI = tempbray.suc[1],        # Bray-Curtis dissimilarity: replacement
    #          sucm1rI = tempbray.suc[2],        # Bray-Curtis dissimilarity: nestedness
    #          sucm1rI = tempbray.suc[3],        # Bray-Curtis dissimilarity: Full index
    #          ammm1rI = tempbray.amm[1],        # Bray-Curtis dissimilarity: replacement
    #          ammm1rI = tempbray.amm[2],        # Bray-Curtis dissimilarity: nestedness
    #          ammm1rI = tempbray.amm[3]        # Bray-Curtis dissimilarity: Full index
    # )
    # # missingness 0.2
    # tempmatx <- delete_MCAR(matx, 0., 1)
    # tempmatx[!!rowSums(is.na(tempmatx)),] <- NA
    # tempbray.suc <- calc.bray.suc(matx)
    # tempbray.amm <- calc.bray.amm(matx)
    # out <- c(out,
    #          sucm1rI = tempbray.suc[1],        # Bray-Curtis dissimilarity: replacement
    #          sucm1rI = tempbray.suc[2],        # Bray-Curtis dissimilarity: nestedness
    #          sucm1rI = tempbray.suc[3],        # Bray-Curtis dissimilarity: Full index
    #          ammm1rI = tempbray.amm[1],        # Bray-Curtis dissimilarity: replacement
    #          ammm1rI = tempbray.amm[2],        # Bray-Curtis dissimilarity: nestedness
    #          ammm1rI = tempbray.amm[3]        # Bray-Curtis dissimilarity: Full index
    # )


  # this "returnparam" is a debugging option.
  # If TRUE, returns the mean, sd, autocorrelation, and correlation of the input dataset
  if (returnparam){
    tmpacf<-acf(matx,plot=FALSE)
    out <- c(out,
             mean = colMeans(matx),
             sd = apply(matx,2,sd),
             ar = sapply(seq(ncol(matx)),function(x) tmpacf[x,x]$acf),
             cor = (cor(matx)[upper.tri(cor(matx))])
    )
  }
  out
}

# data generation by resampling the Lorenz System
  genswitch <- function(nrep=10, prob=0.5, addmean = 50, flipy = FALSE, outcol = 3, xyonly = FALSE){
    lormat <- rootlorenz
    lormat <- cbind(lormat[,1]+addmean,lormat[,2]+addmean,lormat[,3]-mean(lormat[,3])+addmean)

    # group points into two wings by the symmetrical x-axis
    lormat.c1 <- lormat[lormat[,1]>=addmean,]
    lormat.c2 <- lormat[lormat[,1]<addmean,]


    # randomly determine which observation is there going to be a switching
    u <- runif(nrep, min = 0, max = 1)
    u <- (u <= prob)

    # serially resample n observations
    switchoutput = NULL
    switchoutput = c(lormat.c1[sample(nrow(lormat.c1),1),],
                     lormat.c1[sample(nrow(lormat.c1),1),],
                     lormat.c1[sample(nrow(lormat.c1),1),])
    currentpool <- 1 + round(runif(1),0)
    for (i in 2:length(u)){
      if(u[i]){
        if(currentpool ==1){
          currentpool <- 2
        }else{
          currentpool <-1
        }
      }
      if(currentpool==1){
        temprow <- c(lormat.c1[sample(nrow(lormat.c1),1),],
                     lormat.c1[sample(nrow(lormat.c1),1),],
                     lormat.c1[sample(nrow(lormat.c1),1),])
      }else{
        temprow <- c(lormat.c2[sample(nrow(lormat.c2),1),],
                     lormat.c2[sample(nrow(lormat.c2),1),],
                     lormat.c2[sample(nrow(lormat.c2),1),])
      }
      switchoutput <- rbind(switchoutput,temprow)
    }

    # flip the y-axis (at y=addmean) so that the grand mean of each wing will be the same
    if (flipy){
      switchoutput[,2] <- -switchoutput[,2]+addmean*2
      switchoutput[,5] <- -switchoutput[,5]+addmean*2
      switchoutput[,8] <- -switchoutput[,8]+addmean*2
    }
    zoutput <- switchoutput[,c(3,6,9)]
    xyoutput <- switchoutput[,-c(3,6,9)]
    # return
    if (xyonly) {
      xyoutput[,1:outcol]
    }else{
      switchoutput[,1:outcol]
    }
  }

# simulation with different probability of switching
genresult <- function( nobs = 100, simrep = 1, ERn = 6, xyonly = FALSE, seed = 1999){
  set.seed(seed)
  twooutput <- NULL
  for (i in 0:4){
    prob <- 0.1+i*0.2
    # (1) generate data, (2) calculate dissimilarity, (3) repeat for simrep instances
    tempoutput <- t(replicate(simrep,calcdis(genswitch(nrep = nobs,
                                                       prob = prob,
                                                       flipy=TRUE,
                                                       outcol = ERn,
                                                       xyonly = xyonly),returnparam = FALSE)))
    tempoutput <- cbind(prob=rep(prob,nrow(tempoutput)),
                        ERn=rep(ERn,nrow(tempoutput)),
                        nobs=rep(nobs,nrow(tempoutput)),
                        tempoutput)
    # rbind results generated with different probability of switching
    twooutput <- rbind(twooutput,tempoutput)
  }
  twooutput
}

# a wrapper for serial simulation (i.e., not parallel)
  lsimwrapper <- function (siminput){
    output <- NULL
    for (i in 1:nrow(siminput)){
      tmpres <- genresult(simrep = siminput$rep[i],
                          nobs = siminput$n[i],
                          ERn = siminput$ERn[i])
      tmpres <- as.data.frame(tmpres)
      output <- rbindlist(list(output,tmpres), fill = TRUE)
    }
    output
  }

#======================
# Simulation: Data generation
#======================
# prepare parallel processing
rootlorenz <- lorenz(do.plot = FALSE)
rootlorenz <- as.matrix(cbind(rootlorenz$x, rootlorenz$y, rootlorenz$z))
nclust <- parallel::detectCores()
cl <- makeCluster(nclust)
registerDoSNOW(cl)

# run simulation
time_start <- Sys.time()
tab <- foreach(rownum = 1:nrow(siminput), .packages = c( "betapart", "vegan","missMethods"), .combine = rbind) %dopar% {
  # Set seed
  suppressMessages(attach(siminput[rownum, ]))
  set.seed(seed)
  # Simulate data
  df <- genresult(simrep = rep,
                  nobs = n,
                  ERn = ERn)

  df
}
end_time <- Sys.time()
#Close cluster
stopCluster(cl)

time_per_row <- as.numeric(end_time - time_start) / nrow(siminput)
dfcombine <- as.data.frame(tab)

# serial processing - might be faster when the number of repetition is low
# out.tmp <- lsimwrapper(siminput = siminput)
# dfcombine <- as.data.frame(out.tmp)

#======================
# Simulation: Evaluate indices' performance
#======================

respcor <- NULL
for (i in 1:length(list_metrics)){
  tmpres<- pcor(dfcombine[dfcombine$prob>0.0,c("prob","ERn","nobs",colnames(dfcombine)[i+3])])$estimate[,4]
  tmpres[4] <- list_metrics[i]
  respcor<- rbind(respcor,tmpres)
}
respcor <- as.data.frame(respcor)
rescor <- NULL
for (i in 1:length(list_metrics)){
  tmpres<- cor(dfcombine[,c("prob","ERn","nobs",colnames(dfcombine)[i+3])])[,4]
  tmpres[4] <- list_metrics[i]
  rescor<- rbind(rescor,tmpres)
}
rescor <- as.data.frame(rescor)

rownames(respcor) <-  1:length(list_metrics)
rownames(rescor) <-  1:length(list_metrics)
respcor
rescor
# write results as .csv files
write.csv(respcor,paste0("sim2_measurement_respcor_rep",siminput$rep[1]," ",Sys.Date(),".csv")) # partial correlation
write.csv(rescor,paste0("sim2_measurement_rescor_rep",siminput$rep[1]," ",Sys.Date(),".csv")) #  correlation
