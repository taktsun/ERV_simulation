suppressMessages(library('nonlinearTseries'))
library("betapart")
library("vegan")
library(data.table)
setDTthreads(threads = 1) # so it won't take up the parallel processing for running simulation
library(ppcor) # partial correlation

source("metric_functions.R")
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

set.seed(1999)
siminput <- expand.grid(
  # reps
  rep = 1000,
  # nobs of time series
  n = c(30,70,100),
  # Number of ER strategies
  ERn = c(3,6,9)
)
siminput$seed <- sample(1:.Machine$integer.max, nrow(siminput))
list_metrics <- c("sd.ws",
                  "sd.bs",
                  "sd.mom.all",
                  "sd.suc.second",
                  "bray.mom.all",
                  "bray.suc.second",
                  "dBCs.mom.all",
                  "dBCs.suc.second",
                  "dBCe.mom.all",
                  "dBCe.suc.second",
                  "jaccard.mom.all",
                  "jaccard.suc.second",
                  "djs.mom.all",
                  "djs.suc.second",
                  "dje.mom.all",
                  "dje.suc.second",
                  "chord.mom.all",
                  "chord.suc.second",
                  "chisq.mom.all",
                  "chisq.suc.second"
)

calcdis <- function(matx, returnparam = FALSE){

  out <- c(sd.ws = metric_person_within_SD(matx),
           sd.bs = metric_person_between_SD(matx),
           sd = metric_person_SD(matx),
           bray = metric_person_beta(matx, "bray"),
           dBCs = metric_person_beta(matx, "bray",".bal"),
           dBCe = metric_person_beta(matx, "bray",".gra"),
           jaccard = metric_person_beta(matx, "ruz"),
           djs = metric_person_beta(matx, "ruz",".bal"),
           dje = metric_person_beta(matx, "ruz",".gra"),
           chord = metric_person_vegan(matx, "chord"),
           chisq = metric_person_vegan(matx, "chisq")
  )
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


genswitch <- function(nrep=10, prob=0.5, addmean = 50, flipy = FALSE, outcol = 3, xyonly = FALSE){
  lormat <- rootlorenz
  lormat <- cbind(lormat[,1]+addmean,lormat[,2]+addmean,lormat[,3]-mean(lormat[,3])+addmean)
  lormat.c1 <- lormat[lormat[,1]>=addmean,]
  lormat.c2 <- lormat[lormat[,1]<addmean,]

  switchoutput = NULL
  u <- runif(nrep, min = 0, max = 1)
  u <- (u <= prob)
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

  if (flipy){
    switchoutput[,2] <- -switchoutput[,2]+addmean*2
    switchoutput[,5] <- -switchoutput[,5]+addmean*2
    switchoutput[,8] <- -switchoutput[,8]+addmean*2
  }
  zoutput <- switchoutput[,c(3,6,9)]
  xyoutput <- switchoutput[,-c(3,6,9)]
  if (xyonly) {
    xyoutput[,1:outcol]
    # switchoutput <- switchoutput[,-9]
    # switchoutput <- switchoutput[,-6]
    # switchoutput <- switchoutput[,-3]
  }else{
    switchoutput[,1:outcol]
  }
}

genresult <- function( nobs = 100, simrep = 1, ERn = 6, xyonly = FALSE, seed = 1999){
  set.seed(seed)
  twooutput <- NULL
  for (i in 0:4){ #1:100
    prob <- 0.1+i*0.2
    tempoutput <- t(replicate(simrep,calcdis(genswitch(nrep = nobs,
                                                       prob = prob,
                                                       flipy=TRUE,
                                                       outcol = ERn,
                                                       xyonly = xyonly),returnparam = FALSE)))
    tempoutput <- cbind(prob=rep(prob,nrow(tempoutput)),
                        ERn=rep(ERn,nrow(tempoutput)),
                        nobs=rep(nobs,nrow(tempoutput)),
                        tempoutput)
    twooutput <- rbind(twooutput,tempoutput)
  }
  twooutput
}


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

# prepare parallel processing
rootlorenz <- lorenz() # lorenz(time = seq(0, 50, by = 0.02))
rootlorenz <- as.matrix(cbind(rootlorenz$x, rootlorenz$y, rootlorenz$z))
library(doSNOW)
library(parallel)
nclust <- parallel::detectCores()
cl <- makeCluster(nclust)
registerDoSNOW(cl)

# run simulation
time_start <- Sys.time()
# EL: not sure why adding "paths = paths," bring error in my machine. Temporarily remove
tab <- foreach(rownum = 1:nrow(siminput), .packages = c( "betapart", "vegan"), .combine = rbind) %dopar% {
  # Set seed
  suppressMessages(attach(siminput[rownum, ]))
  set.seed(seed)
  # Simulate data
  df <- genresult(simrep = rep,
                  nobs = n,
                  ERn = ERn)

  # CJ: If the sim takes very long or uses a lot of memory, print to text files.
  # CJ: If not, just return the output.
  # CJ: We have to see which works best when we start running a larger chunk!
  # write.table(x = t(out), file = sprintf("results_%d.txt" , Sys.getpid()), sep = "\t", append = TRUE, row.names = FALSE, col.names = FALSE)
  # NULL
  df
}
end_time <- Sys.time()
#Close cluster
stopCluster(cl)

time_per_row <- as.numeric(end_time - time_start) / nrow(siminput)
dfcombine <- as.data.frame(tab)


# serial processing
# out.tmp <- lsimwrapper(siminput = siminput)
# dfcombine <- as.data.frame(out.tmp)

respcor <- NULL
for (i in 1:length(list_metrics)){
  tmpres<- pcor(dfcombine[dfcombine$prob>0.0,c("prob","ERn","nobs",list_metrics[i])])$estimate[,4]
  tmpres[4] <- list_metrics[i]
  respcor<- rbind(respcor,tmpres)
}
respcor <- as.data.frame(respcor)
rescor <- NULL
for (i in 1:length(list_metrics)){
  tmpres<- cor(dfcombine[,c("prob","ERn","nobs",list_metrics[i])])[,4]
  tmpres[4] <- list_metrics[i]
  rescor<- rbind(rescor,tmpres)
}
rescor <- as.data.frame(rescor)

rownames(respcor) <-  list_metrics
rownames(rescor) <-  list_metrics
write.csv(respcor,paste0("sim2_respcor_rep",siminput$rep[1]," ",Sys.Date(),".csv"))
write.csv(rescor,paste0("sim2_rescor_rep ",siminput$rep[1]," ",Sys.Date(),".csv"))
