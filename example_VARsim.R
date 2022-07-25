library("pbapply") # percentage bar for replicate function
library(tsDyn) # package for VAR.sim
withinSD <- 1
crosslag <- 0 # assumes zero cross-lag relationships for now
obs <- 100 # not sure why obs = 2 is not working
nrep <- 500
nvar <- 2 # number of variables (ER strategies) to generate
autocorr <- 0
correlation <- -0.7

CV <- diag(withinSD, nvar,nvar)
CV[CV == 0] <- correlation
B1 <- diag(autocorr, nvar,nvar)
B1[B1 == 0] <- crosslag


testnstgy <- function(){
  var1 <- VAR.sim(B=B1, n=obs, include="none",varcov=CV)
  # # plot scatterplots with x=y line to see how observations are distributed
  # plot(var1)
  # abline(a=0, b=1)
  resact <- acf(var1, plot = FALSE)
  tmp <- list(resact$acf, # AR estimation
              cor(var1), #correlation matrix
              apply(var1,2,sd), # within-variable SD
              sum(var1[,1]>var1[,2])/obs, # % S1>S2
              abs(diff(as.numeric(var1[,1]>var1[,2])))/obs # % rank change over obs
              )
  return(tmp)
}

testres <- pbreplicate(nrep,testnstgy())
acfmean <- apply(simplify2array(testres[1,]),1:3,mean)
acfsd <- apply(simplify2array(testres[1,]),1:3,sd)
extractAR1 <- function(matinput){
  outlist<- NULL
  if (nrow(matinput) >1){
  for (i in 1:nvar){
      outlist<- append(outlist, matinput[2,i,i])
  }
  }else{
    # when obs is too low acf returns only 1 row of estimation
  outlist<- append(outlist, "obs is too low, no acf AR(1) estimates")
  }
  return(outlist)
}

funoutput <- function(){
#  AR(1) mean and SD
print(paste0("AR(1) mean and SD; input autocorrelation = ",autocorr))
print(extractAR1(acfmean))
print(extractAR1(acfsd))
# Correlation matrix: Mean
print(paste0("Correlation matrix: Mean & SD; input correlation = ",correlation))
print(apply(simplify2array(testres[2,]), 1:2, mean))
# Correlation matrix: SD
print(apply(simplify2array(testres[2,]), 1:2, sd))
# mean within-variable ("within-strategy") SD
print("Within-variable  SD")
print(apply((do.call("rbind",testres[3,])),2,mean))
# mean % of S1>S2
paste0("% Var1>Var2: [mean] ",mean(unlist(testres[4,])))
# mean  % of rank change
paste0("% Rank Change: [mean] ",mean(unlist(testres[5,])), " [SD] ",sd(unlist(testres[5,])))
}
funoutput()
