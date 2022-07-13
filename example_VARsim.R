library(tsDyn) # package for VAR.sim
library("pbapply") # percentage bar for replicate function

autocorr <- 0.2
correlation <- 0.5
crosslag <- 0 # assumes zero cross-lag relationships for now
nvar <- 4 # number of variables (ER strategies) to generate

# 2-strategy testing

# test2stgy <- function(){
# var1 <- VAR.sim(B=B1, n=100, include="none",varcov=CV)
# resact <- acf(var1, plot = FALSE)
# tmp <- list(resact$acf[2,1,1],resact$acf[2,2,2],cor(var1)[2,1])
# return(tmp)
# }
#
# CV<-matrix(c(1,correlation,correlation,1),2)
# B1<-matrix(c(autocorr, 0.0, 0.0, autocorr), 2)
# testres<-pbreplicate(500,test2stgy())
# #mean of autocorrelation of 1st variable, autocorrelation of 2nd variable, and correlation between 2
# apply((apply(testres,1,unlist)),2,mean)
# #SD of autocorrelation of 1st variable, autocorrelation of 2nd variable, and correlation between 2
# apply((apply(testres,1,unlist)),2,sd)
#

# n-strategy testing


CV <- diag(1, nvar,nvar)
CV[CV == 0] <- correlation

B1 <- diag(autocorr, nvar,nvar)
B1[B1 == 0] <- crosslag

testnstgy <- function(){
  var1 <- VAR.sim(B=B1, n=100, include="none",varcov=CV)
  resact <- acf(var1, plot = FALSE)
  tmp <- list(resact$acf,cor(var1))
  return(tmp)
}

testres <- pbreplicate(500,testnstgy())
acfmean <- apply(simplify2array(testres[1,]),1:3,mean)
acfsd <- apply(simplify2array(testres[1,]),1:3,sd)
extractAR1 <- function(matinput){
  outlist<- NULL
  for (i in 1:nvar){
    outlist<- append(outlist, matinput[2,i,i])
  }
  return(outlist)
}
#  AR(1) mean and SD
extractAR1(acfmean)
extractAR1(acfsd)
# mean correlation matrix
apply(simplify2array(testres[2,]), 1:2, mean)
apply(simplify2array(testres[2,]), 1:2, sd)

