# remotes::install_github("wviechtb/esmpack")
# remove.packages(c("StanHeaders", "rstan"))
# remove.packages(c("mixedup"))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

library(betapart)
library(vegan)
library(crqa)
library(esmpack)
library(misty) # for calculating ICC in descriptive statistics
library(nlme) # MLM estimation
library(brms) # MLM bayesian method
library(rstan) # required for brms
options(mc.cores = parallel::detectCores()-1)

# ===== load functions from other .R ================

# relative SD function
# from: https://ppw.kuleuven.be/okp/software/relative_variability/
source("RSD/checkInput.R")
source("RSD/maximumVAR.R")
source("RSD/relativeSD.R")
source("RSD/checkOutput.R")

source("metric_functions.R")
source("reanalysis_calculateERV.R")
source("reanalysis_desStat.R")
source("reanalysis_MLM.R")

source
# =====  ================

dfERV1<-loadESMdata(1)
dfERV2<-loadESMdata(2)
dfERV3<-loadESMdata(3)
#warnings are from calculating relative SD: all-zero ratings return NaN due to division of zero

desstat1 <- summarydesstat(dfERV1)
desstat2 <- summarydesstat(dfERV2)
desstat3 <- summarydesstat(dfERV3)

resMLM1 <- MLMresults(dfERV1,1)
resMLM2 <- MLMresults(dfERV2,2)
resMLM3 <- MLMresults(dfERV3,3)
desstat3
resMLM_all<-rbind(resMLM1,resMLM2,resMLM3)

# output .csv if needed:
# write.csv(resMLM1,"desstat1.csv")
# write.csv(resMLM2,"desstat2.csv")
# write.csv(resMLM3,"desstat3.csv")
# write.csv(resMLM_all,"reanalysis_results.csv")
