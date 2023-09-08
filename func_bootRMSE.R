
# Creating Function to obtain R-Squared from the data
boot.rmseWithinRSD <- function(data, indices) {
  val <- data[indices,] # selecting sample with boot
  fit <- lme(fixed=moment_meanNA ~ moment_withinRSDcb + timecw,
             data=val,
             random=~1 | ppnr, correlation = corAR1(),
             control =lmeControl(msMaxIter = 1000, msMaxEval = 1000,opt='optim'),na.action = na.omit)
  return(performance_rmse(fit, normalized = FALSE))
}

boot.rmseBetweenRSD <- function(data, indices) {
  val <- data[indices,] # selecting sample with boot
  fit <- lme(fixed=moment_meanNA ~ moment_betweenRSD.singlecw+moment_betweenRSD.singlecb + timecw,
             data=val,
             random=~1+ moment_betweenRSD.singlecw | ppnr, correlation = corAR1(),
             control =lmeControl(msMaxIter = 1000, msMaxEval = 1000,opt='optim'),na.action = na.omit)
  return(performance_rmse(fit, normalized = FALSE))
}
boot.rmseBetweenRSD.suc <- function(data, indices) {
  val <- data[indices,] # selecting sample with boot
  fit <- lme(fixed=moment_meanNA ~ moment_betweenRSD.succw +moment_betweenRSD.succb + timecw,
             data=val,
             random=~1+ moment_betweenRSD.succw | ppnr, correlation = corAR1(),
             control =lmeControl(msMaxIter = 1000, msMaxEval = 1000,opt='optim'),na.action = na.omit)
  return(performance_rmse(fit, normalized = FALSE))
}


boot.rmseDBC <- function(data, indices) {
  val <- data[indices,] # selecting sample with boot
  fit <- lme(fixed=moment_meanNA ~ moment_bray.all.succw+moment_bray.all.succb+ timecw,
             data=val,
             random=~1+ moment_bray.all.succw | ppnr, correlation = corAR1(),
             control =lmeControl(msMaxIter = 1000, msMaxEval = 1000,opt='optim'),na.action = na.omit)
  return(performance_rmse(fit, normalized = FALSE))
}

boot.rmseDBCsub <- function(data, indices) {
  val <- data[indices,] # selecting sample with boot
  fit <- lme(moment_meanNA ~ moment_bray.bal.succw+ moment_bray.gra.succw+
               moment_bray.bal.succb+ moment_bray.gra.succb+ timecw,
             data=val,
             random=~1+ moment_bray.bal.succw+ moment_bray.gra.succw | ppnr, correlation = corAR1(),
             control =lmeControl(msMaxIter = 1000, msMaxEval = 1000,opt='optim'),na.action = na.omit)
  return(performance_rmse(fit, normalized = FALSE))
}


# Performing 1000 replications with boot
bootRMSEresults <- function (bootrep = 5, saveRdata = FALSE){
outputrmseWithinRSD.1 <- boot(data=dfERV1[dfERV1$b_completeER,], statistic=boot.rmseWithinRSD,
                               R=bootrep, parallel = "snow")
outputrmseWithinRSD.2 <- boot(data=dfERV2[dfERV2$b_completeER,], statistic=boot.rmseWithinRSD,
                               R=bootrep, parallel = "snow")
outputrmseWithinRSD.3 <- boot(data=dfERV3[dfERV3$b_completeER,], statistic=boot.rmseWithinRSD,
                               R=bootrep, parallel = "snow")
outputrmseBetweenRSD.1 <- boot(data=dfERV1[dfERV1$b_completeER,], statistic=boot.rmseBetweenRSD,
                               R=bootrep, parallel = "snow")
outputrmseBetweenRSD.2 <- boot(data=dfERV2[dfERV2$b_completeER,], statistic=boot.rmseBetweenRSD,
                               R=bootrep, parallel = "snow")
outputrmseBetweenRSD.3 <- boot(data=dfERV3[dfERV3$b_completeER,], statistic=boot.rmseBetweenRSD,
                               R=bootrep, parallel = "snow")
outputrmseBetweenRSD.suc.1 <- boot(data=dfERV1[dfERV1$b_completeER,], statistic=boot.rmseBetweenRSD.suc,
                               R=bootrep, parallel = "snow")
outputrmseBetweenRSD.suc.2 <- boot(data=dfERV2[dfERV2$b_completeER,], statistic=boot.rmseBetweenRSD.suc,
                               R=bootrep, parallel = "snow")
outputrmseBetweenRSD.suc.3 <- boot(data=dfERV3[dfERV3$b_completeER,], statistic=boot.rmseBetweenRSD.suc,
                               R=bootrep, parallel = "snow")
outputrmseDBC.1 <- boot(data=dfERV1[dfERV1$b_completeER,], statistic=boot.rmseDBC,
                        R=bootrep, parallel = "snow")
outputrmseDBC.2 <- boot(data=dfERV2[dfERV2$b_completeER,], statistic=boot.rmseDBC,
                        R=bootrep, parallel = "snow")
outputrmseDBC.3 <- boot(data=dfERV3[dfERV3$b_completeER,], statistic=boot.rmseDBC,
                        R=bootrep, parallel = "snow")
outputrmseDBCsub.1 <- boot(data=dfERV1[dfERV1$b_completeER,], statistic=boot.rmseDBCsub,
                        R=bootrep, parallel = "snow")
outputrmseDBCsub.2 <- boot(data=dfERV2[dfERV2$b_completeER,], statistic=boot.rmseDBCsub,
                        R=bootrep, parallel = "snow")
outputrmseDBCsub.3 <- boot(data=dfERV3[dfERV3$b_completeER,], statistic=boot.rmseDBCsub,
                        R=bootrep, parallel = "snow")
# table of mean RMSEs

resbootRMSE <- matrix(c(mean(outputrmseWithinRSD.1$t),
         mean(outputrmseWithinRSD.2$t),
         mean(outputrmseWithinRSD.3$t),
         mean(outputrmseBetweenRSD.1$t),
         mean(outputrmseBetweenRSD.2$t),
         mean(outputrmseBetweenRSD.3$t),
         mean(outputrmseBetweenRSD.suc.1$t),
         mean(outputrmseBetweenRSD.suc.2$t),
         mean(outputrmseBetweenRSD.suc.3$t),
         mean(outputrmseDBC.1$t),
         mean(outputrmseDBC.2$t),
         mean(outputrmseDBC.3$t),
         mean(outputrmseDBCsub.1$t),
         mean(outputrmseDBCsub.2$t),
         mean(outputrmseDBCsub.3$t)),5, byrow = TRUE)

if(saveRdata){
saveRDS(outputrmseWithinRSD.1, paste0("boot",bootrep,"_withinRSD1.RData"))
saveRDS(outputrmseWithinRSD.2, paste0("boot",bootrep,"_withinRSD2.RData"))
saveRDS(outputrmseWithinRSD.3, paste0("boot",bootrep,"_withinRSD3.RData"))
saveRDS(outputrmseBetweenRSD.1, paste0("boot",bootrep,"_btwRSD1.RData"))
saveRDS(outputrmseBetweenRSD.2, paste0("boot",bootrep,"_btwRSD2.RData"))
saveRDS(outputrmseBetweenRSD.3, paste0("boot",bootrep,"_btwRSD3.RData"))
saveRDS(outputrmseBetweenRSD.suc.1, paste0("boot",bootrep,"_btwsucRSD1.RData"))
saveRDS(outputrmseBetweenRSD.suc.2, paste0("boot",bootrep,"_btwsucRSD2.RData"))
saveRDS(outputrmseBetweenRSD.suc.3, paste0("boot",bootrep,"_btwsucRSD3.RData"))
saveRDS(outputrmseDBC.1, paste0("boot",bootrep,"_dBC1.RData"))
saveRDS(outputrmseDBC.2, paste0("boot",bootrep,"_dBC2.RData"))
saveRDS(outputrmseDBC.3, paste0("boot",bootrep,"_dBC3.RData"))
saveRDS(outputrmseDBCsub.1, paste0("boot",bootrep,"_dBCsub1.RData"))
saveRDS(outputrmseDBCsub.2, paste0("boot",bootrep,"_dBCsub2.RData"))
saveRDS(outputrmseDBCsub.3, paste0("boot",bootrep,"_dBCsub3.RData"))
}

colnames(resbootRMSE) <- c(paste0("dataset",1:3))
rownames(resbootRMSE) <- c("within RSD",
                           "between RSD",
                           "between RSD suc dif",
                           "Bray-Curtis (full index)",
                           "Bray-Curtis (subcomponents)")
round(resbootRMSE,3)
}

