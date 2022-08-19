library(ppcor) # partial correlation

n_singleoutput_metric <- 5
list_metrics <- c("multibray",
                  "mbray_gra",
                  "mbray_bal",
                  "betweenSD",
                  "withinSD",
                  # all person level measures go above
                  # remember to amend n_singleoutput_metric
                  "SD",
                  "euclidean",
                  "manhattan",
                  "chord",
                  "logchord",
                  "chisq",
                  "hellinger",
                  "jaccard",
                  "kulczynski",
                  "bray",
                  "bray(beta)",
                  "bray.bal",
                  "bray.gra",
                  "KLdiv"
)


print_result <- function(timemode="successive"){
  if (n_singleoutput_metric > 0){
    temptab <- (tab[,-c(1:n_singleoutput_metric)])
  }else{
    temptab <- tab
  }

# 0 = suc.second, 1 = mom.all, 2 = mom.second
if (timemode == "composite"){
  temptab <- (temptab[, seq_len(ncol(temptab)) %% 3 == 0] + temptab[, seq_len(ncol(temptab)) %% 3 == 2])/2
} else if (timemode == "successive"){
  temptab <- temptab[, seq_len(ncol(temptab)) %% 3 == 0]
} else { #"allmoment"
  temptab <- temptab[, seq_len(ncol(temptab)) %% 3 == 1]
}
temptab <- cbind((tab[,  (n_singleoutput_metric>0):n_singleoutput_metric]),temptab)
colnames(temptab) <- list_metrics
output <-  cbind(siminput,temptab)
output[["seed"]] <- NULL
output[["ER_withinSD"]] <- NULL # replace by "adjSD" if needed

res_pcor <- data.frame()
res_cor <- data.frame()
for (i in 1:length(list_metrics)){
  dftemp <- cbind(output[,list_metrics[i]],output[,1:(ncol(output)-length(list_metrics))])
  dftemp <- cbind("dissimilarity" = dftemp[,1],Filter(function(x)(length(unique(x))>1), dftemp[2:ncol(dftemp)]))
  cortemp <- cor(dftemp)[1,]
  tryCatch({
    pcortemp <- pcor(dftemp)$estimate[1,]
  }, error = function(e){
    # assign NA to pcor results as pcor throws error upon NA input
    pcortemp <<- cortemp
  })

  if(i==1){
    res_pcor <- data.frame(as.list(pcortemp))
    res_cor <- data.frame(as.list(cortemp))
  }else{
    res_pcor <- rbind.data.frame(res_pcor,pcortemp)
    res_cor <- rbind.data.frame(res_cor,cortemp)
  }
  res_pcor[i,1] <- list_metrics[i]
  res_cor[i,1] <- list_metrics[i]
}

list(partial_correlation=res_pcor,correlation=res_cor)
}

# allmoment, successive, or composite
print_result("allmoment")
