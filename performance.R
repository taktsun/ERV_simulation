library(ppcor) # partial correlation

list_metrics <- c("mssd",
                  "suc_euclidean",
                  "betweenSD",
                  "withinSD",
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
                  "KLdiv",
                  "multibray"
)
colnames(tab) <- list_metrics
output <-  cbind(siminput,tab)
output[["seed"]] <- NULL
output[["ER_withinSD"]] <- NULL

res_pcor <- data.frame()
res_cor <- data.frame()
for (i in 1:length(list_metrics)){
  dftemp <- cbind(output[,list_metrics[i]],output[,1:(ncol(output)-length(list_metrics))])
  dftemp <- Filter(function(x)(length(unique(x))>1), dftemp)
  pcortemp <- pcor(dftemp)$estimate[1,]
  cortemp <- cor(dftemp)[1,]

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

res_pcor
res_cor
