metric_mssd <- function(x){
  mean(rowSums((x[-1, ] - x[-nrow(x),])^2))
}

metric_mean_euclidean <- function(x){
  mean(sqrt(rowSums((x[-1, ] - x[-nrow(x),])^2)))
}
