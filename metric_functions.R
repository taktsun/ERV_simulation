#======================
# common functions
#======================

dis_suc_vector <- function(distmat){
  return (c(0,(distmat[row(distmat) == col(distmat) + 1])))
}

# person-level mean dissimilarity across all time points
# 1. the *nrow(x)/(nrow(x)-1) adjustment is needed because there is one 0 in each row/column
# 2. mom is a vector of one-to-all momentary comparisons
# 3. suc is a vector of successive dissimilarity
# 4. only 2nd observation onwards has has successive dissimilarity thus 2:n
composite_mean <- function (matx, successive){
  mom <- apply(matx,1,mean)*nrow(matx)/(nrow(matx)-1)
  suc <- dis_suc_vector(matx)
  com <- (mom + suc*successive)/(1+successive)
  mean(com[(1+successive):nrow(matx)])
}

#======================
# metric functions
#======================


metric_mssd <- function(x){
  tryCatch({
  if(!is.matrix(x)) stop()
    mean(rowSums((x[-1, ] - x[-nrow(x),])^2))
  }, error = function(e){
    mean((x[1:(length(x)-1)] - x[2:(length(x))])^2)
  })
}
metric_mean_euclidean <- function(x){
  tryCatch({
  if(!is.matrix(x)) stop()
    mean(sqrt(rowSums((x[-1, ] - x[-nrow(x),])^2)))
  }, error = function(e){
    mean(sqrt((x[1:(length(x)-1)] - x[2:(length(x))])^2))
  })
}

# between-strategy SD
metric_person_between_SD <- function(x){
  tryCatch({
    if(!is.matrix(x)) stop()
    mean(apply(x,1,sd))
  }, error = function(e){
    NA
  })
}
# within-strategy SD
metric_person_within_SD <- function(x){
  tryCatch({
    if(!is.matrix(x)) stop()
    mean(apply(x,2,sd))
  }, error = function(e){
    sd(x)
  })
}

# composite person mean dissimilarity by between-strategy SD
metric_person_SD <- function(x, successive = TRUE){
  tryCatch({
  if(!is.matrix(x)) stop()
  matx <- abs(outer(apply(x,1,sd),apply(x,1,sd), '-'))
  mom <- apply(matx,1,mean)*nrow(matx)/(nrow(matx)-1)
  suc <- dis_suc_vector(matx)
  com <- (mom + suc*successive)/(1+successive)
  mean(com[(1+successive):nrow(matx)])
  }, error = function(e){
    NA
  })
}

# composite person mean dissimilarity with various methods (see vegdist)
metric_person_vegan <- function(x, method, successive = TRUE, deco = ""){
  if (deco !="") {
    x <- decostand(x,deco)
  }
  matx <- as.matrix(vegdist(x,method = method))
  composite_mean(matx,successive)
}

# composite person mean KL divergence
metric_person_KLdiv <- function(x, successive = TRUE){
  x <- decostand(x,"total")
  tempdist <- c()
  for (k in 1:(nrow(x)-1)){
    for (j in (k+1):nrow(x)){
      tempdist <- append(tempdist,(KL.plugin(x[k,],x[j,])))
    }
  }
  tempdist <- as.matrix(tempdist)
  mat.KLdiv <- matrix(0, nrow = nrow(x), ncol = nrow(x))
  mat.KLdiv[lower.tri(mat.KLdiv, diag = FALSE)] <- tempdist
  mat.KLdiv[upper.tri(mat.KLdiv)] <- t(mat.KLdiv)[upper.tri(mat.KLdiv)]
  composite_mean(mat.KLdiv,successive)
}
