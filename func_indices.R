#==================================================================================
# functions to be used in calculating different variability indices in simulation
#==================================================================================

# return a vector of dissimilarity between t and t-1 in all observations (except t=1, where 0 is assigned)
dis_suc_vector <- function(distmat){
  return (c(0,(distmat[row(distmat) == col(distmat) + 1])))
}

# return all-moment comparison and successive difference results per simulated dataset
composite_mean <- function (matx){
  # all-moment comparison approach
    #  *nrow(x)/(nrow(x)-1) to adjust for the one occasion of 0 dissimilarity
    # e.g., in all-moment comparison, moment 1 is compared with moment 1,2,3,...n
    #       moment 1 vs moment 1 will give 0 dissimilarity
  mom <- apply(matx,1,mean)*nrow(matx)/(nrow(matx)-1)
  # successive difference approach
  suc <- dis_suc_vector(matx)
  # only 2nd observation onwards has successive dissimilarity thus 2:n
  c(mom.all= mean(mom), mom.second = mean(mom[2:nrow(matx)]), suc.second = mean(suc[2:nrow(matx)]))
}
#==================================================================================
# functions that calculate specific variability indices
#==================================================================================

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

# Mimicked Temporal Comparisons with between-strategy SD: all-moment comparison & successive difference approach
metric_person_SD <- function(x){
  tryCatch({
  if(!is.matrix(x)) stop()
  matx <- abs(outer(apply(x,1,sd),apply(x,1,sd), '-'))
  mom <- apply(matx,1,mean)*nrow(matx)/(nrow(matx)-1)
  suc <- dis_suc_vector(matx)
  c(mom.all= mean(mom), mom.second = mean(mom[2:nrow(matx)]), suc.second = mean(suc[2:nrow(matx)]))
  }, error = function(e){
    NA
  })
}

# return person mean dissimilarity with various comparison approaches
# vegdist is a function to calculate dissimilarity in the vegan package
# for chi-squared and chord distance
metric_person_vegan <- function(x, method){
  matx <- as.matrix(vegdist(x,method = method))
  composite_mean(matx)
}
# beta.pair.abund is a function to calculate bray-curtis (or jaccard) dissimilarity
# in full indices or in subcomponents
metric_person_beta <- function(x, index.family = "bray" , extract= ""){
  matx <- (beta.pair.abund(x, index.family = index.family))
  matx <- as.matrix(matx[[paste0("beta.",index.family,extract)]])
  composite_mean(matx)
}

#==================================================================================
# The below are added for the review & resubmit in Sep 2023
#==================================================================================

# Between-strategy-difference SD
metric_person_sdSD <- function(x){
  tryCatch({
    if(!is.matrix(x)) stop()
    mean(apply((x[-1, ] - x[-nrow(x),])^2,1,sd))
  }, error = function(e){
    mean(apply((x[-1, ] - x[-nrow(x),])^2,1,sd))
  })
}

# wrappers to calculate Bray-Curtis dissimilarity in successive difference
# return NA for comparisons with NA
calc.bray.pair <- function(row1, row2) {
  if(is.na(sum(row1, row2))){
    c(beta.bray.bal = NA, beta.bray.gra = NA, beta.bray = NA)
  }else{
    (beta.pair.abund(rbind(row1,row2),"bray"))
  }
}
calc.bray.suc <- function(x) {
  suc.temp <- as.data.frame(do.call(rbind, lapply(1:(nrow(x) - 1), function(i) calc.bray.pair(x[i, ], x[i + 1, ]))))
  suc.temp[] <- as.numeric(unlist(suc.temp))
  suc.temp <- colMeans(suc.temp, na.rm = TRUE)
  suc.temp
}

# wrappers to calculate Bray-Curtis dissimilarity in all-moment comparisons
calc.bray.amm <- function(x) {
  x <- x[complete.cases(x), ]
  unlist(lapply(beta.pair.abund(x,"bray"),mean))
}

