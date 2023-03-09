# correction factor for SD
correctSD <- function(a,n){
  # see section 1.1 of Beran (1994)
  a<- abs(a)
  d <- 1+2*a/(1-a) # eq. 1.14
  cf <- d*(1-(1/(1-a)/n+ a^n/(1-a)/n)) # eq 1.12
  cf^0.25
}

# measurement correction 1: range bound

correct_range_bound <- function (df){
df[df < 0] <- 0.0001
df[df > 1] <- 1
df
}

# tbc: scale up to scalemax/scale down everything to range 0-1?
# tbc: rounding to integer

get_sign_vector <- function (n,ERn) {
# create a vector
signvector <- c(0,sign(rnorm(1))*as.vector(rbind(rep(1,ERn), rep(-1,ERn))))
signvector <- signvector[1:ERn]
# center to 0; repeat n times to fit # of observations
signvector <- rep(signvector-sum(signvector)/ERn, each = n)
signvector
}
