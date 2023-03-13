loadESMdata <- function (datasource){
#------------------------------
masterscalemax <- 6

# Load data

if(datasource == 1){
currentstudy <- "BLANKE 1"
df <- read.csv("dfraw1.csv")
inputEmotions <- c("M_nerv","M_nieder","M_bekuem")
inputER <- c("M_rum2","M_rum1","M_distr2","M_distr1","M_refl2","M_refl1")
scalemin <- 0
scalemax <- 6
CESDmax <- 4
}else if(datasource ==2){
currentstudy <- "BLANKE 2"
df <- read.csv("dfraw2.csv")
inputEmotions <- c("ANGRY","SAD","ANX","DEPRE")
inputER <- c("RUMI.","DIST.","RFLCT.","REAPP.","SUPP.","SOCSH.")
scalemin <- 1
scalemax <- 100
CESDmax <- 3
}else if(datasource == 3){
currentstudy <- "BLANKE 3"
df <- read.csv("dfraw3.csv")
inputEmotions <- c("kwaad","droevig","angstig","depressief")
inputER <- c("Rumi_past","Rumi_future","Dist","Reap","Supp","Sosu")
scalemin <- 0
scalemax <- 100
CESDmax <- 3
}
# ============RUN BELOW=================

#identify the first beep within each person
df$firstbeep <- df$triggerid==ave(df$triggerid, df$ppnr, FUN = min)
df$b_completeER <- complete.cases(df[,inputER])

# Harmonization to make everything between 0 to 6 (or masterscalemax value)
df$person_CESD <- df$person_CESD/CESDmax*masterscalemax
if(scalemin != 0 || scalemax != masterscalemax){
  if(scalemin != 0){
  df[,append(inputEmotions,inputER)] <- sapply(df[,append(inputEmotions,inputER)],
                                               function(x) x-scalemin)
  }
  if((scalemax-scalemin) != masterscalemax){
    df[,append(inputEmotions,inputER)] <- sapply(df[,append(inputEmotions,inputER)],
                                                 function(x) x*masterscalemax/(scalemax-scalemin))
  }
  scalemin = 0
  scalemax = masterscalemax
}


df[,append(inputEmotions,inputER)] <- sapply(df[,append(inputEmotions,inputER)], function(x) as.numeric(x))

#moment-level ER-between-SD without temporal comparison and its person-level
df$moment_betweenSD.single <- apply(df[,inputER],1,sd)
df$person_betweenSD.single <- calc.mean(moment_betweenSD.single, ppnr, data=df,expand=TRUE)

#moment-level ER-between-relative-SD and its person-level
# the relativeSD function cannot handle a vector of NA; also MIN and MAX are hardcoded so wrong values have to be removed
dfTemp <- df[rowSums(is.na(df[,inputER])) < (ncol(df[,inputER])-1), #exclude rows with all NA or only 1 ER score
              append(c("ppnr","triggerid"),inputER)] #include only ER items and binding info
dfTemp$moment_betweenRSD.single <- apply(dfTemp[,inputER],1,relativeSD,MIN = scalemin, MAX=scalemax)
dfTemp$person_betweenRSD.single <- calc.mean(moment_betweenRSD.single, ppnr, data=dfTemp,expand=TRUE)
dfTemp<- dfTemp[,setdiff(names(dfTemp), inputER)]
df<- merge(df,dfTemp, by=c("ppnr","triggerid"),all=TRUE)

#moment-level NA and ER endorsement, and their person-level
df$moment_meanNA <- rowMeans(df[,inputEmotions])
df$moment_meanER <- rowMeans(df[,inputER])
df$moment_meanER.suc <- (lagvar(x = moment_meanER, id = ppnr,
                              obs=triggerid, data=df, day = day) + df$moment_meanER)/2
df$moment_meanNAL1D <-lagvar(x = moment_meanNA, id = ppnr,
                             obs=triggerid, data=df, day = day)

df$person_meanNA <- calc.mean(moment_meanNA, ppnr, data=df,expand=TRUE)
df$person_meanER <- calc.mean(moment_meanER, ppnr, data=df,expand=TRUE)
df$moment_meanER.amm <- (df$person_meanER + df$moment_meanER)/2


#person-mean-center the moment-level measures
df$moment_betweenSD.singlecw <- calc.mcent(moment_betweenSD.single, ppnr, data=df)
df$moment_betweenRSD.singlecw <- calc.mcent(moment_betweenRSD.single, ppnr, data=df)
df$moment_meanNAcw <- calc.mcent(moment_meanNA, ppnr, data=df)
df$moment_meanERcw <- calc.mcent(moment_meanER, ppnr, data=df)
df$moment_meanER.succw <- calc.mcent(moment_meanER.suc, ppnr, data=df)
df$moment_meanER.ammcw <- calc.mcent(moment_meanER.amm, ppnr, data=df)
df$moment_meanNAcwL1D <-lagvar(x = moment_meanNAcw, id = ppnr,
                             obs=triggerid, data=df, day = day)

#person-mean-center time
df$timecw <- calc.mcent(triggerid, ppnr, data=df)

# Create boolean of which rows suitable for succesive comparisons

df$b_successiveL1 <- lagvar(b_completeER, id=ppnr, obs=triggerid, data=df)
df$b_successiveL1[is.na(df$b_successiveL1)] <- FALSE
df$b_successive <- !df$firstbeep & df$b_completeER & df$b_successiveL1
df$b_successiveL1 <- NULL

# create successive dissimilarity

df$moment_betweenSD.suc <- NA  # accepts incomplete ER scores
df$moment_betweenRSD.suc <- NA # accepts incomplete ER scores
df$moment_bray.all.suc <- NA   # if one of two rows are all-zero, return 1
df$moment_bray.bal.suc <- NA   # if one of two rows are all-zero, return NA
df$moment_bray.gra.suc <- NA   # if one of two rows are all-zero, return NA

for (i in 1:nrow(df)){
  if(df$b_successive[i]){
    tempmat <- rbind(df[i,inputER],setNames(df[i-1,inputER],names(df[i,inputER])))
    tempres <- bray.part(tempmat)
    # betapart does not accept na.rm = TRUE but vegan does.
    df$moment_bray.all.suc[i] <- tempres$bray[1]
    df$moment_bray.bal.suc[i] <- tempres$bray.bal[1]
    df$moment_bray.gra.suc[i] <- tempres$bray.gra[1]

  }
  if(!df$firstbeep[i]){
    df$moment_betweenSD.suc[i] <- abs(df$moment_betweenSD.single[i] - df$moment_betweenSD.single[i-1])
    df$moment_betweenRSD.suc[i] <- abs(df$moment_betweenRSD.single[i] - df$moment_betweenRSD.single[i-1])
  }
}

df$moment_betweenSD.succw <- calc.mcent(moment_betweenSD.suc, ppnr, data=df)
df$moment_betweenRSD.succw <- calc.mcent(moment_betweenRSD.suc, ppnr, data=df)
df$moment_bray.all.succw <-calc.mcent(moment_bray.all.suc, ppnr, data=df)
df$moment_bray.bal.succw <-calc.mcent(moment_bray.bal.suc, ppnr, data=df)
df$moment_bray.gra.succw <-calc.mcent(moment_bray.gra.suc, ppnr, data=df)


df$person_betweenSD.suc <- calc.mean(moment_betweenSD.suc, ppnr, data=df,expand=TRUE)
df$person_betweenRSD.suc <- calc.mean(moment_betweenRSD.suc, ppnr, data=df,expand=TRUE)
df$person_bray.all.suc <- calc.mean(moment_bray.all.suc, ppnr, data=df,expand=TRUE)
df$person_bray.bal.suc <- calc.mean(moment_bray.bal.suc, ppnr, data=df,expand=TRUE)
df$person_bray.gra.suc <- calc.mean(moment_bray.gra.suc, ppnr, data=df,expand=TRUE)

# create all-moment dissimilarity
dfTemp <- df[df$b_completeER & df$moment_meanER != 0, #exclude rows with all NA or only 1 ER score, or all zero
             append(c("ppnr","triggerid"),inputER)] #include only ER items and binding info

dfTemp$moment_bray.all.amm <- NA
dfTemp$moment_bray.bal.amm <- NA
dfTemp$moment_bray.gra.amm <- NA
dfTemp$moment_betweenSD.amm <- NA
dfTemp$moment_betweenRSD.amm <- NA

vbray.all <- 0
vbray.bal <- 0
vbray.gra <- 0

for (i in 1:length(unique(dfTemp$ppnr))){
    dfPerson <- dfTemp[dfTemp$ppnr==unique(dfTemp$ppnr)[i],]
    matx <- dfPerson[,inputER]
    nobs <- nrow(dfPerson)
    resbray <- bray.part(matx)
    vbray.all <- append(vbray.all,colMeans(as.matrix(resbray$bray))*nobs/(nobs-1))
    vbray.bal <- append(vbray.bal,colMeans(as.matrix(resbray$bray.bal))*nobs/(nobs-1))
    vbray.gra <- append(vbray.gra,colMeans(as.matrix(resbray$bray.gra))*nobs/(nobs-1))
}

dfTemp$moment_bray.all.amm <- vbray.all[2:length(vbray.all)]
dfTemp$moment_bray.bal.amm <- vbray.bal[2:length(vbray.bal)]
dfTemp$moment_bray.gra.amm <- vbray.gra[2:length(vbray.gra)]

dfTemp<- dfTemp[,setdiff(names(dfTemp), inputER)] # take away the inputER cols
df<- merge(df,dfTemp, by=c("ppnr","triggerid"),all=TRUE)

df$moment_betweenSD.ammcw <- calc.mcent(moment_betweenSD.amm, ppnr, data=df)
df$moment_betweenRSD.ammcw <- calc.mcent(moment_betweenRSD.amm, ppnr, data=df)
df$moment_bray.all.ammcw <-calc.mcent(moment_bray.all.amm, ppnr, data=df)
df$moment_bray.bal.ammcw <-calc.mcent(moment_bray.bal.amm, ppnr, data=df)
df$moment_bray.gra.ammcw <-calc.mcent(moment_bray.gra.amm, ppnr, data=df)

df$person_betweenSD.amm <- calc.mean(moment_betweenSD.amm, ppnr, data=df,expand=TRUE)
df$person_betweenRSD.amm <- calc.mean(moment_betweenRSD.amm, ppnr, data=df,expand=TRUE)
df$person_bray.all.amm <- calc.mean(moment_bray.all.amm, ppnr, data=df,expand=TRUE)
df$person_bray.bal.amm <- calc.mean(moment_bray.bal.amm, ppnr, data=df,expand=TRUE)
df$person_bray.gra.amm <- calc.mean(moment_bray.gra.amm, ppnr, data=df,expand=TRUE)

# calculate within strategy SD per strategy, then person level within-SD
df<-cbind(df, df[,c(inputER)])
names(df)[(length(names(df))+1-length(inputER)):length(names(df))]<-paste0(inputER,"_SD")
df[,c(paste0(inputER,"_SD"))] <- lapply(df[,c(paste0(inputER,"_SD"))],
                                        function(x) ave(x, df$ppnr, FUN = function(x) sd(x, na.rm = TRUE)))
df$person_withinSD <- rowMeans((df[,c(paste0(inputER,"_SD"))]))

# calculate RELATIVE within strategy SD per strategy, then person level within-RSD
df<-cbind(df, df[,c(inputER)])
names(df)[(length(names(df))+1-length(inputER)):length(names(df))]<-paste0(inputER,"_RSD")
df[,c(paste0(inputER,"_RSD"))] <- lapply(df[,c(paste0(inputER,"_RSD"))],
                                        function(x) ave(x, df$ppnr, FUN = function(x) relativeSD(x,MIN = scalemin, MAX=scalemax)))
df$person_withinRSD <- rowMeans((df[,c(paste0(inputER,"_RSD"))]),na.rm = TRUE)

# calculate grand means and the between-person component

df$grand_moment_meanER <- mean(calc.mean(moment_meanER, ppnr, data=df),na.rm=TRUE)
df$grand_moment_meanNA <- mean(calc.mean(moment_meanNA, ppnr, data=df),na.rm=TRUE)
df$grand_moment_betweenSD.single <- mean(calc.mean(moment_betweenSD.single, ppnr, data=df),na.rm=TRUE)
df$grand_moment_betweenRSD.single <- mean(calc.mean(moment_betweenRSD.single, ppnr, data=df),na.rm=TRUE)
df$grand_moment_withinRSD <- mean(calc.mean(person_withinRSD, ppnr, data=df),na.rm=TRUE)

df$grand_moment_betweenSD.amm <- mean(calc.mean(moment_betweenSD.amm, ppnr, data=df),na.rm=TRUE)
df$grand_moment_betweenRSD.amm <- mean(calc.mean(moment_betweenRSD.amm, ppnr, data=df),na.rm=TRUE)
df$grand_moment_bray.all.amm <- mean(calc.mean(moment_bray.all.amm, ppnr, data=df),na.rm=TRUE)
df$grand_moment_bray.bal.amm <- mean(calc.mean(moment_bray.bal.amm, ppnr, data=df),na.rm=TRUE)
df$grand_moment_bray.gra.amm <- mean(calc.mean(moment_bray.gra.amm, ppnr, data=df),na.rm=TRUE)

df$grand_moment_betweenSD.suc <- mean(calc.mean(moment_betweenSD.suc, ppnr, data=df),na.rm=TRUE)
df$grand_moment_betweenRSD.suc <- mean(calc.mean(moment_betweenRSD.suc, ppnr, data=df),na.rm=TRUE)
df$grand_moment_bray.all.suc <- mean(calc.mean(moment_bray.all.suc, ppnr, data=df),na.rm=TRUE)
df$grand_moment_bray.bal.suc <- mean(calc.mean(moment_bray.bal.suc, ppnr, data=df),na.rm=TRUE)
df$grand_moment_bray.gra.suc <- mean(calc.mean(moment_bray.gra.suc, ppnr, data=df),na.rm=TRUE)

df$moment_meanERcb <- df$person_meanER - df$grand_moment_meanER
df$moment_meanNAcb <- df$person_meanNA - df$grand_moment_meanNA
df$moment_betweenSD.singlecb <- df$person_betweenSD.single - df$grand_moment_betweenSD.single
df$moment_betweenRSD.singlecb <- df$person_betweenRSD.single - df$grand_moment_betweenRSD.single
df$moment_withinRSDcb <- df$person_withinRSD - df$grand_moment_withinRSD

df$moment_betweenSD.ammcb <- df$person_betweenSD.amm - df$grand_moment_betweenSD.amm
df$moment_betweenRSD.ammcb <- df$person_betweenRSD.amm - df$grand_moment_betweenRSD.amm
df$moment_bray.all.ammcb <- df$person_bray.all.amm - df$grand_moment_bray.all.amm
df$moment_bray.bal.ammcb <- df$person_bray.bal.amm - df$grand_moment_bray.bal.amm
df$moment_bray.gra.ammcb <- df$person_bray.gra.amm - df$grand_moment_bray.gra.amm

df$moment_betweenSD.succb <- df$person_betweenSD.suc - df$grand_moment_betweenSD.suc
df$moment_betweenRSD.succb <- df$person_betweenRSD.suc - df$grand_moment_betweenRSD.suc
df$moment_bray.all.succb <- df$person_bray.all.suc - df$grand_moment_bray.all.suc
df$moment_bray.bal.succb <- df$person_bray.bal.suc - df$grand_moment_bray.bal.suc
df$moment_bray.gra.succb <- df$person_bray.gra.suc - df$grand_moment_bray.gra.suc

df
}
