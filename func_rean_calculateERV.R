# Variable names - meaning of prefixes:
#   moment = a variable of that particular observation (moment)
#   person = person-mean, i.e., the mean of a certain variable within a person
#   grand  = grand-mean, i.e., the mean across all persons
#
# Variable names - meaning of suffixes:
#   Approach of temporal comparisons
#     .single = indices without temporal comparison (between-strategy SD of a certain moment)
#     .suc = successive difference
#     .amm = all-moment comparison
#   Components of moment-level variables (see Bolger & Niall, 2013 for details)
#     cw = within-person component, given by person-mean-centering
#     cb = between-person component, given by person-mean - grand-mean.
#        this component does not change over time.

loadESMcalculateERV <- function (datasource){
#------------------------------
masterscalemax <- 6

# ====================
# Load data
# ====================

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

#identify the first beep within each person
df$firstbeep <- df$triggerid==ave(df$triggerid, df$ppnr, FUN = min)
# boolean: does the observation has complete ER strategy ratings (i.e., no missingness)?
df$b_completeER <- complete.cases(df[,inputER])
if(min(df$triggerid,na.rm = TRUE)==0){
  df$triggerid <- df$triggerid +1
}
# Create boolean of which rows suitable for successive comparisons. Three conditions:
# 1. if the observation is NOT the first beep of a person
# 2. if the moment of interest has complete ER strategy ratings (i.e., no missingness)
# 3. if the moment before the moment of interest has complete ER strategy ratings

df$b_completeERL1 <- lagvar(b_completeER, id=ppnr, obs=triggerid, data=df)
df$b_completeERL1[is.na(df$b_b_completeERL1)] <- FALSE
df$b_successive <- !df$firstbeep & df$b_completeER & df$b_completeERL1

# Harmonization to make everything between 0 to masterscalemax (set at 6 for this manuscript)

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

# ========================================
# Calculate/Center NA, ER and time
# ========================================

#person-mean-center time
df$timecw <- calc.mcent(triggerid, ppnr, data=df)

#moment-level NA and ER endorsement, and their person-level
df$moment_meanNA <- rowMeans(df[,inputEmotions])
df$moment_meanER <- rowMeans(df[,inputER])
df$person_meanNA <- calc.mean(moment_meanNA, ppnr, data=df,expand=TRUE)
df$person_meanER <- calc.mean(moment_meanER, ppnr, data=df,expand=TRUE)

# ========================================
# Calculate SD-based indices
# ========================================

# moment-level ER-between-SD without temporal comparison and its person-level
df$moment_betweenSD.single <- apply(df[,inputER],1,sd)
df$person_betweenSD.single <- calc.mean(moment_betweenSD.single, ppnr, data=df,expand=TRUE)

# moment-level ER-between-relative-SD and its person-level
# the relativeSD function cannot handle a vector of NA;
dfTemp <- df[rowSums(is.na(df[,inputER])) < (ncol(df[,inputER])-1), #exclude rows with all NA or only 1 ER score
              append(c("ppnr","triggerid"),inputER)] #include only ER items and binding info
dfTemp$moment_betweenRSD.single <- apply(dfTemp[,inputER],1,relativeSD,MIN = scalemin, MAX=scalemax)
dfTemp$person_betweenRSD.single <- calc.mean(moment_betweenRSD.single, ppnr, data=dfTemp,expand=TRUE)
dfTemp<- dfTemp[,setdiff(names(dfTemp), inputER)]
df<- merge(df,dfTemp, by=c("ppnr","triggerid"),all=TRUE)


# person-mean-center the moment-level measures
df$moment_betweenSD.singlecw <- calc.mcent(moment_betweenSD.single, ppnr, data=df)
df$moment_betweenRSD.singlecw <- calc.mcent(moment_betweenRSD.single, ppnr, data=df)
df$moment_meanNAcw <- calc.mcent(moment_meanNA, ppnr, data=df)

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


# ======================================
# create successive dissimilarity
# ======================================

# assign NA values for variability indices calculated by successive difference
# which will be overwritten when such indices can be calculated.
df$moment_betweenSD.suc <- NA  # accepts incomplete ER scores
df$moment_betweenRSD.suc <- NA # accepts incomplete ER scores
df$moment_bray.all.suc <- NA   # if one of two rows are all-zero, return 1
df$moment_bray.bal.suc <- NA   # if one of two rows are all-zero, return NA
df$moment_bray.gra.suc <- NA   # if one of two rows are all-zero, return NA

for (i in 1:nrow(df)){
  if(df$b_successive[i]){
    # create a temporary 2-observation matrix
    tempmat <- rbind(df[i,inputER],setNames(df[i-1,inputER],names(df[i,inputER])))
    tempres <- bray.part(tempmat)
    df$moment_bray.all.suc[i] <- tempres$bray[1] #Bray-Curtis dissimilarity full index
    df$moment_bray.bal.suc[i] <- tempres$bray.bal[1]  #Bray-Curtis dissimilarity replacement subcomponent
    df$moment_bray.gra.suc[i] <- tempres$bray.gra[1] #Bray-Curtis dissimilarity nestedness subcomponent
  }
  if(!df$firstbeep[i]){
    # mimicked temporal comparisons for between-strategy SD (successive difference)
    df$moment_betweenSD.suc[i] <- abs(df$moment_betweenSD.single[i] - df$moment_betweenSD.single[i-1])
    df$moment_betweenRSD.suc[i] <- abs(df$moment_betweenRSD.single[i] - df$moment_betweenRSD.single[i-1])
  }
}
inputIndices <- c("moment_betweenSD.suc",
                  "moment_betweenRSD.suc",
                  "moment_bray.all.suc",
                  "moment_bray.bal.suc",
                  "moment_bray.gra.suc")
df$b_completeIndices <- complete.cases(df[,inputIndices])

# calculate moment-level variability indices centered within a person
df$moment_betweenSD.succw <- calc.mcent(moment_betweenSD.suc, ppnr, data=df)
df$moment_betweenRSD.succw <- calc.mcent(moment_betweenRSD.suc, ppnr, data=df)
df$moment_bray.all.succw <-calc.mcent(moment_bray.all.suc, ppnr, data=df)
df$moment_bray.bal.succw <-calc.mcent(moment_bray.bal.suc, ppnr, data=df)
df$moment_bray.gra.succw <-calc.mcent(moment_bray.gra.suc, ppnr, data=df)

# calculate person-level variability indices,
# i.e., value is the same across all time points within a person
df$person_betweenSD.suc <- calc.mean(moment_betweenSD.suc, ppnr, data=df,expand=TRUE)
df$person_betweenRSD.suc <- calc.mean(moment_betweenRSD.suc, ppnr, data=df,expand=TRUE)
df$person_bray.all.suc <- calc.mean(moment_bray.all.suc, ppnr, data=df,expand=TRUE)
df$person_bray.bal.suc <- calc.mean(moment_bray.bal.suc, ppnr, data=df,expand=TRUE)
df$person_bray.gra.suc <- calc.mean(moment_bray.gra.suc, ppnr, data=df,expand=TRUE)

# ======================================
# create all-moment dissimilarity
# ======================================
dfTemp <- df[df$b_completeER & df$moment_meanER != 0, #exclude rows with all NA or only 1 ER score, or all zero
             append(c("ppnr","triggerid"),inputER)] #include only ER items and binding info

# assign NA values for all-moment comparison variability indices
# which will be overwritten when such indices can be calculated.
dfTemp$moment_bray.all.amm <- NA
dfTemp$moment_bray.bal.amm <- NA
dfTemp$moment_bray.gra.amm <- NA
dfTemp$moment_betweenSD.amm <- NA
dfTemp$moment_betweenRSD.amm <- NA

# create empty (temporary) variable
# to hold a vector of Bray-Curtis dissimilarity calculated in all-moment comparison
vbray.all <- NULL
vbray.bal <- NULL
vbray.gra <- NULL

for (i in 1:length(unique(dfTemp$ppnr))){
    # extract a matrix of ER strategies use for each person
    dfPerson <- dfTemp[dfTemp$ppnr==unique(dfTemp$ppnr)[i],]
    matx <- dfPerson[,inputER]
    nobs <- nrow(dfPerson)
    # then, calculate a matrix of Bray-Curtis dissimilarity
    # based on that person's ER strategy use across time
    resbray <- bray.part(matx)
    # obtain a vector of column (with length = nobs) means of the dissimilarity matrix
    # append it to the temporary vector holder
    vbray.all <- append(vbray.all,colMeans(as.matrix(resbray$bray))*nobs/(nobs-1))
    vbray.bal <- append(vbray.bal,colMeans(as.matrix(resbray$bray.bal))*nobs/(nobs-1))
    vbray.gra <- append(vbray.gra,colMeans(as.matrix(resbray$bray.gra))*nobs/(nobs-1))
}

# after looping all persons, the vbray vector will have a length matching all observations in all persons
# assign it back to the working dataframe
dfTemp$moment_bray.all.amm <- vbray.all
dfTemp$moment_bray.bal.amm <- vbray.bal
dfTemp$moment_bray.gra.amm <- vbray.gra

dfTemp<- dfTemp[,setdiff(names(dfTemp), inputER)] # take away the inputER cols
df<- merge(df,dfTemp, by=c("ppnr","triggerid"),all=TRUE)

# calculate moment-level variability indices centered within a person
df$moment_betweenSD.ammcw <- calc.mcent(moment_betweenSD.amm, ppnr, data=df)
df$moment_betweenRSD.ammcw <- calc.mcent(moment_betweenRSD.amm, ppnr, data=df)
df$moment_bray.all.ammcw <-calc.mcent(moment_bray.all.amm, ppnr, data=df)
df$moment_bray.bal.ammcw <-calc.mcent(moment_bray.bal.amm, ppnr, data=df)
df$moment_bray.gra.ammcw <-calc.mcent(moment_bray.gra.amm, ppnr, data=df)

# calculate person-level variability indices,
# i.e., value is the same across all time points within a person
df$person_betweenSD.amm <- calc.mean(moment_betweenSD.amm, ppnr, data=df,expand=TRUE)
df$person_betweenRSD.amm <- calc.mean(moment_betweenRSD.amm, ppnr, data=df,expand=TRUE)
df$person_bray.all.amm <- calc.mean(moment_bray.all.amm, ppnr, data=df,expand=TRUE)
df$person_bray.bal.amm <- calc.mean(moment_bray.bal.amm, ppnr, data=df,expand=TRUE)
df$person_bray.gra.amm <- calc.mean(moment_bray.gra.amm, ppnr, data=df,expand=TRUE)

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
