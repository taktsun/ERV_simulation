# self defined functions

# a function to extract n characters from the right of a string
# adapted from: https://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# a wrapper to extract fixed effect and random effect of an MLM (nlme)
preparemmresult <- function (m){
  cbind(model = deparse(substitute(m)),
        FEest = summary(m)$tTable[,1],
        SE = summary(m)$tTable[,2],
        DF = summary(m)$tTable[,3],
        pvalue = summary(m)$tTable[,5],
        "ranef" = VarCorr(m)[1:(length(attributes(m$terms)$term.labels)+1)],
        residual = m$sigma^2,
        phi = coef(m$modelStruct$corStruct, unconstrained = FALSE), # needs coef to make this work
        nobs = nobs(m))
}

# a wrapper to extract fixed effect and random effect of an MLM (brms)
preparemmresult.brm <- function (m){
  cbind(model = deparse(substitute(m)),
        FEest = fixef(m)[,1],
        SE = fixef(m)[,2],
        DF = ifelse(substrRight(rownames(fixef(m)),2)=="cb",summary(m)$ngrps$ppnr,summary(m)$nobs),
        pvalue = rep("brms",length(fixef(m)[,1])),
        ranef = VarCorr(m)$ppnr$sd[1:nrow(fixef(m))],
        residual = rep(VarCorr(m)$residual__$sd[1]^2,nrow(fixef(m))),
        phi = rep(summary(m)$cor_pars[1],nrow(fixef(m))),
        nobs = rep(summary(m)$nobs,nrow(fixef(m)))
  )

}

# the main function to be used
MLMresults <- function(df, datasource,completeIndices = FALSE){
modelmoment.intercept<-NULL
modelmoment.withinRSD<-NULL
modelmoment.betweenRSD<-NULL
modelmoment.betweenRSD.suc<-NULL
modelmoment.bray.all.suc<-NULL
modelmoment.bray.part.suc<-NULL

if(completeIndices){
  df<- df[df$b_completeIndices,]
}else{
  df<- df[df$b_completeER,]
}

modelmoment.intercept <- lme(fixed=moment_meanNA ~1,
                        data=df,
                        random=~1 | ppnr, correlation = corAR1(),
                        control =list(msMaxIter = 1000, msMaxEval = 1000),na.action = na.omit)

modelmoment.withinRSD <- lme(fixed=moment_meanNA ~ moment_withinRSDcb + timecw,
                           data=df,
                           random=~1 | ppnr, correlation = corAR1(),
                           control =list(msMaxIter = 1000, msMaxEval = 1000),na.action = na.omit)
modelmoment.betweenRSD <- lme(fixed=moment_meanNA ~ moment_betweenRSD.singlecw+moment_betweenRSD.singlecb + timecw,
                         data=df,
                         random=~1+ moment_betweenRSD.singlecw | ppnr, correlation = corAR1(),
                         control =list(msMaxIter = 1000, msMaxEval = 1000),na.action = na.omit)
modelmoment.betweenRSD.suc <- lme(fixed=moment_meanNA ~ moment_betweenRSD.succw +moment_betweenRSD.succb + timecw,
                             data=df,
                             random=~1+ moment_betweenRSD.succw | ppnr, correlation = corAR1(),
                             control =list(msMaxIter = 1000, msMaxEval = 1000),na.action = na.omit)
modelmoment.bray.all.suc <- lme(fixed=moment_meanNA ~ moment_bray.all.succw+moment_bray.all.succb+ timecw,
                                     data=df,
                                     random=~1+ moment_bray.all.succw | ppnr, correlation = corAR1(),
                                     control =list(msMaxIter = 1000, msMaxEval = 1000),na.action = na.omit)
if (datasource==2){
  modelmoment.bray.part.suc <- brm((moment_meanNA) ~    (moment_bray.bal.succw) + (moment_bray.gra.succw)+
                                          moment_bray.bal.succb+ moment_bray.gra.succb+
                                 (1+moment_bray.bal.succw + moment_bray.gra.succw|ppnr) + ar(),
                               data=df)
}else{
  # nlme cannot converge (singular convergence) for dataset2
  modelmoment.bray.part.suc <- lme(fixed=moment_meanNA ~ moment_bray.bal.succw+ moment_bray.gra.succw+
                                          moment_bray.bal.succb+ moment_bray.gra.succb+ timecw,
                                        data=df,
                                        random=~1+ moment_bray.bal.succw+ moment_bray.gra.succw | ppnr, correlation = corAR1(),
                                        control =list(msMaxIter = 1000, msMaxEval = 1000),na.action = na.omit)
}

resmodelmomentest <- rbind(
  preparemmresult(modelmoment.intercept),
  preparemmresult(modelmoment.withinRSD),
  preparemmresult(modelmoment.betweenRSD),
  preparemmresult(modelmoment.betweenRSD.suc),
  preparemmresult(modelmoment.bray.all.suc)
)
if (datasource==2){
  resbrm<- (preparemmresult.brm(modelmoment.bray.part.suc))
  names(resbrm)[names(resbrm) == 'Estimate'] <- 'phi'
  resmodelmomentest<- rbind(
    resmodelmomentest,
    resbrm)
}else{
  resmodelmomentest<- rbind(
    resmodelmomentest,
    preparemmresult(modelmoment.bray.part.suc)) # error for dataset2
}

cbind( dataset = rep(datasource,nrow(resmodelmomentest)),resmodelmomentest)
}
