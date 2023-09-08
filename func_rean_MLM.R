# extract n characters from the right of a string
#   adapted from: https://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r
#   useful for preparing summary of results later in preparemmresult.brm
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
        nobs = nobs(m),
        RMSE = performance_rmse(m, normalized = FALSE)
        )

}

# a wrapper to extract fixed effect and random effect of an MLM (brms)
preparemmresult.brm <- function (m){
  cbind(model = deparse(substitute(m)),
        FEest = fixef(m)[,1],
        SE = fixef(m)[,2],
        DF = ifelse(substrRight(rownames(fixef(m)),2)=="cb",summary(m)$ngrps$ppnr,summary(m)$nobs),
        pvalue = pt(q = abs(fixef(m)[,1]/fixef(m)[,2]),
                    df = ifelse(substrRight(rownames(fixef(m)),2)=="cb",summary(m)$ngrps$ppnr,summary(m)$nobs),
                    lower.tail = FALSE)*2,
        ranef = VarCorr(m)$ppnr$sd[1:nrow(fixef(m))],
        residual = rep(VarCorr(m)$residual__$sd[1]^2,nrow(fixef(m))),
        phi = rep(summary(m)$cor_pars[1],nrow(fixef(m))),
        nobs = rep(summary(m)$nobs,nrow(fixef(m))),
        RMSE = sqrt(sum(residuals(m)[,1]^2)/length(residuals(m)[,1]))
  )

}

# the main function to be called in reanalysis_main.R
MLMresults <- function(df, datasource,completeIndices = FALSE){
modelmoment.intercept<-NULL
modelmoment.withinRSD<-NULL
modelmoment.betweenRSD<-NULL
modelmoment.betweenRSD.suc<-NULL
modelmoment.bray.all.suc<-NULL
modelmoment.bray.part.suc<-NULL

# completeIndices is default to be FALSE:
# which corresponds to the manuscript, where we estimated multilevel models with available observations for each index
# if completeIndices=TRUE, we exclude an observation if there is ANY missing/NA ER variability indices
if(completeIndices){
  df<- df[df$b_completeIndices,]
}else{
  df<- df[df$b_completeER,]
}

# Multilevel modeling (MLM) goes below
modelmoment.intercept <- lme(fixed=moment_meanNA ~1,
                        data=df,
                        random=~1 | ppnr, correlation = corAR1(),
                        control =lmeControl(msMaxIter = 1000, msMaxEval = 1000,opt='optim'),na.action = na.omit)

modelmoment.withinRSD <- lme(fixed=moment_meanNA ~ moment_withinRSDcb + timecw,
                           data=df,
                           random=~1 | ppnr, correlation = corAR1(),
                           control =lmeControl(msMaxIter = 1000, msMaxEval = 1000,opt='optim'),na.action = na.omit)
modelmoment.betweenRSD <- lme(fixed=moment_meanNA ~ moment_betweenRSD.singlecw+moment_betweenRSD.singlecb + timecw,
                         data=df,
                         random=~1+ moment_betweenRSD.singlecw | ppnr, correlation = corAR1(),
                         control =lmeControl(msMaxIter = 1000, msMaxEval = 1000,opt='optim'),na.action = na.omit)
modelmoment.betweenRSD.suc <- lme(fixed=moment_meanNA ~ moment_betweenRSD.succw +moment_betweenRSD.succb + timecw,
                             data=df,
                             random=~1+ moment_betweenRSD.succw | ppnr, correlation = corAR1(),
                             control =lmeControl(msMaxIter = 1000, msMaxEval = 1000,opt='optim'),na.action = na.omit)
modelmoment.bray.all.suc <- lme(fixed=moment_meanNA ~ moment_bray.all.succw+moment_bray.all.succb+ timecw,
                                     data=df,
                                     random=~1+ moment_bray.all.succw | ppnr, correlation = corAR1(),
                                     control =lmeControl(msMaxIter = 1000, msMaxEval = 1000,opt='optim'),na.action = na.omit)
modelmoment.bray.part.suc <- lme(fixed=moment_meanNA ~ moment_bray.bal.succw+ moment_bray.gra.succw+
                                          moment_bray.bal.succb+ moment_bray.gra.succb+ timecw,
                                        data=df,
                                        random=~1+ moment_bray.bal.succw+ moment_bray.gra.succw | ppnr, correlation = corAR1(),
                                        control =lmeControl(msMaxIter = 1000, msMaxEval = 1000,opt='optim'),na.action = na.omit)
modelmoment.withinSD <- lme(fixed=moment_meanNA ~ moment_withinSDcb + timecw,
                             data=df,
                             random=~1 | ppnr, correlation = corAR1(),
                             control =lmeControl(msMaxIter = 1000, msMaxEval = 1000,opt='optim'),na.action = na.omit)
modelmoment.betweenSD <- lme(fixed=moment_meanNA ~ moment_betweenSD.singlecw+moment_betweenSD.singlecb + timecw,
                              data=df,
                              random=~1+ moment_betweenSD.singlecw | ppnr, correlation = corAR1(),
                              control =lmeControl(msMaxIter = 1000, msMaxEval = 1000,opt='optim'),na.action = na.omit)
modelmoment.betweenSD.suc <- lme(fixed=moment_meanNA ~ moment_betweenSD.succw +moment_betweenSD.succb + timecw,
                                  data=df,
                                  random=~1+ moment_betweenSD.succw | ppnr, correlation = corAR1(),
                                  control =lmeControl(msMaxIter = 1000, msMaxEval = 1000,opt='optim'),na.action = na.omit)


resmodelmomentest <- rbind(
  preparemmresult(modelmoment.intercept),
  preparemmresult(modelmoment.withinSD),
  preparemmresult(modelmoment.betweenSD),
  preparemmresult(modelmoment.betweenSD.suc),
  preparemmresult(modelmoment.withinRSD),
  preparemmresult(modelmoment.betweenRSD),
  preparemmresult(modelmoment.betweenRSD.suc),
  preparemmresult(modelmoment.bray.all.suc),
  preparemmresult(modelmoment.bray.part.suc)
)

cbind( dataset = rep(datasource,nrow(resmodelmomentest)),resmodelmomentest)
}
