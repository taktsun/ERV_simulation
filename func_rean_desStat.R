desstat <- function(x,group){



  c(
    # variable name
    var = deparse(substitute(x)),
    # grand mean
    mean = round(mean(x, na.rm = TRUE),2),

    # grand SD
    # grandSD = round(sd(x, na.rm =TRUE),3)

    # within person SD
    wSD = round(mean(unlist(aggregate(x, by=list(group), FUN=sd, na.rm = TRUE)[2]),na.rm=TRUE),2),

    # between person SD
    bSD = round(sd(unlist(aggregate(x, by=list(group), FUN=mean, na.rm = TRUE)[2]),na.rm=TRUE),2),

    # ICC: how much of the variance is between-person
    ICC = round(multilevel.icc(x, cluster = group),2),

    # n obs
    nobs = sum(!is.na(x))
    )
}


summarydesstat <- function(df){

  rbind(
  desstat(df$moment_meanNA,df$ppnr),
  desstat(df$moment_meanER,df$ppnr),
  c("df$withinRSD",round(mean(df[df$firstbeep,"person_withinRSD"]),2),
    0, round(sd(df[df$firstbeep,"person_withinRSD"]),2), 1,nrow(df[df$firstbeep,])),
  desstat(df$moment_betweenRSD.single,df$ppnr),
  desstat(df$moment_betweenRSD.suc,df$ppnr),
  desstat(df$moment_bray.all.suc,df$ppnr),
  desstat(df$moment_bray.bal.suc,df$ppnr),
  desstat(df$moment_bray.gra.suc,df$ppnr)
)
}

