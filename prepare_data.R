# In this file, write the R-code necessary to load your original data file
# (e.g., an SPSS, Excel, or SAS-file), and convert it to a data.frame. Then,
# use the function open_data(your_data_frame) or closed_data(your_data_frame)
# to store the data.

library(worcs)
library(haven) # for reading sav files
library(lubridate) #for function day

# define function for handling 3 datasets

prepRawData <- function(rawdata, sID, sBeep, sTS, sDateFormat, cEmotion, cER, cCESD){
d <- data.frame(rawdata)
colnames(d)[which(names(d) == sID)] <- "ppnr"
colnames(d)[which(names(d) == sBeep)] <- "triggerid"
colnames(d)[which(names(d) == sTS)] <- "timestamp"
dateformat<-sDateFormat
d$day <-day(parse_date_time(d$timestamp,orders=dateformat))
inputEmotions <- cEmotion
inputER <- cER
inputCESD <- cCESD
inputNeeded <- append(c("ppnr","triggerid","day"),inputEmotions)
inputNeeded <- append(inputNeeded,inputER)
d <- d[,inputNeeded]
d$person_CESD <- rowMeans(rawdata[,inputCESD],na.rm = TRUE)
d
}


# Using Data1.sav, Data2.sav and Data3.sav from https://osf.io/mxjfh/ (Blanke et al., 2020)

raw1 <- read_sav('https://osf.io/download/w8y33/')
raw2 <- read_sav('https://osf.io/download/gm52c/')
raw3 <- read_sav('https://osf.io/download/uvqjh/')

dfraw1 <- prepRawData(rawdata = raw1,
                      sID = "ID_anonym",
                      sBeep = "a_ftl_0",
                      sTS = "A_time",
                      sDateFormat = "ymd HMS",
                      cEmotion = c("M_nerv","M_nieder","M_bekuem"),
                      cER = c("M_rum2","M_rum1","M_distr2","M_distr1","M_refl2","M_refl1"),
                      cCESD = paste0("T1_CESD",seq(1:20)))

dfraw2 <- prepRawData(rawdata = raw2,
                      sID = "PpID",
                      sBeep = "BeepNo",
                      sTS = "date_time",
                      sDateFormat = "ymd HMS",
                      cEmotion = c("ANGRY","SAD","ANX","DEPRE"),
                      cER = c("RUMI.","DIST.","RFLCT.","REAPP.","SUPP.","SOCSH."), # they end with # originally but parsed as . upon turning into data.frame
                      cCESD = paste0("cesd",seq(1:20)))


dfraw3 <- prepRawData(rawdata = raw3,
                      sID = "PpID",
                      sBeep = "BeepNr",
                      sTS = "date_time",
                      sDateFormat = "ymd HMS",
                      cEmotion = c("kwaad","droevig","angstig","depressief"),
                      cER = c("Rumi_past","Rumi_future","Dist","Reap","Supp","Sosu"),
                      cCESD = paste0("cesd",seq(1:20),"_w1"))

open_data(dfraw1)
open_data(dfraw2)
open_data(dfraw3)

