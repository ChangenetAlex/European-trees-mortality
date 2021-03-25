# script on the 13/08/2020
# put together the dfplot final with the new cvalculated mortality and recruitment count at the hectare level

rm(list = ls())
gc()
Dir <- c("~/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
tfinal <- readRDS("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/tfinal.biotic.August2020.rds")
library(parallel)

i <- 1
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir = c(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/", CODE, "/CLIMAP/Recrut.Mortality.2020/"))
  setwd(Dir)
  list.files(Dir, pattern = paste0(".rds"))
  dfplot <- readRDS(paste0("dfplotFinalV2",CODE,seuil,"R.M.rds")) #Base de données plot
  
  test <- tfinal[tfinal$code==CODE,]        ## all trees but only the species we want 
  test <- test[!duplicated(test$plotcode),] ## only unique plot
  print(nrow(dfplot))          #nummber of plot including french plots 
  print(nrow(dfplot[dfplot$country!="FR",])) # expected number oif macthing plots between the two database
  print(nrow(test)) # number of plots (without the french) but before all removing from the script 10.recruitment. 
  length(which(test$plotcode%in%dfplot$plotcode)) 
  # print(length(match(test$plotcode, dfplot$plotcode,incomparables = NA)))
  # print(length(match(dfplot$plotcode, test$plotcode,incomparables = NA)))
  dfplot$sp.mortality.ha.weight <- test$sp.mortality.ha.weight[match(dfplot$plotcode, test$plotcode,incomparables = NA)]
  dfplot$sp.SUM.ha.weight1 <- test$sp.SUM.ha.weight1[match(dfplot$plotcode, test$plotcode,incomparables = NA)]
  dfplot$sp.recruitment.ha.weight <- test$sp.recruitment.ha.weight[match(dfplot$plotcode, test$plotcode,incomparables = NA)]
  dfplot$sp.SUM.ha.weight2 <- test$sp.SUM.ha.weight2[match(dfplot$plotcode, test$plotcode,incomparables = NA)]
  print(nrow(dfplot[!is.na(dfplot$sp.recruitment.ha.weight),])) # number of plots without na values = number of plots in dfplot that are not french 
  print(nrow(dfplot[!is.na(dfplot$sp.mortality.ha.weight),]))
  print(nrow(dfplot[!is.na(dfplot$sp.SUM.ha.weight2),]))
  print(nrow(dfplot[!is.na(dfplot$sp.SUM.ha.weight1),]))
  saveRDS(dfplot, paste0(Dir, "dfplotFinalV3",CODE,seuil,"R.M.rds"))
  rm(list = c("dfplot","test"))
}

dfplot3 <- readRDS(paste0("dfplotFinalV3",CODE,seuil,"R.M.rds")) #Base de données plot




test <- tfinal.biotic[tfinal.biotic$plotcode=="ES100017A1",]
sum(test[test$treestatus_th%in%c("4","2") & !is.na(test$weight1) & test$code=="ACECAM","weight1"])
sum(test$weight1,na.rm=T)
sum(test[test$treestatus_th%in%c("4") & !is.na(test$weight1),"weight1"])


