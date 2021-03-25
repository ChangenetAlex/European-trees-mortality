rm(list = ls())
gc()
library(glmmTMB)
library(car)
library(emmeans)
library(effects)
library(multcomp)
library(MuMIn)
library(DHARMa)
library(broom)
library(broom.mixed)
library(dotwhisker)
library(ggplot2); theme_set(theme_bw())
library(texreg)
library(xtable)
library(huxtable)
library(plyr)
library(beanplot)

Dir <- c("/home/alexandrechangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data")
setwd("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/")

mod_ALLbin <- list()
mod_ALLneg <- list()
dfplot <- list()
dftree <- list()

# Then try to obtain predicted mortality and observed plot mortality 
i = 1
Allcode <- c("ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINPINA","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
Allmod <- c("Mbin14A.19","Mbin13B.26","Mbin13A.22","M2bin13B.22","M2bin13A21","M2bin15B23","Mbin3C.26","M2bin7B.17","Mbin13A.18","Mbin11B.19","Mbin3B.20","Mbin15B.21","Mbin5B.17","Mbin3B.31","Mbin13A.27","M2bin1C.20","M2bin1B.23","Mbin13B.26","Mbin7A.26","Mbin1B.23")

for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir = c("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/") # Directory 
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/binomial/",Allmod[i],"/"))
  setwd(Dir)
  mod_ALLbin[[CODE]] <- get(load(file = list.files(pattern=".rda")))
  
  Dir = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Recrut.Mortality.2020/")
  setwd(Dir)
  dfplot[[CODE]] <- readRDS(paste0("dfplotV2", CODE, seuil, "R.M.rds")) #Base de données plot full No scale !!! sans zone de transition
  dftree[[CODE]] <- readRDS(paste0("Mydf3_",CODE,"_",seuil,"_",seuilC,"R.M.rds"))
}


i = 1
Allcode <- c("ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINPINA","PINHAL","PINNIG","PINPIN","PINSYL","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.8,0.7,0.7,0.7,0.8,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
Allmod <- c("MnbZT14A.20","MnbZT13B.23","MnbZT13A.21","M2nbZT13B.29","M2nbZT13A.22","M2nbZT15B.24","MnbZT3C.24","M2nbZT7B.26","MnbZT13A.19","MnbZT11B.22","MnbZT3B.27","MnbZT15B.28","MnbZT5B.21","MnbZT13A.33","M2nbZT1C.27","M2nbZT1B.24","MnbZT13B.27","MnbZT7A.23","MnbZT1B.27")
for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/Negbin/",Allmod[i],"/"))
  setwd(Dir)
  mod_ALLneg[[CODE]] <- get(load(file = list.files(pattern=".rda")))
}

  
# Predictions for binomial model 
mod_ALLbin <- lapply(mod_ALLbin,function(x) {
  y <- predict(x)
  paste0(round(mean(y,na.rm = T),2),
         "\n(",round(min(y,na.rm = T),2)," - ",
         round(max(y,na.rm = T),2),")")
})
mod_ALLbin <- do.call(rbind,mod_ALLbin)


# Predictions for ZT model where mortality occurs only !
mod_ALLneg <- lapply(mod_ALLneg,function(x) {
  y <- predict(x)
  paste0(round(mean(y,na.rm = T)/10,2),
         "\n(",round(min(y,na.rm = T)/10,2)," - ",
         round(max(y,na.rm = T)/10,2),")")
})
mod_ALLneg <- do.call(rbind,mod_ALLneg)

# individual mortality 
dftreeTot <- lapply(dftree,function(x) {
  paste0(round(nrow(x[x$treestatus_th=="4",])/nrow(x)*100,2))
  })
dftreeTot <- do.call(rbind,dftreeTot)
  


# plot mortality for all ! 
dfplotTotall <- lapply(dfplot,function(x) {
  paste0(round(mean(x$sp.mortality.plot.rate.yr,na.rm = T)/10,2),
         "\n(",round(min(x$sp.mortality.plot.rate.yr,na.rm = T)/10,2)," - ",
  round(max(x$sp.mortality.plot.rate.yr,na.rm = T)/10,2),")")
})
dfplotTotall <- do.call(rbind,dfplotTotall)

nrow(x[x$sp.mortality.plot.rate.yr!="0",])
# need to plot mortality for plots where it occurs only !!!
dfplotTot <- lapply(dfplot,function(x) {
  paste0(round(mean(x[x$sp.mortality.plot.rate.yr!="0","sp.mortality.plot.rate.yr"],na.rm = T)/10,2),
         "\n(",round(min(x[x$sp.mortality.plot.rate.yr!="0","sp.mortality.plot.rate.yr"],na.rm = T)/10,2)," - ",
         round(max(x[x$sp.mortality.plot.rate.yr!="0","sp.mortality.plot.rate.yr"],na.rm = T)/10,2),")")
})
dfplotTot <- do.call(rbind,dfplotTot)
dfplotTotall
dfplotTot
dftreeTot
mod_ALLbin
mod_ALLneg

test <- c(mod_ALLneg[1:13],"POPNIG" = NA ,mod_ALLneg[14:19])
test2 <- cbind(dftreeTot,dfplotTotall,dfplotTot,test,mod_ALLbin)
colnames(test2) <- c("Observed % of dead trees\n(individual scale)",
                     "Observed mortality in all plots \n(mean% dead trees/yr/plot)",
                     "Observed mortality in plots where mortality occurs \n(mean% dead trees/yr/plot)",
                     "Predicted mortality in plots where mortality occurs \n(mean% dead trees/yr/plot)",
                     "Predicted mortality occurrence between census interval (probability)")


write.csv(test2,file = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/mortality.summary.csv"))


  