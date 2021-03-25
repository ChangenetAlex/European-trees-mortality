### Script sylvain 
rm(list = ls())
gc()
require(ade4)
library(sjPlot) 
library(lattice)
library(latticeExtra)
library(stringr)
library(data.table)
library(pscl)
library(MASS)
library(lme4)
library("glmmTMB")
library("bbmle")
library(piecewiseSEM)
library(pgirmess)
library(MuMIn)
library(spaMM)
library(dplyr)
library(plyr)


Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
Allcode <- c("FAGSYL","PINHAL","PINSYL","QUEILE","ABIALB","QUEROB")
Allseuil <- c(0.7,0.7,0.8,0.8,0.7,0.8)
i = 1

# Here we kept only the variables we wanted 
for (i in 1:length(Allcode)){
    Dir = (paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP")) # Directory 
    setwd(Dir)
    assign(paste0("dfplot",Allcode[i]),readRDS(paste0("dfplot",Allcode[i],Allseuil[i],".rds"))) #Base de données plot
    assign(paste0("dfplot2",Allcode[i]),get(paste0("dfplot",Allcode[i]))[,c("plotcode","latitude","longitude","speciesrichness","binomial","country","sp.mortality.plot.rate","mortality.plot.rate","sp.mortality.plot.rate.yr","sp.mortality.plot.count.yr","yearsbetweensurveys","mean_spei12")]) #Base de données plot
    assign(paste0("dfplot5",Allcode[i]),get(paste0("dfplot2",Allcode[i]))[get(paste0("dfplot2",Allcode[i]))$sp.mortality.plot.rate.yr>0 & !is.na(get(paste0("dfplot2",Allcode[i]))$sp.mortality.plot.rate.yr),]) #Base de données plot
    assign(paste0("dfplot3",Allcode[i]),get(paste0("dfplot2",Allcode[i]))[!is.na(get(paste0("dfplot2",Allcode[i]))$sp.mortality.plot.rate.yr),])
    }

# For each species, transform continuous mortality and latitude as class variables made of three levels. 

### ABIALB ###
dfplot3ABIALB$cut.latitude <- cut(dfplot3ABIALB$latitude, quantile(dfplot3ABIALB$latitude,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3ABIALB$cut.latitude) <- c(1,2,3)
# Mortality without excluding zero
dfplot3ABIALB[dfplot3ABIALB$sp.mortality.plot.rate.yr==0,"cut.sp.mortality.plot.rate.yr"] <- 0 # All alive plots
# Next line to split all dead plots in three quantiles 
dfplot3ABIALB[dfplot3ABIALB$sp.mortality.plot.rate.yr>0,"cut.sp.mortality.plot.rate.yr"] <- cut(dfplot3ABIALB[dfplot3ABIALB$sp.mortality.plot.rate.yr>0,"sp.mortality.plot.rate.yr"], quantile(dfplot3ABIALB[dfplot3ABIALB$sp.mortality.plot.rate.yr>0,"sp.mortality.plot.rate.yr"],na.rm=T,probs = c(0,0.333,0.666,1)))
#levels(as.factor(dfplot3ABIALB$cut.sp.mortality.plot.rate.yr)) <- c(1,2,3)
summary(as.factor(dfplot3ABIALB$cut.sp.mortality.plot.rate.yr))
dfplot3ABIALB[which(dfplot3ABIALB$cut.sp.mortality.plot.rate.yr==0),"sp.mortality.plot.rate.yr"]
#Idem SPEI
dfplot3ABIALB$cut.mean_spei12 <- cut(dfplot3ABIALB$mean_spei12, quantile(dfplot3ABIALB$mean_spei12,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3ABIALB$cut.mean_spei12) <- c(1,2,3)


### FAGSYL ###
dfplot3FAGSYL$cut.latitude <- cut(dfplot3FAGSYL$latitude, quantile(dfplot3FAGSYL$latitude,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3FAGSYL$cut.latitude) <- c(1,2,3)
# Mortality without excluding zero
dfplot3FAGSYL[dfplot3FAGSYL$sp.mortality.plot.rate.yr==0,"cut.sp.mortality.plot.rate.yr"] <- 0 # All alive plots
# Next line to split all dead plots in three quantiles 
dfplot3FAGSYL[dfplot3FAGSYL$sp.mortality.plot.rate.yr>0,"cut.sp.mortality.plot.rate.yr"] <- cut(dfplot3FAGSYL[dfplot3FAGSYL$sp.mortality.plot.rate.yr>0,"sp.mortality.plot.rate.yr"], quantile(dfplot3FAGSYL[dfplot3FAGSYL$sp.mortality.plot.rate.yr>0,"sp.mortality.plot.rate.yr"],na.rm=T,probs = c(0,0.333,0.666,1)))
#levels(as.factor(dfplot3FAGSYL$cut.sp.mortality.plot.rate.yr)) <- c(1,2,3)
summary(as.factor(dfplot3FAGSYL$cut.sp.mortality.plot.rate.yr))
dfplot3FAGSYL[which(dfplot3FAGSYL$cut.sp.mortality.plot.rate.yr==0),"sp.mortality.plot.rate.yr"]
#Idem SPEI
dfplot3FAGSYL$cut.mean_spei12 <- cut(dfplot3FAGSYL$mean_spei12, quantile(dfplot3FAGSYL$mean_spei12,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3FAGSYL$cut.mean_spei12) <- c(1,2,3)



### PINHAL ###
dfplot3PINHAL$cut.latitude <- cut(dfplot3PINHAL$latitude, quantile(dfplot3PINHAL$latitude,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3PINHAL$cut.latitude) <- c(1,2,3)
# Mortality without excluding zero
dfplot3PINHAL[dfplot3PINHAL$sp.mortality.plot.rate.yr==0,"cut.sp.mortality.plot.rate.yr"] <- 0 # All alive plots
# Next line to split all dead plots in three quantiles 
dfplot3PINHAL[dfplot3PINHAL$sp.mortality.plot.rate.yr>0,"cut.sp.mortality.plot.rate.yr"] <- cut(dfplot3PINHAL[dfplot3PINHAL$sp.mortality.plot.rate.yr>0,"sp.mortality.plot.rate.yr"], quantile(dfplot3PINHAL[dfplot3PINHAL$sp.mortality.plot.rate.yr>0,"sp.mortality.plot.rate.yr"],na.rm=T,probs = c(0,0.333,0.666,1)))
#levels(as.factor(dfplot3PINHAL$cut.sp.mortality.plot.rate.yr)) <- c(1,2,3)
summary(as.factor(dfplot3PINHAL$cut.sp.mortality.plot.rate.yr))
dfplot3PINHAL[which(dfplot3PINHAL$cut.sp.mortality.plot.rate.yr==0),"sp.mortality.plot.rate.yr"]
#Idem SPEI
dfplot3PINHAL$cut.mean_spei12 <- cut(dfplot3PINHAL$mean_spei12, quantile(dfplot3PINHAL$mean_spei12,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3PINHAL$cut.mean_spei12) <- c(1,2,3)


### PINSYL ###
dfplot3PINSYL$cut.latitude <- cut(dfplot3PINSYL$latitude, quantile(dfplot3PINSYL$latitude,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3PINSYL$cut.latitude) <- c(1,2,3)
# Mortality without excluding zero
dfplot3PINSYL[dfplot3PINSYL$sp.mortality.plot.rate.yr==0,"cut.sp.mortality.plot.rate.yr"] <- 0 # All alive plots
# Next line to split all dead plots in three quantiles 
dfplot3PINSYL[dfplot3PINSYL$sp.mortality.plot.rate.yr>0,"cut.sp.mortality.plot.rate.yr"] <- cut(dfplot3PINSYL[dfplot3PINSYL$sp.mortality.plot.rate.yr>0,"sp.mortality.plot.rate.yr"], quantile(dfplot3PINSYL[dfplot3PINSYL$sp.mortality.plot.rate.yr>0,"sp.mortality.plot.rate.yr"],na.rm=T,probs = c(0,0.333,0.666,1)))
#levels(as.factor(dfplot3PINSYL$cut.sp.mortality.plot.rate.yr)) <- c(1,2,3)
summary(as.factor(dfplot3PINSYL$cut.sp.mortality.plot.rate.yr))
dfplot3PINSYL[which(dfplot3PINSYL$cut.sp.mortality.plot.rate.yr==0),"sp.mortality.plot.rate.yr"]
#Idem SPEI
dfplot3PINSYL$cut.mean_spei12 <- cut(dfplot3PINSYL$mean_spei12, quantile(dfplot3PINSYL$mean_spei12,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3PINSYL$cut.mean_spei12) <- c(1,2,3)


### QUEROB ###
dfplot3QUEROB$cut.latitude <- cut(dfplot3QUEROB$latitude, quantile(dfplot3QUEROB$latitude,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3QUEROB$cut.latitude) <- c(1,2,3)
# Mortality without excluding zero
dfplot3QUEROB[dfplot3QUEROB$sp.mortality.plot.rate.yr==0,"cut.sp.mortality.plot.rate.yr"] <- 0 # All alive plots
# Next line to split all dead plots in three quantiles 
dfplot3QUEROB[dfplot3QUEROB$sp.mortality.plot.rate.yr>0,"cut.sp.mortality.plot.rate.yr"] <- cut(dfplot3QUEROB[dfplot3QUEROB$sp.mortality.plot.rate.yr>0,"sp.mortality.plot.rate.yr"], quantile(dfplot3QUEROB[dfplot3QUEROB$sp.mortality.plot.rate.yr>0,"sp.mortality.plot.rate.yr"],na.rm=T,probs = c(0,0.333,0.666,1)))
#levels(as.factor(dfplot3QUEROB$cut.sp.mortality.plot.rate.yr)) <- c(1,2,3)
summary(as.factor(dfplot3QUEROB$cut.sp.mortality.plot.rate.yr))
dfplot3QUEROB[which(dfplot3QUEROB$cut.sp.mortality.plot.rate.yr==0),"sp.mortality.plot.rate.yr"]
#Idem SPEI
dfplot3QUEROB$cut.mean_spei12 <- cut(dfplot3QUEROB$mean_spei12, quantile(dfplot3QUEROB$mean_spei12,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3QUEROB$cut.mean_spei12) <- c(1,2,3)


### QUEILE ###
dfplot3QUEILE$cut.latitude <- cut(dfplot3QUEILE$latitude, quantile(dfplot3QUEILE$latitude,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3QUEILE$cut.latitude) <- c(1,2,3)
# Mortality without excluding zero
dfplot3QUEILE[dfplot3QUEILE$sp.mortality.plot.rate.yr==0,"cut.sp.mortality.plot.rate.yr"] <- 0 # All alive plots
# Next line to split all dead plots in three quantiles 
dfplot3QUEILE[dfplot3QUEILE$sp.mortality.plot.rate.yr>0,"cut.sp.mortality.plot.rate.yr"] <- cut(dfplot3QUEILE[dfplot3QUEILE$sp.mortality.plot.rate.yr>0,"sp.mortality.plot.rate.yr"], quantile(dfplot3QUEILE[dfplot3QUEILE$sp.mortality.plot.rate.yr>0,"sp.mortality.plot.rate.yr"],na.rm=T,probs = c(0,0.333,0.666,1)))
#levels(as.factor(dfplot3QUEILE$cut.sp.mortality.plot.rate.yr)) <- c(1,2,3)
summary(as.factor(dfplot3QUEILE$cut.sp.mortality.plot.rate.yr))
dfplot3QUEILE[which(dfplot3QUEILE$cut.sp.mortality.plot.rate.yr==0),"sp.mortality.plot.rate.yr"]
#Idem SPEI
dfplot3QUEILE$cut.mean_spei12 <- cut(dfplot3QUEILE$mean_spei12, quantile(dfplot3QUEILE$mean_spei12,na.rm=T,probs = c(0,0.333,0.666,1)))
levels(dfplot3QUEILE$cut.mean_spei12) <- c(1,2,3)


## For each species, sample randomly three plots for each combination of modality. 
# (3 levels for latitude, 3 levels for mortality => 9 modality * 3 plots = 27 rows.)


table(dfplot3ABIALB$cut.latitude,dfplot3ABIALB$cut.sp.mortality.plot.rate.yr,dfplot3ABIALB$cut.mean_spei12)
l=1
for (i in 1:3){
  for (j in 0:3){
    for (k in 1:3){
    assign(paste0("test",l),sample_n(dfplot3ABIALB[dfplot3ABIALB$cut.latitude==i & dfplot3ABIALB$cut.sp.mortality.plot.rate.yr==j & dfplot3ABIALB$cut.mean_spei12==k,],3,replace = T))
    l = l+1
    k=k+1
    }
    j=j+1
  }
  i=i+1
}
dfplot4ABIALB <- rbind(test1, test2, test3, test4, test5, test6, test7, 
                       test8, test9, test10, test11, test12, test13, test14, 
                       test15, test16, test17, test18, test19, test20, test21, 
                       test22, test23, test24, test25, test26, test27, test28, 
                       test29, test30, test31, test32, test33, test34, test35, 
                       test36)



table(dfplot3FAGSYL$cut.latitude,dfplot3FAGSYL$cut.sp.mortality.plot.rate.yr,dfplot3FAGSYL$cut.mean_spei12)
l=1
for (i in 1:3){
  for (j in 0:3){
    for (k in 1:3){
      assign(paste0("test",l),sample_n(dfplot3FAGSYL[dfplot3FAGSYL$cut.latitude==i & dfplot3FAGSYL$cut.sp.mortality.plot.rate.yr==j & dfplot3FAGSYL$cut.mean_spei12==k,],3,replace = T))
      l = l+1
      k=k+1
    }
    j=j+1
  }
  i=i+1
}
dfplot4FAGSYL <- rbind(test1, test2, test3, test4, test5, test6, test7, 
                       test8, test9, test10, test11, test12, test13, test14, 
                       test15, test16, test17, test18, test19, test20, test21, 
                       test22, test23, test24, test25, test26, test27, test28, 
                       test29, test30, test31, test32, test33, test34, test35, 
                       test36)


table(dfplot3PINHAL$cut.latitude,dfplot3PINHAL$cut.sp.mortality.plot.rate.yr,dfplot3PINHAL$cut.mean_spei12)
l=1
for (i in 1:3){
  for (j in 0:3){
    for (k in 1:3){
      assign(paste0("test",l),sample_n(dfplot3PINHAL[dfplot3PINHAL$cut.latitude==i & dfplot3PINHAL$cut.sp.mortality.plot.rate.yr==j & dfplot3PINHAL$cut.mean_spei12==k,],3,replace = T))
      l = l+1
      k=k+1
    }
    j=j+1
  }
  i=i+1
}
dfplot4PINHAL <- rbind(test1, test2, test3, test4, test5, test6, test7, 
                       test8, test9, test10, test11, test12, test13, test14, 
                       test15, test16, test17, test18, test19, test20, test21, 
                       test22, test23, test24, test25, test26, test27, test28, 
                       test29, test30, test31, test32, test33, test34, test35, 
                       test36)


table(dfplot3PINSYL$cut.latitude,dfplot3PINSYL$cut.sp.mortality.plot.rate.yr,dfplot3PINSYL$cut.mean_spei12)
l=1
for (i in 1:3){
  for (j in 0:3){
    for (k in 1:3){
      assign(paste0("test",l),sample_n(dfplot3PINSYL[dfplot3PINSYL$cut.latitude==i & dfplot3PINSYL$cut.sp.mortality.plot.rate.yr==j & dfplot3PINSYL$cut.mean_spei12==k,],3,replace = T))
      l = l+1
      k=k+1
    }
    j=j+1
  }
  i=i+1
}
dfplot4PINSYL <- rbind(test1, test2, test3, test4, test5, test6, test7, 
                       test8, test9, test10, test11, test12, test13, test14, 
                       test15, test16, test17, test18, test19, test20, test21, 
                       test22, test23, test24, test25, test26, test27, test28, 
                       test29, test30, test31, test32, test33, test34, test35, 
                       test36)


table(dfplot3QUEILE$cut.latitude,dfplot3QUEILE$cut.sp.mortality.plot.rate.yr,dfplot3QUEILE$cut.mean_spei12)
l=1
for (i in 1:3){
  for (j in 0:3){
    for (k in 1:3){
      assign(paste0("test",l),sample_n(dfplot3QUEILE[dfplot3QUEILE$cut.latitude==i & dfplot3QUEILE$cut.sp.mortality.plot.rate.yr==j & dfplot3QUEILE$cut.mean_spei12==k,],3,replace = T))
      l = l+1
      k=k+1
    }
    j=j+1
  }
  i=i+1
}
dfplot4QUEILE <- rbind(test1, test2, test3, test4, test5, test6, test7, 
                       test8, test9, test10, test11, test12, test13, test14, 
                       test15, test16, test17, test18, test19, test20, test21, 
                       test22, test23, test24, test25, test26, test27, test28, 
                       test29, test30, test31, test32, test33, test34, test35, 
                       test36)


table(dfplot3QUEROB$cut.latitude,dfplot3QUEROB$cut.sp.mortality.plot.rate.yr,dfplot3QUEROB$cut.mean_spei12)
l=1
for (i in 1:3){
  for (j in 0:3){
    for (k in 1:3){
      assign(paste0("test",l),sample_n(dfplot3QUEROB[dfplot3QUEROB$cut.latitude==i & dfplot3QUEROB$cut.sp.mortality.plot.rate.yr==j & dfplot3QUEROB$cut.mean_spei12==k,],3,replace = T))
      l = l+1
      k=k+1
    }
    j=j+1
  }
  i=i+1
}
dfplot4QUEROB <- rbind(test1, test2, test3, test4, test5, test6, test7, 
                       test8, test9, test10, test11, test12, test13, test14, 
                       test15, test16, test17, test18, test19, test20, test21, 
                       test22, test23, test24, test25, test26, test27, test28, 
                       test29, test30, test31, test32, test33, test34, test35, 
                       test36)

saveRDS(dfplot4ABIALB, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/dfplot4ABIALB.rds"))
saveRDS(dfplot4FAGSYL, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/dfplot4FAGSYL.rds"))
saveRDS(dfplot4PINHAL, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/dfplot4PINHAL.rds"))
saveRDS(dfplot4PINSYL, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/dfplot4PINSYL.rds"))
saveRDS(dfplot4QUEILE, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/dfplot4QUEILE.rds"))
saveRDS(dfplot4QUEROB, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/dfplot4QUEROB.rds"))

####################################
### Partie 2 multispecies sites ! ##
####################################


## 4 species 
df4sp1 <- join_all(list(dfplot3ABIALB,dfplot3FAGSYL,dfplot3PINSYL,dfplot3QUEROB), by='plotcode', type='inner') # 26 
colnames(df4sp1)
df4sp1 <- df4sp1[df4sp1[,9]!=0 | df4sp1[,23]!=0 | df4sp1[,37]!=0 | df4sp1[,51]!=0 & !is.na(df4sp1[,9]) & !is.na(df4sp1[,23]) & !is.na(df4sp1[,37]) & !is.na(df4sp1[,51]),] #8 plots only
df4sp1[,c(9,23,37,51)]
# 8 plots
saveRDS(df4sp1, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/df4sp.rds"))

## 3 species 
library(plyr)
elem <- c("dfplot3ABIALB","dfplot3FAGSYL","dfplot3PINHAL","dfplot3PINSYL","dfplot3QUEILE","dfplot3QUEROB")
t(combn(elem,3))

# These lines are all combinations for which plots containing all three species were found. 
AAAA <- join_all(list(dfplot3ABIALB,dfplot3FAGSYL,dfplot3PINSYL),by='plotcode', type='inner')
colnames(AAAA)
df3sp1 <- AAAA[AAAA[,9]!=0 | AAAA[,23]!=0 | AAAA[,37]!=0 & !is.na(AAAA[,9]) & !is.na(AAAA[,23]) & !is.na(AAAA[,37]),]
df3sp1[,c(9,23,37)] #125

AAAA <- join_all(list(dfplot3ABIALB,dfplot3FAGSYL,dfplot3QUEROB),by='plotcode', type='inner')
colnames(AAAA)
df3sp2 <- AAAA[AAAA[,9]!=0 | AAAA[,23]!=0 | AAAA[,37]!=0 & !is.na(AAAA[,9]) & !is.na(AAAA[,23]) & !is.na(AAAA[,37]),] #8 plots only
df3sp2[,c(9,23,37)] #41


AAAA <- join_all(list(dfplot3ABIALB,dfplot3PINSYL,dfplot3QUEROB),by='plotcode', type='inner')
colnames(AAAA)
df3sp3 <- AAAA[AAAA[,9]!=0 | AAAA[,23]!=0 | AAAA[,37]!=0 & !is.na(AAAA[,9]) & !is.na(AAAA[,23]) & !is.na(AAAA[,37]),] #8 plots only
df3sp3[,c(9,23,37)] #23


AAAA <- join_all(list(dfplot3FAGSYL,dfplot3PINSYL,dfplot3QUEILE),by='plotcode', type='inner')
colnames(AAAA)
df3sp4 <- AAAA[AAAA[,9]!=0 | AAAA[,23]!=0 | AAAA[,37]!=0 & !is.na(AAAA[,9]) & !is.na(AAAA[,23]) & !is.na(AAAA[,37]),] #8 plots only
df3sp4[,c(9,23,37)] #7

AAAA <- join_all(list(dfplot3FAGSYL,dfplot3PINSYL,dfplot3QUEROB),by='plotcode', type='inner')
colnames(AAAA)
df3sp5 <- AAAA[AAAA[,9]!=0 | AAAA[,23]!=0 | AAAA[,37]!=0 & !is.na(AAAA[,9]) & !is.na(AAAA[,23]) & !is.na(AAAA[,37]),] #8 plots only
df3sp5[,c(9,23,37)] #117


AAAA <- join_all(list(dfplot3FAGSYL,dfplot3QUEILE,dfplot3QUEROB),by='plotcode', type='inner')
colnames(AAAA)
df3sp6 <- AAAA[AAAA[,9]!=0 | AAAA[,23]!=0 | AAAA[,37]!=0 & !is.na(AAAA[,9]) & !is.na(AAAA[,23]) & !is.na(AAAA[,37]),] #8 plots only
df3sp6[,c(9,23,37)] #1


AAAA <- join_all(list(dfplot3PINHAL,dfplot3PINSYL,dfplot3QUEILE),by='plotcode', type='inner')
colnames(AAAA)
df3sp7 <- AAAA[AAAA[,9]!=0 | AAAA[,23]!=0 | AAAA[,37]!=0 & !is.na(AAAA[,9]) & !is.na(AAAA[,23]) & !is.na(AAAA[,37]),] #8 plots only
df3sp7[,c(9,23,37)] #25

AAAA <- join_all(list(dfplot3PINHAL,dfplot3QUEILE,dfplot3QUEROB),by='plotcode', type='inner')
colnames(AAAA)
df3sp8 <- AAAA[AAAA[,9]!=0 | AAAA[,23]!=0 | AAAA[,37]!=0 & !is.na(AAAA[,9]) & !is.na(AAAA[,23]) & !is.na(AAAA[,37]),] #8 plots only
df3sp8[,c(9,23,37)] #2


df3sp <- rbind(df3sp1,df3sp2,df3sp3,df3sp4,df3sp5,df3sp6,df3sp7,df3sp8)
saveRDS(df3sp, paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Fundiv.Drive/Species/Syl/df3sp.rds"))







    