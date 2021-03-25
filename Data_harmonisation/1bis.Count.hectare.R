rm(list = ls())
gc()

# Script to calculate mortality and recruitment as a count per hectare and not as a proportion compare to the total weight of the plot 


Dir <- c("~/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
library(ggplot2)
require(scales)
library(dplyr)
library(raster)
library(rgdal)
library(rworldmap)
library(parallel)
#These are the database with which I work 
tfinal.biotic2 <- readRDS("our-data/tfinal.biotic.Sep2019.RDS")
tfinal.biotic <- tfinal.biotic2[tfinal.biotic2$country!="FR",] 
tfinal.biotic <- tfinal.biotic2[tfinal.biotic2$plotcode=="ES100017A1",] 

dat.fundiv <- split(tfinal.biotic[,],as.character(tfinal.biotic$plotcode))

# split my list into several list 
dataSplit <- split(dat.fundiv, rep(1:1001, each = 93))
NamesDF <- paste0("dat.fundiv",1:length(names(dataSplit)))
dir.create(path = "~/SynologyDrive/FUNDIV - NFI - Europe/our-data/data.fundiv")
Ncore <- 24
mcmapply(function(x,y){
  saveRDS(x,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/data.fundiv/",y,".rds"))
  },x=dataSplit,y=NamesDF,mc.cores=Ncore)
rm(list = c("dat.fundiv","dataSplit"))
gc()

x <- 0

lapply(NamesDF, function(p){
  x <<- x+1
  dat.fundiv <- readRDS(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/data.fundiv/",p,".rds"))
  data <- do.call(rbind,mclapply(dat.fundiv,function(DF){
  # idp
  i=DF[1,"plotcode"]
  dat=matrix(NA,nrow(tfinal.biotic[tfinal.biotic$plotcode==i,]),5)
  dat[,5]=c(1:nrow(tfinal.biotic))[tfinal.biotic$plotcode==i]
  if(nrow(DF)>0) {
    S=unique(DF$speciesid)
    S=S[!is.na(S)]
    if(length(S)>0)
      for(s in S) {
        df=DF[!is.na(DF$speciesid) & DF$speciesid==s,]
        # dead trees
        df.m=df[df$treestatus=="4",]
        df.r=df[df$treestatus=="1",]
        dat[tfinal.biotic$speciesid[tfinal.biotic$plotcode==i]==s,2]<-sum(df[df$treestatus_th%in%c("4","2") & !is.na(df$weight1),"weight1"])
        dat[tfinal.biotic$speciesid[tfinal.biotic$plotcode==i]==s,4]<-sum(df[df$treestatus_th%in%c("1","2") & !is.na(df$weight2),"weight2"])
        if(nrow(df.m)==0) dat[tfinal.biotic$speciesid[tfinal.biotic$plotcode==i]==s,1]=0
        if(nrow(df.m)>0) dat[tfinal.biotic$speciesid[tfinal.biotic$plotcode==i]==s,1]=sum(df.m$weight1)
        if(nrow(df.r)==0) dat[tfinal.biotic$speciesid[tfinal.biotic$plotcode==i]==s,3]=0
        if(nrow(df.r)>0) dat[tfinal.biotic$speciesid[tfinal.biotic$plotcode==i]==s,3]=sum(df.r$weight2)
      }
  }
  as.data.frame(dat)},mc.cores=Ncore,mc.silent=F))
  rm(dat.fundiv)
  saveRDS(data,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/data.fundiv/DF.",p,".rds"))
  rm(data)
  print(x)
  gc()
})

data <- list()
for (i in 1:length(NamesDF)){
  data[[i]] <- readRDS(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/data.fundiv/DF.",NamesDF[i],".rds"))
}
data <- do.call(rbind,data)

sp.mortality.ha.weight=as.numeric(data[order(data[,5]),1]) # weight of all dead trees in the plot
sp.SUM.ha.weight1=as.numeric(data[order(data[,5]),2])      # sum of weight1 of all ingrowth and dead trees
sp.recruitment.ha.weight=as.numeric(data[order(data[,5]),3]) # idem recrut 
sp.SUM.ha.weight2=as.numeric(data[order(data[,5]),4])      # sum weight 2
tfinal.biotic[,"sp.mortality.ha.weight"]=sp.mortality.ha.weight
tfinal.biotic[,"sp.SUM.ha.weight1"]=sp.SUM.ha.weight1
tfinal.biotic[,"sp.recruitment.ha.weight"]=sp.recruitment.ha.weight
tfinal.biotic[,"sp.SUM.ha.weight2"]=sp.SUM.ha.weight2
#tfinal.biotic[,"sp.SUM.ha.mean.weight"] <- (tfinal.biotic[,"sp.SUM.ha.weight2"]+tfinal.biotic[,"sp.SUM.ha.weight1"])/2

saveRDS(tfinal.biotic,file="~/SynologyDrive/FUNDIV - NFI - Europe/our-data/tfinal.biotic.August2020.rds")


