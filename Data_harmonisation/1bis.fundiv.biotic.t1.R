### Edit 08/03/2018 : add the calculations based on the first inventory for biotic variables
Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
load("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/tfinal.biotic.RData") # Base de données avec nouvelles variables biotiques
library(ggplot2)
require(scales)
library(dplyr)
library(raster)
library(rgdal)
library(rworldmap)
library(parallel)

ftp.short <- tfinal.biotic[tfinal.biotic$country!="FR",]
treefinalbis <- tfinal.biotic[tfinal.biotic$country=="FR",]
dat.fundiv=split(ftp.short[,],as.character(ftp.short$plotcode))
####BA.ha.plot##### and BA.PLOT 1
data <- do.call(rbind, mclapply(dat.fundiv,function(df){
  
  # idp
  i=df[1,"plotcode"]
  dat=matrix(NA,nrow(ftp.short[ftp.short$plotcode==i,]),2)
  dat[,2]=c(1:nrow(ftp.short))[ftp.short$plotcode==i]
  
  if(nrow(df)>0) {
    
    # set to NA dead trees BA
    #df$ba_ha1[df$treestatus=="4"]=NA Ne pas enlever les arbres mort du calcul, car dead au second inventaires (148 individus au total)
    ba.ha=df$ba_ha1
    res=sum(ba.ha,na.rm=T)
    if(length(res)==0 | res==0) res=NA
    dat[,1]=as.numeric(res)
  }
  as.data.frame(dat)
},mc.cores=12,mc.silent=T))
BA.ha.plot.1=as.numeric(data[order(data[,2]),1])

# 08/03/2018 : The new one is based at t1 => account for the dead before they died. But not for the recrut
BAnei.1=BA.ha.plot.1-ftp.short$ba_ha1


#########################################
### BAI.J.plot (ba.sp.plot) m²/ha/year ##  bis
#########################################
# plot specific BA (m2/ha)
data <- do.call(rbind, mclapply(dat.fundiv,function(DF){
  
  # idp
  i=DF[1,"plotcode"]
  dat=matrix(NA,nrow(ftp.short[ftp.short$plotcode==i,]),2)
  dat[,2]=c(1:nrow(ftp.short))[ftp.short$plotcode==i]
  
  if(nrow(DF)>0) {
    
    S=unique(DF$speciesid)
    S=S[!is.na(S)]
    if(length(S)>0)
      for(s in S) {
        
        df=DF[!is.na(DF$speciesid) & DF$speciesid==s,]
        BAIpj.s=sum(df$bachange_ha_yr,na.rm=T)
        if(length(BAIpj.s)==0 | BAIpj.s==0) BAIpj.s=0 #rajout march 0 instead of na
        dat[ftp.short$speciesid[ftp.short$plotcode==i]==s,1]=BAIpj.s
      }
  }
  
  as.data.frame(dat)
  
},mc.cores=12,mc.silent=T))
BAIj.plot.bis=as.numeric(data[order(data[,2]),1])
summary(BAIj.plot.bis)

##### neighbouring conspecific tree BAI (m2/ha)
BAIneicon.bis=BAIj.plot.bis-ftp.short$bachange_ha_yr


####################################
### BA.J.plot (ba.sp.plot) m²/ha  ## t1
####################################
# plot specific BA (m2/ha) au plot
data <- do.call(rbind, mclapply(dat.fundiv,function(DF){
  
  # idp
  i=DF[1,"plotcode"]
  dat=matrix(NA,nrow(ftp.short[ftp.short$plotcode==i,]),2)
  dat[,2]=c(1:nrow(ftp.short))[ftp.short$plotcode==i]
  
  if(nrow(DF)>0) {
    
    S=unique(DF$speciesid)
    S=S[!is.na(S)]
    if(length(S)>0)
      for(s in S) {
        
        df=DF[!is.na(DF$speciesid) & DF$speciesid==s,]
        BApj.s=sum(df$ba_ha1,na.rm=T)
        dat[ftp.short$speciesid[ftp.short$plotcode==i]==s,1]=BApj.s
        
      }
  }
  
  as.data.frame(dat)
  
},mc.cores=12,mc.silent=T))
BAj.plot.1=as.numeric(data[order(data[,2]),1])
summary(BAj.plot.1)

# neighbouring conspecific tree BA (m2/ha)
BAneicon.1=BAj.plot.1-ftp.short$ba_ha1


ftp.short[,"BA.ha.plot.1"]=BA.ha.plot.1 # To run 
ftp.short[,"BAnei.1"]=BAnei.1 # To run 
ftp.short[,"BAIj.plot.bis"]=BAIj.plot.bis
ftp.short[,"BAIneicon.bis"]=BAIneicon.bis
ftp.short[,"BAj.plot.1"]=BAj.plot.1
ftp.short[,"BAneicon.1"]=BAneicon.1

# Edit on the 18/03/2018 : For the 6 first new columns, 
# it has to be the same as the old ones for french inventory since there is just one. Otherwise it is NA. 

treefinalbis$BA.ha.plot.1 <- treefinalbis$BA.plot.ha
treefinalbis$BAnei.1 <- treefinalbis$BAnei
treefinalbis$BAIj.plot.bis <- treefinalbis$BAIj.plot
treefinalbis$BAIneicon.bis <- treefinalbis$BAIneicon
treefinalbis$BAj.plot.1 <- treefinalbis$BAj.plot
treefinalbis$BAneicon.1 <- treefinalbis$BAneicon

tfinal.biotic <- rbind(ftp.short,treefinalbis)
#Check 23900 plot with code gest in dfplot
summary(as.factor(tfinal.biotic[,"management1"]))
summary(as.factor(tfinal.biotic[,"management2"]))
summary(as.factor(tfinal.biotic[,"gest"]))


# Ajout metrics competition and mortality.count
tfinal.biotic$BA.O.plot <- tfinal.biotic$BA.plot.ha-tfinal.biotic$BAj.plot # Ajout metrics competition
tfinal.biotic$BAI.O.plot <- tfinal.biotic$BAI.plot-tfinal.biotic$BAIj.plot #Idem croissance
tfinal.biotic$BA.O.plot.1 <- tfinal.biotic$BA.ha.plot.1 - tfinal.biotic$BAj.plot.1 # Ajout metrics competition inv1
tfinal.biotic$BAI.O.plot.1 <- tfinal.biotic$BAI.plot-tfinal.biotic$BAIj.plot.bis # Moins de zeros but idem 1
tfinal.biotic$mortality.plot.count <- round((tfinal.biotic$mortality.plot.rate*100),0)
tfinal.biotic$sp.mortality.plot.count <- round((tfinal.biotic$sp.mortality.plot.rate*100),0)
summary(as.factor(tfinal.biotic$mortality.plot.count))
summary(as.factor(tfinal.biotic$sp.mortality.plot.count))
# Check this is okay (none with inferior than 0, -0. or -1 values)
test <- tfinal.biotic[!is.na(tfinal.biotic$BA.O.plot.1)&tfinal.biotic$BA.O.plot.1<0,] 
test2 <- tfinal.biotic[!is.na(tfinal.biotic$BA.O.plot)&tfinal.biotic$BA.O.plot<=-0.00000000000000001,] #16887 < 0.0000000000000001

# Save it as the final database

save(tfinal.biotic,file="/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/tfinal.biotic.RData")

##########################
### Edit 03/06/2018 : add the number of trees and the absolute dead mass (need to add information for the french base)
##########################
rm(list = ls())
gc()
library(ggplot2)
require(scales)
library(dplyr)
library(raster)
library(rgdal)
library(rworldmap)
library(parallel)

Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
load("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/tfinal.biotic.June2018.RData") # Base de données avec nouvelles variables biotiques
load(paste0(Dir,"our-data/ifn.fr.metric.2018.RData")) #IFN with 3500 names corrected and all variables


# Modif noms des ifn.fr.2018 
#ifn.fr$treecode <- paste0("FR",ifn.fr$idp,"_",ifn.fr$a) # Give the right name to the french treecodes
#test <- semi_join(ifn.fr,tfinal.biotic,by="treecode") # What is well attributed = ALL = 637830

#tfinal.biotic$alive <- ifn.fr$alive[match(tfinal.biotic$treecode,ifn.fr$treecode,incomparables = NA)] # The two columns we need to evaluate mortality in french inventory
#tfinal.biotic$veget <- ifn.fr$veget[match(tfinal.biotic$treecode,ifn.fr$treecode,incomparables = NA)]
#tfinal.biotic <- tfinal.biotic[tfinal.biotic$plotcode=="FR2261",]

#tfinal.biotic <- rbind(tfinal.biotic[tfinal.biotic$plotcode=="DE12974_2",],tfinal.biotic[tfinal.biotic$plotcode=="FR300089",]) test

datDF=split(tfinal.biotic[,],as.character(tfinal.biotic$plotcode)) # split mon df

# Nbr of trees in each plot 

data <- do.call(rbind, mclapply(datDF,function(DF){
  i=DF[1,"plotcode"]
  dat=matrix(NA,nrow(tfinal.biotic[tfinal.biotic$plotcode==i,]),2)
  dat[,2]=c(1:nrow(tfinal.biotic))[tfinal.biotic$plotcode==i]
  if(nrow(DF)>0){
    dat[,1]=as.numeric(nrow(DF))} # The number of trees is including the number of line so that include the number of recruit, number of dead and ingrowth 
  as.data.frame(dat)},mc.cores=10,mc.silent=T))
treeNbr=as.numeric(data[order(data[,2]),1])


# Mortality (dead mass/ha/plot)

#### Calculate sum of basal area of dead 
####################################
### plot SP.mortality (% tree /ha)## ba
####################################
# plot specific mortality (% tree /ha): dens dead trees IFN3 / dens live trees IFN2
data <- do.call(rbind, mclapply(datDF,function(DF){
  # idp
  i=DF[1,"plotcode"]
  dat=matrix(NA,nrow(tfinal.biotic[tfinal.biotic$plotcode==i,]),2)
  dat[,2]=c(1:nrow(tfinal.biotic))[tfinal.biotic$plotcode==i]
  # remove 'arbres chablis' (dead or alive)
  if (grepl("FR",DF[,"plotcode"],fixed=TRUE)==T){
    chab=c(1:nrow(DF))[!is.na(DF$veget) & DF$veget=="A"]
    chab=chab[!is.na(chab)]
    if(length(chab)>0) DF=DF[-c(chab),]}
  if(nrow(DF)>0) {
    S=unique(DF$speciesid)
    S=S[!is.na(S)]
    if(length(S)>0)
      for(s in S) {
        df=DF[!is.na(DF$speciesid) & DF$speciesid==s,]
        # dead trees
        if (grepl("FR",DF[,"plotcode"],fixed=TRUE)==T){df.m=df[df$alive==F,]
        } else df.m=df[df$treestatus=="4",]
        
        if(nrow(df.m)==0) dat[tfinal.biotic$speciesid[tfinal.biotic$plotcode==i]==s,1]=0
        if(nrow(df.m)>0) {
          
          #dens.tot=sum(df$ba_ha1[!is.na(df$ba_ha1)]) # No need for it beca  use we want the absolute lost
          if (grepl("FR",DF[,"plotcode"],fixed=TRUE)==T){dens.m=sum(df.m$ba_ha2,na.rm = T)
          } else dens.m=sum(df.m$ba_ha1,na.rm = T)
          #res=dens.m/dens.tot	
          dat[tfinal.biotic$speciesid[tfinal.biotic$plotcode==i]==s,1]=as.numeric(dens.m)
        }
      }
  }
  as.data.frame(dat)
},mc.cores=10,mc.silent=T))
sp.mortality.plot.ba.ABS2=as.numeric(data[order(data[,2]),1])
summary(sp.mortality.plot.ba.ABS2)  

#tfinal.biotic[,"treeNbr"]=treeNbr
tfinal.biotic[,"sp.mortality.plot.ba.ABS2"]=summary(as.factor(sp.mortality.plot.ba.ABS2))
summary(as.factor(sp.mortality.plot.ba.ABS2))
summary(as.factor(tfinal.biotic$sp.mortality.plot.count))


save(tfinal.biotic,file="/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/tfinal.biotic.June2018.RData")
