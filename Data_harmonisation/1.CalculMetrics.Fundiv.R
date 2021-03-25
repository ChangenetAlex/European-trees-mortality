#Script to do the calculation again for all the variable that interest us with the fundivData
Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
library(ggplot2)
require(scales)
library(dplyr)
library(raster)
library(rgdal)
library(rworldmap)
library(parallel)
#These are the database with which I work 
#load(paste0(Dir,"our-data/tree_data_harmonise_fr.fundiv.indexes.RData"))
#load(paste0(Dir,"our-data/tree_data_IFN_france_with_competition-indices.bai_ha_forFUNDIVcomp.RData"))
load(paste0(Dir,"FunDivEUROPE_Inventory_data_7Oct/FunDivEUROPE_all_trees_mortality_rates_10cm.RData"))
#I need to check if they are the same as the one created (to detet any problemi the script)
ifn_plot <- read.csv("/home/achangenet/Documents/FUNDIV\ -\ NFI\ -\ Europe/FunDivEUROPE_Inventory_data_7Oct/FunDivEUROPE_plot_data.csv", sep=",")
ifn_aggregation_75 <- read.csv("/home/achangenet/Documents/FUNDIV\ -\ NFI\ -\ Europe/FunDivEUROPE_Inventory_data_7Oct/FunDivEUROPE_aggregation_75.csv", sep=",")
ifn_trees_10 <- read.csv("/home/achangenet/Documents/FUNDIV\ -\ NFI\ -\ Europe/FunDivEUROPE_Inventory_data_7Oct/FunDivEUROPE_all_trees_10cm.csv", sep=",")


load(paste0(Dir,"our-data/Fundiv_alltree_allmetric.RData")) #last version


####Remove white spaces: 
#Conversion as character
trim.trailing <- function (x) sub("\\s+$", "", x)
ifn_trees_10$plotcode <- as.character(ifn_trees_10$plotcode)
ifn_trees_10$plotcode <- trim.trailing(ifn_trees_10$plotcode)
ifn_aggregation_75$plotcode <- as.character(ifn_aggregation_75$plotcode)
ifn_aggregation_75$plotcode <- trim.trailing(ifn_aggregation_75$plotcode)
ifn_plot$plotcode <- as.character(ifn_plot$plotcode)
ifn_plot$plotcode <- trim.trailing(ifn_plot$plotcode)

all1<- merge(ifn_trees_10,ifn_plot, by.x="plotcode", by.y="plotcode" )
#Here are merged different info

#Still the same number of individual as before = 1384360
#Just a few infomations more
fundiv.tree.prev<- merge(all1,ifn_aggregation_75, by.x="plotcode", by.y="plotcode" ) #here are merged info based at different level of trees. (10 VS 75) This is why ba_ha at the plot level does not fit.
#We lost a number 33032 individual => 1351328
ifn_not_aggregated_75 <- anti_join(all1,fundiv.tree.prev,by="treecode")
ifn_not_aggregated_75 <- anti_join(all1,fundiv.tree,by="treecode")



###Estimate mortality.plot and sp.mortality.plot rates in the same way of IFN-France: 
#test on a party or all

#new line 30/01/2018 to correct the problem with NA values and mod data

ftp.short <- test
ftp.short <- fundiv.tree
ftp.short$plotcode <- as.character(ftp.short$plotcode)
ftp.short$treecode <- as.character(ftp.short$treecode)
ftp.short$plotcode = paste0(ftp.short$country.x,ftp.short$plotcode) #new line 30/01/2018 unique plotcode name
ftp.short$treecode = paste0(ftp.short$country.x,ftp.short$treecode) #new line 30/01/2018 unique treecode name
#ft.short.test <- ftp.short[ftp.short$plotcode=="DE10001_2",c(1:6,18,19,32:35)]
ftp.short <- ftp.short[,c(1:31)] #here select the lines we want to test


#these data are not going to be processed and used to calculate anything apart from the recruitment
#Here we delete these information that we don't need. 

ftp.short[which(complete.cases(ftp.short$ba1_mod)==TRUE|complete.cases(ftp.short$dbh1_mod)==TRUE),"ba1"] <- NA
ftp.short[which(complete.cases(ftp.short$ba1_mod)==TRUE|complete.cases(ftp.short$dbh1_mod)==TRUE),"dbh1"] <- NA
ftp.short[which(complete.cases(ftp.short$ba1_mod)==TRUE|complete.cases(ftp.short$dbh1_mod)==TRUE),"ba_ha1"] <- NA
ftp.short[!is.na(ftp.short$ba_ha2) & ftp.short$ba_ha2==0,"ba_ha2"] <- NA #Those trees that have 0 as Ba

ftp.short[ftp.short$country!="DE","weight1"] <- ftp.short[ftp.short$country!="DE","weight1"]*(10000/(ftp.short[ftp.short$country!="DE","weight1"]^3*3.14159265))
ftp.short[ftp.short$country!="DE","weight2"] <- ftp.short[ftp.short$country!="DE","weight2"]*(10000/(ftp.short[ftp.short$country!="DE","weight2"]^3*3.14159265))

test2 = ftp.short$ba1*ftp.short$weight1
cor.test(test2,ftp.short$ba_ha1)
test2 = ftp.short$ba2*ftp.short$weight2
cor.test(test2,ftp.short$ba_ha2)

dat.fundiv=split(ftp.short[ ,c("plotcode","treecode","speciesid","treestatus_th","dbh1","dbh2","height1","height2", "ba1", "ba_ha1", "ba2", "ba_ha2", "bachange_ha", "bachange_ha_yr", "dbh1_mod", "ba1_mod", "bachange_mod", "weight1", "weight2", "country.x", "cluster", "country.y", "longitude", "latitude", "yearsbetweensurveys", "surveydate1", "surveydate2", "speciesrichness", "ba_ha", "n_ha","country")],as.character(ftp.short$plotcode))
ftp.short[,c("plotcode","treecode","speciesid","treestatus_th","dbh1","dbh2","height1","height2", "ba1", "ba_ha1", "ba2", "ba_ha2", "bachange_ha", "bachange_ha_yr", "dbh1_mod", "ba1_mod", "bachange_mod", "weight1", "weight2", "country.x", "cluster", "country.y", "longitude", "latitude", "yearsbetweensurveys", "surveydate1", "surveydate2", "speciesrichness", "ba_ha", "n_ha","country")]
colnames(ftp.short)
colnames(treefinal)




#################################
### plot mortality (% tree /ha)##
#################################
data <- do.call(rbind, mclapply(dat.fundiv,function(df){
  # idp
  i=df[1,"plotcode"]
  dat=matrix(NA,nrow(ftp.short[ftp.short$plotcode==i,]),2)
  dat[,2]=c(1:nrow(ftp.short))[ftp.short$plotcode==i]
  if(nrow(df)>0) {
    
    # dead trees
    df.m=df[df$treestatus=="4",]
    
    if(nrow(df.m)==0) dat[,1]=0
    
    if(nrow(df.m)>0) {
      
      #dens.tot=sum(df$weight1)#Probleme quand on a un recrutement et weight 1 = NA
      dens.tot=sum(df$weight1[!is.na(df$weight1)])
      dens.m=sum(df.m$weight1) 
      res=dens.m/dens.tot
      dat[,1]=as.numeric(res)
      
    }
  }
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
mortality.plot.rate=as.numeric(data[order(data[,2]),1])
summary(mortality.plot.rate)

####################################
### plot SP.mortality (% tree /ha)##
####################################
# plot specific mortality (% tree /ha): dens dead trees IFN3 / dens live trees IFN2
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
        # dead trees
        df.m=df[df$treestatus=="4",]
        
        if(nrow(df.m)==0) dat[ftp.short$speciesid[ftp.short$plotcode==i]==s,1]=0
        
        if(nrow(df.m)>0) {
          
          #dens.tot=sum(df$weight1) pbm des NA
          dens.tot=sum(df$weight1[!is.na(df$weight1)])
          dens.m=sum(df.m$weight1)
          res=dens.m/dens.tot	
          dat[ftp.short$speciesid[ftp.short$plotcode==i]==s,1]=as.numeric(res)
          #dat[,1]=as.numeric(res)
          #      dat[fundiv.tree$plotcode[fundiv.tree$plotcode==i]==s,1]=as.numeric(res)
          
        }
      }
    # }
  }
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
sp.mortality.plot.rate=as.numeric(data[order(data[,2]),1])
summary(sp.mortality.plot.rate)

#################################
### plot mortality (% tree /ha)## # ba_ha
#################################
data <- do.call(rbind, mclapply(dat.fundiv,function(df){
  # idp
  i=df[1,"plotcode"]
  dat=matrix(NA,nrow(ftp.short[ftp.short$plotcode==i,]),2)
  dat[,2]=c(1:nrow(ftp.short))[ftp.short$plotcode==i]
  if(nrow(df)>0) {
    
    # dead trees
    df.m=df[df$treestatus=="4",]
    
    if(nrow(df.m)==0) dat[,1]=0
    
    if(nrow(df.m)>0) {
      
      #dens.tot=sum(df$weight1)#Probleme quand on a un recrutement et weight 1 = NA
      dens.tot=sum(df$ba_ha1[!is.na(df$ba_ha1)])
      dens.m=sum(df.m$ba_ha1) 
      res=dens.m/dens.tot
      dat[,1]=as.numeric(res)
      
    }
  }
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
mortality.plot.ba=as.numeric(data[order(data[,2]),1])
summary(mortality.plot.ba)

####################################
### plot SP.mortality (% tree /ha)## ba
####################################
# plot specific mortality (% tree /ha): dens dead trees IFN3 / dens live trees IFN2
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
        # dead trees
        df.m=df[df$treestatus=="4",]
        
        if(nrow(df.m)==0) dat[ftp.short$speciesid[ftp.short$plotcode==i]==s,1]=0
        
        if(nrow(df.m)>0) {
          
          #dens.tot=sum(df$weight1) pbm des NA
          dens.tot=sum(df$ba_ha1[!is.na(df$ba_ha1)])
          dens.m=sum(df.m$ba_ha1)
          res=dens.m/dens.tot	
          dat[ftp.short$speciesid[ftp.short$plotcode==i]==s,1]=as.numeric(res)
          #dat[,1]=as.numeric(res)
          #      dat[fundiv.tree$plotcode[fundiv.tree$plotcode==i]==s,1]=as.numeric(res)
          
        }
      }
    # }
  }
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
sp.mortality.plot.ba=as.numeric(data[order(data[,2]),1])
summary(sp.mortality.plot.ba)


#################################
### plot recrutment (% tree /ha)#
#################################
data <- do.call(rbind, mclapply(dat.fundiv,function(df){
  # idp
  i=df[1,"plotcode"]
  dat=matrix(NA,nrow(ftp.short[ftp.short$plotcode==i,]),2)
  dat[,2]=c(1:nrow(ftp.short))[ftp.short$plotcode==i]
  if(nrow(df)>0) {
    
    # dead trees
    df.m=df[df$treestatus=="1",]
    
    if(nrow(df.m)==0) dat[,1]=0
    
    if(nrow(df.m)>0) {
      
      #dens.tot=sum(df$weight1)#Probleme quand on a un recrutement et weight 1 = NA
      dens.tot=sum(df$weight2[!is.na(df$weight2)]) #somme de tous ceux qui ne sont pas mort
      dens.m=sum(df.m$weight2) 
      res=dens.m/dens.tot
      dat[,1]=as.numeric(res)
      
    }
  }
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
recruitment.plot.rate=as.numeric(data[order(data[,2]),1])
summary(recruitment.plot.rate)

####################################
### Sp.plot recrutment (% tree /ha)#
####################################
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
        # dead trees
        df.m=df[df$treestatus=="1",]
        
        if(nrow(df.m)==0) dat[ftp.short$speciesid[ftp.short$plotcode==i]==s,1]=0
        
        if(nrow(df.m)>0) {
          
          #dens.tot=sum(df$weight2) probleme des NA
          dens.tot=sum(df$weight2, na.rm = T)
          dens.tot=sum(df$weight2[!is.na(df$weight2)])
          dens.m=sum(df.m$weight2)
          res=dens.m/dens.tot	
          dat[ftp.short$speciesid[ftp.short$plotcode==i]==s,1]=as.numeric(res)
          #dat[,1]=as.numeric(res)
          #      dat[fundiv.tree$plotcode[fundiv.tree$plotcode==i]==s,1]=as.numeric(res)
          
        }
      }
    # }
  }
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
sp.recruitment.plot.rate=as.numeric(data[order(data[,2]),1])
summary(sp.recruitment.plot.rate)

#################################
### plot recrutment (% tree /ha)# ba
#################################
data <- do.call(rbind, mclapply(dat.fundiv,function(df){
  # idp
  i=df[1,"plotcode"]
  dat=matrix(NA,nrow(ftp.short[ftp.short$plotcode==i,]),2)
  dat[,2]=c(1:nrow(ftp.short))[ftp.short$plotcode==i]
  if(nrow(df)>0) {
    
    # dead trees
    df.m=df[df$treestatus=="1",]
    
    if(nrow(df.m)==0) dat[,1]=0
    
    if(nrow(df.m)>0) {
      
      #dens.tot=sum(df$weight1)#Probleme quand on a un recrutement et weight 1 = NA
      dens.tot=sum(df$ba_ha2[!is.na(df$ba_ha2)]) #somme de tous ceux qui ne sont pas mort = survivant ou recru
      dens.m=sum(df.m$ba_ha2) #juste recru 
      res=dens.m/dens.tot
      dat[,1]=as.numeric(res)
      
    }
  }
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
recruitment.plot.ba=as.numeric(data[order(data[,2]),1])
summary(recruitment.plot.ba)

####################################
### Sp.plot recrutment (% tree /ha)# ba
####################################
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
        # dead trees
        df.m=df[df$treestatus=="1",]
        
        if(nrow(df.m)==0) dat[ftp.short$speciesid[ftp.short$plotcode==i]==s,1]=0
        
        if(nrow(df.m)>0) {
          
          #dens.tot=sum(df$weight2) probleme des NA
          #dens.tot=sum(df$weight2, na.rm = T)
          dens.tot=sum(df$ba_ha2[!is.na(df$ba_ha2)]) #Vivant (recru ou survivors)
          dens.m=sum(df.m$ba_ha2) #juste les recru 
          res=dens.m/dens.tot	
          dat[ftp.short$speciesid[ftp.short$plotcode==i]==s,1]=as.numeric(res)
          #dat[,1]=as.numeric(res)
          #      dat[fundiv.tree$plotcode[fundiv.tree$plotcode==i]==s,1]=as.numeric(res)
          
        }
      }
    # }
  }
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
sp.recruitment.plot.ba=as.numeric(data[order(data[,2]),1])
summary(sp.recruitment.plot.ba)


ftp.short[,"mortality.plot.rate"]=mortality.plot.rate
ftp.short[,"sp.mortality.plot.rate"]=sp.mortality.plot.rate
ftp.short[,"mortality.plot.ba"]=mortality.plot.ba
ftp.short[,"sp.mortality.plot.ba"]=sp.mortality.plot.ba
ftp.short[,"recruitment.plot.rate"]=recruitment.plot.rate
ftp.short[,"sp.recruitment.plot.rate"]=sp.recruitment.plot.rate
ftp.short[,"recruitment.plot.ba"]=recruitment.plot.ba
ftp.short[,"sp.recruitment.plot.ba"]=sp.recruitment.plot.ba


####BA.ha.plot##### 
data <- do.call(rbind, mclapply(dat.fundiv,function(df){
  
  # idp
  i=df[1,"plotcode"]
  dat=matrix(NA,nrow(ftp.short[ftp.short$plotcode==i,]),2)
  dat[,2]=c(1:nrow(ftp.short))[ftp.short$plotcode==i]
  
  if(nrow(df)>0) {
    
    # set to NA dead trees BA
    df$ba_ha2[df$treestatus=="4"]=NA # Na pour les dead qui ont une valeur au second = 148 individus 
    
    # set mean dbh for alive trees
    # if known dbh data >= r%
    #   r.dbh=nrow(df[!is.na(df$bachange_ha_yr),])/nrow(df)# ratio of known dbh data
    #   if(thres==F | r.dbh>=th) {
    
    #if(any(is.na(df$dbh))) df$dbh[is.na(df$dbh)]=mean(df$dbh,na.rm=T)
    ba.ha=df$ba_ha2
    #ba.ha=df$ba2*df$weight2
    res=sum(ba.ha,na.rm=T)
    if(length(res)==0 | res==0) res=NA
    dat[,1]=as.numeric(res)
  }
  #}
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
BA.ha.plot=as.numeric(data[order(data[,2]),1])

#Should be the same as the Ba_ha (FOR THE SECOND ONE)

# neighbouring tree BA (m2/ha). Based at T2. Ba-ha/plot measure the sum of BAha a T2. 
#BAnei = BAha - poids de l'arbe en question x ba.
#Cette mesure compte les recru mais pas les morts.
#BAnei=(BA.ha.plot-(ftp.short$ba2*ftp.short$weight2)) # !!!!!!!!! Mistake test correlation ba_ha et nuvelkle metric
#Avant Banei = ba_ha - ba2*weight2

BAnei=BA.ha.plot-ftp.short$ba_ha2

####################################################
######BAI.Plot = total basal area per plot; ########
######BAI.nei = neighbors basal area increment;#####
######BAnei = neighbors basal area sum at T2  ######
####################################################

# plot BA (m2/ha) This is the sum of the basal area increment / ha / year of each tree for a plot (without recrut and dead)
data <- do.call(rbind, mclapply(dat.fundiv,function(df){
  

  # idp
  i=df[1,"plotcode"]
  dat=matrix(NA,nrow(ftp.short[ftp.short$plotcode==i,]),2)
  dat[,2]=c(1:nrow(ftp.short))[ftp.short$plotcode==i]
  
  if(nrow(df)>0) {
    
    # set to NA dead trees BA
    df$bachange_ha_yr[df$treestatus=="4"]=NA
    # set mean dbh for alive trees
    # if known dbh data >= r%
    #   r.dbh=nrow(df[!is.na(df$bachange_ha_yr),])/nrow(df)# ratio of known dbh data
    #   if(thres==F | r.dbh>=th) {
    
    #if(any(is.na(df$dbh))) df$dbh[is.na(df$dbh)]=mean(df$dbh,na.rm=T)
   bai.ha=df$bachange_ha_yr
   res=sum(bai.ha,na.rm=T)
    if(length(res)==0 | res==0) res=NA
    dat[,1]=as.numeric(res)
  }
  #}
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
BAI.plot=as.numeric(data[order(data[,2]),1])

# neighbouring tree BAI (m2/ha) = basal area increment au plot à l'année de tous les arbres sauf l'individus e question
#This measure is an individual oe. herefore, dead and new recrut trees are NA's
BAInei=BAI.plot-ftp.short$bachange_ha_yr


#########################################
### BAI.J.plot (ba.sp.plot) m²/ha/year ##
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
        # dead trees
        #df.m=df[df$treestatus=="4",]
        #if(nrow(df.m)==0) dat[ftp.short$speciesid[ftp.short$plotcode==i]==s,1]=0
        #if(nrow(df.m)>0) {
        BAIpj.s=sum(df$bachange_ha_yr,na.rm=T)
        if(length(BAIpj.s)==0 | BAIpj.s==0) BAIpj.s=NA #rajout february
        dat[ftp.short$speciesid[ftp.short$plotcode==i]==s,1]=BAIpj.s
        
        #}
      }
  }
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
BAIj.plot=as.numeric(data[order(data[,2]),1])
summary(BAIj.plot)

##### neighbouring conspecific tree BAI (m2/ha)
BAIneicon=BAIj.plot-ftp.short$bachange_ha_yr


####################################
### BA.J.plot (ba.sp.plot) m²/ha  ##
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
        
        df=DF[!is.na(DF$speciesid) & DF$speciesid==s,] # On comptabilise les arbres morts ??? 148 en tout
        # dead trees
        #df.m=df[df$treestatus=="4",]
        #if(nrow(df.m)==0) dat[ftp.short$speciesid[ftp.short$plotcode==i]==s,1]=0
        #if(nrow(df.m)>0) {
        BApj.s=sum(df$ba_ha2,na.rm=T)
        #if(length(BAIpj.s)==0 | BAIpj.s==0) BAIpj.s=NA
        dat[ftp.short$speciesid[ftp.short$plotcode==i]==s,1]=BApj.s
        
        #}
      }
  }
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
BAj.plot=as.numeric(data[order(data[,2]),1])
summary(BAj.plot)

# neighbouring conspecific tree BAI (m2/ha)
BAneicon=BAj.plot-ftp.short$ba_ha2


ftp.short[,"BA.ha.plot"]=BA.ha.plot
ftp.short[,"BA.ha.plot2"]=BA.ha.plot2

ftp.short[,"BAnei"]=BAnei

ftp.short[,"BAI.plot"]=BAI.plot
ftp.short[,"BAInei"]=BAInei
ftp.short[,"BAnei"]=BAnei
ftp.short[,"BAIj.plot"]=BAIj.plot
ftp.short[,"BAIneicon"]=BAIneicon
ftp.short[,"BAj.plot"]=BAj.plot
ftp.short[,"BAneicon"]=BAneicon
#Reste plusiurs metrics


#Tester 100% corrélation ba.ha.plot 1 et 2 et ba_ha de Sophie. Limit
#Test correlation deux taux de mortalité et taux de recruitment 
#On obtient presque 100% dasn les trois cas. 
#

save(ftp.short,file="/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/Fundiv_alltree_allmetric.2018.RData")





