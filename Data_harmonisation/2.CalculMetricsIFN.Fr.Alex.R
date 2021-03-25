#Script to do the calculation again for all the variable that interest us with the IfnFr
library(parallel)
Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
load(paste0(Dir,"our-data/ifn.fr.renamed.RData")) #IFN with 3500 names corrected
load(paste0(Dir,"our-data/ifn.fr.metric.2018.RData")) #IFN with 3500 names corrected

#All the estimations .plot are calculated per ha  
load("/home/achangenet/Documents/FUNDIV\ -\ NFI\ -\ Europe/our-data/tree_data_IFN_france_with_competition-indices.bai_ha_forFUNDIVcomp.RData")
#Need to compare the rate in both databases in the end
#Correct a few lines
ifn.fr.renamed$espar <- as.character(ifn.fr.renamed$espar)
ifn.fr.renamed$binomial <- as.character(ifn.fr.renamed$binomial)
ifn.fr.renamed[!is.na(ifn.fr.renamed$binomial)&ifn.fr.renamed$binomial=="Salix_pentandra",6] <- "2500000"
ifn.fr.renamed[!is.na(ifn.fr.renamed$binomial)&ifn.fr.renamed$binomial=="Salix_trianda",6] <- "25000"


# 07/02/2018 add line to remove data above 7.5 cm dbh but keep dead. 
# 13/02/2018 add blocks to calculate mortality and recruitment differently 
# 19/02/18 corrections mistakes

ifn.fr.renamed <- ifn.fr.renamed[- which(ifn.fr.renamed$dbh<=9.9 & ifn.fr.renamed$dbh>=7.5),] #780802 to 637830
ifn.fr <- ifn.fr.renamed[,-c(17:22,40:44)]


#################################
###   BA.ha.yr (baI/ha/yr)     ##
###      m²/hecatre/year     #### 
#################################
colnames(ifn.fr)
BAI.ha=ifn.fr$BAI*ifn.fr$w
ifn.fr[,"BAI.ha"]=BAI.ha #bachange_ha_year de fundiv
BA.ha=ifn.fr$BA*ifn.fr$w
ifn.fr[,"BA.ha"]=BA.ha #equal to ba_ha2 of fundiv. The sum of this is BA.plot
str(ifn.fr)

###check in a small sample: 
#ifn.short <- ifn.fr[1:100,]
#Here is divided by plot: as many lists as number of plots. 
dat.fr=split(ifn.fr[,c("idp","alive","dbh","BA","BAI","veget","binomial","w","BAI.ha")],as.character(ifn.fr$idp))

# use unknown dbh data threshold ?
thres=F# no because of the use of w to weight
th=.66# if yes (What does it mean ?)

#################################
###     Number of species      ##
#################################
data <- do.call(rbind, mclapply(dat.fr,function(df){
  
  # idp
  i=df[1,"idp"]
  dat=matrix(NA,nrow(ifn.fr[ifn.fr$idp==i,]),2)
  dat[,2]=c(1:nrow(ifn.fr))[ifn.fr$idp==i]
  
  # remove 'arbres chablis' (dead or alive)
  chab=c(1:nrow(df))[!is.na(df$veget) & df$veget=="A"] ; chab=chab[!is.na(chab)]
  if(length(chab)>0) df=df[-c(chab),]
  
  if(nrow(df)>0) {
    
    # set to NA dead trees BA
    df$BA[df$alive==F]=NA
    # set mean dbh for alive trees
    # if known dbh data >= r%
    r.dbh=nrow(df[!is.na(df$dbh),])/nrow(df)# ratio of known dbh data
    if(thres==F | r.dbh>=th) {
      
      #if(any(is.na(df$dbh))) df$dbh[is.na(df$dbh)]=mean(df$dbh,na.rm=T)
      total=length(unique(df$binomial))
      res=sum(total,na.rm=T)
      #    if(length(res)==0 | res==0) res=NA
      dat[,1]=as.numeric(res)
      
    }
  }
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
richness.plot=as.numeric(data[order(data[,2]),1])
summary(richness.plot)

####################################
###   BAI.yr.plot (baI/yr) m²/ha  ##
####################################
# plot BAI (m2/ha) => Par an ou pas ? a priori oui
data <- do.call(rbind, mclapply(dat.fr,function(df){
  
  # idp
  i=df[1,"idp"]
  dat=matrix(NA,nrow(ifn.fr[ifn.fr$idp==i,]),2)
  dat[,2]=c(1:nrow(ifn.fr))[ifn.fr$idp==i]
  
  # remove 'arbres chablis' (dead or alive)
  #chab=c(1:nrow(df))[!is.na(df$veget) & df$veget=="A"] ; chab=chab[!is.na(chab)]
  #if(length(chab)>0) df=df[-c(chab),]
  
  if(nrow(df)>0) {
    
    # set to NA dead trees BA
    df$BAI.ha[df$alive==F]=NA
    # set mean dbh for alive trees
    # if known dbh data >= r%
    #r.dbh=nrow(df[!is.na(df$dbh),])/nrow(df)# ratio of known dbh data
    #if(thres==F | r.dbh>=th) {
      
      #if(any(is.na(df$dbh))) df$dbh[is.na(df$dbh)]=mean(df$dbh,na.rm=T)
      #bai.ha=df$BAI*df$w #Garder cette ligne car équivalent au BAI de fundiv
      #res=sum(bai.ha,na.rm=T)
      res=sum(df$BAI.ha,na.rm=T)
      if(length(res)==0 | res==0) res=NA
      dat[,1]=as.numeric(res)
      
    #}
  }
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
BAI.plot=as.numeric(data[order(data[,2]),1])
summary(BAI.plot)

####################################
###     BAInei (baI/yr/ha)      ###
####################################
# neighbouring tree BA (m2/ha)

#BAInei=BAI.plot-ifn.fr$BAI*ifn.fr$w #avant
BAInei=BAI.plot-BAI.ha #now


####################################
### BAI.J.plot (ba.sp.plot) m²/ha ##
####################################
# plot specific BA (m2/ha)
data <- do.call(rbind, mclapply(dat.fr,function(DF){
  
  # idp
  i=DF[1,"idp"]
  dat=matrix(NA,nrow(ifn.fr[ifn.fr$idp==i,]),2)
  dat[,2]=c(1:nrow(ifn.fr))[ifn.fr$idp==i]
  
  # remove 'arbres chablis' (dead or alive)
  chab=c(1:nrow(DF))[!is.na(DF$veget) & DF$veget=="A"] ; chab=chab[!is.na(chab)]
  if(length(chab)>0) DF=DF[-c(chab),]
  
  if(nrow(DF)>0) {
    
    S=unique(DF$binomial)
    S=S[!is.na(S)]
    if(length(S)>0)
      for(s in S) {
        
        df=DF[!is.na(DF$binomial) & DF$binomial==s,]
        # set to NA dead trees BA
        df$BAI[df$alive==F]=NA
        # set mean dbh for alive trees if needed
        # if known dbh data >= r%
        r.dbh=nrow(df[!is.na(df$dbh),])/nrow(df)# ratio of known dbh data
        if(thres==F | r.dbh>=th) {
          
          #if(any(is.na(df$dbh))) df$dbh[is.na(df$dbh)]=mean(df$dbh,na.rm=T)
          
          BAIpj.s=sum(df$BAI*df$w,na.rm=T) #ou bien sum de BAI.ha
          if(length(BAIpj.s)==0 | BAIpj.s==0) BAIpj.s=NA
          dat[ifn.fr$binomial[ifn.fr$idp==i]==s,1]=BAIpj.s
          
        }
      }
  }
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
BAIj.plot=as.numeric(data[order(data[,2]),1])
summary(BAIj.plot)


####################################################
#  neighbouring conspecific tree BAI (m2/ha).    ###
# Croissance des conspecifique moins nous au plot  #
####################################################

#BAIneicon=BAIj.plot-ifn.fr$BAI*ifn.fr$w #avant
BAIneicon=BAIj.plot-ifn.fr$BAI.ha #now


ifn.fr[,"BAInei"]=BAInei
ifn.fr[,"BAIneicon"]=BAIneicon
ifn.fr[,"BAI.plot"]=BAI.plot
ifn.fr[,"BAIj.plot"]=BAIj.plot
ifn.fr[,"richness.plot"]=richness.plot

dat.fr=split(ifn.fr[,c("idp","alive","dbh","BA.ha","veget","binomial","w")],as.character(ifn.fr$idp)) #correction BA.ha instead of BA !

# use unknown dbh data threshold ?
thres=F# no because of the use of w to weight
th=.66# if yes

####################################
### BA.plot (ba.plot) m²/ha !!!!  ##
####################################
data <- do.call(rbind, mclapply(dat.fr,function(df){
  
  # idp
  i=df[1,"idp"]
  dat=matrix(NA,nrow(ifn.fr[ifn.fr$idp==i,]),2)
  dat[,2]=c(1:nrow(ifn.fr))[ifn.fr$idp==i]
  
  
  
  if(nrow(df)>0) {
    
    # set to NA dead trees BA
    df$BA.ha[df$alive==F]=NA
    # set mean dbh for alive trees
    # if known dbh data >= r%
    #r.dbh=nrow(df[!is.na(df$dbh),])/nrow(df)# ratio of known dbh data
    #if(thres==F | r.dbh>=th) {
      
      #if(any(is.na(df$dbh))) df$dbh[is.na(df$dbh)]=mean(df$dbh,na.rm=T)
      #ba.ha=df$BA*df$w
      res=sum(df$BA.ha,na.rm=T)
      if(length(res)==0 | res==0) res=NA
      dat[,1]=as.numeric(res)
      
    #}
  }
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
BA.plot=as.numeric(data[order(data[,2]),1])
summary(BA.plot)
# neighbouring tree BA (m2/ha)

#BAnei=BA.plot-ifn.fr$BA*ifn.fr$w #avant
BAnei=BA.plot-ifn.fr$BA.ha #now


####################################
### BA.J.plot (ba.sp.plot) m²/ha  ##
####################################
# plot specific BA (m2/ha) au plot
data <- do.call(rbind, mclapply(dat.fr,function(DF){
  
  # idp
  i=DF[1,"idp"]
  dat=matrix(NA,nrow(ifn.fr[ifn.fr$idp==i,]),2)
  dat[,2]=c(1:nrow(ifn.fr))[ifn.fr$idp==i]
  
  # remove 'arbres chablis' (dead or alive)
  chab=c(1:nrow(DF))[!is.na(DF$veget) & DF$veget=="A"] ; chab=chab[!is.na(chab)]
  if(length(chab)>0) DF=DF[-c(chab),]
  
  if(nrow(DF)>0) {
    
    S=unique(DF$binomial)
    S=S[!is.na(S)]
    if(length(S)>0)
      for(s in S) {
        
        df=DF[!is.na(DF$binomial) & DF$binomial==s,]
        # set to NA dead trees BA
        df$BA.ha[df$alive==F]=NA
        # set mean dbh for alive trees if needed
        # if known dbh data >= r%
        r.dbh=nrow(df[!is.na(df$dbh),])/nrow(df)# ratio of known dbh data
        if(thres==F | r.dbh>=th) {
          
          #if(any(is.na(df$dbh))) df$dbh[is.na(df$dbh)]=mean(df$dbh,na.rm=T)
          
          BApj.s=sum(df$BA.ha,na.rm=T) #ramené a l'hecatre
          if(length(BApj.s)==0 | BApj.s==0) BApj.s=NA
          dat[ifn.fr$binomial[ifn.fr$idp==i]==s,1]=BApj.s
          
        }
      }
  }
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
BAj.plot=as.numeric(data[order(data[,2]),1])
summary(BAj.plot)

# neighbouring conspecific tree BA (m2/ha)

#BAneicon=BAj.plot-ifn.fr$BA*ifn.fr$w #avant
BAneicon=BAj.plot-ifn.fr$BA.ha #now

####################################
###  plot mortality (% tree /ha) ###
####################################
data <- do.call(rbind, mclapply(dat.fr,function(df){
  # idp
  i=df[1,"idp"]
  dat=matrix(NA,nrow(ifn.fr[ifn.fr$idp==i,]),2)
  dat[,2]=c(1:nrow(ifn.fr))[ifn.fr$idp==i]
  # remove 'arbres chablis' (dead or alive)
  chab=c(1:nrow(df))[!is.na(df$veget) & df$veget=="A"] ; chab=chab[!is.na(chab)]
  if(length(chab)>0) df=df[-c(chab),]
  if(nrow(df)>0) {
    
    # dead trees
    df.m=df[df$alive==F,]
    
    if(nrow(df.m)==0) dat[,1]=0
    
    if(nrow(df.m)>0) {
      
      # not accounting for sampling effort
      #res=nrow(df.m)/nrow(df)
      # if known dbh data >= r% -> set dbh (when NA) to mean tree DBH to account for sampling effort
      r.dbh=nrow(df[!is.na(df$dbh),])/nrow(df)# ratio of known dbh data
      if(thres==F | r.dbh>=th) {
        dens.tot=sum(df$w) #enlever les NA ???
        dens.m=sum(df.m$w)
        res=dens.m/dens.tot
        dat[,1]=as.numeric(res)
        
      }
    }
  }
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
mortality.plot.rate=as.numeric(data[order(data[,2]),1])
summary(mortality.plot.rate)

####################################
### plot SP.mortality (% tree /ha)##
####################################: dens dead trees IFN3 / dens live trees IFN2
data <- do.call(rbind, mclapply(dat.fr,function(DF){
  
  # idp
  i=DF[1,"idp"]
  dat=matrix(NA,nrow(ifn.fr[ifn.fr$idp==i,]),2)
  dat[,2]=c(1:nrow(ifn.fr))[ifn.fr$idp==i]
  # remove 'arbres chablis' (dead or alive)
  chab=c(1:nrow(DF))[!is.na(DF$veget) & DF$veget=="A"] ; chab=chab[!is.na(chab)]
  if(length(chab)>0) DF=DF[-c(chab),]
  
  if(nrow(DF)>0) {
    
    S=unique(DF$binomial)
    S=S[!is.na(S)]
    if(length(S)>0)
      for(s in S) {
        
        df=DF[!is.na(DF$binomial) & DF$binomial==s,]
        # dead trees
        df.m=df[df$alive==F,]
        
        if(nrow(df.m)==0) dat[ifn.fr$binomial[ifn.fr$idp==i]==s,1]=0
        
        if(nrow(df.m)>0) {
          
          # not accounting for sampling effort
          #res=nrow(df.m)/nrow(df)
          
          # if known dbh data >= r% -> set dbh (when NA) to mean tree DBH to account for sampling effort
          r.dbh=nrow(df[!is.na(df$dbh),])/nrow(df)# ratio of known dbh data
          if(thres==F | r.dbh>=th) {
            
            # old fashion
            #if(any(is.na(df$dbh))) {df[is.na(df$dbh),"dbh"]=mean(df$dbh,na.rm=T) ; df.m[is.na(df.m$dbh),"dbh"]=mean(df$dbh,na.rm=T)}
            #dens.tot=sum(c(toha(nrow(df[df$dbh<22.5,]),10),toha(nrow(df[df$dbh>=22.5 & df$dbh<37.5,]),30),toha(nrow(df[df$dbh>=37.5,]),50)),na.rm=T)
            #dens.m=sum(c(toha(nrow(df.m[df.m$dbh<22.5,]),10),toha(nrow(df.m[df.m$dbh>=22.5 & df.m$dbh<37.5,]),30),toha(nrow(df.m[df.m$dbh>=37.5,]),50)),na.rm=T)
            
            dens.tot=sum(df$w) #Enlever les NA ? 
            dens.m=sum(df.m$w)
            res=dens.m/dens.tot
            dat[ifn.fr$binomial[ifn.fr$idp==i]==s,1]=as.numeric(res)
            
          }
        }
      }
  }
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
sp.mortality.plot.rate=as.numeric(data[order(data[,2]),1])
summary(sp.mortality.plot.rate)


####################################
###  plot mortality (% tree /ha) ###  BA
####################################
data <- do.call(rbind, mclapply(dat.fr,function(df){
  
  # idp
  i=df[1,"idp"]
  dat=matrix(NA,nrow(ifn.fr[ifn.fr$idp==i,]),2)
  dat[,2]=c(1:nrow(ifn.fr))[ifn.fr$idp==i]
  # remove 'arbres chablis' (dead or alive)
  chab=c(1:nrow(df))[!is.na(df$veget) & df$veget=="A"] ; chab=chab[!is.na(chab)]
  if(length(chab)>0) df=df[-c(chab),]
  
  if(nrow(df)>0) {
    
    # dead trees
    df.m=df[df$alive==F,]
    
    if(nrow(df.m)==0) dat[,1]=0
    
    if(nrow(df.m)>0) {
      
      # not accounting for sampling effort
      #res=nrow(df.m)/nrow(df)
      
      # if known dbh data >= r% -> set dbh (when NA) to mean tree DBH to account for sampling effort
      r.dbh=nrow(df[!is.na(df$dbh),])/nrow(df)# ratio of known dbh data
      if(thres==F | r.dbh>=th) {
        
        # old fashion
        #if(any(is.na(df$dbh))) {df[is.na(df$dbh),"dbh"]=mean(df$dbh,na.rm=T) ; df.m[is.na(df.m$dbh),"dbh"]=mean(df$dbh,na.rm=T)}
        #dens.tot=sum(c(toha(nrow(df[df$dbh<22.5,]),10),toha(nrow(df[df$dbh>=22.5 & df$dbh<37.5,]),30),toha(nrow(df[df$dbh>=37.5,]),50)),na.rm=T)
        #dens.m=sum(c(toha(nrow(df.m[df.m$dbh<22.5,]),10),toha(nrow(df.m[df.m$dbh>=22.5 & df.m$dbh<37.5,]),30),toha(nrow(df.m[df.m$dbh>=37.5,]),50)),na.rm=T)
        
        dens.tot=sum(df$BA.ha[!is.na(df$BA.ha)]) #enlever les NA ???
        dens.m=sum(df.m$BA.ha)
        res=dens.m/dens.tot
        dat[,1]=as.numeric(res)
        
      }
    }
  }
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
mortality.plot.ba=as.numeric(data[order(data[,2]),1])
summary(mortality.plot.ba)

####################################
### plot SP.mortality (% tree /ha)##  BA
####################################: dens dead trees IFN3 / dens live trees IFN2
data <- do.call(rbind, mclapply(dat.fr,function(DF){
  
  # idp
  i=DF[1,"idp"]
  dat=matrix(NA,nrow(ifn.fr[ifn.fr$idp==i,]),2)
  dat[,2]=c(1:nrow(ifn.fr))[ifn.fr$idp==i]
  # remove 'arbres chablis' (dead or alive)
  chab=c(1:nrow(DF))[!is.na(DF$veget) & DF$veget=="A"] ; chab=chab[!is.na(chab)]
  if(length(chab)>0) DF=DF[-c(chab),]
  
  if(nrow(DF)>0) {
    
    S=unique(DF$binomial)
    S=S[!is.na(S)]
    if(length(S)>0)
      for(s in S) {
        
        df=DF[!is.na(DF$binomial) & DF$binomial==s,]
        # dead trees
        df.m=df[df$alive==F,]
        
        if(nrow(df.m)==0) dat[ifn.fr$binomial[ifn.fr$idp==i]==s,1]=0
        
        if(nrow(df.m)>0) {
          
          # not accounting for sampling effort
          #res=nrow(df.m)/nrow(df)
          
          # if known dbh data >= r% -> set dbh (when NA) to mean tree DBH to account for sampling effort
          r.dbh=nrow(df[!is.na(df$dbh),])/nrow(df)# ratio of known dbh data
          if(thres==F | r.dbh>=th) {
            
            # old fashion
            #if(any(is.na(df$dbh))) {df[is.na(df$dbh),"dbh"]=mean(df$dbh,na.rm=T) ; df.m[is.na(df.m$dbh),"dbh"]=mean(df$dbh,na.rm=T)}
            #dens.tot=sum(c(toha(nrow(df[df$dbh<22.5,]),10),toha(nrow(df[df$dbh>=22.5 & df$dbh<37.5,]),30),toha(nrow(df[df$dbh>=37.5,]),50)),na.rm=T)
            #dens.m=sum(c(toha(nrow(df.m[df.m$dbh<22.5,]),10),toha(nrow(df.m[df.m$dbh>=22.5 & df.m$dbh<37.5,]),30),toha(nrow(df.m[df.m$dbh>=37.5,]),50)),na.rm=T)
            
            dens.tot=sum(df$BA.ha[!is.na(df$BA.ha)]) #Enlever les NA ? 
            dens.m=sum(df.m$BA.ha)
            res=dens.m/dens.tot
            dat[ifn.fr$binomial[ifn.fr$idp==i]==s,1]=as.numeric(res)
            
          }
        }
      }
  }
  
  as.data.frame(dat)
  
},mc.cores=15,mc.silent=T))
sp.mortality.plot.ba=as.numeric(data[order(data[,2]),1])
summary(sp.mortality.plot.ba)


ifn.fr[,"BAnei"]=BAnei
ifn.fr[,"BAneicon"]=BAneicon
ifn.fr[,"BA.plot"]=BA.plot
ifn.fr[,"BAj.plot"]=BAj.plot
ifn.fr[,"mortality.plot.rate"]=mortality.plot.rate
ifn.fr[,"sp.mortality.plot.rate"]=sp.mortality.plot.rate
ifn.fr[,"mortality.plot.ba"]=mortality.plot.ba
ifn.fr[,"sp.mortality.plot.ba"]=sp.mortality.plot.ba

#Verif mortality vs mortality BA

save(ifn.fr,file="/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/ifn.fr.metric.2018.RData")



save(ifn.fr,file="/home/marta/Dropbox/IFN/tree_data_IFN_france_with_competition-indices.bai_ha_forFUNDIVcomp.RData")
save(ifn.fr,file="/home/marta/gdata/FUNDIV/data/tree_data_IFN_france_with_competition-indices.bai_ha_forFUNDIVcomp.RData")



