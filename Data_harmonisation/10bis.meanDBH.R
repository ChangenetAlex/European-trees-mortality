# Alex on the 02/07/2018
# Script to add a calulation of the mean dbh and mean BAIj.plot. Also an absolute dbh sum for the plot  

require(parallel)

SavingInfo = "MeanDBH(dir='bureau' or 'home',
CODE = 'BETPEN',
seuil = 0.8
seuilC = 0.6)
This fonction calculate the mean dbh and BAIj for each plot and also the cumulated dbh"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))

MeanDBH <- function(dir="bureau",
                  CODE,
                  seuil = 0.7,
                  seuilC = 0.6) {
  if (dir == "home") {
    Dir = c("/Users/alexandrechangenet/Dropbox/FUNDIV/")
  } else if (dir == "bureau") {
    Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
  }
  setwd(Dir)
  Dir = c(paste0(Dir, "our-data/species/", CODE, "/CLIMAP"))
  setwd(Dir)
  list.files(Dir, pattern = paste0(".rds"))
  Mydf <- readRDS(paste0("Mydf_", CODE, "_", seuil, "_", seuilC, ".rds")) #Base de données plot full
  
  # Once loaded, we add the transformed dbh (*10 in the french inventory)
  Mydf[Mydf$country=="FR","dbh1"] <- Mydf[Mydf$country=="FR","dbh"]*10 # dbh for french inventory on the same scale (in the column dbh1)
  Mydf$bachange_ha_yr.1 <- Mydf$bachange_ha_yr # dbh for french inventory on the same scale (in the column dbh1)
  Mydf[which(Mydf$bachange_ha_yr.1<0),"bachange_ha_yr.1"] <- 0 # The ones under trhe threshold
  
  # Then we apply a fonction that calculate the sum dbh for each plot
  dat.fundiv=split(Mydf[],as.character(Mydf$plotcode))
  data <- do.call(rbind, mclapply(dat.fundiv,function(df){
    i=df[1,"plotcode"]
    dat=matrix(NA,nrow(Mydf[Mydf$plotcode==i,]),2)
    dat[,2]=c(1:nrow(Mydf))[Mydf$plotcode==i]
    if(nrow(df)>0) {
      dbh=df$dbh1
      res=sum(dbh,na.rm=T)
      if(length(res)==0 | res==0) res=NA
      dat[,1]=as.numeric(res)
    }
    as.data.frame(dat)
  },mc.cores=10,mc.silent=T))
  dbh.plot=as.numeric(data[order(data[,2]),1])
  
  data <- do.call(rbind, mclapply(dat.fundiv,function(df){
    i=df[1,"plotcode"]
    dat=matrix(NA,nrow(Mydf[Mydf$plotcode==i,]),2)
    dat[,2]=c(1:nrow(Mydf))[Mydf$plotcode==i]
    if(nrow(df)>0) {
      bai.1=df$bachange_ha_yr.1
      res=sum(bai.1,na.rm=T)
      if(length(res)==0 | res==0) res=NA
      dat[,1]=as.numeric(res)
    }
    as.data.frame(dat)
  },mc.cores=10,mc.silent=T))
  BAIj.plot.1=as.numeric(data[order(data[,2]),1])
  
  
  Mydf[,"BAIj.plot.1"]=BAIj.plot.1 # Here is the calculation 
  Mydf[,"dbh.plot"]=dbh.plot # Here is the calculation 
  Mydf[,"dbh.plot.mean"]=Mydf$dbh.plot/Mydf$treeNbr # Here is the average dbh of my plot (proxy of the age)
  Mydf[,"BAIj.plot.mean"]=Mydf$BAIj.plot.bis/Mydf$treeNbr # Also a proxy of the age
  Mydf[,"BAIj.plot.1.mean"]=Mydf$BAIj.plot.1/Mydf$treeNbr # Also a proxy of the age
  
  #cor.test(Mydf$BAIj.plot.mean,Mydf$BAIj.plot.1.mean)
  #tapply(Mydf$BAIj.plot.1.mean, Mydf$country, mean,na.rm=T)
  #tapply(Mydf$BAIj.plot.mean, Mydf$country, mean,na.rm=T)
  #tapply(Mydf$dbh.plot.mean, Mydf$country, mean,na.rm=T)
  #tapply(Mydf$dbh.plot, Mydf$country, mean,na.rm=T)
  
  saveRDS(get("Mydf"), paste0(Dir, "/Mydf2_", CODE, "_", seuil, "_", seuilC, ".rds")) # Work at the plot scale
  
}



