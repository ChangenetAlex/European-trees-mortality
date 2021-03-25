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


Dir <- c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data")
setwd("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/")

mod_ALLbin <- list()
mod_ALLneg <- list()
dfplot <- list()
dfplot2 <- list()
dfplot3 <- list()
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
  Dir = c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/") # Directory 
  Dir =c(paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/binomial/",Allmod[i],"/"))
  setwd(Dir)
  mod_ALLbin[[CODE]] <- get(load(file = list.files(pattern=".rda")))$data
  
  Dir = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/CLIMAP2019/")
  setwd(Dir)
  dfplot[[CODE]] <- readRDS(paste0("dfplot", CODE, seuil, ".rds")) #Base de données plot full No scale !!! sans zone de transition
  dfplot2[[CODE]] <- readRDS(paste0("dfplot2", CODE, seuil, ".rds")) #Base de données plot full No scale !!! sans zone de transition
  dfplot3[[CODE]] <- readRDS(paste0("dfplot3", CODE, seuil, ".rds")) #Base de données plot full No scale !!! sans zone de transition
  dftree[[CODE]] <- readRDS(paste0("Mydf_",CODE,"_",seuil,"_",seuilC,".rds"))
}

# Nbr individus total 
Mytest <- lapply(dftree, function(x){
  nrow(x)
})
sum(do.call(rbind,Mytest)[,1])
test <- do.call(rbind,dftree)
test[test$country=="FR" & test$treestatus_th==0,"treestatus_th"] <- 4  ### Dead french trees noted 3 instead of 0. (Like other inventories) (for calculation of treeNBR)
nrow(test[test$treestatus_th=="4",])
nrow(Mydf[Mydf$treestatus_th=="4",])

# 1989158 # Avant tri based on tfinal.biotic (in the table)

# 1595968 # Après tri (Qu'est ce qui a bougé ?)

# 153892 from tfinal.biotic (in the table)
C <- do.call(rbind,dfplot)
length(unique(C$plotcode))

# 134590 after removed managed

B <- do.call(rbind,dfplot2)
length(unique(B$plotcode))
# 110680 after removed transition zone !

A <- do.call(rbind,dfplot3)
length(unique(A$plotcode))
# 51466



## Here we summarize all the variables we used in the models 
C <- mapply(function (x,xx){
  apply(x[which(rownames(x)%in%rownames(xx)),c("sp.mortality.plot.rate.yr",
                                               "dbhJ.IMall.plot.mean",
                                               "BAIj.plot.1.mean.J.I",
                                               "treeNbr",
                                               "mean_spei12",
                                               "BA.O.plot.1", 
                                               "BA.ha.plot.1",
                                               "BAj.plot.1",
                                               "tmean.djf_climate_mean.30",
                                               "bio1_climate_mean.30",
                                               "bio5_climate_mean.30",
                                               "bio12_climate_mean.30",
                                               "bio13_climate_mean.30",
                                               "bio14_climate_mean.30",
                                               "ppet.mean_climate_mean.30")],2,function(y){paste0(round(mean(y),2),"\n(",round(sd(y),2),")\n",round(min(y),2),"-",round(max(y),2))})
},xx=mod_ALLbin,x=dfplot,SIMPLIFY = F)


C <- do.call(rbind,C)

write.csv2(C,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/QuantitableMort",".csv"))


## Here make a graph with distrib number of event par plot par country 
##
# individual mortality 
Dir=c("~/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
tfinal.biotic2 <- readRDS("our-data/tfinal.biotic.Sep2019.RDS")

length(unique(tfinal.biotic2$plotcode))
length(Mydf[-c(Mydf$treestatus_th%in%c("3","5")),])
nrow(Mydf[Mydf$treestatus_th%in%c("3","5"),])

# nbr ligne par plot. 
# nbr dead par plot moi,ns chablis et moins manegemed
# Idem * weight 
# Idem * weight (check autre script ou je l'ai déjà fait )
Mydf <- tfinal.biotic2
Mydf[Mydf$country=="FR" & Mydf$treestatus_th==0,"treestatus_th"] <- 4  ### Dead french trees noted 3 instead of 0. (Like other inventories) (for calculation of treeNBR)
Mydf[Mydf$country=="FR" & Mydf$treestatus_th==1,"treestatus_th"] <- 2  ### Ingrowth french trees noted 2 instead of 1. (Like other inventories) (fo calculation of treeNBR)
Mydf[Mydf$treestatus_th==6,"treestatus_th"] <- 5  ### For PICABI and PINSYL some trees (42 and 142 were noted as 6 and are effectively dead which false the calculations later on)
Mydf[Mydf$country=="FR","weight1"] <- Mydf[Mydf$country=="FR","weight2"]  ### Ingrowth french trees noted 2 instead of 1. (Like other inventories) (fo calculation of treeNBR)

# Test
Mydf <- Mydf[1:10000,] 

dat.fundiv=split(Mydf[],as.character(Mydf$plotcode))
# split my list into several list 
dataSplit <- split(dat.fundiv, rep(1:1000, each = 154))
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
data <- do.call(rbind, mclapply(dat.fundiv,function(df){
  i=df[1,"plotcode"]
  dat=matrix(NA,nrow(Mydf[Mydf$plotcode==i,]),7)
  dat[,7]=c(1:nrow(Mydf))[Mydf$plotcode==i]
  if(df$country!="FR"){
  if(nrow(df)>0) {
      dat[,1] <- as.numeric(nrow(df)) # Number of lines = total number of trees of the species (both time)
      dat[,2] <- sum(df[!is.na(df$weight1),"weight1"])
      dat[,3] <- as.numeric(nrow(df[df$treestatus_th=="4" & !is.na(df$weight1),])) # Real dead (used in calculation of mortality rate)
      dat[,4] <- sum(df[df$treestatus_th=="4" & !is.na(df$weight1),"weight1"]) # Weight Real dead (used in calculation of mortality rate)
      dat[,5] <- as.numeric(nrow(df[df$treestatus_th%in%c("3","5")  & !is.na(df$weight1),])) #dead trees with values of dbh (excluding some french trees for the caclulation of meandbh)
      dat[,6] <- sum(df[df$treestatus_th%in%c("3","5")  & !is.na(df$weight1),"weight1"]) #dead trees with values of dbh (excluding some french trees for the caclulation of meandbh)
    }}
  
  if(df$country=="FR"){
    # remove 'arbres chablis' (dead or alive)
    chab=c(1:nrow(df))[!is.na(df$veget) & df$veget=="A"] ; chab=chab[!is.na(chab)]
    if(length(chab)>0) df=df[-c(chab),]
    if(nrow(df)>0) {
      dat[,1] <- as.numeric(nrow(df)) # Number of lines = total number of trees of the species (both time)
      dat[,2] <- sum(df[!is.na(df$weight1),"weight1"])
      dat[,3] <- as.numeric(nrow(df[df$treestatus_th=="4" & !is.na(df$weight1),])) # Real dead (used in calculation of mortality rate)
      dat[,4] <- sum(df[df$treestatus_th=="4" & !is.na(df$weight1),"weight1"]) # Weight Real dead (used in calculation of mortality rate)
      dat[,5] <- NA #dead trees with values of dbh (excluding some french trees for the caclulation of meandbh)
      dat[,6] <- NA #dead trees with values of dbh (excluding some french trees for the caclulation of meandbh)
    }}
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
data

NBRtree.plot=as.numeric(data[order(data[,7]),1]) # weight of all dead trees in the plot
NBRtree.ha=as.numeric(data[order(data[,7]),2])      # sum of weight1 of all ingrowth and dead trees
NBRtreedead.plot=as.numeric(data[order(data[,7]),3]) # idem recrut 
NBRtreedead.ha=as.numeric(data[order(data[,7]),4]) # idem recrut 
NBRtreeFdead.plot=as.numeric(data[order(data[,7]),5]) # idem recrut 
NBRtreeFdead.ha=as.numeric(data[order(data[,7]),6]) # idem recrut

tfinal.biotic <- tfinal.biotic2

tfinal.biotic[,"NBRtree.plot"]=NBRtree.plot
tfinal.biotic[,"NBRtree.ha"]=NBRtree.ha
tfinal.biotic[,"NBRtreedead.plot"]=NBRtreedead.plot
tfinal.biotic[,"NBRtreedead.ha"]=NBRtreedead.ha
tfinal.biotic[,"NBRtreeFdead.plot"]=NBRtreeFdead.plot
tfinal.biotic[,"NBRtreeFdead.ha"]=NBRtreeFdead.ha

tfinal.biotic[,"logNBRtreedead.plot"]=log(1+tfinal.biotic[,"NBRtreedead.plot"])
tfinal.biotic[,"logNBRtreedead.ha"]=log(1+tfinal.biotic[,"NBRtreedead.ha"])

saveRDS(tfinal.biotic,file="~/SynologyDrive/FUNDIV - NFI - Europe/our-data/tfinal.biotic.Novembre2020.rds")
tfinal.biotic2 <- readRDS("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/tfinal.biotic.Novembre2020.rds")
tfinal.biotic <-tfinal.biotic2

by(tfinal.biotic$country,unique(tfinal.biotic$plotcode))
summary(as.factor(tfinal.biotic[tfinal.biotic$treestatus_th=="4","country"]))

### Now draw the plot by country for the 6 variables (= 4 plots with 6 lines on it)
dplot1 <- ggplot(tfinal.biotic) + 
  stat_density(aes(x=NBRtree.plot,color=country),size=1.25,key_glyph = "abline",geom = "line",trim=T,position="identity",adjust=1)+
  #geom_line(aes(x=NBRtree.plot,color=country,y=..density..),size=1.5,key_glyph = "abline",stat = "density")+
  labs(y=NULL, x="trees/plot",las=1)+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.85,0.55),legend.direction ="vertical",
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.9, "cm"),
        panel.border = element_rect(colour = "black",fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
dplot1


dplot2 <- ggplot(tfinal.biotic) + 
  stat_density(aes(x=NBRtree.ha,color=country),size=1.25,key_glyph = "abline",geom = "line",trim = T,adjust=1)+
  #geom_line(aes(x=NBRtree.ha,color=country,y=..count..),size=1.25,key_glyph = "abline",stat = "density")+
  labs(y=NULL, x="trees/ha",las=1)+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.85,0.55),legend.direction ="vertical",
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.9, "cm"),
        panel.border = element_rect(colour = "black",fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
dplot2


dplot3 <- ggplot(tfinal.biotic) + 
  stat_density(aes(x=log(1+NBRtreedead.plot),color=country),size=1.25,key_glyph = "abline",geom = "line",trim = T,adjust=1)+
  labs(y=NULL,x="dead trees/plot (log)",las=1)+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.85,0.55),legend.direction ="vertical",
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.9, "cm"),
        panel.border = element_rect(colour = "black",fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
dplot3


dplot4 <- ggplot(tfinal.biotic) + 
  stat_density(aes(x=log(1+NBRtreedead.ha),color=country),size=1.25,key_glyph = "abline",geom = "line",trim = T,adjust=1)+
  labs(y=NULL,x="dead trees/ha (log)",las=1)+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.85,0.55),legend.direction ="vertical",
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.9, "cm"),
        panel.border = element_rect(colour = "black",fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
dplot4


library(cowplot)
library(grid)
library(gridExtra)
pall<-plot_grid(dplot1,dplot2,dplot3,dplot4,labels = c('a)', 'b)','c)', 'd)'),
            align="hv", label_size = 14,nrow = 2,ncol=2)
pall

y.grob <- textGrob("Frequency (%)", 
                   gp=gpar(fontface="bold", col="black", fontsize=18), rot=90)
#add to plot
g <- arrangeGrob(pall, left = y.grob) #generates g

save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/trees_density.pdf"),plot = g, base_width = 12.37, dpi=800,units = "in",nrow=2)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/trees_densityV2.png"),plot = g, base_width = 12.37, dpi=800,units = "in",nrow=3,ncol=2)


?plot_grid



#### Table simple and interaction effects summary 
BINsimple <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/BIN.Table.Simple.Review.csv")
NBsimple <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/ZT.Table.Simple.Review.csv")



sumNEGBIN <- apply(BINsimple, 2, function(x) sum(grepl("^-.*\\*$",x,fixed=F)==T))
sumPOSBIN <- apply(BINsimple, 2, function(x) sum(grepl("^[0-9].*\\*$",x,fixed=F)==T))

sumNEGNB <- apply(NBsimple, 2, function(x) sum(grepl("^-.*\\*$",x,fixed=F)==T))
sumPOSNB <- apply(NBsimple, 2, function(x) sum(grepl("^[0-9].*\\*$",x,fixed=F)==T))


A <- cbind(sumNEGBIN,sumPOSBIN,sumNEGNB,sumPOSNB)



BINsimple <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/BIN.Table.Inter.Review.csv")
NBsimple <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/ZT.Table.Inter.Review.csv")



sumNEGBIN <- apply(BINsimple, 2, function(x) sum(grepl("^-.*\\*$",x,fixed=F)==T))
sumPOSBIN <- apply(BINsimple, 2, function(x) sum(grepl("^[0-9].*\\*$",x,fixed=F)==T))

sumNEGNB <- apply(NBsimple, 2, function(x) sum(grepl("^-.*\\*$",x,fixed=F)==T))
sumPOSNB <- apply(NBsimple, 2, function(x) sum(grepl("^[0-9].*\\*$",x,fixed=F)==T))


NBsimple[,"Precipitation.X.BA.O"] <- NA

A <- cbind(sumNEGBIN,sumPOSBIN)
B <- cbind(sumNEGNB,sumPOSNB)



