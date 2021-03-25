library(knitr)
library(spaMM)
library(fields)
library(piecewiseSEM)
library(raster)
library(rworldmap)
library(rgdal)
library(rasterVis)
library(rgeos)
library(ggplot2)


i = 1
Allcode <- c("ABIALB","ACEPSE","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINSYL","POPNIG","QUEILE")
Allmod <- c("Mbin14A.19","Mbin13B.26","M2bin13B.22","M2bin13A21","M2bin15B23","Mbin3C.26","M2bin7B.17","Mbin11B.19","Mbin5B.17","Mbin3B.31","M2bin1C.20")
for (i in 1:length(Allcode)){
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/binomial/",Allmod[i],"/"))
  setwd(Dir)
  assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
  #rm("x")
  # New section 
  r1 <- data.frame(x$data$longitude, x$data$latitude, x$data$sp.mortality.plot.rate.yr,x$data$sp.mort.bin,x$data$Plotcat)
  colnames(r1) <- c("X","Y","rate.yr","mort.bin","Plotcat")
  #r1$rate.yr <- log(r1$rate.yr+1)
  print(summary(r1$rate.yr))
  r2 <- r1[r1$rate.yr>0,]
  p <- ggplot(r1, aes(x=Plotcat, y=rate.yr,fill=Plotcat)) + geom_violin() +
    geom_boxplot(outlier.colour="black", outlier.shape=8,
                 outlier.size=1)+
    stat_summary(fun.y = "mean",geom="point",shape=23,size=7,fill="white")
  
  p <-p + labs(title=paste0('Observed mortality by region of ',Allcode[i]," ALL plot"), y="Proportion of dead by year", x="Marginality")+
    scale_fill_manual(values=c("yellow", "red", "blue"))+
    scale_x_discrete(labels=c("Core", "RE", "LE"))+
    coord_cartesian(ylim = c(0,12))
  p
  ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Figures.all/Synthesis/Boxplot.margins/",Allcode[i],"_.Allplot.png"),plot = p, width = 12.37, height = 7.04, dpi=150,units = "in")

  p <- ggplot(r2, aes(x=Plotcat, y=rate.yr,fill=Plotcat)) + geom_violin() +
    geom_boxplot(outlier.colour="black", outlier.shape=8,
                 outlier.size=1)+
    stat_summary(fun.y = "mean",geom="point",shape=23,size=7,fill="white")
  
  p <-p + labs(title=paste0('Observed mortality by region of ',Allcode[i]," with dead plots"), y="Proportion of dead by year", x="Marginality")+
    scale_fill_manual(values=c("yellow", "red", "blue"))+
    scale_x_discrete(labels=c("Core", "RE", "LE"))
  p
  ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Figures.all/Synthesis/Boxplot.margins/",Allcode[i],"_.Withdeadplot.png"),plot = p, width = 12.37, height = 7.04, dpi=150,units = "in")
}  
  
  
  
i = 1
Allcode <- c("PINPINA","ALNGLU","PINNIG","PINPIN","POPTRE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allmod <- c("Mbin13A.18","Mbin13A.22","Mbin3B.20","Mbin15B.21","Mbin13A.27","M2bin1B.23","Mbin13B.26","Mbin7A.26","Mbin1B.23")
for (i in 1:length(Allcode)){
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/binomial/",Allmod[i],"/"))
  setwd(Dir)
  assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
  #rm("x")
  # New section 
  r1 <- data.frame(x$data$longitude, x$data$latitude, x$data$sp.mortality.plot.rate.yr,x$data$sp.mort.bin,x$data$Plotcat)
  colnames(r1) <- c("X","Y","rate.yr","mort.bin","Plotcat")
  #r1$rate.yr <- log(r1$rate.yr+1)
  print(summary(r1$rate.yr))
  r2 <- r1[r1$rate.yr>0,]
  p <- ggplot(r1, aes(x=Plotcat, y=rate.yr,fill=Plotcat)) + geom_violin() +
    geom_boxplot(outlier.colour="black", outlier.shape=8,
                 outlier.size=1)+
    stat_summary(fun.y = "mean",geom="point",shape=23,size=7,fill="white")
  p <-p + labs(title=paste0('Observed mortality by region of ',Allcode[i]," ALL plot"), y="Proportion of dead by year", x="Marginality")+
    scale_fill_manual(values=c("yellow", "blue", "red"))+
    scale_x_discrete(labels=c("Core", "LE", "RE"))+
    coord_cartesian(ylim = c(0,12))
  p
  ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Figures.all/Synthesis/Boxplot.margins/",Allcode[i],"_.Allplot.png"),plot = p, width = 12.37, height = 7.04, dpi=150,units = "in")
  p <- ggplot(r2, aes(x=Plotcat, y=rate.yr,fill=Plotcat)) + geom_violin() +
    geom_boxplot(outlier.colour="black", outlier.shape=8,
                 outlier.size=1)+
    stat_summary(fun.y = "mean",geom="point",shape=23,size=7,fill="white")
  p <-p + labs(title=paste0('Observed mortality by region of ',Allcode[i]," with dead plots"), y="Proportion of dead by year", x="Marginality")+
    scale_fill_manual(values=c("yellow", "blue", "red"))+
    scale_x_discrete(labels=c("Core", "LE", "RE"))
  p
  ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Figures.all/Synthesis/Boxplot.margins/",Allcode[i],"_.Withdeadplot.png"),plot = p, width = 12.37, height = 7.04, dpi=150,units = "in")
}    
