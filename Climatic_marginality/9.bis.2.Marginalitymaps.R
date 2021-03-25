# Alex on the 20 December 2018
# Scipt to do all the maps without plotting the points on it (just the marginality)
rm(list = ls())
gc()
require(rgdal)
require(adegenet)
require(ade4)
require(raster)
require(parallel)
require(fields)
library(raster)
library(rgdal) 
library(rworldmap)
library(lattice)
require(spatial.tools)
library(maptools)
require(rworldxtra)
library(rgeos)
library(RStoolbox)
require(stringr)
library(data.table)
library(ggplot2) 
library(rangeBuilder)
Dir = c("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/") # Directory 
i = 1
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
Allcode <- c("PINPINA","ALNGLU","PINNIG","POPTRE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.6,0.7,0.7,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.5,0.6,0.6,0.6,0.6,0.6,0.6)


Allcode= c("QUESUB")
Allseuil <- c(0.7)
AllseuilC <- c(0.6)

size = 2400

for (i in 1:length(Allcode)){
  test <- raster(paste0(Dir,"species/",Allcode[i],"/CLIMAP/Climap_",Allseuil[i],"_",AllseuilC[i],".tif"))
  extent.map=c(-15, 45, 35, 70)
  test <- extend(test,extent.map) # to keep the same area
  test <- ratify(test, count=T)
  rat <- levels(test)[[1]]
  rat$landcover <- c('Core', 'Edge1', 'Edge2',"Transition")
  levels(test) <- rat
  png(file=paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/Fig.Margins/",Allcode[i],"_Marginality_",Allseuil[i],"_",AllseuilC[i],"2400.png"),width=size,height=size*0.9,res=size/12)
  plot(test,"landcover",col=c("yellow",'blue',"red","gray"),main = paste0(Allcode[i]),asp=NA,cex.main=2.5,cex.axis = 1.8, cex.sub = 1.8, cex.lab = 1.8,las=1)
  worldmap <- getMap(resolution = "high")
  europe <- worldmap[which(worldmap$REGION=="Europe"),]
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & as.character(worldmap$NAME) != "Russia"),] 
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3)),] 
  europe <-spTransform(europe,CRS(proj4string(test))) # Convert to the right system
  europe = crop(europe,extent.map)
  plot(europe,add=T)
  dev.off()
}
