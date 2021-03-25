#Alex on the 08/01/2018
# Script to load and extract species disribution range from euforgen for any species 
SavingInfo = "ARGUMENTS
# 'CODE' character giving the species we want to extract. See AllspFinal for the full list of species. 
# 'SaveALL' Logical to save the maps
# 'MyPointsAll' Logical to add my plots on the maps
# Exemple : Spatialisation(CODE=species,SaveALL=Y,MyPointsAll=F)

# This scirpt also loads the function Plotmaps which is called in the function spatialiation to record al the maps
# This function needs a few parameters that are loaded with this script 
1) Dir = c(/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/) # Directory in which we are located before running the script
2) Allspfinal <- read.csv2(paste0(Dir,Allspfinal.csv)) # This database is the one containing all the species 
3) path =  c(paste0(Dir,climate/)) # Path of the climatic data... Not useful here

4) Para to extract among the 21 variables. 
bioclim=paste0(bio,c(1,2,5,6,12,13,14))# WorldClim type parameters
T.seas=paste0(tmean.,c(djf,mam,jja,son))# seasonal mean temperature
P.seas=paste0(prec.,c(djf,mam,jja,son))# seasonal precipitation
pet=paste0(pet.,c(mean,min,max))# potential evapotranspiration
ppet=paste0(ppet.,c(mean,min,max))# water balance
eumedclim.vars=c(bioclim,T.seas,P.seas,pet,ppet)# all

=> eumedclim.vars = eumedclim.vars[c(1:21)] #Here we can specify the files we want to extract. By defaut, all 21 variables are extracted 
varclim <- eumedclim.vars

5) Para to plot the map
type <- c(_T_Dist,_T,_T_Tst,..) # This corresponds to the kind of map we want to extract among The full, the cut map or the cut ans scaled map

The second function : 
Plot.map(CODE = PINSYL,                    The species we want
          Climvar = eumedclim.vars[1],      The parameter to plot 
          Maptype=type[1],                  What map we want to plot (Let 1 By defaut)
          top = T,                          Should we added europe on the top of the map
          Save,                             Should we save the map
          myPoints                          Should we add the points)

"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))

#Load libraries
require(rgdal)
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

#########################################
##### Load and dl Sp Distribution #######
#########################################

Dir = c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/") # Directory 
Allspfinal <- read.csv2(paste0(Dir,"Allspfinal.csv"))
path =  c(paste0(Dir,"climate/")) # Path of the climatic data

# Para to extract
bioclim=paste0("bio",c(1,2,5,6,12,13,14))# WorldClim type parameters
T.seas=paste0("tmean.",c("djf","mam","jja","son"))# seasonal mean temperature
P.seas=paste0("prec.",c("djf","mam","jja","son"))# seasonal precipitation
pet=paste0("pet.",c("mean","min","max"))# potential evapotranspiration
ppet=paste0("ppet.",c("mean","min","max"))# water balance
eumedclim.vars=c(bioclim,T.seas,P.seas,pet,ppet)# all
eumedclim.vars = eumedclim.vars[c(1:21)] #Here we can specify the files we want to extract

# Para to plot the map
type <- c("_T_Dist","_T","_T_Tst","")
newext <<- c(-15, 45, 35, 70) #deined in global variable 

################################################################
####   Function to plot any of the previously loaded map : #####
################################################################

Plot.map = function(CODE,Climvar = eumedclim.vars[1], Maptype=type[1], top = T, Save, myPoints){
  
  if(Save==T){ 
    dir.create(path=paste0(Dir,"species/",CODE,"/"))
    jpeg(file=paste0(Dir,"species/",CODE,"/",CODE,"_",Climvar,"_",Maptype,".jpeg"),width=750)}

  Mymap <- paste0(Climvar,Maptype)
  plot(get(Mymap),col=tim.colors(),bigplot=c(.1,.75,.1,.95),legend=F,axes=F,box=F,asp=NA)
  
  # box
  if (str_detect(Mymap,"_T")==T) extent.map=c(-15, 45, 35, 70) else extent.map=c(-20,60,20,72)
  axis(1,extent.map[1:2],labels=F,lwd.ticks=0,pos=extent.map[3])
  axis(2,extent.map[3:4],labels=F,lwd.ticks=0,pos=extent.map[1])
  axis(3,extent.map[1:2],labels=F,lwd.ticks=0,pos=extent.map[4])
  axis(4,extent.map[3:4],labels=F,lwd.ticks=0,pos=extent.map[2])
  # axes
  xmn.lab=seq(-15,45,5)[(seq(-15,45,5)-extent.map[1])>0][1]
  ymn.lab=seq(35,70,5)[(seq(35,70,5)-extent.map[3])>0][1]
  axis(1,seq(xmn.lab,extent.map[2],5),tck=-.02,lwd.ticks=1,lwd=0,cex.axis=.8,pos=extent.map[3])
  axis(2,seq(ymn.lab,extent.map[4],5),tck=-.02,lwd.ticks=1,lwd=0,cex.axis=.8,pos=extent.map[1])
  # legend
  if(any(substr(Mymap,1,5)==c("bio12","bio13","bio14","pet.m","prec.","ppet."))) unit.var="(mm)" else unit.var='(°C)'
  plot(get(Mymap),col=tim.colors(),legend.only=T,smallplot=c(.80,.82,.25,.75),legend.args=list(text=paste(Climvar,unit.var),side=3,adj=.2,font=2,line=1.5,cex=1))
  
  if (top){
    worldmap <- getMap(resolution = "high")
    europe <- worldmap[which(worldmap$REGION=="Europe"),]
    europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & as.character(worldmap$NAME) != "Russia"),]
    europe <-spTransform(europe,CRS(proj4string(get(Mymap)))) # Convert to the right system
    europe = crop(europe,newext)
    plot(europe,add=T)}
  
  # Ajout des points de mon espèce for a giving map for instance
  if (myPoints==T){
    A=summary_all[summary_all$Clim_Var==Climvar,c(2:4)]
    points(c(A[,1]),A[,2],pch=3,cex=0.1)
  }
  
  if(Save==T){ 
    #dir.create(path=paste0(Dir,"species/",CODE,"/"))
    #dev.print(file=paste0(CODE,"_",Climvar,"_",Maptype,".jpeg"),device=jpeg,width=1100)
    dev.off()
  }
}


Spatialisation <- function(CODE,SaveAll = T,myPointsAll = F){
  
  if (myPointsAll==T){load(paste0(Dir,"species/",CODE,"/",CODE,"_summary.RData"),envir = globalenv())}
  
  Binom = as.character(Allspfinal[Allspfinal$code==CODE,"binomial"]) #Find the corresponding binomial name from the code of the species
  if (length(Binom)>1){
    Binom <- Binom[1]
    Binom <- unlist(strsplit(Binom,"_"))
    Binom <- paste0(Binom[1],"_",Binom[2])
  }
  BinSp <- sub("_"," ", Binom, ignore.case = FALSE,fixed = T)
  n <- length(list.files(paste0(Dir,"species/",CODE,"/"),pattern = ".shp$",full.names = T))
  if (n==0){
    if (dir.exists(paste0(Dir,"species/species.maps/",BinSp))==T){
      file.copy(from=list.files(paste0(Dir,"species/species.maps/",BinSp,"/shapefiles/"),pattern="_plg.shp$",full.names = T), #Problem if subspecies
                to = paste0(Dir,"species/",CODE), overwrite = FALSE,recursive = TRUE,copy.mode = TRUE)
      file.copy(from=list.files(paste0(Dir,"species/species.maps/",BinSp,"/shapefiles/"),pattern="_plg.shx$",full.names = T), #Problem if subspecies
                to = paste0(Dir,"species/",CODE), overwrite = FALSE,recursive = TRUE,copy.mode = TRUE)
      file.copy(from=list.files(paste0(Dir,"species/species.maps/",BinSp,"/shapefiles/"),pattern="_plg.dbf$",full.names = T), #Problem if subspecies
                to = paste0(Dir,"species/",CODE), overwrite = FALSE,recursive = TRUE,copy.mode = TRUE)
      file.copy(from=list.files(paste0(Dir,"species/species.maps/",BinSp,"/shapefiles/"),pattern="_plg.prj$",full.names = T), #Problem if subspecies
                to = paste0(Dir,"species/",CODE), overwrite = FALSE,recursive = TRUE,copy.mode = TRUE)
      message("The Species distribution has been downloaded from the paper of the italian Guy")
    } else {
      A <- try(download.file(url=paste0('http://www.euforgen.org/fileadmin/templates/euforgen.org/upload/Documents/Maps/Shapefile/',Binom,'.zip'),
                             destfile=paste0(Dir,"species/",CODE,"/",Binom,".zip"),method='auto',cacheOK=T),silent=T)
      if (A==0) {unzip(paste0(Dir,"species/",CODE,"/",Binom,".zip"), exdir = paste0(Dir,"species/",CODE,"/"))
        message(paste0("The Species distribution of",BinSp," has been downloaded from Euforgen because not available in the database"))
      } else stop(paste("The species ",BinSp,'is not provided by Euforgen nor in the other database'))}
  } else {message(paste0("This species ",BinSp," was already used gros"))}
  setwd(paste0(Dir,"species/",CODE,"/"))
  if (n>1){
    A = list(NULL)
    for (i in 1:n) {A[[i]] <- shapefile(list.files(path=paste0(Dir,"species/",CODE),pattern = ".shp$",full.names = T)[i])}
    SpDistr <- do.call(bind,A)
  }else SpDistr <- shapefile(list.files(path=paste0(Dir,"species/",CODE),pattern = ".shp$",full.names = T))
  
  ###############################################################
  ###   Load all our rasters and give them the right name   ##### 
  ###############################################################
  
  # Extract the rasters we want and give the right name.Also give it a new extend and scale. The output is three rasters for each variables.
  i = 1
  for (i in 1:length(eumedclim.vars)) {
    assign(paste0(eumedclim.vars[i]),raster(list.files(path=path,pattern=paste0("map_",eumedclim.vars[i],"_1901"), full.names=T)),envir = globalenv()) #Map normal
    assign(paste0(eumedclim.vars[i],"_T"),crop(get(eumedclim.vars[i]),newext),envir = globalenv())                                                 #Map croped
    assign(paste0(eumedclim.vars[i],"_T_Tst"),scale(get(paste0(eumedclim.vars[i],"_T")),center = TRUE, scale = TRUE),envir = globalenv())          #Map crop and scaled
    i = i+1
  }
  ####Mask for the 21 variables. THIS OPERATION CAN TAKE LONG. 
  if (compareCRS(SpDistr,get(paste0(eumedclim.vars[1],"_T")), unknown=FALSE, verbatim=FALSE, verbose=TRUE)==F)
  {SpDistr <-spTransform(SpDistr,CRS(proj4string(get(paste0(eumedclim.vars[1],"_T")))))}
  try(writeOGR(SpDistr,"SpDistr",layer=CODE,driver="ESRI Shapefile"),silent=T)
  i = 1
  for (i in 1:length(eumedclim.vars)) {
  assign(paste0(eumedclim.vars[i],"_T_Dist"),mask(get(paste0(eumedclim.vars[i],"_T")),SpDistr),envir = globalenv()) #Mask with the species distribution. Maybe shoud do this later (after the climatic niche characterisation). Problem is we did not extratc the cliamte for all the repartition of the species
    i = i + 1
  }
  
  #Once loaded on the memory, we load it all on the disk
  #We could do that for any of the other rasters created with the previous function
  i = 1
  for (i in 1:length(eumedclim.vars)) {
    writeRaster(get(paste0(eumedclim.vars[i],"_T_Dist")),filename=paste0(eumedclim.vars[i],"_T_Dist.tif"),options="INTERLEAVE=BAND",overwrite=T)
    i = i + 1
  }
  #Save all climatic variables map
  i = 1
  if (SaveAll == T & myPointsAll == T){
    for (i in 1:length(eumedclim.vars)){
      Plot.map(CODE = CODE, Climvar = eumedclim.vars[i], Maptype = type[1], top = T, Save = T, myPoints = T)
      i = i+1}
  }else if (SaveAll == T & myPointsAll == F){
    for (i in 1:length(eumedclim.vars)){
      Plot.map(CODE = CODE, Climvar = eumedclim.vars[i], Maptype = type[1], top = T, Save = T, myPoints = F)
      i = i+1}
  }else if (SaveAll == F & myPointsAll == T){
    for (i in 1:length(eumedclim.vars)){
      Plot.map(CODE = CODE, Climvar = eumedclim.vars[i], Maptype = type[1], top = T, Save = F, myPoints = T)
      i = i+1}
  }else if (SaveAll == F & myPointsAll == F){
    for (i in 1:length(eumedclim.vars)){
      Plot.map(CODE = CODE, Climvar = eumedclim.vars[i], Maptype = type[1], top = T, Save = F, myPoints = F)
      i = i+1}}
  gc()
  dev.off()
}

