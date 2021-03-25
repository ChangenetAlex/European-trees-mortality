#Alex on the 15/01/2018
# Script to run and analyse CPA with the distribution of our species
#Load libraries
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


SavingInfo = "ARGUMENTS

# 'Dir is a character. This is the path in which different folders and tables are stored. Already defined before to run the function
# 'CODE' character giving the species we want to extract. See AllspFinal for the full list of species. 
# 'nsample' number of sample to take to run the PCA in a first time
# 'NF' Number of axes to be kept after the PCA 
# 'seuil' above this threshold, you are a marginal individual ! 
# 'seuilC' below this threshold, you are a core individual ! 
# 'save' logical. If TRUE, lot of table and plots are going to be saved.

# Exemple : 
Acp1000('PINSYL',nsample=10000,NF=2,seuil=0.8,seuilC=0.6,save=T)
Marginality_Levels(seuil=0.8,seuilC=0.6,save=T)

"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))



#########################################
##### Load and dl Sp Distribution #######
#########################################


#CODE = "PINPINA" # Sp we want to extract
Dir = c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/") # Directory 

# Para
bioclim=paste0("bio",c(1,2,5,6,12,13,14))# WorldClim type parameters
T.seas=paste0("tmean.",c("djf","mam","jja","son"))# seasonal mean temperature
P.seas=paste0("prec.",c("djf","mam","jja","son"))# seasonal precipitation
pet=paste0("pet.",c("mean","min","max"))# potential evapotranspiration
ppet=paste0("ppet.",c("mean","min","max"))# water balance
eumedclim.vars=c(bioclim,T.seas,P.seas,pet,ppet)# all
eumedclim.vars = eumedclim.vars[c(1:21)] #Here we can specify the files we want to extract

################################################################
####   Perform a PCA with dudiPCA to obtain good figures   ####
################################################################

Acp1000 <- function(CODE,nsample=10000,NF=2,seuil=0.8,seuilC=0.6,save=T){
  dir.create(paste0(Dir,"species/",CODE,"/PCA/")) #Create a specific file 
  #Select all variables and perform a PCA : 
  mesrasters = list.files(paste0(Dir,"species/",CODE,""),pattern=".tif", full.names=T)
  map <- stack(list.files(path=paste0(Dir,"species/",CODE,"/"),pattern=".tif", full.names=T)) #map made of 21 layers
  set.seed(5) #repeatability
  smap <- sampleRandom(map, nsample) # sample 10000 random grid cells # Mettre aussi en argument 
  # write PCA model to file if we want to 
  #Run the PCA
  acp=dudi.pca(smap,scannf=FALSE,nf=NF,center=TRUE,scale=TRUE) # Ajout de nf en argument 
  # Inertie :
  inertie<-acp$eig/sum(acp$eig)*100
  if (save==T) jpeg(file=paste0(Dir,"species/",CODE,"/PCA/Eigenvalues.jpeg"),width=750)
  barplot(inertie, names.arg=round(inertie,1), ann=FALSE)
  title(ylab="% d'inertie")
  title("Eboulis des valeurs propres")
  acp$eig[c(1,2)]/sum(acp$eig)*100 #var explained
  if (save==T) dev.off()
  sum(acp$eig[c(1,2)]/sum(acp$eig))*100 #Cumul 
  acp$co #resume des composantes
  acp$li #nouvelles donnee selon les 2 axes
  
  # Plot du cumul du % d'inertie
  par(cex=1)
  cumul<-(cumsum(acp$eig/sum(acp$eig)*100))
  if (save==T) jpeg(file=paste0(Dir,"species/",CODE,"/PCA/Eigen_Cumul.jpeg"),width=750)
  plot(cumul, type="o", ann=FALSE)
  title(ylab="% d'inertie")
  title("Cumul du % d'inertie")
  if (save==T) dev.off()
  
  # Cercle des corrélations :
  if (save==T) jpeg(file=paste0(Dir,"species/",CODE,"/PCA/Corcircle.jpeg"),width=750)
  s.corcircle(acp$co)
  title(xlab= "Composante 1")
  title(ylab= "Composante 2")
  title("Cercle des corrélations des variables")
  if (save==T) dev.off() ; jpeg(file=paste0(Dir,"species/",CODE,"/PCA/Variables.jpeg"),width=750)
  plot(acp$co,cex=0.01)
  text(acp$co[,1],acp$co[,2],row.names(acp$co),cex=0.6)
  if (save==T) dev.off()
  
  #Individual representation and margins : Using normed or on normed data is in fact the same (li or l1)
  
  Marginality=apply(acp$li,1,function(x) sum(abs(x)*c(acp$eig[1:2]/sum(acp$eig)*100)))
  Margin = ifelse(Marginality>quantile(Marginality,seuil,type=7),T,F) # Par défaut 0.8 mais a mettre en parametre
  Core = ifelse(Marginality<quantile(Marginality,seuilC,type=7),T,F) ## Par defaut 0.5 mais a mettre en parametre 
  res=data.frame(acp$li,Marginality,Margin,Core)
  
  # Marginalité, scores and margins 
  if (save==T) jpeg(file=paste0(Dir,"species/",CODE,"/PCA/Indiv.Margin_Core.jpeg"),width=750)
  s.label(acp$li,label=NULL)
  points(acp$li[,1:2],pch=19,col=gray(.7),cex=.5,xpd=T)
  points(res[,1:2],pch=19,col=ifelse(res$Margin,3,1),cex=.5,xpd=T)
  points(res[,1:2],pch=19,col=ifelse(res$Margin==F & res$Core==T,gray(.7),ifelse(res$Core==F & res$Margin==F,gray(.5),gray(.1))),cex=.5,xpd=T)
  s.arrow(acp$co*23,lab=names(acp$tab),boxes=F,clabel=.8,add.plot=T)
  iner=acp$eig[1:4]/sum(acp$eig)
  var.exp=paste0('axis-',1:4,': ',round(iner*100),'%')
  legend('bottomright',bty='n',inset=c(0,0),legend=var.exp,title='inertia',xpd=T,cex=.8)
  if (save==T) dev.off()
  saveRDS(acp, file=paste0(Dir,"species/",CODE,"/PCA/acp10000.rds")) # Save the results of the ACP 
}

#######################################
#####                           #######
#####   Function to obtain      #######
#####     Core and Climatic     #######
#####         Margins           #######
#####                           #######
#######################################

Marginality_Levels <- function(CODE,seuil=0.8,seuilC=0.6,save=T){
  mesrasters = list.files(paste0(Dir,"species/",CODE,""),pattern=".tif", full.names=T)
  map <- stack(list.files(path=paste0(Dir,"species/",CODE,"/"),pattern=".tif", full.names=T)) #map made of 21 layers
  A=as.data.frame(x = map,xy=TRUE,na.rm=TRUE) #Concertir mes rasters en df
  if (save==T) saveRDS(A, file=paste0(Dir,"species/",CODE,"/PCA/RasterAsDf_",CODE,".rds")) # Save the raster as a dataframe  
  acp <- readRDS(file=paste0(Dir,"species/",CODE,"/PCA/acp10000.rds")) # Load the acp10000 results (must check if this works)
  envsup.acp <- suprow(acp, A[,3:23]) #Projeter individus supplémentaires
  par(mfrow=c(1,1))
  Marginality=apply(envsup.acp$lisup,1,function(x) sum(abs(x)*c(acp$eig[1:2]/sum(acp$eig)*100)))
  Margin = ifelse(Marginality>quantile(Marginality,seuil,type=7),T,F) # Margins VS Core
  Core = ifelse(Marginality<quantile(Marginality,seuilC,type=7),T,F)
  res=data.frame(envsup.acp$lisup,Marginality,Margin,Core) #Put together individuals, and marginality
  if (save==T){
    capture.output(print(summary(as.factor(Margin))), file=paste0(Dir,"species/",CODE,"/PCA/MarginCore.Table.txt")) # Output as a latex wrapped in a txt file
    capture.output(print(summary(as.factor(Core))), file=paste0(Dir,"species/",CODE,"/PCA/MarginCore.Table.txt"),append = T) # Output as a latex wrapped in a txt file
  }
  dir.create(paste0(Dir,"species/",CODE,"/CLUSTER/")) #Create a specific file 
  if (save==T) jpeg(file=paste0(Dir,"species/",CODE,"/CLUSTER/Indiv.Margin.Core.AllPoints_",seuil,"_",seuilC,".jpeg"),width=750)
  s.label(envsup.acp$lisup,label=NULL)
  points(res[,1:2],pch=19,col=ifelse(res$Margin==F & res$Core==T,gray(.7),ifelse(res$Core==F & res$Margin==F,gray(.5),gray(.1))),cex=.5,xpd=T)
  s.arrow(acp$co*23,lab=names(acp$tab),boxes=F,clabel=.8,add.plot=T) #Plot these 2 categories of individuals
  iner=acp$eig[1:4]/sum(acp$eig)
  var.exp=paste0('axis-',1:4,': ',round(iner*100),'%')
  legend('bottomright',bty='n',inset=c(0,0),legend=var.exp,title='inertia',xpd=T,cex=.8)
  if (save==T) dev.off()
  
  #################################
  #### Analyse discriminante ######
  #################################
  
  Final <- cbind(A[,1:23],res) #coordonées géo + scores + all variables + marginalité 
  FinalMargin <- Final[Final$Margin=="TRUE",] #Select only margin individuals 
  GrMarg = find.clusters.data.frame(FinalMargin[,-c(1:2,24:28)],
                                    n.clust=2,stat = c("BIC", "AIC", "WSS"),
                                    scale=TRUE,n.pca=50)
  #Scale, Normed clustering. Force to 2 groups
  
  dapc1<-dapc(FinalMargin[,-c(1:2,24:28)], GrMarg$grp,n.pca=21,n.da=1) #Keep one function and all axes
  # Save as rds HERE 
  if (save==T){capture.output(print(summary(dapc1)),file=paste0(Dir,"species/",CODE,"/CLUSTER/DAPC_",seuil,"_",seuilC,".summary.txt"))
    jpeg(file=paste0(Dir,"species/",CODE,"/CLUSTER/Discrimine_",seuil,"_",seuilC,".jpeg"),width=750)}
  scatter(dapc1, scree.da=TRUE, scree.pca =TRUE, bg="white", pch=20, cstar=0, col=c("blue","red"), solid=.4,
          cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:2),cell=2)
  if(save==T) dev.off()
  
  #Put together the core individuals and the new splitted margins 
  
  test=FinalMargin #All our previous table (marginality)
  test$groupes <- dapc1$grp[match(row.names(test), names(dapc1$grp))] #Give them their attributed group
  test2=Final[Final$Margin=="FALSE"&Final$Core=="TRUE",] #All our core individuals
  test2bis=Final[Final$Margin=="FALSE"&Final$Core=="FALSE",] #All our transition zone
  test3=c(rep(0,nrow(test2))) # New column whose length is equal to the number of core individuals and filled with 0
  test3bis=c(rep(10,nrow(test2bis))) # New column whose length is equal to the number of transition individuals and filled with 10
  test4=cbind(test2,test3) #New table made of these two previous tables
  test4bis=cbind(test2bis,test3bis) #New table made of these two previous tables
  colnames(test4)[29] <- "groupes" # Renamed last column
  colnames(test4bis)[29] <- "groupes" # Renamed last column
  test5=rbind(test4,test4bis,test) # Put together core and margins 
  if (save==T) capture.output(print(summary(as.factor(test5$groupes))),file=paste0(Dir,"species/",CODE,"/CLUSTER/Groupes_",seuil,"_",seuilC,
                                                                                   ".summary.txt"))
  #Save the summary 
  test6 <- test5[order(as.numeric(row.names(test5))),] #Order our cells to trasnform it to a raster
  test6[,29] <- as.numeric(test6[,29]) # Put it numeric (raster's requirement)
  ras <- rasterFromXYZ(test6[,c("x","y","groupes")],crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ") #New raster with our three categories 
  
  extent.map=c(-15, 45, 35, 70)
  ras <- extend(ras,extent.map) # to keep the same area
  # test 
  ras <- ratify(ras, count=T)
  rat <- levels(ras)[[1]]
  rat$landcover <- c('Core', 'Edge1', 'Edge2',"Transition")
  levels(ras) <- rat
  assign(paste0("CLIMAP_",seuil,"_",seuilC),ras,envir = .GlobalEnv)
  
  if (save == T) jpeg(file=paste0(Dir,"species/",CODE,"/CLUSTER/Landscape_",seuil,"_",seuilC,".jpeg"),width=1100,height=1044)
  plot(get(paste0("CLIMAP_",seuil,"_",seuilC)),"landcover",col=c("yellow",'red',"blue","gray"),main = paste0("Climatic Core, Edges and transition zones",seuil,"_",seuilC),asp=NA)
  worldmap <- getMap(resolution = "high")
  europe <- worldmap[which(worldmap$REGION=="Europe"),]
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & as.character(worldmap$NAME) != "Russia"),] 
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3)),] 
  europe <-spTransform(europe,CRS(proj4string(get(paste0("CLIMAP_",seuil,"_",seuilC))))) # Convert to the right system
  europe = crop(europe,extent.map)
  plot(europe,add=T)
  if (save == T) dev.off()
  
  # Save the table and the raster before cleaning memory 
  
  dir.create(paste0(Dir,"species/",CODE,"/CLIMAP/"))
  save(test6,file=paste0(Dir,"species/",CODE,"/CLIMAP/Climap_",CODE,"_",seuil,"_",seuilC,".RData"))
  writeRaster(get(paste0("CLIMAP_",seuil,"_",seuilC)),filename=paste0(Dir,"species/",CODE,"/CLIMAP/Climap_",seuil,"_",seuilC,".tif"),options="INTERLEAVE=BAND",overwrite=T)
  
  ###################################   !!!!!!!
  # ADDED on THE 29th OF JANUARY ####   !!!!!!!
  ###################################   !!!!!!!
  
  # In order to keep the attributes table, it is necessary to save the file as a .grd with ' bandorder="bil" ' instead of options = ...
  # This generates two different files : .grd and .gri

  
  # line to remove the test 6 object
  # Ajout des points de mon espèce for a giving map for instance
  load(paste0(Dir,"species/",CODE,"/",CODE,"_summary.RData"),envir = globalenv())
  load(paste0(Dir,"species/",CODE,"/CLIMAP/Climap_",CODE,"_",seuil,"_",seuilC,".RData"),envir = globalenv())
  
  #Determine the color according to the number
  A=summary_all[summary_all$Clim_Var=="bio1",c(2:4)] # my plots
  Plotcat <- raster::extract(get(paste0("CLIMAP_",seuil,"_",seuilC)),A[,1:2])
  capture.output(print(summary(as.factor(Plotcat))), file=paste0(Dir,"species/",CODE,"/CLIMAP/MarginCore.Table.txt"),append=T) # Summary number of plots !!! 
  #capture.output(print(summary(as.factor(Plotcat))), file=paste0(Dir,"species/MarginCore.Table.txt"),append = T) # Summary number of plots !!! 
  Myplot = cbind(A,Plotcat)
  if (min(Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),"latitude"])<min(Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),"latitude"])
      & max(Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),"latitude"])>max(Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),"latitude"]))
  {Mycol = c("yellow","red","blue","gray")
  } else if (min(Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),"latitude"])<min(Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),"latitude"])
             & max(Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),"latitude"])>max(Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),"latitude"]))
  {Mycol = c("yellow","blue","red","gray")
  } else if (mean(Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),"latitude"])<mean(Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),"latitude"]))
  {Mycol = c("yellow","blue","red","gray")
  } else if (mean(Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),"latitude"])<mean(Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),"latitude"]))
  {Mycol = c("yellow","red","blue","gray")
  message(paste0("The Species distribution of",CODE," was confusing. We based our definition on the average latitude for the two zones."))}
  
  # Plot a map with all we want : 
  par(mfrow=c(1,1))
  if (save == T) jpeg(file=paste0(Dir,"species/",CODE,"/CLIMAP/Myplots_",seuil,"_",seuilC,".jpeg"),width=800)
  plot(get(paste0("CLIMAP_",seuil,"_",seuilC)),"landcover",col=Mycol,bigplot=c(.1,.75,.1,.95),legend=F,axes=F,box=F,main=c("My plots"),asp=NA)
  # box
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
  plot(get(paste0("CLIMAP_",seuil,"_",seuilC)),"landcover",col=Mycol,legend.only=T,smallplot=c(.80,.82,.25,.75),legend.args=list(text="Climatic Areas",side=3,adj=.2,font=2,line=1.5,cex=1),asp=NA)
  worldmap <- getMap(resolution = "high")
  europe <- worldmap[which(worldmap$REGION=="Europe"),]
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & as.character(worldmap$NAME) != "Russia"),] 
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3)),] 
  europe <-spTransform(europe,CRS(proj4string(get(paste0("CLIMAP_",seuil,"_",seuilC))))) # Convert to the right system
  europe = crop(europe,extent.map)
  plot(europe,add=T)
  #Add all our points 
  points(c(A[,1]),A[,2],pch=3,cex=0.01)
  #Save the map 
  dev.off()
  
  #Extraction de mes plots par catégories
  Core = Myplot[Myplot$Plotcat=="0"&!is.na(Myplot$Plotcat),]
  Transit = Myplot[Myplot$Plotcat=="10"&!is.na(Myplot$Plotcat),]
  if (min(Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),"latitude"])<min(Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),"latitude"])
      & max(Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),"latitude"])>max(Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),"latitude"]))
    {LE = Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),]
    RE = Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),]
    } else if (min(Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),"latitude"])<min(Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),"latitude"])
              & max(Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),"latitude"])>max(Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),"latitude"]))
    {LE = Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),]
    RE = Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),]
    } else if (mean(Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),"latitude"])<mean(Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),"latitude"]))
      {LE = Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),]
      RE = Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),]
    } else if (mean(Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),"latitude"])<mean(Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),"latitude"]))
    {LE = Myplot[Myplot$Plotcat=="2"&!is.na(Myplot$Plotcat),]
    RE = Myplot[Myplot$Plotcat=="1"&!is.na(Myplot$Plotcat),]}
  Nodata = Myplot[is.na(Myplot$Plotcat),]
  
  #Plot all three and save
  if (save==T) jpeg(file=paste0(Dir,"species/",CODE,"/CLIMAP/Myplots_Cat_",seuil,"_",seuilC,".jpeg"),width=1100,height=1044)
  par(mfrow=c(2,2))
  
  plot(get(paste0("CLIMAP_",seuil,"_",seuilC)),"landcover",col=Mycol,main=c(paste0("Core ",seuil,"_",seuilC," : ",nrow(Core)," plots")),cex.main=2,cex.lab=1,cex.axis=1,legend=F,las=1,asp=NA)
  points(c(Core[,1]),Core[,2],pch=3,cex=0.01)
  plot(get(paste0("CLIMAP_",seuil,"_",seuilC)),"landcover",col=Mycol,legend.only=T,legend.width=1,legend.shrink=1,legend.args=list(text="Area",side=4,adj=0.5,font=1,line=2,cex=1),axis.args=list(cex.axis=1),asp=NA)
  plot(europe,add=T)
  
  plot(get(paste0("CLIMAP_",seuil,"_",seuilC)),"landcover",col=Mycol,main=c(paste0("Lead edge ",seuil,"_",seuilC," : ",nrow(LE)," plots")),cex.main=2,cex.lab=1,cex.axis=1,legend=F,las=1,asp=NA)
  points(c(LE[,1]),LE[,2],pch=3,cex=1)
  plot(get(paste0("CLIMAP_",seuil,"_",seuilC)),"landcover",col=Mycol,legend.only=T,legend.width=1,legend.shrink=1,legend.args=list(text="Area",side=4,adj=0.5,font=1,line=2,cex=1),axis.args=list(cex.axis=1),asp=NA)
  plot(europe,add=T)
  
  plot(get(paste0("CLIMAP_",seuil,"_",seuilC)),"landcover",col=Mycol,main=c(paste0("Rear edge ",seuil,"_",seuilC," : ",nrow(RE)," plots")),cex.main=2,cex.lab=1,cex.axis=1,legend=F,las=1,asp=NA)
  points(c(RE[,1]),RE[,2],pch=3,cex=0.01)
  plot(get(paste0("CLIMAP_",seuil,"_",seuilC)),"landcover",col=Mycol,legend.only=T,legend.width=1,legend.shrink=1,legend.args=list(text="Area",side=4,adj=0.5,font=1,line=2,cex=1),axis.args=list(cex.axis=1),asp=NA)
  plot(europe,add=T)
  
  plot(get(paste0("CLIMAP_",seuil,"_",seuilC)),"landcover",col=Mycol,main=c(paste0("Transition zone ",seuil,"_",seuilC," : ",nrow(Transit)," plots")),cex.main=2,cex.lab=1,cex.axis=1,legend=F,las=1,asp=NA)
  points(c(Transit[,1]),Transit[,2],pch=3,cex=0.01)
  plot(get(paste0("CLIMAP_",seuil,"_",seuilC)),"landcover",col=Mycol,legend.only=T,legend.width=1,legend.shrink=1,legend.args=list(text="Area",side=4,adj=0.5,font=1,line=2,cex=1),axis.args=list(cex.axis=1),asp=NA)
  plot(europe,add=T)
  dev.off()
  
  if (save==T) jpeg(file=paste0(Dir,"species/",CODE,"/CLIMAP/Myplots_Cat_NO.DATA.",seuil,"_",seuilC,".jpeg"),width=1750,height=830)
  par(mfrow=c(1,2))
  plot(get(paste0("CLIMAP_",seuil,"_",seuilC)),"landcover",col=Mycol,main=c(paste0("Core, LE and RE ",seuil,"_",seuilC," : ",nrow(Core)+nrow(LE)+nrow(RE)," plots")),cex.main=2,cex.lab=1,cex.axis=1,legend=F,las=1,asp=NA)
  points(c(Core[,1],LE[,1],RE[,1]),c(Core[,2],LE[,2],RE[,2]),pch=3,cex=0.01)
  plot(get(paste0("CLIMAP_",seuil,"_",seuilC)),"landcover",col=Mycol,legend.only=T,legend.width=1,legend.shrink=1,legend.args=list(text="Area",side=4,adj=0.5,font=1,line=2,cex=1),axis.args=list(cex.axis=1),asp=NA)
  plot(europe,add=T)
  
  plot(get(paste0("CLIMAP_",seuil,"_",seuilC)),"landcover",col=Mycol,main=c(paste0("No Cat ",seuil," : ",nrow(Nodata)," plots")),cex.main=2,cex.lab=1,cex.axis=1,legend=F,las=1,asp=NA)
  points(c(Nodata[,1]),Nodata[,2],pch=3,cex=0.01)
  plot(get(paste0("CLIMAP_",seuil,"_",seuilC)),"landcover",col=Mycol,legend.only=T,legend.width=1,legend.shrink=1,legend.args=list(text="Area",side=4,adj=0.5,font=1,line=2,cex=1),axis.args=list(cex.axis=1),asp=NA)
  plot(europe,add=T)
  dev.off()
  
  
  # Merge these cateagories with all my informations and save it 
  rownames(Myplot) <- Myplot$plotcode
  colnames(Myplot)[1:2] <- c("longitude","latitude")
  
  #Load the two dataset. The one we want is made of every single sumarized variable (the second)
  load(paste0(Dir,"species/",CODE,"/",CODE,"_allVariable.RData"),envir = globalenv())
  
  #Merge climate data with inventory individual data
  assign(paste0("Mydf_",CODE,"_",seuil,"_",seuilC),merge(Myplot,Sp2,by=c("plotcode")),envir = globalenv())
  saveRDS(get(paste0("Mydf_",CODE,"_",seuil,"_",seuilC)), file = paste0(Dir,"species/",CODE,"/CLIMAP/Mydf_",CODE,"_",seuil,"_",seuilC,".rds"))
  rm(list = ls())
  gc()
}


#Check these are the good points with correlation between latitudes and longitudes 







