library(raster)
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
Dir = c("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/") # Directory 
nsample=10000
NF=2
seuil=0.8
seuilC=0.6
# Para
bioclim=paste0("bio",c(1,2,5,6,12,13,14))# WorldClim type parameters
T.seas=paste0("tmean.",c("djf","mam","jja","son"))# seasonal mean temperature
P.seas=paste0("prec.",c("djf","mam","jja","son"))# seasonal precipitation
pet=paste0("pet.",c("mean","min","max"))# potential evapotranspiration
ppet=paste0("ppet.",c("mean","min","max"))# water balance
eumedclim.vars=c(bioclim,T.seas,P.seas,pet,ppet)# all
eumedclim.vars = eumedclim.vars[c(1:21)] #Here we can specify the files we want to extract
CODE <- c("PINSYL")
acptest <- readRDS(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/PCA/acp10000.rds"))

Marginality=apply(acptest$li,1,function(x) sum(abs(x)*c(acptest$eig[1:2]/sum(acptest$eig)*100)))
Margin = ifelse(Marginality>quantile(Marginality,seuil,type=7),T,F) # Par défaut 0.8 mais a mettre en parametre
Core = ifelse(Marginality<quantile(Marginality,seuilC,type=7),T,F) ## Par defaut 0.5 mais a mettre en parametre 
res=data.frame(acptest$li,Marginality,Margin,Core)

size = 1200
png(file=paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/Indiv.Margin_Core",size,".png"),width=size,height = size*0.9,res=size/9)
s.label(acptest$li,label=NULL)
points(acptest$li[,1:2],pch=19,col=gray(.7),cex=.5,xpd=T)
#points(res[,1:2],pch=19,col=ifelse(res$Margin,3,1),cex=.5,xpd=T)
points(res[,1:2],pch=19,col=ifelse(res$Margin==F & res$Core==T,"yellow",ifelse(res$Core==F & res$Margin==F,"gray75","gray55")),cex=.5,xpd=T)
#s.arrow(acptest$co*23,lab=names(acptest$tab),boxes=F,clabel=1.2,add.plot=T)
iner=acptest$eig[1:2]/sum(acptest$eig)
var.exp=paste0('axis-',1:2,': ',round(iner*100),'%')
legend('bottomright',bty='n',inset=c(0,0),legend=var.exp,title='inertia',xpd=T,cex=1.5)
dev.off()


# Test pour refaire les figures en meilleur résolution 


seuil=0.8
seuilC=0.6

mesrasters = list.files(paste0(Dir,"species/",CODE,""),pattern=".tif", full.names=T)
map <- stack(list.files(path=paste0(Dir,"species/",CODE,"/"),pattern=".tif", full.names=T)) #map made of 21 layers
A=as.data.frame(x = map,xy=TRUE,na.rm=TRUE) #Concertir mes rasters en df
acp <- readRDS(file=paste0(Dir,"species/",CODE,"/PCA/acp10000.rds")) # Load the acp10000 results (must check if this works)
envsup.acp <- suprow(acp, A[,3:23]) #Projeter individus supplémentaires
par(mfrow=c(1,1))
Marginality=apply(envsup.acp$lisup,1,function(x) sum(abs(x)*c(acp$eig[1:2]/sum(acp$eig)*100)))
Margin = ifelse(Marginality>quantile(Marginality,seuil,type=7),T,F) # Margins VS Core
Core = ifelse(Marginality<quantile(Marginality,seuilC,type=7),T,F)
res=data.frame(envsup.acp$lisup,Marginality,Margin,Core) #Put together individuals, and marginality
Final <- cbind(A[,1:23],res) #coordonées géo + scores + all variables + marginalité 
FinalMargin <- Final[Final$Margin=="TRUE",] #Select only margin individuals 
GrMarg = find.clusters.data.frame(FinalMargin[,-c(1:2,24:28)],
                                  n.clust=2,stat = c("BIC", "AIC", "WSS"),
                                  scale=TRUE,n.pca=50)
#Scale, Normed clustering. Force to 2 groups

dapc1<-dapc(FinalMargin[,-c(1:2,24:28)], GrMarg$grp,n.pca=21,n.da=1) #Keep one function and all axes
# Save as rds HERE 
png(file=paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/2Discrimine_",seuil,"_",seuilC,".png"),width=size,height = size*0.9,res=size/9)
scatter(dapc1, scree.da=F, scree.pca =F, bg="white", pch=20, cstar=0, col=c("blue","red"), solid=.8,
        cex=5,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:2),cell=2,cex.axis=2,cex.names=2)
dev.off()

