rm(list = ls())
gc()
require(rgdal)
require(adegenet)
require(ade4)
require(parallel)
require(fields)
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
library(knitr)
library(spaMM)
library(fields)
library(piecewiseSEM)
library(raster)
library(rasterVis)
library(cowplot)
library(sf)
library(rnaturalearth)
library(grid)
library(svglite)
library("ggplotify")
library(parallel)

code <- list("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
seuil <- list(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.8,0.7,0.7,0.8,0.7)
seuilC <- list(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
mod <- list("MnbZT13A.19","MnbZT14A.20","MnbZT13B.23","MnbZT13A.21","M2nbZT13B.29","M2nbZT13A.22","M2nbZT15B.24","MnbZT3C.24","M2nbZT7B.26","MnbZT11B.22","MnbZT3B.27","MnbZT15B.28","MnbZT5B.21","MnbZT13A.33","M2nbZT1C.27","M2nbZT1B.24","MnbZT13B.27","MnbZT7A.23","MnbZT1B.27")
i <- 1

for (i in c(1:19)){
  if (code[i]%in%c("PINPIN","PINNIG","QUEPET","QUEROB","PINPINA","POPTRE","ALNGLU","QUEPYR","QUESUB")){
    Mycol <- c("0"="yellow","2"="red","1"="blue","10"="gray")
    Mylab <- c("0"="Core","2"="Trailing edge","1"="Leading edge","10"="Transition area")
  }else{
    Mycol <- c("0"="yellow","1"="red","2"="blue","10"="gray")
    Mylab <- c("0"="Core","1"="Trailing edge","2"="Leading edge","10"="Transition area")
  } 
  
  Dir = c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/") # Directory 
  test6 <- get(load(file=paste0(Dir,"species/",code[i],"/CLIMAP/CLIMAP2019/Climap_",code[i],"_",seuil[i],"_",seuilC[i],".RData")))
  Dir =c(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",code[i],"/CLIMAP/Models/Negbin/",mod[i],"/"))
  setwd(Dir)
  x <- get(load(file = paste0(mod[i],".rda")))
  r1 <- data.frame(x$data$longitude, x$data$latitude, predict(x))
  colnames(r1) <- c("X","Y","Z")
  r1$bin <- ifelse(r1$Z>quantile(r1$Z,0.75,type=7),"0",ifelse(r1$Z<quantile(r1$Z,0.25,type=7),"2","1")) # mortality as classes 
  
  par(mfrow=c(1,1))
  worldmap <- getMap(resolution = "high")
  #europe <- worldmap[which(worldmap$REGION=="Europe"),]
  #europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & as.character(worldmap$NAME) != "Russia"),] 
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3)|grepl("North Africa", worldmap$GEO3)),] 
  europe <-spTransform(europe,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ") # Convert to the right system
  europe <- gSimplify(europe, tol = 0.00001)
  #europe = crop(europe,c(-15, 45, 35, 70))
  europe1 = crop(europe,c(-10, 32, 36, 70))
  A <- rasterFromXYZ(test6[,c(1,2,29)])
  A2 <- aggregate(A, 10,fun=modal) # this takes some time to reduced 10 times de resol
  A3 <- as.data.frame(rasterToPoints(A2))
  
  pobs <- ggplot() + theme(panel.background = element_rect(fill="white", colour="black", size=3, 
                                                           linetype=1, color="black"),legend.key.size = unit(3, "cm"),legend.position=c(0.5,0.5)) +
    theme(panel.grid.major = element_line(colour="blue", size=4, linetype=0,lineend="square", color="red"),
          panel.grid.minor = element_line(colour="blue", size=1, linetype=0,lineend="round", color="blue"))+
    theme(legend.background = element_rect(fill="white",size=1, linetype="solid", 
                                           colour ="black"))+
    #theme(legend.position = "none")+ ## Allow us to disable the full legend
    #geom_tile(data=A3[,c(1:2)],aes(x=x, y=y),fill="gray90",na.rm = T) +
    geom_point(data=A3[,c(1:2)], aes(x=x, y = y),col="gray90",size=0.001)+ #instead of Tile 
    #geom_point(data=A, aes(x=x, y = y),col="gray90",size=0.001)+ #instead of Tile 
    scale_y_continuous(expand = c(0.02,0),limits = c(36,70)) +
    scale_x_continuous(expand = c(0.02,0),limits = c(-10,32)) +
    coord_fixed(1.3,expand=T) +  
    geom_polygon(data = europe1, aes(x=long, y = lat,group=group),fill=NA,col="black") + #gray 60 & gray 60
    geom_point(data = r1, aes(x = X, y = Y, group=bin, fill = bin, shape = bin, size=bin),alpha=0.55,color="black",stroke=0.1)+
    #scale_colour_manual(values = c("red2","springgreen","blue3"),name="Mortality occurence \n probability",labels = c("High (>Q3)", "Medium \n(Q1 - Q3)","Low (<Q1)")) +
    scale_fill_manual(values = c("darkred","chocolate1","olivedrab1"),name="Mortality amount \n (proportion/plot/year)",labels = c("High (>Q3)", "Medium \n(Q1 - Q3)","Low (<Q1)")) +
    scale_shape_manual(values = c(21,21,21),name="Mortality amount \n (proportion/plot/year)",labels = c("High (>Q3)", "Medium \n(Q1 - Q3)","Low (<Q1)")) +
    scale_size_manual(values= c(2.5,2.5,2.5),name="Mortality amount \n (proportion/plot/year)",labels = c("High (>Q3)", "Medium \n(Q1 - Q3)","Low (<Q1)"))+ # 0.7 0.5 0.7 de base
    guides(shape = guide_legend(override.aes = list(size = 15)))+
    #xlim(-10,32)+
    labs(title=paste0(code[i]," - NB"), y=paste0("Latitude"), x="Longitude")+
    theme(text = element_text(face="bold"),legend.direction ="vertical",legend.position = "top",
          axis.text.x = element_text(size=18,color="black"),axis.text.y = element_text(size=18,color="black"),
          axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
          legend.key = element_rect(fill = "white", colour = "black",size = 1),
          legend.background=element_rect(fill="white",colour="black",size=0),
          legend.key.heigh = unit(4,"line"),
          legend.key.width = unit(0.8,"line"),
          legend.text=element_text(size=16),
          legend.title=element_text(size=17),
          panel.border = element_rect(colour = "black", fill=NA, size=0),
          axis.line = element_line(colour="black"),
          plot.margin = margin(0,0,0,0, "cm"),
          plot.title = element_text(size=12,hjust = 0.02,vjust = -1),
          plot.caption = element_text(face="bold.italic"))
  #pobs
  #legend <- get_legend(pobs)
  #saveRDS(legend,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Legend.Quantiles.NB.rds"))
  ## plot et inset 
  save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Ag.pNB.",code[i],".png"),plot = pobs,base_width = 5, base_height = 5, dpi = 300 ,units = "in",nrow=1,ncol=1)
  saveRDS(pobs,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Ag.pNB.",code[i],".rds"))
  # assign(paste0("pall",Allcode[i]),pallV2,envir = .GlobalEnv)
  rm("x")
}



A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3"),pattern = "Ag.pNB.*rds$",full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("pall",y),readRDS(x),envir = .GlobalEnv)
},x=A[1:12],y=code[c(2:12,1)]) ## be carreful with the order 


library(patchwork)
pall <- pallABIALB+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        #axis.title.y=element_blank(),
        plot.margin = unit(c(0,0,-1,0), "cm"))+
  pallACEPSE+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(0,0,-1,0), "cm"))+
  pallALNGLU+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,-1,0), "cm"))+
  pallBETPEN+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   #axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   #axis.title.y=element_blank(),
                   plot.margin = unit(c(-1,0,-1,0), "cm"))+
  pallCASSAT+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(-1,0,-1,0), "cm"))+
  pallFAGSYL+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(-1,0,-1,0), "cm"))+
  pallFRAEXC+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   #axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   #axis.title.y=element_blank(),
                   plot.margin = unit(c(-1,0,-1,0), "cm"))+
  pallPICABI+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(-1,0,-1,0), "cm"))+
  pallPINHAL+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(-1,0,-1,0), "cm"))+
  pallPINNIG+theme(legend.position = "none",
                   #axis.text.x=element_blank(),
                   #axis.text.y=element_blank(),
                   #axis.title.x=element_blank(),
                   #axis.title.y=element_blank(),
                   plot.margin = unit(c(-1,0,0,0), "cm"))+
  pallPINPIN+theme(legend.position = "none",
                   #axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   #axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(-1,0,0,0), "cm"))+
  pallPINPINA+theme(legend.position = "none",
                    #axis.text.x=element_blank(),
                    axis.text.y=element_blank(),
                    #axis.title.x=element_blank(),
                    axis.title.y=element_blank(),
                    plot.margin = unit(c(-1,0,0,0), "cm"))+
  plot_layout(ncol=3,nrow=4)


save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Ag.Prediction.Mort.NB.Part1.png"),
          plot = pall,base_width = 5, base_height = 5, dpi = 300 ,units = "in",nrow=4,ncol=3)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Ag.Prediction.Mort.NB.Part1.pdf"),
          plot = pall,base_width = 5, base_height = 5, dpi = 300 ,units = "in",nrow=3,ncol=3)

### make some space 
rm(list = ls())
gc()

### Seconde moitié 
Allcode <- list("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")

A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3"),pattern = "Ag.pNB.*rds$",full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("pall",y),readRDS(x),envir = .GlobalEnv)
},x=A[13:19],y=Allcode[13:19])

Marginality <- readRDS(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Legend.Marginality.rds"))
Quantile <- readRDS(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Legend.Quantiles.NB.rds"))
Marginality<-plot_grid(Marginality,align="hv")
Quantile<-plot_grid(Quantile,align="hv")


plegend <-plot_spacer()+
  Quantile+theme(plot.margin = unit(c(-2,1.1,-2,0), "cm"))+
  plot_spacer()+plot_spacer()+
  plot_layout(ncol=4)

pall <- pallPINSYL+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        #axis.title.y=element_blank(),
        plot.margin = unit(c(0,0,-2,0), "cm"))+
  pallPOPTRE+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(0,0,-2,0), "cm"))+
  pallQUEILE+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,-2,0), "cm"))+
  pallQUEPET+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   #axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   #axis.title.y=element_blank(),
                   plot.margin = unit(c(-2,0,-2,0), "cm"))+
  pallQUEPYR+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(-2,0,-2,0), "cm"))+
  pallQUEROB+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(-2,0,-2,0), "cm"))+
  pallQUESUB+theme(legend.position = "none",
                   #axis.text.x=element_blank(),
                   #axis.text.y=element_blank(),
                   #axis.title.x=element_blank(),
                   #axis.title.y=element_blank(),
                   plot.margin = unit(c(-2,0,-2,0), "cm"))+
  Quantile+theme(plot.margin = unit(c(-2,0,-2,0), "cm"))+
  plot_spacer()+
  plot_layout(ncol=3,nrow=3)
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Ag.Prediction.Mort.NB.Part2bis.png"),
          plot = pall,base_width = 5, base_height = 5, dpi = 300 ,units = "in",nrow=3,ncol=3)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Ag.Prediction.Mort.NB.Part2bis.pdf"),
          plot = pall,base_width = 5, base_height = 5, dpi = 300 ,units = "in",nrow=3,ncol=3)

