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

QUANT.all <- data.frame(matrix(nrow=20,ncol=4,NA))
colnames(QUANT.all) <- c("Q1 Occurence","Q3 Occurence","Q1 Amount","Q3 Amount")
i = 1
Allcode <- list("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
                "QUEPYR","QUEROB","QUESUB")
Allseuil <- list(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- list(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
rownames(QUANT.all) <-Allcode
Allmod <- list("Mbin13A.18","Mbin14A.19","Mbin13B.26","Mbin13A.22","M2bin13B.22","M2bin13A21","M2bin15B23","Mbin3C.26","M2bin7B.17","Mbin11B.19","Mbin3B.20","Mbin15B.21","Mbin5B.17","Mbin3B.31","Mbin13A.27","M2bin1C.20","M2bin1B.23","Mbin13B.26","Mbin7A.26","Mbin1B.23")

mcmapply(function(code,seuil,seuilC,mod){
  if (code%in%c("PINPIN","PINNIG","QUEPET","QUEROB","PINPINA","POPTRE","ALNGLU","QUEPYR","QUESUB")){
    Mycol <- c("0"="yellow","2"="red","1"="blue","10"="gray")
    Mylab <- c("0"="Core","2"="Trailing edge","1"="Leading edge","10"="Transition area")
  }else{
    Mycol <- c("0"="yellow","1"="red","2"="blue","10"="gray")
    Mylab <- c("0"="Core","1"="Trailing edge","2"="Leading edge","10"="Transition area")
  } 
  
  Dir = c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/") # Directory 
  test6 <- get(load(file=paste0(Dir,"species/",code,"/CLIMAP/CLIMAP2019/Climap_",code,"_",seuil,"_",seuilC,".RData")))
  Dir =c(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",code,"/CLIMAP/Models/binomial/",mod,"/"))
  setwd(Dir)
  x <- get(load(file = paste0(mod,".rda")))
  r1 <- data.frame(x$data$longitude, x$data$latitude, predict(x))
  colnames(r1) <- c("X","Y","Z")
  r1$bin <- ifelse(r1$Z>quantile(r1$Z,0.75,type=7),"0",ifelse(r1$Z<quantile(r1$Z,0.25,type=7),"2","1")) # mortality as classes 
  
  ## Inset map
  par(mfrow=c(1,1))
  worldmap <- getMap(resolution = "high")
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3)|grepl("North Africa", worldmap$GEO3)),] 
  europe <-spTransform(europe,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ") # Convert to the right system
  europe <- gSimplify(europe, tol = 0.00001)
  europe1 = crop(europe,c(-10, 32, 36, 70))
  europe = crop(europe,c(min(r1$X), max(r1$X), min(r1$Y), max(r1$Y)))
  plot(europe1)
  # Here is the inset map in which we will zoom
  insetMap <- ggplot() + 
    labs(x = NULL, y = NULL)+
    geom_raster(data=test6[,c(1:2,29)],aes(x=x, y=y,fill=as.factor(groupes)))+#fill="gray80",color="gray80",interpolate=T) +
    # xlim(-10,32)+
    # ylim(36,70)+
    scale_y_continuous(expand = c(0.02,0),limits = c(36,70)) +
    scale_x_continuous(expand = c(0.02,0),limits = c(-10,32)) +
    coord_fixed(1.3,expand=T) + 
    geom_polygon(data = europe1, aes(x=long, y = lat,group=group),fill=NA,col="black") + #gray 60 & gray 60
    #geom_point(data = r1, aes(x = X, y = Y,col=Z,fill=Z),shape=21,size=1)+
    scale_fill_manual(values = Mycol,name="Climatic areas:",labels=Mylab)+
    geom_rect(data = data.frame(),aes(xmin=min(r1$X), xmax=max(r1$X), ymin=min(r1$Y), ymax=max(r1$Y)),colour = "red", fill = NA,size=1)+
    theme(text = element_text(face="bold"),
          legend.direction ="vertical",legend.position = "none",
          legend.key = element_rect(fill = "white",colour = "black"),
          legend.key.heigh = unit(3,"line"),
          legend.key.width = unit(2.5,"line"),
          legend.text=element_text(size=16),
          legend.title=element_text(size=17),
          legend.justification = "center",
          legend.margin = margin(0,0,0,0),
          legend.background=element_rect(fill="white",colour="black",size=0,linetype="solid"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3),
          # axis.text.x = element_text(size=18,color="black"),axis.text.y = element_text(size=18,color="black"),
          # axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(0,0,0,0, "cm"),#margin(0,0,0,0,unit="pt"),
          panel.background = element_rect(fill="white", colour="black", size=0.4,linetype=1))
  myinset <- ggplotGrob(insetMap) #transform to a grob object 
  # legend <- get_legend(insetMap)
  # saveRDS(legend,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Legend.Marginality.rds"))
  Inset<-plot_grid(insetMap,align="hv") #Need to trasnform it to a grid (for the inset)
  myinset <- ggplotGrob(Inset) #transform to a grob 
  #plot(myinset)
  
  
  par(mfrow=c(1,1))
  worldmap <- getMap(resolution = "high")
  #europe <- worldmap[which(worldmap$REGION=="Europe"),]
  #europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & as.character(worldmap$NAME) != "Russia"),] 
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3)|grepl("North Africa", worldmap$GEO3)),] 
  europe <-spTransform(europe,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ") # Convert to the right system
  europe = crop(europe,c(-15, 45, 35, 70))
  pobs <- ggplot() + theme(panel.background = element_rect(fill="white", colour="black", size=3, 
                                                           linetype=1, color="black"),legend.key.size = unit(3, "cm"),legend.position=c(0.5,0.5)) +
    theme(panel.grid.major = element_line(colour="blue", size=4, linetype=0,lineend="square", color="red"),
          panel.grid.minor = element_line(colour="blue", size=1, linetype=0,lineend="round", color="blue"))+
    theme(legend.background = element_rect(fill="white",size=1, linetype="solid", 
                                           colour ="black"))+
    #theme(legend.position = "none")+ ## Allow us to disable the full legend
    #geom_tile(data=test6[,c(1:2)],aes(x=x, y=y),fill="gray80",na.rm = T) +
    geom_tile(data=test6[,c(1:2)],aes(x=x, y=y),fill="gray90",na.rm = T) +
    #geom_polygon(data=test6, aes(x=x, y = y,fill=Marginality),col="gray80")+
    scale_y_continuous(expand = c(0.02,0),limits = c(36,70)) +
    scale_x_continuous(expand = c(0.02,0),limits = c(-10,32)) +
    coord_fixed(1.3,expand=T) +  
    geom_polygon(data = europe1, aes(x=long, y = lat,group=group),fill=NA,col="black") + #gray 60 & gray 60
    geom_point(data = r1, aes(x = X, y = Y, group=bin, fill = bin, shape = bin, size=bin),alpha=0.55,color="black",stroke=0.1)+
    #scale_colour_manual(values = c("red2","springgreen","blue3"),name="Mortality occurence \n probability",labels = c("High (>Q3)", "Medium \n(Q1 - Q3)","Low (<Q1)")) +
    scale_fill_manual(values = c("darkred","chocolate1","olivedrab1"),name="Mortality occurence \n probability",labels = c("High (>Q3)", "Medium \n(Q1 - Q3)","Low (<Q1)")) +
    scale_shape_manual(values = c(21,21,21),name="Mortality occurence \n probability",labels = c("High (>Q3)", "Medium \n(Q1 - Q3)","Low (<Q1)")) +
    scale_size_manual(values= c(2.5,2.5,2.5),name="Mortality occurence \n probability",labels = c("High (>Q3)", "Medium \n(Q1 - Q3)","Low (<Q1)"))+ # 0.7 0.5 0.7 de base
    guides(shape = guide_legend(override.aes = list(size = 15)))+
    #xlim(-10,32)+
    labs(title=paste0(code," - BIN"), y=paste0("Latitude"), x="Longitude")+
    theme(text = element_text(face="bold"),legend.direction ="vertical",legend.position = "none",
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
  #saveRDS(legend,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Legend.Quantiles.rds"))
  
  ## recup coord
  xrange <- unlist(ggplot_build(pobs)$layout$panel_params[[1]][1])
  yrange <- unlist(ggplot_build(pobs)$layout$panel_params[[1]][8])
  xMin = min(xrange)
  xMax = max(xrange)
  xdelta = xMax-xMin
  yMin = min(yrange)
  yMax = max(yrange)
  ydelta = yMax-yMin 
  
  ## plot et inset 
  pallV2 <- pobs + annotation_custom(myinset, xmin=xMin+(xdelta*0.43)*0.03, xmax=xMin+(xdelta*0.43)+(xdelta*0.43)*0.03, 
                                     ymin = yMax-(ydelta*0.43)*0.02-ydelta*0.43, ymax=yMax-(ydelta*0.43)*0.02)  
  
  save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Inset",code,".png"),plot = pallV2,base_width = 5, base_height = 5, dpi = 300 ,units = "in",nrow=1,ncol=1)
  saveRDS(pallV2,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Inset",code,".rds"))
  # assign(paste0("pall",Allcode[i]),pallV2,envir = .GlobalEnv)
  rm("x")
},code=Allcode,seuil=Allseuil,seuilC=AllseuilC,mod=Allmod,mc.cores = 20)


A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3"),pattern = "Inset.*rds$",full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("pall",y),readRDS(x),envir = .GlobalEnv)
},x=A[1:12],y=Allcode[c(2:11,1,12)]) ## be carreful with the order 


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

save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Prediction.Mort.bin.Part1.png"),
          plot = pall,base_width = 5, base_height = 5, dpi = 300 ,units = "in",nrow=4,ncol=3)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Prediction.Mort.bin.Part1.pdf"),
          plot = pall,base_width = 5, base_height = 5, dpi = 300 ,units = "in",nrow=3,ncol=3)

### make some space 
rm(list = ls())
gc()

### Seconde moitié 


Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")

A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3"),pattern = "Inset.*rds$",full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("pall",y),readRDS(x),envir = .GlobalEnv)
},x=A[13:20],y=Allcode[13:20])

Marginality <- readRDS(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Legend.Marginality.rds"))
Quantile <- readRDS(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Legend.Quantiles.rds"))
Marginality<-plot_grid(Marginality,align="hv")
Quantile<-plot_grid(Quantile,align="hv")


plegend <-plot_spacer()+
  Quantile+theme(plot.margin = unit(c(-2,1.1,-2,0), "cm"))+
  plot_spacer()+Marginality+theme(plot.margin = unit(c(-2,0,-2,1.1), "cm"))+
  plot_layout(ncol=4)

pall <- pallPINSYL+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        #axis.title.y=element_blank(),
        plot.margin = unit(c(0,0,-2,0), "cm"))+
  pallPOPNIG+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(0,0,-2,0), "cm"))+
  pallPOPTRE+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(0,0,-2,0), "cm"))+
  pallQUEILE+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   #axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   #axis.title.y=element_blank(),
                   plot.margin = unit(c(-2,0,-2,0), "cm"))+
  pallQUEPET+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(-2,0,-2,0), "cm"))+
  pallQUEPYR+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(-2,0,-2,0), "cm"))+
  pallQUEROB+theme(legend.position = "none",
                   #axis.text.x=element_blank(),
                   #axis.text.y=element_blank(),
                   #axis.title.x=element_blank(),
                   #axis.title.y=element_blank(),
                   plot.margin = unit(c(-2,0,-2,0), "cm"))+
  pallQUESUB+theme(legend.position = "none",
                   #axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   #axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   plot.margin = unit(c(-2,0,-2,0), "cm"))+
  plegend+theme(plot.margin = unit(c(-0.5,-0.5,-1,0), "cm"))+
  plot_layout(ncol=3,nrow=3)
#pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Prediction.Mort.bin.Part2.png"),
          plot = pall,base_width = 5, base_height = 5, dpi = 300 ,units = "in",nrow=3,ncol=3)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Prediction.Mort.bin.Part2.pdf"),
          plot = pall,base_width = 5, base_height = 5, dpi = 300 ,units = "in",nrow=3,ncol=3)

