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

code <- list("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
seuil <- list(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
seuilC <- list(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
mod <- list("Mbin13A.18","Mbin14A.19","Mbin13B.26","Mbin13A.22","M2bin13B.22","M2bin13A21","M2bin15B23","Mbin3C.26","M2bin7B.17","Mbin11B.19","Mbin3B.20","Mbin15B.21","Mbin5B.17","Mbin3B.31","Mbin13A.27","M2bin1C.20","M2bin1B.23","Mbin13B.26","Mbin7A.26","Mbin1B.23")
Dir = c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/") # Directory 
i <- 13

for (i in c(1:20)){
  if (code[i]%in%c("PINPIN","PINNIG","QUEPET","QUEROB","PINPINA","POPTRE","ALNGLU","QUEPYR","QUESUB")){
    Mycol <- c("0"="yellow","2"="red","1"="blue","10"="gray")
    Mylab <- c("0"="Core","2"="Trailing edge","1"="Leading edge","10"="Transition area")
  }else{
    Mycol <- c("0"="yellow","1"="red","2"="blue","10"="gray")
    Mylab <- c("0"="Core","1"="Trailing edge","2"="Leading edge","10"="Transition area")
  } 
  test6 <- get(load(file=paste0(Dir,"species/",code[i],"/CLIMAP/CLIMAP2019/Climap_",code[i],"_",seuil[i],"_",seuilC[i],".RData")))
  x <- get(load(file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",code[i],"/CLIMAP/Models/binomial/",mod[i],"/",mod[i],".rda")))
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
  
  A <- rasterFromXYZ(test6[,c(1,2,29)])
  A2 <- aggregate(A, 5,fun=modal) # this takes some time to reduced 10 times de resol
  A3 <- as.data.frame(rasterToPoints(A2))
  # Here is the inset map in which we will zoom
  insetMap <- ggplot() + 
    labs(x = NULL, y = NULL)+
    geom_tile(data=A3[,c(1:2,3)],aes(x=x, y=y,fill=as.factor(groupes)))+#fill="gray80",color="gray80",interpolate=T) +
    #geom_point(data=A, aes(x=x, y = y,color=as.factor(groupes)),size=0.001)+
    #geom_tile(data=test6[,c(1:2,29)],aes(x=x, y=y,fill=as.factor(groupes)))+#fill="gray80",color="gray80",interpolate=T) +
    # xlim(-10,32)+
    # ylim(36,70)+
    scale_y_continuous(expand = c(0.02,0),limits = c(36,70)) +
    scale_x_continuous(expand = c(0.02,0),limits = c(-10,32)) +
    coord_fixed(1.3,expand=T) + 
    geom_polygon(data = europe1, aes(x=long, y = lat,group=group),fill=NA,col="black") + #gray 60 & gray 60
    #geom_point(data = r1, aes(x = X, y = Y,col=Z,fill=Z),shape=21,size=1)+
    #scale_color_manual(values = Mycol,name="Climatic areas:",labels=Mylab)+ # for points instead of tile
    scale_fill_manual(values = Mycol,name="Climatic areas:",labels=Mylab)+
    #geom_rect(data = data.frame(),aes(xmin=min(r1$X), xmax=max(r1$X), ymin=min(r1$Y), ymax=max(r1$Y)),colour = "red", fill = NA,size=1)+
    labs(title=paste0(code[i]), y=paste0("Latitude"), x="Longitude")+
    theme(text = element_text(face="bold"),
          legend.direction ="horizontal",legend.position = c(0.5, 0.5),
          legend.key = element_rect(fill = "white",colour = "black"),
          legend.key.heigh = unit(1,"line"),
          legend.key.width = unit(10,"line"),
          legend.text=element_text(size=16),
          legend.title=element_text(size=17),
          legend.justification = "center",
          legend.margin = margin(0,0,0,0),
          legend.background=element_rect(fill="white",colour="black",size=0,linetype="solid"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3),
          axis.text.x = element_text(size=16,color="black"),axis.text.y = element_text(size=16,color="black"),
          axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
          plot.title = element_text(vjust = -1.25,hjust = 0.5),
          # axis.text = element_blank(),
          # axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(0,0,0,0, "cm"),#margin(0,0,0,0,unit="pt"),
          panel.background = element_rect(fill="white", colour="black", size=0.4,linetype=1))+
    guides(fill = guide_legend(title.position = "top", 
                                # hjust = 0.5 centres the title horizontally
                                title.hjust = 0.5,
                                label.position = "bottom"))
  insetMap #transform to a grob object
  saveRDS(insetMap,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/CM_Figures/",code[i],".rds"))
  save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/CM_Figures/",code[i],".pdf"),plot = insetMap,base_width = 5, base_height = 5, dpi = 300 ,units = "in",nrow=1,ncol=1)
  # legend <- get_legend(insetMap)
  # plot(legend)
  # saveRDS(legend,paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/CM_Figures/Legend.Marginality.rds"))
}
  
code <- c("ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPINA","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")

A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/CM_Figures"),pattern = ".*rds$",full.names = T)
A <- A[c(1:7,9:21)] 
# Load only the first half 
mapply(function(x,y){
  assign(paste0("pall",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=code) ## be carreful with the order 

## Put the legend right of the figure to be then placed at the center
TEST <- pallABIALB+theme(legend.position = c(0.5, 0.5))
Legend2 <- get_legend(TEST)
plot(Legend2)
library(patchwork)
pall <- 
  pallABIALB+theme(legend.position = "none",
        axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        #axis.title.y=element_blank(),
        axis.ticks.x = element_blank(),
        #axis.ticks.y = element_blank(),
        plot.title = element_text(vjust = -8,hjust = 0.1),
        plot.margin = unit(c(0,0,0,0), "cm"))+
  
  pallACEPSE+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(vjust = -8,hjust = 0.1),
        plot.margin = unit(c(0,0,0,0), "cm"))+
  
  pallALNGLU+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.ticks.y = element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   plot.margin = unit(c(0,0,0,0), "cm"))+
  
  pallBETPEN+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.ticks.y = element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   plot.margin = unit(c(0,0,0,0), "cm"))+
  
  pallCASSAT+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   #axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   #axis.title.y=element_blank(),
                   axis.ticks.x = element_blank(),
                   #axis.ticks.y = element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   plot.margin = unit(c(0,0,0,0), "cm"))+
  
  pallFAGSYL+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.ticks.y = element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   plot.margin = unit(c(0,0,0,0), "cm"))+
  
  pallFRAEXC+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.ticks.y = element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   plot.margin = unit(c(0,0,0,0), "cm"))+
  
  pallPICABI+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.x = element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   axis.ticks.y = element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm"))+
  
  pallPINHAL+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   #axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   #axis.title.y=element_blank(),
                   axis.ticks.x = element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   #axis.ticks.y = element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm"))+
  
  pallPINNIG+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.x = element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   axis.ticks.y = element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm"))+
  
  pallPINPIN+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.x = element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   axis.ticks.y = element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm"))+
  
  pallPINPINA+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.title.x=element_blank(),
                    axis.title.y=element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   axis.ticks.x = element_blank(),
                   axis.ticks.y = element_blank(),
                    plot.margin = unit(c(0,0,0,0), "cm"))+
  
  pallPINSYL+theme(legend.position = "none",
        axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        plot.title = element_text(vjust = -8,hjust = 0.1),
        #axis.title.y=element_blank(),
        axis.ticks.x = element_blank(),
        #axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))+
  
  pallPOPNIG+theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title = element_text(vjust = -8,hjust = 0.1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))+
  
  pallPOPTRE+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   axis.text.y=element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.ticks.y = element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm"))+
                     
  pallQUEILE+theme(legend.position = "none",
                   axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.ticks.y = element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm"))+
    
  pallQUEPET+theme(legend.position = "none",
                   #axis.text.x=element_blank(),
                   #axis.text.y=element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   #axis.title.x=element_blank(),
                   #axis.title.y=element_blank(),
                   #axis.ticks.x = element_blank(),
                   #axis.ticks.y = element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm"))+
    
  pallQUEPYR+theme(legend.position = "none",
                   #axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   #axis.title.x=element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   axis.title.y=element_blank(),
                   #axis.ticks.x = element_blank(),
                   axis.ticks.y = element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm"))+
    
  pallQUEROB+theme(legend.position = "none",
                   #axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   #axis.title.x=element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   axis.title.y=element_blank(),
                   #axis.ticks.x = element_blank(),
                   axis.ticks.y = element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm"))+
    
  pallQUESUB+theme(legend.position = "none",
                   #axis.text.x=element_blank(),
                   axis.text.y=element_blank(),
                   #axis.title.x=element_blank(),
                   plot.title = element_text(vjust = -8,hjust = 0.1),
                   axis.title.y=element_blank(),
                   #axis.ticks.x = element_blank(),
                   axis.ticks.y = element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm"))+
    plot_layout(ncol=4,nrow=5)
  

pall2 <- pall/Legend2+plot_layout(heights = c(1,0.05))

save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/CM_Figures/All.legend.pdf"),
          plot = pall2,base_width = 5, base_height = 5, dpi = 600 ,units = "in",nrow=5,ncol=4)

