# Alex on the 16/08/2019
# Adapted from Sophie Radcliffe and Juliette Archambau 
# Multiple functions to process the coeffficient of my models according to the latitude 

SavingInfo = "
####### function 4 ########
ggeffect allow to plot this parameters
x is my model
y is the kind of info that i want either 'ABS' for absolute values (compared to 0), either 'REL' for comparison between each parameters
effects is the parameters we xould like to plot. it can be either 'indiv' for simple effects, or 'sum' for synthetic effect such as copetition or climate

Example : ggEffect <- function(x,y='REL',effect='sum')
"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))

#### Load libraries #######
library(reshape2)         #
require(raster)           #
library(ggplot2)          #
library(mgcv)             # 
library(grid)             #
library(cowplot)          #
#library(Cairo)           #
###########################


ggEffect <- function(x,y="REL",effect="sum",band=T){
  if (y=="ABS"){clim.bio <- get(load(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/",Mod,"/clim_bio_abs_",Mod,".RData")))  # One model 
  }else if(y=="REL"){clim.bio <- get(load(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/",Mod,"/clim_bio_rel_",Mod,".RData")))}  # One model 
  ## Effect names 
  EffectCol <- grep(as.character(unique(clim.bio$variable)),pattern = "effect",fixed=T,value=T,invert=F) # Simple effect columns for which coef were extracted
  
  ## Individual effects vs sum of the effects 
  if (effect=="indiv"){clim.bio <- clim.bio[clim.bio$variable%in%EffectCol,] # If individual, keep these individual effects 
  }else if (effect=="sum"){
    clim.bio <- clim.bio[!(clim.bio$variable%in%EffectCol),] # If not, keep all the orther effects who have not "effect in their names"
    EffectCol <- as.character(unique(clim.bio$variable))}    # Exatrct these sum effects and put it in the initial vector
  clim.bio[,"variable"] <- as.character(clim.bio[,"variable"]) # Transfrom the factor as a charcater vector
  
  # Find my missing values 
  clim.bio <- clim.bio[order(clim.bio$latitude),] # Put the latitude in the right order (added on the 30th july)
  A <- unique(clim.bio[is.na(clim.bio$mean)&is.na(clim.bio$se),"latitude"]) # Need to be 0 both in the mean and in the SE
  if (length(A)!=0){
    missing <- data.frame(xmin=A-0.5,xmax=A+0.5, ymin=rep(Inf,length=length(A)), ymax=rep(-Inf,length=length(A))) # Identify the missing values automatically 
    # A-0.5
    i = 1
    while (is.na(missing$xmin[i+1])==F){
      if (missing$xmax[i+1]-missing$xmax[i]==0.5){
        while (missing$xmax[i+1]-missing$xmax[i]==0.5 & is.na(missing$xmax[i+1])==F){
          missing$xmax[i] <- missing$xmax[i+1]
          missing <- missing[-c(i+1),]
        }} else i=i+1
    }
    missing[missing$xmin<min(clim.bio$latitude),"xmin"] <- min(clim.bio$latitude) # If missing go under or above the max and min latitude
    missing[missing$xmax>max(clim.bio$latitude),"xmax"] <- max(clim.bio$latitude)
  }else missing <- NULL
  Mycol <- c("red","dark green", "dodgerblue3","orange","yellow","gray","black","green","lightblue","gray20","gray50") #### The colors that I will use 
  Myline <- c(1,3,6,1,3,6,1,3,6,1,3,6)
  
  # The plot theme 
  p<-ggplot(clim.bio) + 
    theme_bw()
  if (y=="ABS"){p <- p + ylab("Predicted absolute importance") + xlab("Latitude") +
    scale_y_continuous(limits=c(-0.2,max(clim.bio$mean)))
  }else if(y=='REL'){p <- p + ylab("Predicted relative importance") + xlab("Latitude") +
    scale_y_continuous(limits=c(-0.2,1))}
  p <- p + theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
  ) +
    theme(text = element_text(face="bold"),#
          legend.background=element_rect(fill="white",colour="black",size=0.2),#
          panel.border = element_rect(colour = "black", fill=NA, size=0.8),#
          axis.line = element_line(colour="black"),#
          plot.title = element_text(size=17,hjust = 0.5),#
          plot.caption = element_text(face="bold.italic"),#
          axis.text.y = element_text(size=13), #17  
          axis.text.x = element_text(size=13), #17  
          axis.title.x = element_text(size=13), #17
          axis.title.y = element_text(size=13), #17
          legend.position="bottom",
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.key.size = unit(5, 'lines'),
          legend.text = element_text(size=13,face="bold"), #17
          legend.key.height=unit(1.8,"line"), #2
          legend.key.width=unit(4.5,"line") #5
    )
  # The plot 
  par(mar=c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  p <- p + geom_line(aes(latitude, mean, colour = as.character(variable), linetype=as.character(variable)), size=1)+ #size = 0.8
    guides(fill=FALSE) +
    guides(col=guide_legend(ncol=ifelse(length(EffectCol)==4,2,3), byrow=F))
  if (effect=="indiv"){p <- p + scale_colour_manual(values=c(Mycol[1:length(EffectCol)]),labels=paste0("  ",EffectCol,"  ")) +
    scale_linetype_manual(values=c(Myline[1:length(EffectCol)]),labels=paste0("  ",EffectCol,"  "))
  } else if (effect=="sum"){p <- p + scale_colour_manual(values=c(Mycol[1:length(EffectCol)]),labels=paste0("  ",EffectCol,"  ")) +
    scale_linetype_manual(values=c(rep(Myline[1],times=length(EffectCol))),labels=paste0("  ",EffectCol,"  "))}
  if (band==T){p <- p + geom_ribbon(data=clim.bio,aes(latitude, mean, ymin=lwr, ymax=upr, colour=variable, fill=variable),alpha=0.05, linetype=2,size=0.5) # no size argument
  }else p <- p
  if (is.null(missing)==F){p <- p + geom_rect(data=missing, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),colour="light grey", fill="white", inherit.aes=FALSE)
  }else p <- p
  p <- p + 
    #geom_rect(xmin=min(clim.bio$latitude), xmax=max(clim.bio$latitude), ymin=-0.10, ymax=-0.05,colour="black", fill="black")+  # Rectangle noir
    geom_rect(xmin=0,xmax=42.5,ymin=-Inf,ymax=-0.1,colour="black", fill="red",alpha=0.2)+
    geom_rect(xmin=42.5,xmax=58,ymin=-Inf,ymax=-0.1,colour="black", fill="green",alpha=0.5)+
    geom_rect(xmin=58,xmax=Inf,ymin=-Inf,ymax=-0.1,colour="black", fill="blue",alpha=0.8)+
    geom_label(label = "Mediterranean", x=42.5, y=-0.01, size=3.5,vjust=2.2,label.padding = unit(0.15, "lines"), # size=3.5
               label.r = unit(0.15, "lines"),colour="white",fill="red",label.size = 0.2)+ 
    geom_label(label = "Temperate", x=50, y=-0.01, size=3.5,vjust=2.2,label.padding = unit(0.15, "lines"), # size=3.5
               label.r = unit(0.15, "lines"),colour="white",fill="darkgreen",label.size = 0.2)+
    geom_label(label = "Boreal", x=58, y=-0.01, size=3.5,vjust=2.2,label.padding = unit(0.15, "lines"),# size=3.5
               label.r = unit(0.15, "lines"),colour="white",fill="blue",label.size = 0.2)+
    geom_segment(x=min(clim.bio$latitude), y=-0.1, xend=min(clim.bio$latitude), yend=Inf, colour="black", size=0.1,linetype=11) + 
    geom_segment(x=max(clim.bio$latitude), y=-0.1, xend=max(clim.bio$latitude), yend=Inf, colour="black", size=0.1,linetype=11) +
    geom_segment(x=58, y=-0.1, xend=58, yend=Inf, colour="light grey", size=0.1,linetype=2) +
    geom_segment(x=42.5, y=-0.1, xend=42.5,yend=Inf, colour="light grey", size=0.1,linetype=2)
  print(p)
  assign(paste0("p",CODE),p,envir = .GlobalEnv)
  if (band==T){ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",
                                        "Results.all/Species/test",y,"_",effect,"_",CODE,"_",Mod,"_150band.png"),plot = p, width = 12.37, height = 7.04, dpi=150,units = "in") # Save
    #ggsave(filename = paste0(y,"_",effect,"_",CODE,"_",Mod,"_300band.png"),plot = p, width = 12.37, height = 7.04, dpi=300,units = "in")
    ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/test",y,"_",effect,"_",CODE,"_",Mod,"_600band.png"),plot = p, width = 12.37, height = 7.04, dpi=600,units = "in")
    ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/test",y,"_",effect,"_",CODE,"_",Mod,"_1000band.png"),plot = p, width = 12.37, height = 7.04, dpi=600,units = "in")
  }else ggsave(filename = paste0(y,"_",effect,"_",CODE,"_",Mod,"_NO.band.png"),plot = p, width = 12.37, height = 7.04, dpi=150,units = "in")
}



Allcode <- c("ABIALB","ACEPSE","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINSYL","POPNIG","QUEILE",
             "PINPINA","ALNGLU","PINNIG","PINPIN","POPTRE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allmod <- c("Mbin14A","Mbin13B","M2bin13B","M2bin13A","M2bin15B","Mbin3C","M2bin7B","Mbin11B","Mbin5B","Mbin3B","M2bin1C",
            "Mbin13A","Mbin13A","Mbin3B","Mbin15B","Mbin13A","M2bin1B","Mbin13B","Mbin7A","Mbin1B")
i= 1
for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  Mod <- Allmod[i]
  ggEffect(Mod,'REL',"sum",band=T)
}
legend <- get_legend(pABIALB)

###########################
### Figure 3 papier 1 #####
###########################

pall<-plot_grid(
  plot_grid(pQUEILE + theme(legend.position = "none"),pFAGSYL + theme(legend.position = "none"),
          pFRAEXC + theme(legend.position = "none"),pABIALB + theme(legend.position = "none"),
          labels = c('a) QUEILE', 'b) FAGSYL','c) FRAEXC', 'd) ABIALB'),align="hv", label_size = 14,hjust = -0.50,vjust =0.3),
          legend, nrow = 2, rel_heights = c(1, 0.15), scale = c(0.95,1)
  )
save_plot(filename = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure.coef/fig3_effect__800band.png"),plot = pall, base_width = 12.37, dpi=800,units = "in",nrow=2)


#

#

#

###########################
### Figure s8 papier 1 ####
###########################

pall<-plot_grid(
  plot_grid(pCASSAT + theme(legend.position = "none"),pPINPIN + theme(legend.position = "none"),
            pPINNIG + theme(legend.position = "none"),pPINHAL + theme(legend.position = "none"),
            labels = c('a) CASSAT', 'b) PINPIN','c) PINNIG', 'd) PINHAL'),align="hv", label_size = 14,hjust = -0.50,vjust =0.3),
  legend, nrow = 2, rel_heights = c(1, 0.15), scale = c(0.95,1)
)
save_plot(filename = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure.coef/figS8_effect__800band.png"),plot = pall, base_width = 12.37, dpi=800,units = "in",nrow=2)
#

#

#

###########################
### Figure s9 papier 1 ####
###########################

pall<-plot_grid(
  plot_grid(pQUEPET + theme(legend.position = "none"),pACEPSE + theme(legend.position = "none"),
            pPICABI + theme(legend.position = "none"),pPOPTRE + theme(legend.position = "none"),
            pPINSYL + theme(legend.position = "none"),pBETPEN + theme(legend.position = "none"),
            labels = c('a) QUEPET', 'b) ACEPSE','c) PICABI', 'd) POPTRE', 'd) PINSYL', 'd) BETPEN'),align="hv",ncol=2, label_size = 14,hjust = -0.50,vjust =0.3),
  legend, nrow = 2, rel_heights = c(1, 0.15), scale = c(0.95,1)
)

save_plot(filename = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure.coef/figS9_effect__800band.png"),plot = pall, base_width = 12.37, dpi = 800 ,units = "in",nrow=3)

#

#

#

###########################
### Figure s10 papier 1 ####
###########################

pall<-plot_grid(
  plot_grid(pALNGLU + theme(legend.position = "none"),pPINPINA + theme(legend.position = "none"),
            pPOPNIG + theme(legend.position = "none"),pQUEPYR+ theme(legend.position = "none"),
            pQUEROB + theme(legend.position = "none"),pQUESUB + theme(legend.position = "none"),
            labels = c('a) ALNGLU', 'b) PINPINA','c) POPNIG', 'd) QUEPYR', 'd) QUEROB', 'd) QUESUB'),align="hv",ncol=2, label_size = 14,hjust = -0.50,vjust =0.3),
  legend, nrow = 2, rel_heights = c(1, 0.15), scale = c(0.95,1)
)
save_plot(filename = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure.coef/figS10_effect__800band.png"),plot = pall, base_width = 12.37, dpi=800,units = "in",nrow=3)
