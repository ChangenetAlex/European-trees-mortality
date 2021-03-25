library(reshape2)
library(plyr)
library(ggplot2)
library(cowplot)
rm(list = ls())
gc()

Dir <- c("~/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
#### Table simple and interaction effects summary 
BINsimple <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/BIN.Table.Simple.Review.csv")
BINinter <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/BIN.Table.Inter.Review.csv")



BIN <- cbind(BINsimple,BINinter)
BIN<-BIN[,-c(20)]
BIN2 <- as.data.frame(t(BIN))
colnames(BIN2) <- BIN$species
BIN2[,"variable"] <- rownames(BIN2)
BIN2 <- BIN2[-1,]
BIN2 <- data.frame(lapply(BIN2, as.character), stringsAsFactors=FALSE)
rownames(BIN2) <- BIN2[,"variable"]


BIN2[,"MARGIN"] <- ifelse(grepl("Marginality",BIN2$variable,fixed = T),"M"," ") ## Here add a categorical for Marginality variables

### Rename variable 
BIN2$variable <- as.character(BIN2$variable)
BIN2$variable <- sub("Marginality..TE.","TE",BIN2$variable,fixed = T)
BIN2$variable <- sub("Marginality..LE.","LE",BIN2$variable,fixed = T)
BIN2$variable <- sub("BAIj.plot.1.mean","G.mean",BIN2$variable,fixed = T)
BIN2$variable <- sub("BAIj.plot.1","G.sum",BIN2$variable,fixed = T)
BIN2$variable <- sub(".X."," * ",BIN2$variable,fixed = T)
BIN2$variable <- sub(".2.","",BIN2$variable,fixed = T)
BIN2$variable <- sub("I.","Q.",BIN2$variable,fixed = T)
BIN2$variable <- sub("BAj","C.Inter",BIN2$variable,fixed = T)
BIN2$variable <- sub("BA.O","C.Intra",BIN2$variable,fixed = T)
BIN2$variable <- sub("BA","C.Tot",BIN2$variable,fixed = T)
BIN2$variable <- sub("dbh.plot.mean","DBH",BIN2$variable,fixed = T)
BIN2$variable <- sub("yearsbetweensurveys","CI",BIN2$variable,fixed = T)
BIN2$variable <- sub("treeNbr","D",BIN2$variable,fixed = T)
BIN2$variable <- sub("Temperature","T°",BIN2$variable,fixed = T)
BIN2$variable <- sub("Precipitation","P°",BIN2$variable,fixed = T)

BIN2 <- BIN2[grep("minSPEI",BIN2$variable,fixed=F,invert=T),]

BIN2[,"sumNEGBIN"] <- apply(BIN2, 1, function(x) sum(grepl("^-.*\\*$",x,fixed=F)==T))
BIN2[,"sumPOSBIN"] <- apply(BIN2, 1, function(x) sum(grepl("^[0-9].*\\*$",x,fixed=F)==T))
BIN2[,"sumSIGNIF"] <- BIN2[,"sumPOSBIN"]+BIN2[,"sumNEGBIN"]
BIN2[,"FREQsignif"] <- round(BIN2[,"sumSIGNIF"]/19,2)
BIN2[,"FREQsignifPOS"] <- round(BIN2[,"sumPOSBIN"]/19,2)
BIN2[,"FREQsignifNEG"] <- round(BIN2[,"sumNEGBIN"]/19,2)


library(plyr)
library(ggplot2)
library(cowplot)

### Here obtain barplot of frequency for CON models 
BIN2 <- melt(BIN2, id=c("variable","FREQsignif","MARGIN"),measure.vars = c("FREQsignifPOS", "FREQsignifNEG"),variable.name = "Signe")
BIN2_sorted <- arrange(BIN2, variable,desc(Signe)) 
BIN2_cumsum <- ddply(BIN2_sorted, "variable",
                     transform, label_ypos=cumsum(value)-0.5*value)
p1 <- ggplot(data=BIN2_cumsum, aes(x=reorder(variable,-value), y=value,fill=Signe)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1,
                                        linetype=1, color="black")) +
  theme(panel.grid.major = element_line(size=0, linetype=0,lineend="square", color="red"),
        panel.grid.minor = element_line(size=0, linetype=0,lineend="round", color="blue"))+
  geom_bar(stat="identity",color="black",size=0.5)+
  geom_text(aes(y=label_ypos, label=ifelse(value==0,NA,value*100)),color="black", vjust=-0.01, size=3,fontface='bold')+
  labs(y=paste0("Significiance frequency"), x="Variable",title=paste0('A) Main effects and interaction importance across species (BIN)'))+
  scale_colour_manual(values = c("blue","red"),name="Effect sign",labels=c("Positive","Negative")) +
  scale_fill_manual(values = c("FREQsignifPOS"="steelblue","FREQsignifNEG"="indianred"),name="Effect sign",labels=c("Positive","Negative"))+
  scale_y_continuous(labels = scales::percent,limits = c(0,1), expand = c(0,0.01)) +
  geom_text(aes(label = MARGIN,y=FREQsignif+0.03),colour="black",size=5,fontface = "bold")+
  theme(text = element_text(face="bold"),
        legend.direction ="vertical",
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=0.5),
        legend.box.margin = margin(6, 6, 6, 6),
        axis.text.x = element_text(size=11,color="black",angle = 45,hjust = 1,vjust = 1),axis.text.y = element_text(size=18,color="black"),
        axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
        legend.key = element_rect(fill = "white", colour = "black",size = 1),
        legend.key.heigh = unit(3,"line"),
        legend.key.width = unit(1.5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),
        #legend.justification = "center",
        legend.margin = margin(0,0,0,0),
        #legend.background=element_rect(fill="white",colour="black",size=0.5,linetype="solid"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.3),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=17,hjust = 0,vjust = 0),
        plot.margin = margin(0.2,0.2,0.2,0.2, "cm"),
        plot.caption = element_text(face="bold.italic"))
p1
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Freq.Barplot.BIN.pdf"),plot = p,base_width = 10, base_height = 7, dpi = 1000 ,units = "in",nrow=1,ncol=1)




### Idem with NB part of the model
NBsimple <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/ZT.Table.Simple.Review.csv")
NBinter <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/ZT.Table.Inter.Review.csv")

NB <- cbind(NBsimple,NBinter)
NB<-NB[,-c(20)]
NB2 <- as.data.frame(t(NB))
colnames(NB2) <- NB$species
NB2[,"variable"] <- rownames(NB2)
NB2 <- NB2[-1,]
NB2 <- data.frame(lapply(NB2, as.character), stringsAsFactors=FALSE)
rownames(NB2) <- NB2[,"variable"]

NB2[,"MARGIN"] <- ifelse(grepl("Marginality",NB2$variable,fixed = T),"M",NA) ## Here add a categorical for Marginality variables


### Rename variable 
NB2$variable <- as.character(NB2$variable)
NB2$variable <- sub("Marginality..TE.","TE",NB2$variable,fixed = T)
NB2$variable <- sub("Marginality..LE.","LE",NB2$variable,fixed = T)
NB2$variable <- sub("BAIj.plot.1.mean","G.mean",NB2$variable,fixed = T)
NB2$variable <- sub("BAIj.plot.1","G.sum",NB2$variable,fixed = T)
NB2$variable <- sub(".X."," * ",NB2$variable,fixed = T)
NB2$variable <- sub(".2.","",NB2$variable,fixed = T)
NB2$variable <- sub("I.","Q.",NB2$variable,fixed = T)
NB2$variable <- sub("BAj","C.Inter",NB2$variable,fixed = T)
NB2$variable <- sub("BA.O","C.Intra",NB2$variable,fixed = T)
NB2$variable <- sub("BA","C.Tot",NB2$variable,fixed = T)
NB2$variable <- sub("dbh.plot.mean","DBH",NB2$variable,fixed = T)
NB2$variable <- sub("yearsbetweensurveys","CI",NB2$variable,fixed = T)
NB2$variable <- sub("treeNbr","D",NB2$variable,fixed = T)
NB2$variable <- sub("Temperature","T°",NB2$variable,fixed = T)
NB2$variable <- sub("Precipitation","P°",NB2$variable,fixed = T)

NB2 <- BNB2[grep("minSPEI",BIN2$variable,fixed=F,invert=T),] ### Remove the min SPEI effect or not


NB2[,"sumNEGNB"] <- apply(NB2, 1, function(x) sum(grepl("^-.*\\*$",x,fixed=F)==T))
NB2[,"sumPOSNB"] <- apply(NB2, 1, function(x) sum(grepl("^[0-9].*\\*$",x,fixed=F)==T))
NB2[,"sumSIGNIF"] <- NB2[,"sumPOSNB"]+NB2[,"sumNEGNB"]
NB2[,"FREQsignif"] <- round(NB2[,"sumSIGNIF"]/19,2)
NB2[,"FREQsignifPOS"] <- round(NB2[,"sumPOSNB"]/19,2)
NB2[,"FREQsignifNEG"] <- round(NB2[,"sumNEGNB"]/19,2)



### Here obtain barplot of frequency for CON models 
NB2 <- melt(NB2, id=c("variable","FREQsignif","MARGIN"),measure.vars = c("FREQsignifPOS", "FREQsignifNEG"),variable.name = "Signe")
NB2_sorted <- arrange(NB2, variable, desc(Signe)) 
NB2_cumsum <- ddply(NB2_sorted, "variable",
                     transform, label_ypos=cumsum(value)-0.5*value)
p <- ggplot(data=NB2_cumsum, aes(x=reorder(variable,-value), y=value,fill=Signe)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1,
                                        linetype=1, color="black")) +
  theme(panel.grid.major = element_line(size=0, linetype=0,lineend="square", color="red"),
        panel.grid.minor = element_line(size=0, linetype=0,lineend="round", color="blue"))+
  geom_bar(stat="identity",color="black",size=0.5)+
  geom_text(aes(y=label_ypos, label=ifelse(value==0,NA,value*100)),color="black", vjust=-0.01, size=3,fontface='bold')+
  labs(y=paste0("Significiance frequency"), x="Variable",title=paste0('B) Main effects and interaction importance across species (NB)'))+
  scale_colour_manual(values = c("blue","red"),name="Effect sign",labels=c("Positive","Negative")) +
  scale_fill_manual(values = c("FREQsignifPOS"="steelblue","FREQsignifNEG"="indianred"),name="Effect sign",labels=c("Positive","Negative"))+
  scale_y_continuous(labels = scales::percent,limits = c(0,1), expand = c(0,0.01)) +
  geom_text(aes(label = MARGIN,y=FREQsignif+0.02),colour="black",size=5,fontface = "bold")+
  theme(text = element_text(face="bold"),
        legend.direction ="vertical",
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=0.5),
        legend.box.margin = margin(6, 6, 6, 6),
        axis.text.x = element_text(size=11,color="black",angle = 45,hjust = 1,vjust = 1),axis.text.y = element_text(size=18,color="black"),
        axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
        legend.key = element_rect(fill = "white", colour = "black",size = 1),
        legend.key.heigh = unit(3,"line"),
        legend.key.width = unit(1.5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),
        #legend.justification = "center",
        legend.margin = margin(0,0,0,0),
        #legend.background=element_rect(fill="white",colour="black",size=0.5,linetype="solid"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.3),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=17,hjust = 0,vjust = 0),
        plot.margin = margin(0.2,0.2,0.2,0.2, "cm"),
        plot.caption = element_text(face="bold.italic"))
p
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Freq.Barplot.NB.pdf"),plot = p,base_width = 10, base_height = 7, dpi = 1000 ,units = "in",nrow=1,ncol=1)



library(patchwork)
pall <- p1+theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  p+theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  plot_layout(ncol=1,nrow=2)

save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/FreqBoth.pdf"),
          plot = pall,base_width = 10, dpi = 1000, base_height = 7,units = "in",nrow=2,ncol=1)
# save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure3/Ag.Prediction.Mort.bin.Part1.pdf"),
#           plot = pall,base_width = 5, base_height = 5, dpi = 300 ,units = "in",nrow=3,ncol=3)





