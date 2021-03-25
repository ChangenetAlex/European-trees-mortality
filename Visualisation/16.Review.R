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


## Correlations pred et obs VS RE LE and CORE

#Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
#            "QUEPYR","QUEROB","QUESUB")
#Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
#AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
#Allmod <- c("Mbin13A.18","Mbin14A.19","Mbin13B.26","Mbin13A.22","M2bin13B.22","M2bin13A21","M2bin15B23","Mbin3C.26","M2bin7B.17","Mbin11B.19","Mbin3B.20","Mbin15B.21","Mbin5B.17","Mbin3B.31","Mbin13A.27","M2bin1C.20","M2bin1B.23","Mbin13B.26","Mbin7A.26","Mbin1B.23")



i = 1
Allcode <- c("ABIALB","ACEPSE","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINSYL","POPNIG","QUEILE")
Allseuil <- c(0.7,0.55,0.8,0.7,0.7,0.7,0.8,0.7,0.8,0.7,0.8)
AllseuilC <- c(0.6,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
Allmod <- c("Mbin14A.19","Mbin13B.26","M2bin13B.22","M2bin13A21","M2bin15B23","Mbin3C.26","M2bin7B.17","Mbin11B.19","Mbin5B.17","Mbin3B.31","M2bin1C.20")

for (i in 1:length(Allcode)){
  Dir = c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/") # Directory 
  test6 <- get(load(file=paste0(Dir,"species/",Allcode[i],"/CLIMAP/CLIMAP2019/Climap_",Allcode[i],"_",Allseuil[i],"_",AllseuilC[i],".RData")))
  Dir =c(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/binomial/",Allmod[i],"/"))
  setwd(Dir)
  #assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
  x <- get(load(file = paste0(Allmod[i],".rda")))
  #rm("x")
  # New section 
  x$data$Plotcat
  r1 <- data.frame(x$data$longitude, x$data$latitude, predict(x), x$data$Plotcat)
  colnames(r1) <- c("X","Y","Z","Cat")
  r1$Cat <- as.character(r1$Cat)
  str(r1)
  r1[r1$Cat=="0","Cat"] <- "Core"
  r1[r1$Cat=="2","Cat"] <- "LE" # To change when reverse. 2=LE in normal. 2=RE in reverse
  r1[r1$Cat=="1","Cat"] <- "TE"
  r1$bin <- ifelse(r1$Z>quantile(r1$Z,0.75,type=7),"> Q3",ifelse(r1$Z<quantile(r1$Z,0.25,type=7),"< Q1","Q1 - Q3")) ## Percentage of mortality as classes 
  print(summary(as.factor(r1$Cat)))
  print(Allcode[i])
  r1$bin <- as.factor(r1$bin)
  r1$bin = factor(r1$bin,levels(r1$bin)[c(1,3,2)])
  A <- data.table::melt(prop.table(table(r1$bin,r1$Cat), 2))
  B <- 100*prop.table(table(r1$Cat,r1$bin),1)
  p <- c(0.25,0.50,0.25)
  B1 <- unlist(chisq.test(B[1,],p=p,correct = F)[c(1,3)])
  B2 <- unlist(chisq.test(B[2,],p=p,correct = F)[c(1,3)])
  B3 <- unlist(chisq.test(B[3,],p=p,correct = F)[c(1,3)])
  B <- rbind(B1,B2,B3)
  B <- round(B,3)
  A[,"x2"] <- rep(B[,1],each=3)
  A[,"pvalue"] <- rep(B[,2],each=3)
  A$value <- round(A$value,3)*100
  A[c(2,5,8),"value"] <- A[c(2,5,8),"value"]/2 
  colnames(A) <- c("Quantiles","Zone","value","x2","pvalue")
  print(A)
  A$Zone = factor(A$Zone,levels(A$Zone)[c(3,1,2)])
  A[,"signif"] <- ifelse(A$pvalue<0.05,T,F)
  A[,"Psignif"] <- ""
  A[A$Quantiles=="Q1 - Q3","Psignif"] <- ifelse(A[A$Quantiles=="Q1 - Q3","pvalue"]<0.05,"*","")
  p <- ggplot(data=A,aes(x=Zone,y=value,fill=Quantiles,linetype=signif,alpha=signif,size=signif)) + theme(panel.background = element_rect(fill="white", colour="black", size=3, 
                                                                                                         linetype=1, color="black"),legend.key.size = unit(2, "cm"),legend.position = c(0.5,0.5)) + # 0.955,0.85 verticale. 
    theme(panel.grid.major = element_line(colour="blue", size=4, linetype=0,lineend="square", color="red"),
          panel.grid.minor = element_line(colour="blue", size=1, linetype=0,lineend="round", color="blue"))+
    theme(legend.background = element_rect(fill="white",size=1, linetype="solid", 
                                           colour ="black"))+
    theme(legend.position = "bottom",legend.box = "vertical")+ ## Allow us to disable the full legend
    #ylim(0,max(A$value)+max(A$value)/100*10)+
    geom_bar(stat="identity",position="dodge",width = 0.65,colour="black")+
    geom_text(aes(label = Psignif,y=Inf),colour="black",size=30,fontface = "bold")+
    geom_text(aes(label=value),vjust=-0.2,colour="black",position=position_dodge(0.7),size=8,fontface = "bold")+
    #annotate("text", aes(label=(Psignif),y=max(A$value)+max(A$value)/100*10) +
    coord_cartesian(ylim = c(0,max(A$value)+max(A$value)/100*10), clip = "off")+
    annotate("label", x = 1, y = 0, label = paste0("X-squared = ",A$x2[1],"; p-value = ",A$pvalue[1]),size=4.7,fontface = "bold")+
    annotate("label", x = 2, y = 0, label = paste0("X-squared = ",A$x2[4],"; p-value = ",A$pvalue[4]),size=4.7,fontface = "bold")+
    annotate("label", x = 3, y = 0, label = paste0("X-squared = ",A$x2[7],"; p-value = ",A$pvalue[7]),size=4.7,fontface = "bold")+
    scale_fill_manual(values = c("< Q1"="olivedrab1","Q1 - Q3"="chocolate1","> Q3"="darkred"),name="Quantile",
                      breaks=waiver(),
                      labels=c("Low (<Q1)", "Medium (>Q1 and <Q3)","High (>Q3)"))+
    scale_linetype_manual(name="Significiance",values = c("FALSE"="blank","TRUE"="solid"),
                          guide=F,
                          labels = c("Not significant", "Significant"))+
    scale_alpha_manual(name="Significiance",values = c("FALSE"=0.55,"TRUE"=1),
                       guide=F,
                       labels = c("Not significant", "Significant"))+
    scale_size_manual(name="Significiance",values = c("FALSE"=1,"TRUE"=1.25),
                      guide=F,
                      labels = c("Not significant", "Significant"))+
    labs(title=paste0(Allcode[i],"\nBIN"), y=paste0("Observed quantiles repartition in %"), x="Zone")+
    theme(text = element_text(face="bold"),legend.direction ="horizontal",
          axis.text.x = element_text(size=18,color="black"),axis.text.y = element_text(size=18,color="black"),
          axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
          legend.key = element_rect(fill = "gray60", colour = "black",size = 0.2),
          legend.background=element_rect(fill="white",colour="black",size=0.2),
          legend.text = element_text(size=24),
          legend.title = element_text(size=24),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=22,hjust = 0.01,vjust = -12),
          plot.caption = element_text(face="bold.italic"))
  #ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Test/300_Small_",Allcode[i],"_",Allmod[i],"_","Intersection.Mort.png"),plot = p, width = 12.37, height = 7.04, dpi=300,units = "in")
  #ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Test/300_Big_",Allcode[i],"_",Allmod[i],"_","Intersection.Mort.png"),plot = p, width = 20, height = 11.05, dpi=300,units = "in")
  #print(p)
  saveRDS(p,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure.intersection/pBIN",Allcode[i],".rds"))
  assign(paste0("p",Allcode[i]),p,envir = .GlobalEnv)
  rm("x")
}
#dev.off()
legend <- get_legend(p)
plot(legend)


i = 1
Allcode <- c("PINPIN","PINPINA","ALNGLU","PINNIG","POPTRE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.7,0.8,0.6,0.7,0.7,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.6,0.6,0.6,0.6,0.6,0.6)
Allmod <- c("Mbin15B.21","Mbin13A.18","Mbin13A.22","Mbin3B.20","Mbin13A.27","M2bin1B.23","Mbin13B.26","Mbin7A.26","Mbin1B.23")
for (i in 1:length(Allcode)){
  Dir = c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/") # Directory 
  test6 <- get(load(file=paste0(Dir,"species/",Allcode[i],"/CLIMAP/CLIMAP2019/Climap_",Allcode[i],"_",Allseuil[i],"_",AllseuilC[i],".RData")))
  Dir =c(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/binomial/",Allmod[i],"/"))
  setwd(Dir)
  #assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
  x <- get(load(file = paste0(Allmod[i],".rda")))
  # New section 
  x$data$Plotcat
  r1 <- data.frame(x$data$longitude, x$data$latitude, predict(x), x$data$Plotcat)
  colnames(r1) <- c("X","Y","Z","Cat")
  r1$Cat <- as.character(r1$Cat)
  str(r1)
  r1[r1$Cat=="0","Cat"] <- "Core"
  r1[r1$Cat=="1","Cat"] <- "LE" # To change when reverse. 2=LE in normal. 2=RE in reverse
  r1[r1$Cat=="2","Cat"] <- "TE"
  r1$bin <- ifelse(r1$Z>quantile(r1$Z,0.75,type=7),"> Q3",ifelse(r1$Z<quantile(r1$Z,0.25,type=7),"< Q1","Q1 - Q3")) ## Percentage of mortality as classes 
  print(summary(as.factor(r1$Cat)))
  print(Allcode[i])
  r1$bin <- as.factor(r1$bin)
  r1$bin = factor(r1$bin,levels(r1$bin)[c(1,3,2)])
  A <- melt(prop.table(table(r1$bin,r1$Cat), 2))
  B <- 100*prop.table(table(r1$Cat,r1$bin),1)
  p <- c(0.25,0.50,0.25)
  B1 <- unlist(chisq.test(B[1,],p=p,correct = F)[c(1,3)])
  B2 <- unlist(chisq.test(B[2,],p=p,correct = F)[c(1,3)])
  B3 <- unlist(chisq.test(B[3,],p=p,correct = F)[c(1,3)])
  B <- rbind(B1,B2,B3)
  B <- round(B,3)
  A[,"x2"] <- rep(B[,1],each=3)
  A[,"pvalue"] <- rep(B[,2],each=3)
  A$value <- round(A$value,3)*100
  A[c(2,5,8),"value"] <- A[c(2,5,8),"value"]/2 
  colnames(A) <- c("Quantiles","Zone","value","x2","pvalue")
  print(A)
  str(A$Zone)
  A$Zone = factor(A$Zone,levels(A$Zone)[c(3,1,2)])
  str(A$Zone)
  A[,"signif"] <- ifelse(A$pvalue<0.05,T,F)
  A[,"Psignif"] <- ""
  A[A$Quantiles=="Q1 - Q3","Psignif"] <- ifelse(A[A$Quantiles=="Q1 - Q3","pvalue"]<0.05,"*","")
  
  p <- ggplot(data=A,aes(x=Zone,y=value,fill=Quantiles,linetype=signif,alpha=signif,size=signif)) + theme(panel.background = element_rect(fill="white", colour="black", size=3, 
                                                                                                                                          linetype=1, color="black"),legend.key.size = unit(2, "cm"),legend.position = c(0.5,0.5)) + # 0.955,0.85 verticale. 
    theme(panel.grid.major = element_line(colour="blue", size=4, linetype=0,lineend="square", color="red"),
          panel.grid.minor = element_line(colour="blue", size=1, linetype=0,lineend="round", color="blue"))+
    theme(legend.background = element_rect(fill="white",size=1, linetype="solid", 
                                           colour ="black"))+
    theme(legend.position = "bottom",legend.box = "vertical")+ ## Allow us to disable the full legend
    #coord_fixed(ylim = 0,max(A$value)+max(A$value)/100*10,clip = 'off')+
    #ylim(0,max(A$value)+max(A$value)/100*10)+
    geom_bar(stat="identity",position="dodge",width = 0.65,colour="black")+
    #geom_text(aes(label=Psignif,y=max(A$value)+max(A$value)/100*10),vjust="top",colour="black",size=30,fontface = "bold")+
    geom_text(aes(label = Psignif,y=Inf),colour="black",size=30,fontface = "bold")+
    geom_text(aes(label=value),vjust=-0.2,colour="black",position=position_dodge(0.7),size=8,fontface = "bold")+
    #annotate("text", aes(label=(Psignif),y=max(A$value)+max(A$value)/100*10) +
    coord_cartesian(ylim = c(0,max(A$value)+max(A$value)/100*10), clip = "off")+
    #geom_text(aes(label=Label),vjust=-0.1,colour="black",position=position_dodge(0.7),size=3,fontface = "bold")+
    annotate("label", x = 1, y = 0, label = paste0("X-squared = ",A$x2[1],"; p-value = ",A$pvalue[1]),size=4.7,fontface = "bold")+
    annotate("label", x = 2, y = 0, label = paste0("X-squared = ",A$x2[4],"; p-value = ",A$pvalue[4]),size=4.7,fontface = "bold")+
    annotate("label", x = 3, y = 0, label = paste0("X-squared = ",A$x2[7],"; p-value = ",A$pvalue[7]),size=4.7,fontface = "bold")+
    scale_fill_manual(values = c("< Q1"="olivedrab1","Q1 - Q3"="chocolate1","> Q3"="darkred"),name="Quantile",
                      breaks=waiver(),
                      labels=c("Low (<Q1)", "Medium (>Q1 and <Q3)","High (>Q3)"))+
    scale_linetype_manual(name="Significiance",values = c("FALSE"="blank","TRUE"="solid"),
                          guide=F,
                          labels = c("Not significant", "Significant"))+
    scale_alpha_manual(name="Significiance",values = c("FALSE"=0.55,"TRUE"=1),
                       guide=F,
                       labels = c("Not significant", "Significant"))+
    scale_size_manual(name="Significiance",values = c("FALSE"=1,"TRUE"=1.25),
                      guide=F,
                      labels = c("Not significant", "Significant"))+
    labs(title=paste0(Allcode[i],"\nBIN"), y=paste0("Observed quantiles repartition in %"), x="Zone")+
    theme(text = element_text(face="bold"),legend.direction ="horizontal",
          axis.text.x = element_text(size=18,color="black"),axis.text.y = element_text(size=18,color="black"),
          axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
          legend.key = element_rect(fill = "gray60", colour = "black",size = 0.2),
          legend.background=element_rect(fill="white",colour="black",size=0.2),
          legend.text = element_text(size=24),
          legend.title = element_text(size=24),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=22,hjust = 0.01,vjust = -12),
          plot.caption = element_text(face="bold.italic"))
  #ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Test/300_Small_",Allcode[i],"_",Allmod[i],"_","Intersection.Mort.png"),plot = p, width = 12.37, height = 7.04, dpi=300,units = "in")
  #ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Test/300_Big_",Allcode[i],"_",Allmod[i],"_","Intersection.Mort.png"),plot = p, width = 20, height = 11.05, dpi=300,units = "in")
  #print(p)
  saveRDS(p,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure.intersection/pBIN",Allcode[i],".rds"))
  assign(paste0("p",Allcode[i]),p,envir = .GlobalEnv)
  rm("x")
}

# Cowplot 10 figures sans legend OK 
# pall <- plot_grid(
#   pABIALB+theme(legend.position = "none"),pACEPSE+theme(legend.position = "none"),
#   pALNGLU+theme(legend.position = "none"),pBETPEN+theme(legend.position = "none"),
#   pCASSAT+theme(legend.position = "none"),pFAGSYL+theme(legend.position = "none"), 
#   pFRAEXC+theme(legend.position = "none"), pPICABI+theme(legend.position = "none"),
#   pPINHAL+theme(legend.position = "none"), pPINNIG+theme(legend.position = "none"),align="hv",ncol=2,nrow=5)
# pall
# save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Intersection.Mort.bin.pdf"),plot = pall,base_width = 12.34, base_height = 7.04, dpi = 800 ,units = "in",nrow=6,ncol=2)
# save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Intersection.Mort.binV2.pdf"),plot = pall,base_width = 12.34, base_height = 7.04, dpi = 800 ,units = "in",nrow=5,ncol=2) #test

# 10 premiers + legend a part
ptest <- pABIALB+theme(legend.position=c(0.5,0.65))
legend <- get_legend(ptest)
plot(legend)

pall <- plot_grid(plot_grid(
  pABIALB+theme(legend.position = "none"),pACEPSE+theme(legend.position = "none"),
  pALNGLU+theme(legend.position = "none"),pBETPEN+theme(legend.position = "none"),
  pCASSAT+theme(legend.position = "none"),pFAGSYL+theme(legend.position = "none"), 
  pFRAEXC+theme(legend.position = "none"), pPICABI+theme(legend.position = "none"),
  pPINHAL+theme(legend.position = "none"), pPINNIG+theme(legend.position = "none"),ncol=2,nrow=5),
  legend,nrow = 2,rel_heights = c(1, 0.05), scale = c(1,1))
#pall

save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Intersection.Mort.bin1.legend.FINAL.pdf"),plot = pall,base_width = 12, base_height = 7, dpi = 300 ,units = "in",nrow=6,ncol=2)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Intersection.Mort.bin1.legend.FINAL.png"),plot = pall,base_width = 11, base_height = 7, dpi = 500 ,units = "in",nrow=6,ncol=2)

# 10 derniers + legend a part

pall <- plot_grid(plot_grid(
  pPINPIN+theme(legend.position = "none"),pPINPINA+theme(legend.position = "none"),
  pPINSYL+theme(legend.position = "none"),pPOPNIG+theme(legend.position = "none"),
  pPOPTRE+theme(legend.position = "none"),pQUEILE+theme(legend.position = "none"), 
  pQUEPET+theme(legend.position = "none"), pQUEPYR+theme(legend.position = "none"),
  pQUEROB+theme(legend.position = "none"), pQUESUB+theme(legend.position = "none"),align="hv",ncol=2,nrow=5),
  legend,nrow = 2,rel_heights = c(1, 0.05), scale = c(1,1))
#pall

save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Intersection.Mort.bin2.legend.FINAL.pdf"),plot = pall,base_width = 12, base_height = 7, dpi = 300 ,units = "in",nrow=6,ncol=2)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Intersection.Mort.bin2.legend.FINAL.png"),plot = pall,base_width = 12, base_height = 7, dpi = 500 ,units = "in",nrow=6,ncol=2)


##################
#### Amount ######
##################
i = 1
Allcode <- c("ABIALB","ACEPSE","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINSYL","QUEILE")
Allseuil <- c(0.7,0.55,0.8,0.7,0.7,0.7,0.8,0.7,0.8,0.8)
AllseuilC <- c(0.6,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
Allmod <- c("MnbZT14A.20","MnbZT13B.23","M2nbZT13B.29","M2nbZT13A.22","M2nbZT15B.24","MnbZT3C.24","M2nbZT7B.26","MnbZT11B.22","MnbZT5B.21","M2nbZT1C.27")

for (i in 1:length(Allcode)){
  Dir = c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/") # Directory 
  test6 <- get(load(file=paste0(Dir,"species/",Allcode[i],"/CLIMAP/CLIMAP2019/Climap_",Allcode[i],"_",Allseuil[i],"_",AllseuilC[i],".RData")))
  Dir =c(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/Negbin/",Allmod[i],"/"))
  setwd(Dir)
  #assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
  x <- get(load(file = paste0(Allmod[i],".rda")))
  # New section 
  x$data$Plotcat
  r1 <- data.frame(x$data$longitude, x$data$latitude, predict(x), x$data$Plotcat)
  colnames(r1) <- c("X","Y","Z","Cat")
  r1$Cat <- as.character(r1$Cat)
  str(r1)
  r1[r1$Cat=="0","Cat"] <- "Core"
  r1[r1$Cat=="2","Cat"] <- "LE" # To change when reverse. 2=LE in normal. 2=RE in reverse
  r1[r1$Cat=="1","Cat"] <- "TE"
  r1$bin <- ifelse(r1$Z>quantile(r1$Z,0.75,type=7),"> Q3",ifelse(r1$Z<quantile(r1$Z,0.25,type=7),"< Q1","Q1 - Q3")) ## Percentage of mortality as classes 
  print(summary(as.factor(r1$Cat)))
  print(Allcode[i])
  r1$bin <- as.factor(r1$bin)
  r1$bin = factor(r1$bin,levels(r1$bin)[c(1,3,2)])
  A <- melt(prop.table(table(r1$bin,r1$Cat), 2))
  B <- 100*prop.table(table(r1$Cat,r1$bin),1)
  p <- c(0.25,0.50,0.25)
  B1 <- unlist(chisq.test(B[1,],p=p,correct = F)[c(1,3)])
  B2 <- unlist(chisq.test(B[2,],p=p,correct = F)[c(1,3)])
  B3 <- unlist(chisq.test(B[3,],p=p,correct = F)[c(1,3)])
  B <- rbind(B1,B2,B3)
  B <- round(B,3)
  A[,"x2"] <- rep(B[,1],each=3)
  A[,"pvalue"] <- rep(B[,2],each=3)
  A$value <- round(A$value,3)*100
  A[c(2,5,8),"value"] <- A[c(2,5,8),"value"]/2 
  colnames(A) <- c("Quantiles","Zone","value","x2","pvalue")
  print(A)
  A$Zone = factor(A$Zone,levels(A$Zone)[c(3,1,2)])
  A[,"signif"] <- ifelse(A$pvalue<0.05,T,F)
  A[,"Psignif"] <- ""
  A[A$Quantiles=="Q1 - Q3","Psignif"] <- ifelse(A[A$Quantiles=="Q1 - Q3","pvalue"]<0.05,"*","")
  
  p <- ggplot(data=A,aes(x=Zone,y=value,fill=Quantiles,linetype=signif,alpha=signif,size=signif)) + theme(panel.background = element_rect(fill="white", colour="black", size=3, 
                                                                                                                                          linetype=1, color="black"),legend.key.size = unit(2, "cm"),legend.position = c(0.5,0.5)) + # 0.955,0.85 verticale. 
    theme(panel.grid.major = element_line(colour="blue", size=4, linetype=0,lineend="square", color="red"),
          panel.grid.minor = element_line(colour="blue", size=1, linetype=0,lineend="round", color="blue"))+
    theme(legend.background = element_rect(fill="white",size=1, linetype="solid", 
                                           colour ="black"))+
    theme(legend.position = "bottom",legend.box = "vertical")+ ## Allow us to disable the full legend
    #coord_fixed(ylim = 0,max(A$value)+max(A$value)/100*10,clip = 'off')+
    #ylim(0,max(A$value)+max(A$value)/100*10)+
    geom_bar(stat="identity",position="dodge",width = 0.65,colour="black")+
    #geom_text(aes(label=Psignif,y=max(A$value)+max(A$value)/100*10),vjust="top",colour="black",size=30,fontface = "bold")+
    geom_text(aes(label = Psignif,y=Inf),colour="black",size=30,fontface = "bold")+
    geom_text(aes(label=value),vjust=-0.2,colour="black",position=position_dodge(0.7),size=8,fontface = "bold")+
    #annotate("text", aes(label=(Psignif),y=max(A$value)+max(A$value)/100*10) +
    coord_cartesian(ylim = c(0,max(A$value)+max(A$value)/100*10), clip = "off")+
    #geom_text(aes(label=Label),vjust=-0.1,colour="black",position=position_dodge(0.7),size=3,fontface = "bold")+
    annotate("label", x = 1, y = 0, label = paste0("X-squared = ",A$x2[1],"; p-value = ",A$pvalue[1]),size=4.7,fontface = "bold")+
    annotate("label", x = 2, y = 0, label = paste0("X-squared = ",A$x2[4],"; p-value = ",A$pvalue[4]),size=4.7,fontface = "bold")+
    annotate("label", x = 3, y = 0, label = paste0("X-squared = ",A$x2[7],"; p-value = ",A$pvalue[7]),size=4.7,fontface = "bold")+
    scale_fill_manual(values = c("< Q1"="olivedrab1","Q1 - Q3"="chocolate1","> Q3"="darkred"),name="Quantile",
                      breaks=waiver(),
                      labels=c("Low (<Q1)", "Medium (>Q1 and <Q3)","High (>Q3)"))+
    scale_linetype_manual(name="Significiance",values = c("FALSE"="blank","TRUE"="solid"),
                          guide=F,
                          labels = c("Not significant", "Significant"))+
    scale_alpha_manual(name="Significiance",values = c("FALSE"=0.55,"TRUE"=1),
                       guide=F,
                       labels = c("Not significant", "Significant"))+
    scale_size_manual(name="Significiance",values = c("FALSE"=1,"TRUE"=1.25),
                      guide=F,
                      labels = c("Not significant", "Significant"))+
    labs(title=paste0(Allcode[i],"\nNB"), y=paste0("Observed quantiles repartition in %"), x="Zone")+
    theme(text = element_text(face="bold"),legend.direction ="horizontal",
          axis.text.x = element_text(size=18,color="black"),axis.text.y = element_text(size=18,color="black"),
          axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
          legend.key = element_rect(fill = "gray60", colour = "black",size = 0.2),
          legend.background=element_rect(fill="white",colour="black",size=0.2),
          legend.text = element_text(size=24),
          legend.title = element_text(size=24),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=22,hjust = 0.01,vjust = -12),
          plot.caption = element_text(face="bold.italic"))
  #ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Test/300_Small_",Allcode[i],"_",Allmod[i],"_","Intersection.Mort.png"),plot = p, width = 12.37, height = 7.04, dpi=300,units = "in")
  #ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Test/300_Big_",Allcode[i],"_",Allmod[i],"_","Intersection.Mort.png"),plot = p, width = 20, height = 11.05, dpi=300,units = "in")
  #print(p)
  saveRDS(p,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure.intersection/pNB",Allcode[i],".rds"))
  assign(paste0("p",Allcode[i]),p,envir = .GlobalEnv)
  rm("x")
}


i = 1
Allcode <- c("PINPINA","ALNGLU","PINNIG","PINPIN","POPTRE","QUEPET","QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.6,0.7,0.7,0.7,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
Allmod <- c("MnbZT13A.19","MnbZT13A.21","MnbZT3B.27","MnbZT15B.28","MnbZT13A.33","M2nbZT1B.24","MnbZT13B.27","MnbZT7A.23","MnbZT1B.27")
for (i in 1:length(Allcode)){
  Dir = c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/") # Directory 
  test6 <- get(load(file=paste0(Dir,"species/",Allcode[i],"/CLIMAP/CLIMAP2019/Climap_",Allcode[i],"_",Allseuil[i],"_",AllseuilC[i],".RData")))
  Dir =c(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/Negbin/",Allmod[i],"/"))
  setwd(Dir)
  #assign(sub(".rda","",list.files(pattern=".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern=".rda")))) 
  x <- get(load(file = paste0(Allmod[i],".rda")))
  # New section 
  x$data$Plotcat
  r1 <- data.frame(x$data$longitude, x$data$latitude, predict(x), x$data$Plotcat)
  colnames(r1) <- c("X","Y","Z","Cat")
  r1$Cat <- as.character(r1$Cat)
  str(r1)
  r1[r1$Cat=="0","Cat"] <- "Core"
  r1[r1$Cat=="1","Cat"] <- "LE" # To change when reverse. 2=LE in normal. 2=RE in reverse
  r1[r1$Cat=="2","Cat"] <- "TE"
  r1$bin <- ifelse(r1$Z>quantile(r1$Z,0.75,type=7),"> Q3",ifelse(r1$Z<quantile(r1$Z,0.25,type=7),"< Q1","Q1 - Q3")) ## Percentage of mortality as classes 
  print(summary(as.factor(r1$Cat)))
  print(Allcode[i])
  r1$bin <- as.factor(r1$bin)
  r1$bin = factor(r1$bin,levels(r1$bin)[c(1,3,2)])
  A <- melt(prop.table(table(r1$bin,r1$Cat), 2))
  B <- 100*prop.table(table(r1$Cat,r1$bin),1)
  p <- c(0.25,0.50,0.25)
  B1 <- unlist(chisq.test(B[1,],p=p,correct = F)[c(1,3)])
  B2 <- unlist(chisq.test(B[2,],p=p,correct = F)[c(1,3)])
  B3 <- unlist(chisq.test(B[3,],p=p,correct = F)[c(1,3)])
  B <- rbind(B1,B2,B3)
  B <- round(B,3)
  A[,"x2"] <- rep(B[,1],each=3)
  A[,"pvalue"] <- rep(B[,2],each=3)
  A$value <- round(A$value,3)*100
  A[c(2,5,8),"value"] <- A[c(2,5,8),"value"]/2 
  colnames(A) <- c("Quantiles","Zone","value","x2","pvalue")
  print(A)
  A$Zone = factor(A$Zone,levels(A$Zone)[c(3,1,2)])
  A[,"signif"] <- ifelse(A$pvalue<0.05,T,F)
  A[,"Psignif"] <- ""
  A[A$Quantiles=="Q1 - Q3","Psignif"] <- ifelse(A[A$Quantiles=="Q1 - Q3","pvalue"]<0.05,"*","")
  
  p <- ggplot(data=A,aes(x=Zone,y=value,fill=Quantiles,linetype=signif,alpha=signif,size=signif)) + theme(panel.background = element_rect(fill="white", colour="black", size=3, 
                                                                                                                                          linetype=1, color="black"),legend.key.size = unit(2, "cm"),legend.position = c(0.5,0.5)) + # 0.955,0.85 verticale. 
    theme(panel.grid.major = element_line(colour="blue", size=4, linetype=0,lineend="square", color="red"),
          panel.grid.minor = element_line(colour="blue", size=1, linetype=0,lineend="round", color="blue"))+
    theme(legend.background = element_rect(fill="white",size=1, linetype="solid", 
                                           colour ="black"))+
    theme(legend.position = "bottom",legend.box = "vertical")+ ## Allow us to disable the full legend
    #coord_fixed(ylim = 0,max(A$value)+max(A$value)/100*10,clip = 'off')+
    #ylim(0,max(A$value)+max(A$value)/100*10)+
    geom_bar(stat="identity",position="dodge",width = 0.65,colour="black")+
    #geom_text(aes(label=Psignif,y=max(A$value)+max(A$value)/100*10),vjust="top",colour="black",size=30,fontface = "bold")+
    geom_text(aes(label = Psignif,y=Inf),colour="black",size=30,fontface = "bold")+
    geom_text(aes(label=value),vjust=-0.2,colour="black",position=position_dodge(0.7),size=8,fontface = "bold")+
    #annotate("text", aes(label=(Psignif),y=max(A$value)+max(A$value)/100*10) +
    coord_cartesian(ylim = c(0,max(A$value)+max(A$value)/100*10), clip = "off")+
    #geom_text(aes(label=Label),vjust=-0.1,colour="black",position=position_dodge(0.7),size=3,fontface = "bold")+
    annotate("label", x = 1, y = 0, label = paste0("X-squared = ",A$x2[1],"; p-value = ",A$pvalue[1]),size=4.7,fontface = "bold")+
    annotate("label", x = 2, y = 0, label = paste0("X-squared = ",A$x2[4],"; p-value = ",A$pvalue[4]),size=4.7,fontface = "bold")+
    annotate("label", x = 3, y = 0, label = paste0("X-squared = ",A$x2[7],"; p-value = ",A$pvalue[7]),size=4.7,fontface = "bold")+
    scale_fill_manual(values = c("< Q1"="olivedrab1","Q1 - Q3"="chocolate1","> Q3"="darkred"),name="Quantile",
                      breaks=waiver(),
                      labels=c("Low (<Q1)", "Medium (>Q1 and <Q3)","High (>Q3)"))+
    scale_linetype_manual(name="Significiance",values = c("FALSE"="blank","TRUE"="solid"),
                          guide=F,
                          labels = c("Not significant", "Significant"))+
    scale_alpha_manual(name="Significiance",values = c("FALSE"=0.55,"TRUE"=1),
                       guide=F,
                       labels = c("Not significant", "Significant"))+
    scale_size_manual(name="Significiance",values = c("FALSE"=1,"TRUE"=1.25),
                      guide=F,
                      labels = c("Not significant", "Significant"))+
    labs(title=paste0(Allcode[i],"\nNB"), y=paste0("Observed quantiles repartition in %"), x="Zone")+
    theme(text = element_text(face="bold"),legend.direction ="horizontal",
          axis.text.x = element_text(size=18,color="black"),axis.text.y = element_text(size=18,color="black"),
          axis.title.x = element_text(size=18),axis.title.y = element_text(size=18),
          legend.key = element_rect(fill = "gray60", colour = "black",size = 0.2),
          legend.background=element_rect(fill="white",colour="black",size=0.2),
          legend.text = element_text(size=24),
          legend.title = element_text(size=24),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=22,hjust = 0.01,vjust = -12),
          plot.caption = element_text(face="bold.italic"))
  #ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Test/300_Small_",Allcode[i],"_",Allmod[i],"_","Intersection.Mort.png"),plot = p, width = 12.37, height = 7.04, dpi=300,units = "in")
  #ggsave(filename = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Test/300_Big_",Allcode[i],"_",Allmod[i],"_","Intersection.Mort.png"),plot = p, width = 20, height = 11.05, dpi=300,units = "in")
  #print(p)
  saveRDS(p,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure.intersection/pNB",Allcode[i],".rds"))
  assign(paste0("p",Allcode[i]),p,envir = .GlobalEnv)
  rm("x")
}

pall <- plot_grid(plot_grid(
  pABIALB+theme(legend.position = "none"),pACEPSE+theme(legend.position = "none"),
  pALNGLU+theme(legend.position = "none"),pBETPEN+theme(legend.position = "none"),
  pCASSAT+theme(legend.position = "none"),pFAGSYL+theme(legend.position = "none"), 
  pFRAEXC+theme(legend.position = "none"), pPICABI+theme(legend.position = "none"),
  pPINHAL+theme(legend.position = "none"), pPINNIG+theme(legend.position = "none"),align="hv",ncol=2,nrow=5),
  legend,nrow = 2,rel_heights = c(1, 0.05), scale = c(1,1))
#pall

save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Intersection.Mort.NB1.legend.FINAL.pdf"),plot = pall,base_width = 12, base_height = 7, dpi = 300 ,units = "in",nrow=6,ncol=2)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Intersection.Mort.NB1.legend.FINAL.png"),plot = pall,base_width = 12, base_height = 7, dpi = 500 ,units = "in",nrow=6,ncol=2)

# 10 derniers + legend a part

pall <- plot_grid(plot_grid(
  pPINPIN+theme(legend.position = "none"),pPINPINA+theme(legend.position = "none"),
  pPINSYL+theme(legend.position = "none"),
  pPOPTRE+theme(legend.position = "none"),pQUEILE+theme(legend.position = "none"), 
  pQUEPET+theme(legend.position = "none"), pQUEPYR+theme(legend.position = "none"),
  pQUEROB+theme(legend.position = "none"), pQUESUB+theme(legend.position = "none"),align="hv",ncol=2,nrow=5),
  legend,nrow = 2,rel_heights = c(1, 0.05), scale = c(1,1))
#pall

save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Intersection.Mort.NB2.legend.FINAL.pdf"),plot = pall,base_width = 12, base_height = 7, dpi = 300 ,units = "in",nrow=6,ncol=2)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Intersection.Mort.NB2.legend.FINAL.png"),plot = pall,base_width = 12, base_height = 7, dpi = 500 ,units = "in",nrow=6,ncol=2)


