rm(list = ls())
gc()
library(glmmTMB)
library(car)
library(emmeans)
library(effects)
library(multcomp)
library(MuMIn)
library(DHARMa)
library(broom)
library(broom.mixed)
library(dotwhisker)
library(ggplot2); theme_set(theme_bw())
library(texreg)
library(xtable)
library(huxtable)
library(plyr)
library(beanplot)
library(parallel)


### load all the data 

Dir <- c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data")
setwd("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/")

mod_ALLbin <- list()
mod_ALLneg <- list()
dfplot <- list()
dfplot2 <- list()
dfplot3 <- list()
dftree <- list()

# Then try to obtain predicted mortality and observed plot mortality 
i = 1
Allcode <- c("ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINPINA","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
Allmod <- c("Mbin14A.19","Mbin13B.26","Mbin13A.22","M2bin13B.22","M2bin13A21","M2bin15B23","Mbin3C.26","M2bin7B.17","Mbin13A.18","Mbin11B.19","Mbin3B.20","Mbin15B.21","Mbin5B.17","Mbin3B.31","Mbin13A.27","M2bin1C.20","M2bin1B.23","Mbin13B.26","Mbin7A.26","Mbin1B.23")

for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir = c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/") # Directory 
  Dir =c(paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/binomial/",Allmod[i],"/"))
  setwd(Dir)
  #mod_ALLbin[[CODE]] <- get(load(file = list.files(pattern=".rda")))$data
  mod_ALLbin[[CODE]] <- get(load(file = list.files(pattern=".rda")))
  rm(x)
  Dir = paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/CLIMAP2019/")
  setwd(Dir)
  # dfplot[[CODE]] <- readRDS(paste0("dfplot", CODE, seuil, ".rds")) #Base de données plot full No scale !!! sans zone de transition
  # dfplot2[[CODE]] <- readRDS(paste0("dfplot2", CODE, seuil, ".rds")) #Base de données plot full No scale !!! sans zone de transition
  # dfplot3[[CODE]] <- readRDS(paste0("dfplot3", CODE, seuil, ".rds")) #Base de données plot full No scale !!! sans zone de transition
  # #dftree[[CODE]] <- readRDS(paste0("Mydf_",CODE,"_",seuil,"_",seuilC,".rds"))
}

## load negbinom models 

i = 1
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
Allmod <- c("MnbZT13A.19","MnbZT14A.20","MnbZT13B.23","MnbZT13A.21","M2nbZT13B.29","M2nbZT13A.22","M2nbZT15B.24","MnbZT3C.24","M2nbZT7B.26","MnbZT11B.22","MnbZT3B.27","MnbZT15B.28","MnbZT5B.21","MnbZT13A.33","M2nbZT1C.27","M2nbZT1B.24","MnbZT13B.27","MnbZT7A.23","MnbZT1B.27")
for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  seuilC <- AllseuilC[i]
  seuil <- Allseuil[i]
  Dir = c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/") # Directory 
  Dir =c(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",Allcode[i],"/CLIMAP/Models/Negbin/",Allmod[i],"/"))
  setwd(Dir)
  #mod_ALLneg[[CODE]] <- get(load(file = list.files(pattern=".rda")))$data
  mod_ALLneg[[CODE]] <- get(load(file = list.files(pattern=".rda")))
  rm(x)
}


### keep fit and observed !

mod_ALLbin <- mclapply(mod_ALLbin, function(x) {
  x[["data"]][,"fit"] <- predict(x)
  ;x[["data"]][,c("sp.mortality.plot.rate.yr","fit")]},mc.cores = 18)

mod_ALLneg <- mclapply(mod_ALLneg, function(x) {
  x[["data"]][,"fit"] <- predict(x)
  ;x[["data"]][,c("sp.mortality.plot.rate.yr","fit")]},mc.cores = 18)

## Convert it to a df
mod_ALLbin <- mapply(function(x,y) {
  x[,"species"] <- y
  return(x)
},
x=mod_ALLbin,y=names(mod_ALLbin),SIMPLIFY = F)

mod_ALLneg <- mapply(function(x,y) {
  x[,"species"] <- y
  return(x)
},
x=mod_ALLneg,y=names(mod_ALLneg),SIMPLIFY = F)

mod_ALLbin <- do.call(rbind,mod_ALLbin) # From list to df
mod_ALLneg <- do.call(rbind,mod_ALLneg) # From list to df

### reshape to have a factor for obseved VS fit 
library(reshape)
mod_ALLneg2 <- melt(mod_ALLneg)
mod_ALLbin2 <- melt(mod_ALLbin)

### Boxplot NEGBIN pred and obs
p3 <- ggplot(mod_ALLneg2,aes(x=species, y=value,colour=variable))+
  #geom_point(size=0.3,alpha=0.1)+
  geom_boxplot(outlier.size = 1,outlier.stroke = 0.5,outlier.alpha = 0.1)+
  stat_summary(fun.y=mean, geom="point", shape=23, size=5,position = position_dodge2(width = 0.8,preserve = c("total"),padding = 0.1))+
  scale_colour_manual(NULL,
                     values = c("sp.mortality.plot.rate.yr"="black","fit"="blue"),
                     labels=c("sp.mortality.plot.rate.yr"="Observed mortality rate","fit"="Predicted mortality rate"))+
  labs(y="Mortality rates (‰ of tree/hectare/year)", x="Species")+
  #scale_shape_manual("", values=c("Mean values"="\U25C7"))+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.36,0.85),legend.direction ="vertical",legend.key.size = unit(1.5, "cm"),
        #legend.key = element_rect(fill = "white", colour = "black",size = 0.2),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        axis.text.x = element_text(size=11,color="black",angle = 45,vjust = 0.6),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
p3

## take the legend of the second plot and put it on the new one 
legend <- get_legend(p2)
plot(legend)
pall <- p3+annotation_custom(legend, xmin = 12, xmax = 12, # ABIALB,ACEPSE,ALNGLU, QUEROB
                                ymin = 280, ymax = 280)

pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure_S7.pdf"),plot = pall,base_width = 12, base_height = 7, dpi = 300 ,units = "in",nrow=1,ncol=1)


### Boxplot BIN occurrence
p2 <- ggplot(mod_ALLbin,aes(x=species, y=fit))+
  #geom_point(size=0.3,alpha=0.1)+
  geom_boxplot(col="blue",outlier.size = 1,outlier.stroke = 0.5,outlier.alpha = 0.1)+
  #geom_jitter(size=0.3,alpha=0.05,position=position_jitter(0.4))+
  stat_summary(fun.y=mean,aes(shape="Mean values"),geom="point",col="blue", size=5,position = position_dodge2(width = 0.8,preserve = c("total"),padding = 0.1))+
  scale_shape_manual("", values=c("Mean values"="\U25C7"))+
  # scale_colour_manual(NULL,
  #                     values = c("sp.mortality.plot.rate.yr"="black","fit"="blue"),
  #                     labels=c("sp.mortality.plot.rate.yr"="Observed mortality rate","fit"="Predicted mortality rate"))+
  labs(y="Occurrence of mortality probability", x="Species")+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.5,0.9),legend.direction ="horizontal",legend.key.size = unit(1.5, "cm"),
        legend.key = element_rect(fill = "white", colour = "black",size = 0.2),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        axis.text.x = element_text(size=11,color="black",angle = 45,vjust = 0.6),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
cairo_pdf("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure_S7.pdf", family="Helvetica",width = 12, height = 7)
p2
dev.off()
#save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure_S6.pdf"),plot = p2,base_width = 12, base_height = 7, dpi = 300 ,units = "in",nrow=1,ncol=1)


## Boxplot BIN observed
p4 <- ggplot(mod_ALLbin,aes(x=species, y=sp.mortality.plot.rate.yr))+
  #geom_point(size=0.3,alpha=0.1)+
  geom_boxplot(col="black",outlier.size = 1,outlier.stroke = 0.5,outlier.alpha = 0.1)+
  #geom_jitter(size=0.3,alpha=0.05,position=position_jitter(0.4))+
  stat_summary(fun.y=mean,aes(shape="Mean values"),geom="point",col="black", size=5,position = position_dodge2(width = 0.8,preserve = c("total"),padding = 0.1))+
  scale_shape_manual("", values=c("Mean values"="\U25C7"))+
  labs(y="Observed mortality rates (‰ of tree/hectare/year)", x="Species")+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.5,0.9),legend.direction ="horizontal",legend.key.size = unit(1.5, "cm"),
        legend.key = element_rect(fill = "white", colour = "black",size = 0.2),
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        axis.text.x = element_text(size=11,color="black",angle = 45,vjust = 0.6),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
cairo_pdf("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure_S5.pdf", family="Helvetica",width = 12, height = 7)
p4
dev.off()

#save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure_S5.pdf"),plot = p4,base_width = 12, base_height = 7, dpi = 300 ,units = "in",nrow=1,ncol=1)



## Fit vs obs (intensity only) ## bonus 
dplot <- ggplot(mod_ALLneg) + 
  geom_line(aes(x=sp.mortality.plot.rate.yr,color="observed"),size=1.5,key_glyph = "abline",stat = "density") +
  geom_line(aes(x=fit,color="fit"),size=1.5,key_glyph = "abline", stat = "density") +
  #key_glyph = "rect"
  facet_wrap(vars(species),scales = "free")+
  scale_color_manual("Counts",values = c("observed"="black","fit"="blue"),labels=c("observed"="Observed distribution","fit"="Fitted distribution"))+
  labs(y="Frequency", x="Mortality rates (‰ of tree/hectare/year)",las=1)+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.9,0.1),legend.direction ="vertical",
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.8, "cm"),
        panel.border = element_rect(colour = "black",fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))

cairo_pdf("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure_S4.pdf", family="Helvetica",width = 12, height = 7)
dplot
dev.off()
#save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/Figure_S4.pdf"),plot = dplot,base_width = 12, base_height = 7, dpi = 500 ,units = "in",nrow=1,ncol=1)


## Fit vs obs (intensity only) (residuals) ## bonus too 
p3 <- ggplot(mod_ALLneg, aes(x=sp.mortality.plot.rate.yr, y=fit)) + 
  geom_point(size=0.3,alpha=0.1)+
  geom_abline(slope = 1,intercept = 0,size=0.3)+
  #geom_violin()+
  facet_wrap(vars(species),scales="free")+
  #scale_color_manual("",values = Mycol,labels=Mylabels)+
  geom_jitter(size=0.3,alpha=0.05,position=position_jitter(0.4))+
  #geom_boxplot(width=0.1,outlier.size = 0)+
  labs(y="Fitted values", x="Observed values")+
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(0.8,0.05),legend.direction ="horizontal",
        axis.text.x = element_text(size=11,color="black"),axis.text.y = element_text(size=13,color="black"),
        #axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
p3
