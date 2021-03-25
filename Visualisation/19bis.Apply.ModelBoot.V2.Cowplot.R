rm(list = ls())
gc()
require(ade4)
library(sjPlot) 
library(lattice)
library(latticeExtra)
library(stringr)
library(data.table)
library(pscl)
library(MASS)
library(lme4)
library("glmmTMB")
library("bbmle")
library(piecewiseSEM)
library(pgirmess)
library(MuMIn)
library(spaMM)

Dir <- c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/")
setwd(Dir)
# Load the function 
source(paste0(Dir,"Myscripts/Fundiv.project/function3.ModelBoot.R"))
i = 1
Allcode <- c("ABIALB","ACEPSE","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINSYL","POPNIG","QUEILE",
             "PINPIN","PINPINA","ALNGLU","PINNIG","POPTRE","QUEPET","QUEPYR","QUEROB","QUESUB",
             "ABIALB","ACEPSE","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINSYL","QUEILE",
             "PINPINA","ALNGLU","PINNIG","PINPIN","POPTRE","QUEPET","QUEPYR","QUEROB","QUESUB"
)

Allmod <- c("Mbin14A.19","Mbin13B.26","M2bin13B.22","M2bin13A21","M2bin15B23","Mbin3C.26","M2bin7B.17","Mbin11B.19","Mbin5B.17","Mbin3B.31","M2bin1C.20",
            "Mbin15B.21","Mbin13A.18","Mbin13A.22","Mbin3B.20","Mbin13A.27","M2bin1B.23","Mbin13B.26","Mbin7A.26","Mbin1B.23",
            "MnbZT14A.20","MnbZT13B.23","M2nbZT13B.29","M2nbZT13A.22","M2nbZT15B.24","MnbZT3C.24","M2nbZT7B.26","MnbZT11B.22","MnbZT5B.21","M2nbZT1C.27",
            "MnbZT13A.19","MnbZT13A.21","MnbZT3B.27","MnbZT15B.28","MnbZT13A.33","M2nbZT1B.24","MnbZT13B.27","MnbZT7A.23","MnbZT1B.27"
)
Allseuil <- c(0.7,0.55,0.8,0.7,0.7,0.7,0.8,0.7,0.8,0.7,0.8,
              0.7,0.8,0.6,0.7,0.7,0.7,0.7,0.8,0.7,
              0.7,0.55,0.8,0.7,0.7,0.7,0.8,0.7,0.8,0.8,
              0.8,0.6,0.7,0.7,0.7,0.7,0.7,0.8,0.7
)

## Here extaction of parameters in a document for all models (two by species) to then extract the right numbers to run the function 
test.files5 <- file(paste0(Dir,"our-data/species/test.extract.Rout"), open="wt")
sink(test.files5, type="output")

for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  mod <- Allmod[i]
  seuil <- Allseuil[i]
  Dir = (paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP")) # Directory 
  setwd(Dir)
  #list.files(Dir,pattern = paste0(seuil,".rds"))
  dfplot <- readRDS(paste0("dfplot",CODE,seuil,".rds")) #Base de données plot
  #dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
  dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
  dfplotbis <- readRDS(paste0("dfplotbis",CODE,seuil,".rds")) #Base de données plot
  # This one is the destination where only the model of interest is (binom and negbin)
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/",mod,"/"))
  setwd(Dir)
  assign(eval(mod),get(load(file = paste0(mod,".rda"))))
  rm("x")
  #summary(get(mod))
  print(paste0(CODE,".",eval(mod)))
  print(paste0(Extraction(get(mod))))
  print(paste0("\n"))
}
close(test.files5) 



#### Run the normal function bootstrap for all three regions and not reverse species 
i = 1
Allcode <- c("ABIALB","PICABI","PICABI","PICABI","PICABI",
             "PICABI","PINSYL","PINSYL","PINSYL","ABIALB",
             "ABIALB","BETPEN","BETPEN","FAGSYL","FAGSYL",
             "FRAEXC","FAGSYL","BETPEN","BETPEN","PINHAL",
             "PICABI","PICABI","PINSYL","PINSYL","QUEILE",
             "QUEILE","PINHAL","ACEPSE","ACEPSE","QUEILE")

Allmod <- c("Mbin14A.19","M2bin7B.17","M2bin7B.17","M2bin7B.17","M2bin7B.17",
            "M2bin7B.17","Mbin5B.17","Mbin5B.17","Mbin5B.17","MnbZT14A.20",
            "MnbZT14A.20","M2nbZT13B.29","M2nbZT13B.29","M2nbZT15B.24","M2nbZT15B.24",
            "MnbZT3C.24","M2bin15B23","M2bin13B.22","M2bin13B.22","Mbin11B.19",
            "M2nbZT7B.26","M2nbZT7B.26","MnbZT5B.21","MnbZT5B.21","M2nbZT1C.27",
            "M2nbZT1C.27","MnbZT11B.22","Mbin13B.26","Mbin13B.26","M2bin1C.20"
)

Allseuil <- c(0.7,0.8,0.8,0.8,0.8,
              0.8,0.8,0.8,0.8,0.7,
              0.7,0.8,0.8,0.7,0.7,
              0.7,0.7,0.8,0.8,0.7,
              0.8,0.8,0.8,0.8,0.8,
              0.8,0.7,0.55,0.55,0.8
)


ParaInter <- c(11,9,10,7,4,
               5,8,6,5,5,
               6,6,7,6,7,
               8,9,7,5,8,
               9,5,11,5,8,
               4,5,9,4,5
)
ParaPlotcat <- c(7,8,8,8,8,
                 8,7,7,7,11,
                 11,5,5,10,10,
                 10,8,10,10,11,
                 7,7,8,8,9,
                 9,6,10,10,10
)
BootstrapFigureErrors <- file(paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/BootstrapFigureErrors.txt"), open="wt")
sink(BootstrapFigureErrors, type="message")
for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  mod <- Allmod[i]
  seuil <- Allseuil[i]
  PInter <- ParaInter[i]
  PPlotcat <- ParaPlotcat[i]
  Dir = (paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP")) # Directory 
  setwd(Dir)
  #list.files(Dir,pattern = paste0(seuil,".rds"))
  dfplot <- readRDS(paste0("dfplot",CODE,seuil,".rds")) #Base de données plot
  dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
  dfplotbis <- readRDS(paste0("dfplotbis",CODE,seuil,".rds")) #Base de données plot
  # This one is the destination where only the model of interest is (binom and negbin)
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/",mod,"/"))
  setwd(Dir)
  assign(eval(mod),get(load(file = paste0(mod,".rda"))))
  rm("x")
  Extraction(get(mod))
  try(ModelBoot.All(get(mod),PInter,PPlotcat,LvL=30,CAT=CAT,nBoot=10,Yportion = 0.66,saveboot=T,nCoeur=13),silent=T)
}
close(BootstrapFigureErrors) 



#### Run the normal function bootstrap for all three regions and reverse species 
i = 1
Allcode <- c("PINPIN","PINPIN","PINPIN","PINNIG",
             "PINNIG","PINNIG","QUEPET","QUEROB",
             "PINPINA","PINPINA","PINPINA","POPTRE",
             "PINPINA","PINPINA","ALNGLU","PINNIG",
             "PINNIG")

Allmod <- c("Mbin15B.21","Mbin15B.21","Mbin15B.21","Mbin3B.20",
            "Mbin3B.20","Mbin3B.20","M2nbZT1B.24","MnbZT7A.23",
            "Mbin13A.18","Mbin13A.18","Mbin13A.18","Mbin13A.27",
            "MnbZT13A.19","MnbZT13A.19","MnbZT13A.21","MnbZT3B.27",
            "MnbZT3B.27")

Allseuil <- c(0.7,0.7,0.7,0.7,
              0.7,0.7,0.7,0.8,
              0.8,0.8,0.8,0.7,
              0.8,0.8,0.6,0.7,
              0.7)

ParaInter <- c(8,3,4,5,
               6,2,9,8,
               8,9,10,6,
               11,12,11,7,
               9)

ParaPlotcat <- c(9,9,9,8,
                 8,8,11,10,
                 6,6,6,10,
                 9,9,6,8,
                 8)

BootstrapFigureErrors <- file(paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/BootstrapFigureErrors.txt"), open="wt")
sink(BootstrapFigureErrors, type="message")
for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  mod <- Allmod[i]
  seuil <- Allseuil[i]
  PInter <- ParaInter[i]
  PPlotcat <- ParaPlotcat[i]
  Dir = (paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP")) # Directory 
  setwd(Dir)
  #list.files(Dir,pattern = paste0(seuil,".rds"))
  dfplot <- readRDS(paste0("dfplot",CODE,seuil,".rds")) #Base de données plot
  dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
  dfplotbis <- readRDS(paste0("dfplotbis",CODE,seuil,".rds")) #Base de données plot
  # This one is the destination where only the model of interest is (binom and negbin)
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/",mod,"/"))
  setwd(Dir)
  assign(eval(mod),get(load(file = paste0(mod,".rda"))))
  rm("x")
  Extraction(get(mod))
  try(ModelBoot.r.All(get(mod),PInter,PPlotcat,LvL=30,CAT=CAT,nBoot=10,Yportion = 0.66,saveboot=T,nCoeur=13),silent=T)
}
close(BootstrapFigureErrors)



### Normal sens -LE


i = 1
Allcode <- c("CASSAT","POPNIG")
Allmod <- c("M2bin13A.21","Mbin3B.31")
Allseuil <- c(0.7,0.7)
ParaInter <- c(8,7)
ParaPlotcat <- c(10,9)
BootstrapFigureErrors <- file(paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/BootstrapFigureErrors.txt"), open="wt")
sink(BootstrapFigureErrors, type="message")
for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  mod <- Allmod[i]
  seuil <- Allseuil[i]
  PInter <- ParaInter[i]
  PPlotcat <- ParaPlotcat[i]
  Dir = (paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP")) # Directory 
  setwd(Dir)
  #list.files(Dir,pattern = paste0(seuil,".rds"))
  dfplot <- readRDS(paste0("dfplot",CODE,seuil,".rds")) #Base de données plot
  dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
  dfplotbis <- readRDS(paste0("dfplotbis",CODE,seuil,".rds")) #Base de données plot
  # This one is the destination where only the model of interest is (binom and negbin)
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/",mod,"/"))
  setwd(Dir)
  assign(eval(mod),get(load(file = paste0(mod,".rda"))))
  rm("x")
  Extraction(get(mod))
  try(ModelBoot.RE(get(mod),PInter,PPlotcat,LvL=30,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=13),silent=T)
}
close(BootstrapFigureErrors)

### Reverse sens - LE

i = 1
Allcode <- c("QUESUB","QUEROB","QUEROB","QUEPYR","QUESUB",
             "QUESUB","ALNGLU","ALNGLU","POPTRE","QUEROB")
Allmod <- c("Mbin1B.23","Mbin7A.26","Mbin7A.26","Mbin13B.26","Mbin1B.23",
            "Mbin1B.23","MnbZT13A.21","MnbZT13A.21","MnbZT13A.33","MnbZT7A.23")
Allseuil <- c(0.7,0.8,0.8,0.7,0.7,
              0.7,0.6,0.6,0.7,0.8)
ParaInter <- c(5,8,6,2,9,
               10,9,10,6,5)
ParaPlotcat <- c(11,12,12,8,11,
                 11,7,7,9,10)
BootstrapFigureErrors <- file(paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/BootstrapFigureErrors.txt"), open="wt")
sink(BootstrapFigureErrors, type="message")
for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  mod <- Allmod[i]
  seuil <- Allseuil[i]
  PInter <- ParaInter[i]
  PPlotcat <- ParaPlotcat[i]
  Dir = (paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP")) # Directory 
  setwd(Dir)
  #list.files(Dir,pattern = paste0(seuil,".rds"))
  dfplot <- readRDS(paste0("dfplot",CODE,seuil,".rds")) #Base de données plot
  dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
  dfplotbis <- readRDS(paste0("dfplotbis",CODE,seuil,".rds")) #Base de données plot
  # This one is the destination where only the model of interest is (binom and negbin)
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/",mod,"/"))
  setwd(Dir)
  assign(eval(mod),get(load(file = paste0(mod,".rda"))))
  rm("x")
  Extraction(get(mod))
  try(ModelBoot.r.RE(get(mod),PInter,PPlotcat,LvL=30,CAT=CAT,nBoot=10,Yportion = 0.66,saveboot=T,nCoeur=13),silent=T)
}
close(BootstrapFigureErrors)


##### Normale - RE 

i = 1
Allcode <- c("ACEPSE","ACEPSE")
Allmod <- c("MnbZT13B.23","MnbZT13B.23")
Allseuil <- c(0.55)
ParaInter <- c(5,9)
ParaPlotcat <- c(6,6)
BootstrapFigureErrors <- file(paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/BootstrapFigureErrors.txt"), open="wt")
sink(BootstrapFigureErrors, type="message")
for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  mod <- Allmod[i]
  seuil <- Allseuil[i]
  PInter <- ParaInter[i]
  PPlotcat <- ParaPlotcat[i]
  Dir = (paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP")) # Directory 
  setwd(Dir)
  #list.files(Dir,pattern = paste0(seuil,".rds"))
  dfplot <- readRDS(paste0("dfplot",CODE,seuil,".rds")) #Base de données plot
  dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
  dfplotbis <- readRDS(paste0("dfplotbis",CODE,seuil,".rds")) #Base de données plot
  # This one is the destination where only the model of interest is (binom and negbin)
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/",mod,"/"))
  setwd(Dir)
  assign(eval(mod),get(load(file = paste0(mod,".rda"))))
  rm("x")
  Extraction(get(mod))
  try(ModelBoot.LE(get(mod),PInter,PPlotcat,LvL=30,CAT=CAT,nBoot=10,Yportion = 0.66,saveboot=T,nCoeur=13),silent=T)
}
close(BootstrapFigureErrors)


#### Reverse -RE 

i = 1
Allcode <- c("QUEPYR","QUEPYR","QUEROB")
Allmod <- c("MnbZT13B.27","MnbZT13B.27","MnbZT7A.23")
Allseuil <- c(0.7,0.7,0.8)
ParaInter <- c(8,4,9)
ParaPlotcat <- c(5,5,10)
BootstrapFigureErrors <- file(paste0("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/BootstrapFigureErrors.txt"), open="wt")
sink(BootstrapFigureErrors, type="message")
for (i in 1:length(Allcode)){
  CODE <- Allcode[i]
  mod <- Allmod[i]
  seuil <- Allseuil[i]
  PInter <- ParaInter[i]
  PPlotcat <- ParaPlotcat[i]
  Dir = (paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP")) # Directory 
  setwd(Dir)
  #list.files(Dir,pattern = paste0(seuil,".rds"))
  dfplot <- readRDS(paste0("dfplot",CODE,seuil,".rds")) #Base de données plot
  dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
  dfplotbis <- readRDS(paste0("dfplotbis",CODE,seuil,".rds")) #Base de données plot
  # This one is the destination where only the model of interest is (binom and negbin)
  Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/",mod,"/"))
  setwd(Dir)
  assign(eval(mod),get(load(file = paste0(mod,".rda"))))
  rm("x")
  Extraction(get(mod))
  try(ModelBoot.r.LE(get(mod),PInter,PPlotcat,LvL=30,CAT=CAT,nBoot=10,Yportion = 0.66,saveboot=T,nCoeur=13),silent=T)
}
close(BootstrapFigureErrors)




############### Now use cowplot to create the figure I want #######################
legend <- get_legend(pABIALB)


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





plot(pABIALB_Mbin14A.19_0.7_mean_spei12_Plotcat)+theme(legend.position = "none")