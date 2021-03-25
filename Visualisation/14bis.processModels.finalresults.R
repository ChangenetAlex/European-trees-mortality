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
Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
# Load the function 
source(paste0(Dir,"Myscripts/Fundiv.project/function1.R.squared.r"))
source(paste0(Dir,"Myscripts/Fundiv.project/function2.Saving.R"))
source(paste0(Dir,"Myscripts/Fundiv.project/function3.ModelBoot.R"))
source(paste0(Dir,"Myscripts/Fundiv.project/function4.Diagnostic.R"))
source(paste0(Dir,"Myscripts/Fundiv.project/function5.Premodel.R"))
source(paste0(Dir,"Myscripts/Fundiv.project/function6.Coefs.R"))

#Extraction of the database of the wanted species : 
CODE = "BETPEN" 
#"ALNGLU" 
#"ABIALB"
#"BETPEN"
#"PICABI"
#"PINPINA"
#"FAGSYL"
#"PINHAL"
#"QUEROB"
#"QUEILE"
#"PINNIG"
#"QUEPET"
#"CASSAT"
#"ABIALB"
#"QUEPUB"
#"QUEPYR"
#"FRAEXC"
#"PINPIN"
#"QUESUB"
#"BETPEN"
#"ALNGLU"
#"POPTRE"
#"ACEPSE"
#"LARDEC"
#"POPNIG"

# Load the df
seuil = 0.8
#Dir = c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models")) # Directory 
Dir = (paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP")) # Directory 
setwd(Dir)
list.files(Dir,pattern = paste0(seuil,".rds"))
dfplot <- readRDS(paste0("dfplot",CODE,seuil,".rds")) #Base de données plot
#dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
dfplotbis <- readRDS(paste0("dfplotbis",CODE,seuil,".rds")) #Base de données plot
# This one is the destination where only the model of interest is (binom and negbin)
Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/"))
setwd(Dir)


###########################
## Premodel analysis ###### 
###########################
# Perform the anaysis on all variable three by three
Resp <- c("sp.mortality.plot.count.yr")
Resp <- c("sp.mort.bin")
Explain <- c("BA.ha.plot.1","BA.O.plot.1","BAj.plot.1","dbh.plot.mean","BAIj.plot.1.mean","BAIj.plot.1",
             "logBA.ha.plot.1","logBA.O.plot.1","logBAj.plot.1","logdbh.plot.mean","logBAIj.plot.1.mean","logBAIj.plot.1",
             "sqrtBA.ha.plot.1","sqrtBA.O.plot.1","sqrtBAj.plot.1","sqrtdbh.plot.mean","sqrtBAIj.plot.1.mean","sqrtBAIj.plot.1",
             "treeNbr","yearsbetweensurveys","bio14_climate_mean.30","bio1_climate_mean.30","min_spei12","mean_spei12","Plotcat")

for (i in 1:((length(Explain)-1)/3)){
  Explain <- c("BA.ha.plot.1","BA.O.plot.1","BAj.plot.1","dbh.plot.mean","BAIj.plot.1.mean","BAIj.plot.1",
               "logBA.ha.plot.1","logBA.O.plot.1","logBAj.plot.1","logdbh.plot.mean","logBAIj.plot.1.mean","logBAIj.plot.1",
               "sqrtBA.ha.plot.1","sqrtBA.O.plot.1","sqrtBAj.plot.1","sqrtdbh.plot.mean","sqrtBAIj.plot.1.mean","sqrtBAIj.plot.1",
               "treeNbr","yearsbetweensurveys","bio14_climate_mean.30","bio1_climate_mean.30","min_spei12","mean_spei12","Plotcat")
  Explain <- Explain[c((3*i-2):(3*i),length(Explain))]
  print(Explain)
  Premodel(z=dfplot2,Resp=Resp,Explain=Explain,size=2,save=T)   # Apply the functictn on both database on all trasnformed and no trasnformed data
} 

# Before this step we need to consider write the variables that are set in the selected model.
# Premodel analysis for 12 variables and to obtain the global VIF indices
Explain <- c("bio14_climate_mean.30","bio1_climate_mean.30","min_spei12","mean_spei12",
             "sqrtBA.ha.plot.1","logBAj.plot.1","sqrtdbh.plot.mean","sqrtBAIj.plot.1.mean","sqrtBAIj.plot.1","treeNbr","yearsbetweensurveys","Plotcat")
Premodel(z=dfplot2,Resp=Resp,Explain=Explain,size=3,save=T) # Obtain the global VIF for all variable (not transformed)




Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/"))
setwd(Dir)
### Load preexisting models 
list.dirs(getwd())
setwd(list.dirs(getwd())[2]) # Chose the model file
list.files()
assign(sub(".rda","",list.files()[18], ignore.case = FALSE,fixed = T),get(load(file = list.files()[18]))) 
rm("x")
summary(get(sub(".rda","",list.files()[12])))


Diagnostic(Mbin13A,0.66,F) # Pbm with the coreelog function 
ggEffect(Mbin13A,'REL',"sum",band=T) # Only two that i need to apply again
ggEffect(Mbin13A,'ABS',"sum",band=T) # # Only two that i need to apply again
require("mgcv")
###########################
#####   Bootstraps   ######
###########################
setwd(Dir) # Set dir in the new affined
list.dirs(getwd())
setwd(list.dirs(getwd())[3]) # Chose the model file
list.files()
assign(sub(".rda","",list.files()[18], ignore.case = FALSE,fixed = T),get(load(file = list.files()[18]))) 
rm("x")
Diagnostic(Mbin13A.22,0.66,F)
summary(Mbin13A.22)
Extraction(Mbin13A.22)
ModelBoot(Mbin13A.22,7,6,LvL=30,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=13)
ModelBoot(Mbin13A.22,11,6,LvL=30,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=13)

# ZT 
Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/"))
setwd(Dir)
list.dirs(getwd())
setwd(list.dirs(getwd())[6]) # Chose the model file ZT
list.files()
assign(sub(".rda","",list.files()[28], ignore.case = FALSE,fixed = T),get(load(file = list.files()[28]))) 
rm("x")
detach("package:mgcv", unload=TRUE)
Diagnostic(MnbZT13A.21,0.66,F)
summary(MnbZT13A.21)
Extraction(MnbZT13A.21)
ModelBoot(MnbZT13A.21,4,7,LvL=30,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=13)
ModelBoot(MnbZT13A.21,9,7,LvL=30,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=13)
ModelBoot(MnbZT13A.21,10,7,LvL=30,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=13)

