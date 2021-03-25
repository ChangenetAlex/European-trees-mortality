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


###################################################################################
###                                                                            ####
###                                                                            #### # Here, 6 scripts are run
###     PART two : how to extract the information we want for all species      ####
###                                                                            ####
###                                                                            ####
###################################################################################
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

CODE = "ACEPSE" 
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

seuil = 0.55
#Dir = c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models")) # Directory 
Dir = (paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP")) # Directory 
setwd(Dir)
list.files(Dir,pattern = paste0(seuil,".rds"))

dfplot <- readRDS(paste0("dfplot",CODE,seuil,".rds")) #Base de données plot
#dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
dfplotbis <- readRDS(paste0("dfplotbis",CODE,seuil,".rds")) #Base de données plot
#dir.create(path=paste0(Dir,"/Models"))
#dir.create(path=paste0(Dir,"/Models/Negbin"))
#dir.create(path=paste0(Dir,"/Models/binomial"))
Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models"))
# This one is the destination where only the model of interest is (binom and negbin)
Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/"))
setwd(Dir)

#summary(as.factor(dfplot2$Plotcat))
#max(dfplot2[dfplot2$Plotcat=="1","latitude"])
#max(dfplot2[dfplot2$Plotcat=="2","latitude"])

#min(dfplot2[dfplot2$Plotcat=="1","latitude"])
#min(dfplot2[dfplot2$Plotcat=="2","latitude"])

#mean(dfplot2[dfplot2$Plotcat=="1","latitude"])
#mean(dfplot2[dfplot2$Plotcat=="2","latitude"])

#nrow(dfplot2[dfplot2$Plotcat=="1",])
#nrow(dfplot2[dfplot2$Plotcat=="2",])


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
             "sqrtBA.ha.plot.1","logBAj.plot.1","logdbh.plot.mean","sqrtBAIj.plot.1.mean","sqrtBAIj.plot.1","treeNbr","yearsbetweensurveys","Plotcat")
Premodel(z=dfplot2,Resp=Resp,Explain=Explain,size=3,save=T) # Obtain the global VIF for all variable (not transformed)


#########################
#####  Models Full  #####
#########################
# This one is the destination where all models are
#Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/binomial"))
#Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/Negbin"))


# This one is the destination where only the model of interest is (binom and negbin)
Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/"))


setwd(Dir)
### Load preexisting models 
list.dirs(getwd())
setwd(list.dirs(getwd())[2]) # Chose the model file
list.files()
assign(sub(".rda","",list.files()[12], ignore.case = FALSE,fixed = T),get(load(file = list.files()[12]))) 
rm("x")
summary(get(sub(".rda","",list.files()[12])))




##########################
#### Coef + valid ########
##########################
Diagnostic(Mbin13B,0.66,F) # Pbm with the coreelog function 
ExtracTest(Mbin13B) # Me donne la liste des paramètres de mon modèles. 
for (i in c(A[c(1:11)])){
  Effect_coef(Mbin14A,i)}     # Para que l'on veut regarder  parmis ceux citer  
Effect_summary(Mbin7A,"biotic") #
ggEffect(Mbin13B,'REL',"sum",band=T) # Only two that i need to apply again
ggEffect(Mbin13B,'ABS',"sum",band=T) # # Only two that i need to apply again

###########################
#####    Affined     ######
###########################
require("mgcv")
Mymod <- "Mbin7A"
num <- ""

num = num + 1
assign(paste0(Mymod,num),fitme(sp.mort.bin ~ sqrtBAIj.plot.1 + sqrtBAIj.plot.1.mean + logdbh.plot.mean + 
                                 treeNbr + yearsbetweensurveys + bio5_climate_mean.30 + bio13_climate_mean.30 + 
                                 min_spei12 + mean_spei12 + sqrtBA.ha.plot.1 + sqrtBA.O.plot.1 + 
                                 Plotcat + I(bio5_climate_mean.30^2) + I(bio13_climate_mean.30^2) + 
                                 I(min_spei12^2) + I(mean_spei12^2) + bio5_climate_mean.30:bio13_climate_mean.30 + 
                                 bio5_climate_mean.30:min_spei12 + bio5_climate_mean.30:mean_spei12 + 
                                 bio5_climate_mean.30:sqrtBA.ha.plot.1 + bio5_climate_mean.30:sqrtBA.O.plot.1 + 
                                 bio13_climate_mean.30:min_spei12 + bio13_climate_mean.30:mean_spei12 + 
                                 bio13_climate_mean.30:sqrtBA.ha.plot.1 + bio13_climate_mean.30:sqrtBA.O.plot.1 + 
                                 min_spei12:mean_spei12 + min_spei12:sqrtBA.ha.plot.1 + min_spei12:sqrtBA.O.plot.1 + 
                                 mean_spei12:sqrtBA.ha.plot.1 + mean_spei12:sqrtBA.O.plot.1 + 
                                 sqrtBA.ha.plot.1:sqrtBA.O.plot.1 + +bio5_climate_mean.30:Plotcat + 
                                 bio13_climate_mean.30:Plotcat + min_spei12:Plotcat + mean_spei12:Plotcat + 
                                 sqrtBA.ha.plot.1:Plotcat + sqrtBA.O.plot.1:Plotcat + (1 |  country), data=dfplot2,family=binomial,method='REML'))
Saving(get(paste0(Mymod,num)))
MyMod <- get(paste0(Mymod,num))
MyMod <- as.data.frame(summary(MyMod)$beta_table)
MyMod2 <- MyMod[order(abs(MyMod$'t-value')),]
MyMod2
MyMod


###########################
#####   Bootstraps   ######
###########################
setwd(Dir) # Set dir in the new affined
### Load preexisting models 
list.dirs(getwd())
setwd(list.dirs(getwd())[3]) # Chose the model file
list.files()
assign(sub(".rda","",list.files()[8], ignore.case = FALSE,fixed = T),get(load(file = list.files()[8]))) 
rm("x")
summary(get(sub(".rda","",list.files()[8])))

Diagnostic(Mbin13B.26,0.66,F)
Extraction(Mbin13B.26)
summary(Mbin13B.26)
ModelBoot(Mbin13B.26,4,10,LvL=30,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=13)
ModelBoot(Mbin13B.26,9,10,LvL=30,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=13)
ModelBoot(Mbin13B.26,9,11,LvL=30,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=13)





#######################
##### AC Spatial ######  # Try to add spatial AC if possible
#######################
spaMM.options(separation_max=6363)
Mbin3B.19.AC <- fitme(sp.mort.bin ~ sqrtBAIj.plot.1 + sqrtBAIj.plot.1.mean + logdbh.plot.mean + 
                         treeNbr + sqrtBA.ha.plot.1 + sqrtBA.O.plot.1 + Plotcat + 
                         I(bio14_climate_mean.30^2) + I(min_spei12^2) + tmean.djf_climate_mean.30:bio14_climate_mean.30 + 
                         tmean.djf_climate_mean.30:min_spei12 + tmean.djf_climate_mean.30:sqrtBA.ha.plot.1 + 
                         min_spei12:sqrtBA.O.plot.1 + mean_spei12:sqrtBA.ha.plot.1 + 
                         mean_spei12:sqrtBA.O.plot.1 + sqrtBA.ha.plot.1:sqrtBA.O.plot.1 + 
                         min_spei12:Plotcat + mean_spei12:Plotcat + (1 | country) + Matern(1|latitude + longitude),
                         data=dfplot2, family = binomial,method='REML')


###########################
##### ZT to affined  ###### # Same process with the ZT model 
###########################
Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/Negbin")) # All models 
# This one is the destination where only the model of interest is (binom and negbin)
Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/"))

setwd(Dir)
### Load preexisting models 
list.dirs(getwd())
setwd(list.dirs(getwd())[5]) # Chose the model file
list.files()
assign(sub(".rda","",list.files()[32], ignore.case = FALSE,fixed = T),get(load(file = list.files()[32]))) 
rm("x")
summary(get(sub(".rda","",list.files()[32])))
detach("package:mgcv", unload=TRUE)

Mymod <- "MnbZT7A"
num <- ""

num = num + 1
assign(paste0(Mymod,num),fitme(sp.mortality.plot.rate.yr ~ sqrtBAIj.plot.1 + sqrtBAIj.plot.1.mean + 
                                 logdbh.plot.mean + treeNbr + yearsbetweensurveys + bio5_climate_mean.30 + 
                                 bio13_climate_mean.30 + min_spei12 + mean_spei12 + sqrtBA.ha.plot.1 + 
                                 sqrtBA.O.plot.1 + Plotcat + I(bio5_climate_mean.30^2) + I(bio13_climate_mean.30^2) + 
                                 I(min_spei12^2) + I(mean_spei12^2) + bio5_climate_mean.30:bio13_climate_mean.30 + 
                                 bio5_climate_mean.30:min_spei12 + bio5_climate_mean.30:mean_spei12 + 
                                 bio5_climate_mean.30:sqrtBA.ha.plot.1 + bio5_climate_mean.30:sqrtBA.O.plot.1 + 
                                 bio13_climate_mean.30:min_spei12 + bio13_climate_mean.30:mean_spei12 + 
                                 bio13_climate_mean.30:sqrtBA.ha.plot.1 + bio13_climate_mean.30:sqrtBA.O.plot.1 + 
                                 min_spei12:mean_spei12 + min_spei12:sqrtBA.ha.plot.1 + min_spei12:sqrtBA.O.plot.1 + 
                                 mean_spei12:sqrtBA.ha.plot.1 + mean_spei12:sqrtBA.O.plot.1 + 
                                 sqrtBA.ha.plot.1:sqrtBA.O.plot.1 + +bio5_climate_mean.30:Plotcat + 
                                 bio13_climate_mean.30:Plotcat + min_spei12:Plotcat + mean_spei12:Plotcat + 
                                 sqrtBA.ha.plot.1:Plotcat + sqrtBA.O.plot.1:Plotcat + (1 |country),
                                 data=subset(dfplot2,sp.mortality.plot.rate.yr>0),family=negbin(),method='REML'))

Saving(get(paste0(Mymod,num)))
MyMod <- get(paste0(Mymod,num))
MyMod <- as.data.frame(summary(MyMod)$beta_table)
MyMod2 <- MyMod[order(abs(MyMod$'t-value')),]
MyMod2
MyMod

###########################
#####   Bootstraps   ######  ## ZT
###########################
Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/Negbin")) # All models 
# This one is the destination where only the model of interest is (binom and negbin)
Dir =c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/"))
setwd(Dir) # Set dir in the new affined
### Load preexisting models 
list.dirs(getwd())
setwd(list.dirs(getwd())[104]) # Chose the model file
list.files()
assign(sub(".rda","",list.files()[8], ignore.case = FALSE,fixed = T),get(load(file = list.files()[8]))) 
rm("x")
summary(get(sub(".rda","",list.files()[2])))
Diagnostic(MnbZT13B.23,0.66,F)
Extraction(MnbZT13B.23)
summary(MnbZT13B.23)
ModelBoot(MnbZT13B.23,3,6,LvL=30,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=13)
ModelBoot(MnbZT13B.23,5,6,LvL=30,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=13)
ModelBoot(MnbZT13B.23,9,6,LvL=30,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=13)
ModelBoot(MnbZT7A.23,9,10,LvL=30,CAT=CAT,nBoot=1000,Yportion = 0.66,saveboot=T,nCoeur=13)


#######################
##### AC Spatial ######  # Try to add spatial AC if possible
#######################
MnbZT3C.24.AC <- fitme(sp.mortality.plot.rate.yr ~ sqrtBAIj.plot.1.mean + logdbh.plot.mean + 
                           yearsbetweensurveys + bio12_climate_mean.30 + min_spei12 + 
                           sqrtBA.O.plot.1 + logBAj.plot.1 + I(min_spei12^2) + bio5_climate_mean.30:sqrtBA.O.plot.1 + 
                           bio12_climate_mean.30:min_spei12 + bio12_climate_mean.30:mean_spei12 + 
                           min_spei12:sqrtBA.O.plot.1 + mean_spei12:sqrtBA.O.plot.1 + 
                           bio5_climate_mean.30:Plotcat + (1 | country) + Matern(1|latitude + longitude),
                        data=subset(dfplot2,sp.mortality.plot.rate.yr>0),family=negbin(),method='REML')
Saving(MnbZT3C.24.AC)
Diagnostic(M2nbZT13B.29.AC,0.66,F)   # A voir si on le fait because too long 
