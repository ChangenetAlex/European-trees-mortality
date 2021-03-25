#########################################
#####   Work and Process data     #######
#########################################




#bioclim=paste("bio",c(1,12,13,14,2,5,6),"_climate_",sep="") #
#T.seas=paste0("tmean.",c("djf","jja","mam","son"),"_climate_")#
#P.seas=paste0("prec.",c("djf","jja","mam","son"),"_climate_")#
#pet=paste0("pet.",c("max","mean","min"),"_climate_") #
#ppet=paste0("ppet.",c("max","mean","min"),"_climate_")#
#eumedclim.vars=(c(bioclim, pet, ppet, P.seas, T.seas))
#Myvariable <- paste0(rep(eumedclim.vars,length(Interv)*length(metric)),rep(metric,length(Interv)*length(eumedclim.vars)),rep(Interv,length(eumedclim.vars)*length(metric)))


############################################
####                                    ####
#### Fonction to generate all the model ####
####       I want to evaluate :         ####
####                                    ####
############################################
metric = c("mean.","min.","max.") #mean
Interv = c(30) #5 10 15 30

Allvariable.axe1 = c(paste0("bio1_climate_mean.",Interv), #Axe 1
                     paste0("tmean.djf_climate_mean.",Interv),
                     paste0("bio5_climate_mean.",Interv),
                     paste0("tmean.jja_climate_mean.",Interv))
Allvariable.axe2 = c(paste0("bio12_climate_mean.",Interv), #Axe 2
                     paste0("bio13_climate_mean.",Interv), #Axe 2
                     paste0("ppet.mean_climate_mean.",Interv),#Axe 2
                     paste0("bio14_climate_mean.",Interv))#Axe2

#Fonction pour écrire mes modèles : 
expand.grid(Allvariable.axe1,Allvariable.axe2)
n <- expand.grid(Allvariable.axe1,Allvariable.axe2,stringsAsFactors = F)
n <- n[n$Var1!=n$Var2,]

VarTrans <- c("sqrtBA.ha.plot.1","sqrtBA.O.plot.1","logBAj.plot.1")
n2 <- t(combn(VarTrans,2))

i=1
j=1
Nom <- LETTERS[seq(from = 1, to = nrow(n2))]

vardep = c("sp.mortality.plot.rate.yr","sp.mort.bin") #keep this one 
varcat = c("Plotcat")
random <- c("country") #This will be kept and modified
for (i in 1:nrow(n)){
  for (j in 1:nrow(n2)){
    varAll = c(n[i,1],n[i,2],"min_spei12","mean_spei12",n2[j,1],n2[j,2])
    vardep1 <- paste0("try(Mbin",i,Nom[j]," <- eval_fork(fitme(",vardep[2]," ~ ",collapse = "")
    vardep11 <- paste0("try(M2bin",i,Nom[j]," <- eval_fork(fitme(",vardep[2]," ~ ",collapse = "")
    vardep2 <- paste0("try(MnbZT",i,Nom[j]," <- eval_fork(fitme(",vardep[1]," ~ ",collapse = "")
    vardep22 <- paste0("try(M2nbZT",i,Nom[j]," <- eval_fork(fitme(",vardep[1]," ~ ",collapse = "")
    varAll1 <- paste0(varAll,collapse=" + ")
    varcat1 <- paste0(" + ",varcat,collapse="")
    carre <- paste0(" + I(",varAll[1:4],"^2)",collapse = "")
    interac1 <- t(combn(varAll,2))
    interac1 = paste0(" + ",interac1[,1],":",interac1[,2],collapse="")
    interac2 = paste0(" + ",expand.grid(varAll,varcat)[,1],":",expand.grid(varAll,varcat)[,2],collapse="")
    random1 <- paste0(" + (1|",random,")",collapse = "")
    finfun <- (", data=dfplot2, family = binomial,method='REML'),timeout = 1200),silent=T)")
    finfun1 <- (", data=subset(dfplot2,sp.mortality.plot.rate.yr>0),family=negbin(),method='REML'),timeout = 1200),silent=T)")
    
    capture.output(cat(paste0("#Mbin",i,Nom[j]," = ",varAll[1]," ~ ",varAll[2])," + ",n2[j,1]," + ",n2[j,2],"\n",
                       paste0("Model <- deparse(substitute(Mbin",i,Nom[j],"))"),"\n",
                       paste0(vardep1,"sqrtBAIj.plot.1 + sqrtBAIj.plot.1.mean + logdbh.plot.mean + treeNbr + yearsbetweensurveys + ",varAll1," + ",varcat,carre,interac1," + ",interac2,random1,finfun,collapse = ""),"\n",sep=""),
                   file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
    
    capture.output(cat(paste0("try(Saving(Mbin",i,Nom[j],"),silent=T)"),"\n"),
                   file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
    
    capture.output(cat(paste0("rm(list = c('Mbin",i,Nom[j],"'))"),"\n","gc()","\n"),
                   file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
    
    capture.output(cat(paste0("#M2bin",i,Nom[j]," = ",varAll[1]," ~ ",varAll[2])," + ",n2[j,1]," + ",n2[j,2],"\n",
                       paste0("Model <- deparse(substitute(M2bin",i,Nom[j],"))"),"\n",
                       paste0(vardep11,"sqrtBAIj.plot.1 + sqrtBAIj.plot.1.mean + sqrtdbh.plot.mean + treeNbr + yearsbetweensurveys + ",varAll1," + ",varcat,carre,interac1," + ",interac2,random1,finfun,collapse = ""),"\n",sep=""),
                   file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
    
    capture.output(cat(paste0("try(Saving(M2bin",i,Nom[j],"),silent=T)"),"\n"),
                   file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
    
    capture.output(cat(paste0("rm(list = c('M2bin",i,Nom[j],"'))"),"\n","gc()","\n"),
                   file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
    
    capture.output(cat(paste0("#MnbZT",i,Nom[j],"nb = ",varAll[1]," ~ ",varAll[2])," + ",n2[j,1]," + ",n2[j,2],"\n",
                       paste0("Model <- deparse(substitute(MnbZT",i,Nom[j],"))"),"\n",
                       paste0(vardep2,"sqrtBAIj.plot.1 + sqrtBAIj.plot.1.mean + logdbh.plot.mean + treeNbr + yearsbetweensurveys + ",varAll1," + ",varcat,carre,interac1," + ",interac2,random1,finfun1,collapse = ""),"\n",sep=""),
                   file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
    
    capture.output(cat(paste0("try(Saving(MnbZT",i,Nom[j],"),silent=T)"),"\n"),
                   file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
    
    capture.output(cat(paste0("rm(list = c('MnbZT",i,Nom[j],"'))"),"\n","gc()","\n"),
                   file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
    
    capture.output(cat(paste0("#M2nbZT",i,Nom[j],"nb = ",varAll[1]," ~ ",varAll[2])," + ",n2[j,1]," + ",n2[j,2],"\n",
                       paste0("Model <- deparse(substitute(M2nbZT",i,Nom[j],"))"),"\n",
                       paste0(vardep22,"sqrtBAIj.plot.1 + sqrtBAIj.plot.1.mean + sqrtdbh.plot.mean + treeNbr + yearsbetweensurveys + ",varAll1," + ",varcat,carre,interac1," + ",interac2,random1,finfun1,collapse = ""),"\n",sep=""),
                   file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
    
    capture.output(cat(paste0("try(Saving(M2nbZT",i,Nom[j],"),silent=T)"),"\n"),
                   file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
    
    capture.output(cat(paste0("rm(list = c('M2nbZT",i,Nom[j],"'))"),"\n","gc()","\n"),
                   file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
    
    }}
capture.output(cat("}} \n close(Errors.files) \n sink(type='message')"),file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Myscripts/Fundiv.project/temporaire.R",append = T,type="output",split=F)
