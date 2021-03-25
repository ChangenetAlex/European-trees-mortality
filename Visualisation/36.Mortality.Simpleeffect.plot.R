rm(list = ls())
gc()
### 24/11/2020
## Simple effects in mortality models 

## Need to load the model first and check if we can make predictions 
### 24/11/2020
## Simple effects in mortality models 
library(spaMM)
library(parallel)
i = 1
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
Allmod <- c("Mbin13A.18","Mbin14A.19","Mbin13B.26","Mbin13A.22","M2bin13B.22","M2bin13A21","M2bin15B23","Mbin3C.26","M2bin7B.17","Mbin11B.19","Mbin3B.20","Mbin15B.21","Mbin5B.17","Mbin3B.31","Mbin13A.27","M2bin1C.20","M2bin1B.23","Mbin13B.26","Mbin7A.26","Mbin1B.23")


## Here we load all the model in BIN ! 
Mortmod.all.bin <- mapply(function(x,y,z){
  Dir =c(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",x,"/CLIMAP/Models/binomial/",y,"/"))
  setwd(Dir)
  x <- get(load(file = paste0(y,".rda")))
},x=Allcode,y=Allmod,SIMPLIFY = F)


## Idem ZT 
i = 1
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
Allmod <- c("MnbZT13A.19","MnbZT14A.20","MnbZT13B.23","MnbZT13A.21","M2nbZT13B.29","M2nbZT13A.22","M2nbZT15B.24","MnbZT3C.24","M2nbZT7B.26","MnbZT11B.22","MnbZT3B.27","MnbZT15B.28","MnbZT5B.21","MnbZT13A.33","M2nbZT1C.27","M2nbZT1B.24","MnbZT13B.27","MnbZT7A.23","MnbZT1B.27")
Mortmod.all.zt <- mapply(function(x,y,z){
  Dir =c(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/",x,"/CLIMAP/Models/Negbin/",y,"/"))
  setwd(Dir)
  x <- get(load(file = paste0(y,".rda")))
},x=Allcode,y=Allmod,SIMPLIFY = F)

## Here we can extract all the simple effect
Extraction <- function(x){
  A <- paste0(getCall(x)[2])
  A <- unlist(strsplit(A, "~",fixed = T))[2]
  A <- unlist(strsplit(A, "+",fixed = T))
  A <- grep(A,pattern = "|",fixed=T,value=T,invert=T)
  A <- unlist(strsplit(A, ":",fixed = T))
  #A <- grep(A,pattern = ":",fixed=T,value=T,invert=T)
  A <- grep(A,pattern = "e)",fixed=T,value=T,invert=T) #Remove the AC term which remains
  A <- grep(A,pattern = "^[ ]$",value=T,invert=T) #Remove empty characters
  A <- sub("I(","", A, ignore.case = FALSE,fixed = T)
  A <- sub("^2)","", A, ignore.case = FALSE,fixed = T)
  A <- sub("offset(log(","", A, ignore.case = FALSE,fixed = T)
  A <- sub("))","", A, ignore.case = FALSE,fixed = T)
  A <- sub("\n ","", A, ignore.case = FALSE,fixed = T) # New line
  A <- sub(" ","", A, ignore.case = FALSE,fixed = T)
  A <- sub(" ","", A, ignore.case = FALSE,fixed = T) #Twice
  A <- sub(" ","", A, ignore.case = FALSE,fixed = T) #Twice
  A <- sub(" ","", A, ignore.case = FALSE,fixed = T) #Twice
  A <- unique(A)
  A
}

MyEffectbin <- lapply(Mortmod.all.bin,Extraction)
MyEffectzt <- lapply(Mortmod.all.zt,Extraction)

MyEffectbin.unique <- unique(unlist(MyEffectbin))
MyEffectzt.unique <- unique(unlist(MyEffectzt))
            
LvL <- 30  

# x <- Mortmod.all.zt$ABIALB
# y <- names(Mortmod.all.zt)[3]  
# i <- 3

for (i in 1:length(MyEffectzt.unique)){
  VARIABLE <- MyEffectzt.unique[i]
Species.effect.ZT  <- mcmapply(function(x,y){
    nameVAR <- Extraction(x)
    if (VARIABLE%in%nameVAR & VARIABLE!="Plotcat"){
      ## Simulation 
      data <- x$data[,colnames(x$data)%in%nameVAR[1:length(nameVAR)]] # Extract by colnames and in the right order
      data <- data[,c(nameVAR[1:length(nameVAR)])] # Reorder my data Need in all cases !!!
      # Creer table avec nouvelles valeurs 
      myDF = as.data.frame(matrix(ncol = length(nameVAR), nrow = LvL, NA)) # Df vierge
      colnames(myDF) <- colnames(data) # Bon
      ## We create a matrix of 30 values that we will predict by variable + 30 by level of plotcat. 

      #####
      # 1 # Remplir tout à la moyenne et tous les Plotcat à 0 par défaut. Et séquence sur le range 
      #####
      myDF[,nameVAR[nameVAR%in%colnames(data) & nameVAR!="Plotcat"]] <- apply(data[,nameVAR[nameVAR%in%colnames(data) & nameVAR!="Plotcat"]],2,
                                                                              function(y) as.numeric(rep(mean(y),nrow(myDF))))
      myDF[,"Plotcat"] <- as.factor(rep("0",nrow(myDF)))
      myDF[,VARIABLE] <- seq(min(data[,VARIABLE]),max(data[,VARIABLE]),length=LvL)
      #####
      # 2 # Prediction 
      #####
      sub_x = update(x, data=x$data) # hack because the formula is absent in the old model (so refit to take the formula for the predictions)
      #predict(sub_x,newdata=numeric(0)) ## TO RUN IF THERE IS A BUG
      pred.boot=predict(sub_x, newdata=myDF,re.form=NA,type="response",binding="fit",list(predVar=TRUE)) # Predict the values for each cases describe in MyDf
      pred.boot[,"species"] <- y
      pred.boot <- cbind(pred.boot,get_intervals(sub_x, newdata=myDF,re.form=NA,intervals="predVar",level = c(0.90),variances=list(linPred=TRUE,disp=FALSE,predVar=T,residVar=FALSE,cov=TRUE)))
      Var <- get_predVar(sub_x, newdata=myDF,re.form=NA,level = c(0.95),which="predVar",variances=list(predVar=TRUE,linPred=TRUE,cov=FALSE))
      pred.boot[,"SEneg"] <- pred.boot[,"fit"]-Var
      pred.boot[,"SEpos"] <- pred.boot[,"fit"]+Var
      return(pred.boot[,c(VARIABLE,"fit","species","predVar_0.05","predVar_0.95","SEneg","SEpos")])
      }else {return(NULL)}
},x=Mortmod.all.zt,y=names(Mortmod.all.zt),SIMPLIFY = F,mc.cores=20,mc.silent=F)
Species.effect.ZT <- do.call(rbind,Species.effect.ZT)
saveRDS(get("Species.effect.ZT"), paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/",VARIABLE,".zt.rds"))
assign(paste0(VARIABLE,".zt"),Species.effect.ZT,envir = .GlobalEnv)}



#### Idem avec les modèles truncated
for (i in 1:length(MyEffectbin.unique)){
  VARIABLE <- MyEffectbin.unique[i]
  Species.effect.bin  <- mcmapply(function(x,y){
    nameVAR <- Extraction(x)
    if (VARIABLE%in%nameVAR & VARIABLE!="Plotcat"){
      ## Simulation 
      data <- x$data[,colnames(x$data)%in%nameVAR[1:length(nameVAR)]] # Extract by colnames and in the right order
      data <- data[,c(nameVAR[1:length(nameVAR)])] # Reorder my data Need in all cases !!!
      # Creer table avec nouvelles valeurs 
      myDF = as.data.frame(matrix(ncol = length(nameVAR), nrow = LvL, NA)) # Df vierge
      colnames(myDF) <- colnames(data) # Bon
      ## We create a matrix of 30 values that we will predict by variable + 30 by level of plotcat. 
      
      #####
      # 1 # Remplir tout à la moyenne et tous les Plotcat à 0 par défaut. Et séquence sur le range 
      #####
      myDF[,nameVAR[nameVAR%in%colnames(data) & nameVAR!="Plotcat"]] <- apply(data[,nameVAR[nameVAR%in%colnames(data) & nameVAR!="Plotcat"]],2,
                                                                              function(y) as.numeric(rep(mean(y),nrow(myDF))))
      myDF[,"Plotcat"] <- as.factor(rep("0",nrow(myDF)))
      myDF[,VARIABLE] <- seq(min(data[,VARIABLE]),max(data[,VARIABLE]),length=LvL)
      #####
      # 2 # Prediction 
      #####
      sub_x = update(x, data=x$data) # hack because the formula is absent in the old model (so refit to take the formula for the predictions)
      #predict(sub_x,newdata=numeric(0)) ## TO RUN IF THERE IS A BUG
      pred.boot=predict(sub_x, newdata=myDF,re.form=NA,type="response",binding="fit") # Predict the values for each cases describe in MyDf
      pred.boot[,"species"] <- y
      pred.boot <- cbind(pred.boot,get_intervals(sub_x, newdata=myDF,re.form=NA,intervals="predVar",level = c(0.90),variances=list(linPred=TRUE,disp=FALSE,predVar=T,residVar=FALSE,cov=TRUE)))
      Var <- get_predVar(sub_x, newdata=myDF,re.form=NA,level = c(0.95),which="predVar",variances=list(predVar=TRUE,linPred=TRUE,cov=FALSE))
      pred.boot[,"SEneg"] <- pred.boot[,"fit"]-Var
      pred.boot[,"SEpos"] <- pred.boot[,"fit"]+Var
      return(pred.boot[,c(VARIABLE,"fit","species","predVar_0.05","predVar_0.95","SEneg","SEpos")])
    }else {return(NULL)}
  },x=Mortmod.all.bin,y=names(Mortmod.all.bin),SIMPLIFY = F,mc.cores=20,mc.silent=F)
  Species.effect.bin <- do.call(rbind,Species.effect.bin)
  saveRDS(get("Species.effect.bin"), paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/",VARIABLE,".bin.rds"))
  assign(paste0(VARIABLE,".bin"),Species.effect.bin,envir = .GlobalEnv)}





