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
  #A <- grep(A,pattern = "|",fixed=T,value=T,invert=T)
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
  A <- gsub("^.*country.*$","country", A)
  A <- unique(A)
  A
}


ExtractionRE <- function(x){
  A <- paste0(getCall(x)[2])
  A <- unlist(strsplit(A, "~",fixed = T))[2]
  A <- unlist(strsplit(A, "+",fixed = T))
  A <- grep(A,pattern = "country",fixed=T,value=T,invert=F)
  A <- grepl(A,pattern = "country")
  A
}

# x <- Mortmod.all.bin$ABIALB
# # x <- Mortmod.all.zt$ABIALB
# y <- names(Mortmod.all.bin)[2]  
# #i <- 3




Species.RE.ZT  <- mcmapply(function(x,y){
    Reffect <- ExtractionRE(x) # check if random effect
    if (Reffect==T){ # if there is 
      nameVAR <- Extraction(x) #extract the other effects 
      ## Simulation 
      data <- x$data[,colnames(x$data)%in%nameVAR[1:length(nameVAR)]] # Extract by colnames and in the right order
      data <- data[,c(nameVAR[1:length(nameVAR)])] # Reorder my data Need in all cases !!!
      # Creer table avec nouvelles valeurs 
      COUNTRY <- unique(data$country)
      myDF = as.data.frame(matrix(ncol = length(nameVAR), nrow = length(COUNTRY), NA)) # Df vierge
      colnames(myDF) <- colnames(data) # Bon
      ## We create a matrix of 30 values that we will predict by variable + 30 by level of plotcat. 
      #####
      # 1 # Remplir tout à la moyenne et tous les Plotcat à 0 par défaut. Et séquence sur le range 
      #####
      myDF[,nameVAR[nameVAR%in%colnames(data) & nameVAR!="Plotcat" & nameVAR!="country"]] <- apply(data[,nameVAR[nameVAR%in%colnames(data) & nameVAR!="Plotcat" & nameVAR!="country"]],2,
                                                                              function(b) as.numeric(rep(mean(b),nrow(myDF))))
      myDF[,"Plotcat"] <- as.factor(rep("0",nrow(myDF)))
      myDF[,"country"] <- COUNTRY
      #####
      # 2 # Prediction 
      #####
      sub_x = update(x, data=x$data) # hack because the formula is absent in the old model (so refit to take the formula for the predictions)
      #predict(sub_x,newdata=numeric(0)) ## TO RUN IF THERE IS A BUG
      pred.boot=predict(sub_x, newdata=myDF,re.form=NULL,type="response",binding="fit",list(predVar=TRUE)) # Predict the values for each cases describe in MyDf
      pred.boot[,"species"] <- y
      pred.boot <- cbind(pred.boot,get_intervals(sub_x, newdata=myDF,re.form=NULL,intervals="predVar",level = c(0.90),type="response",variances=list(linPred=TRUE,disp=FALSE,predVar=T,residVar=FALSE,cov=FALSE)))
      Var <- get_predVar(sub_x, newdata=myDF,re.form=NULL,level = c(0.95),which="residVar",variances=list(predVar=TRUE,linPred=TRUE,disp=FALSE,cov=FALSE))
      pred.boot[,"SEneg"] <- pred.boot[,"fit"]-Var
      pred.boot[,"SEpos"] <- pred.boot[,"fit"]+Var
      return(pred.boot[,c("country","fit","species","predVar_0.05","predVar_0.95","SEneg","SEpos")])
    }else {return(NULL)}
  },x=Mortmod.all.zt,y=names(Mortmod.all.zt),SIMPLIFY = F,mc.cores=20,mc.silent=F)
Species.RE.ZT <- do.call(rbind,Species.RE.ZT)
saveRDS(get("Species.RE.ZT"), paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/Country.zt.rds"))


#### Idem avec les modèles truncated
Species.RE.BIN  <- mcmapply(function(x,y){
  Reffect <- ExtractionRE(x) # check if random effect
  if (Reffect==T){ # if there is 
    nameVAR <- Extraction(x) #extract the other effects 
    ## Simulation 
    data <- x$data[,colnames(x$data)%in%nameVAR[1:length(nameVAR)]] # Extract by colnames and in the right order
    data <- data[,c(nameVAR[1:length(nameVAR)])] # Reorder my data Need in all cases !!!
    # Creer table avec nouvelles valeurs 
    COUNTRY <- unique(data$country)
    myDF = as.data.frame(matrix(ncol = length(nameVAR), nrow = length(COUNTRY), NA)) # Df vierge
    colnames(myDF) <- colnames(data) # Bon
    ## We create a matrix of 30 values that we will predict by variable + 30 by level of plotcat. 
    #####
    # 1 # Remplir tout à la moyenne et tous les Plotcat à 0 par défaut. Et séquence sur le range 
    #####
    myDF[,nameVAR[nameVAR%in%colnames(data) & nameVAR!="Plotcat" & nameVAR!="country"]] <- apply(data[,nameVAR[nameVAR%in%colnames(data) & nameVAR!="Plotcat" & nameVAR!="country"]],2,
                                                                                                 function(b) as.numeric(rep(mean(b),nrow(myDF))))
    myDF[,"Plotcat"] <- as.factor(rep("0",nrow(myDF)))
    myDF[,"country"] <- COUNTRY
    #####
    # 2 # Prediction 
    #####
    sub_x = update(x, data=x$data) # hack because the formula is absent in the old model (so refit to take the formula for the predictions)
    #predict(sub_x,newdata=numeric(0)) ## TO RUN IF THERE IS A BUG
    pred.boot=predict(sub_x, newdata=myDF,re.form=NULL,type="response",binding="fit",list(predVar=TRUE)) # Predict the values for each cases describe in MyDf
    pred.boot[,"species"] <- y
    pred.boot <- cbind(pred.boot,get_intervals(sub_x, newdata=myDF,re.form=NULL,intervals="predVar",level = c(0.90),type="response",variances=list(linPred=TRUE,disp=FALSE,predVar=T,residVar=FALSE,cov=FALSE)))
    Var <- get_predVar(sub_x, newdata=myDF,re.form=NULL,level = c(0.95),which="residVar",variances=list(predVar=TRUE,linPred=TRUE,disp=FALSE,cov=FALSE))
    pred.boot[,"SEneg"] <- pred.boot[,"fit"]-Var
    pred.boot[,"SEpos"] <- pred.boot[,"fit"]+Var
    return(pred.boot[,c("country","fit","species","predVar_0.05","predVar_0.95","SEneg","SEpos")])
  }else {return(NULL)}
},x=Mortmod.all.bin,y=names(Mortmod.all.bin),SIMPLIFY = F,mc.cores=20,mc.silent=F)
Species.RE.BIN <- do.call(rbind,Species.RE.BIN)
saveRDS(get("Species.RE.BIN"), paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/Country.bin.rds"))









### Figure
### And figures
Allcode <- c("PINPINA",
             "ABIALB",
             "ACEPSE",
             "ALNGLU",
             "BETPEN",
             "CASSAT",
             "FAGSYL",
             "FRAEXC",
             "PICABI",
             "PINHAL",
             "PINNIG",
             "PINPIN",
             "PINSYL",
             "POPNIG",
             "POPTRE",
             "QUEILE",
             "QUEPET",
             "QUEPYR",
             "QUEROB",
             "QUESUB")

Myshape20 <- c(0:19)
Mycol20 <- c(          # All colors 
  "#F8766D",
  "#E88526",
  "green",
  "red3",
  "gray70",
  "#5EB300",
  "gray50",
  "#00BF74",
  "#00C19F",
  "black",
  "#00B9E3",
  "gray30",
  "lightskyblue",
  "blue3",
  "#DB72FB",
  "gray15",
  "deeppink3",
  "orange",
  "yellow",
  "steelblue")


###############
PredVAR <- c( #
  "zt",       #
  "bin")      #
######################################################
Predname <- c(                                       #
  "Predicted mortality intensity (number of events)",#
  "Mortality occurrence probability")                #
######################################################
Mymod <- PredVAR[2] # maybe do both automatically 


myDF <- Species.RE.BIN

Myspecies <- unique(myDF$species)
MycolSignif <- Mycol20[match(Allcode,Myspecies)] # Here association of the signif species names with the associated colors (according to number)
MyshapeSignif <- Myshape20[which(Allcode%in%Myspecies)] # Here association of the signif species names with the associated colors (according to number)



### Change color for country !!!
## And inverse speciues on x axis and c


try(rm(p1),silent = T)
p1<-ggplot(myDF,aes(x=species, y=fit,color=species,group=country,shape=country))+
  #geom_line(size=1)+
  geom_point(size=2.5,stroke=1.25,position = position_dodge2(width = 0.85,preserve = c("total")))+
#p1 <- p1+facet_wrap(vars(species),scales = "fixed")+scale_color_manual(values = rep("black",length(Myspecies))) # to delete or not 
  #geom_linerange(aes(ymin = SEneg, ymax = SEpos),alpha=1,size=0.8,position = position_dodge2(width = 0.85,preserve = c("total")))+ #to test
  geom_linerange(aes(ymin = predVar_0.05, ymax = predVar_0.95),alpha=1,size=0.8,position = position_dodge2(width = 0.85,preserve = c("total")))+ #to test
labs(y=Predname[PredVAR==Mymod], x="Species")+
  # title=paste0(Predname[PredVAR==Pred]," VS ",nameReal[nameVAR==VARIABLE])
  scale_color_manual(name=NULL,values = MycolSignif)+
  scale_shape_manual(name="Country",values = MyshapeSignif[1:6])+
  scale_fill_manual(name="",values = MycolSignif)+
  guides(colour = FALSE,shape = guide_legend(ncol = 2))+
  #ylim(0,max(myDF[,Pred])*1.1)+
  #scale_y_continuous(expand = c(0.1,0)) + # for the loop (most of the figure)
  scale_y_continuous(expand = c(0.02,0)) + # after the loop
  #coord_cartesian(ylim=c(0,180))+
  ## Ajout commande pour dépasser 
  
  #scale_x_continuous(expand = c(0.03,0)) +
  theme_classic(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(.01, .99),
        legend.justification = c("left", "top"),legend.direction ="vertical",
        axis.text.x = element_text(size=13,color="black",angle = 45,hjust = 1,vjust = 1),axis.text.y = element_text(size=13,color="black"),
        legend.text=element_text(size=12),
        legend.key.size = unit(1.5,"line"),
        legend.margin = margin(1,1,1,1),
        legend.box.margin = margin(1,1,1,1),
        #legend.title = element_blank(),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))
print(p1)
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/Plot_Country.",Mymod,"_TRUE.png"),
          plot = p1,base_width = 12, base_height = 7, dpi = 400 ,units = "in",1)
saveRDS(p1,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/Plot_Country.",Mymod,"_TRUE.rds"))


### Random effect country 

MyVAR <- "Country"
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.TRUE.rds$"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=PredVAR[c(1:2)]) ## be carreful with the order 

pall <- p.zt+
  labs(title=c("a) Country random effect (BIN models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  
  p.bin+
  labs(title=c("b) Country random effect (NB models)"))+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  
  plot_layout(ncol=2)
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/",MyVAR,"Effect_bothmodel.pdf"),
          plot = pall,base_width = 9, base_height = 7, dpi = 300 ,units = "in",nrow=1,ncol=2)


## Association of species with colors and shape for all figures: 
names(Mycol20) <- levels(as.factor(Allcode))
names(Myshape20) <- levels(as.factor(Allcode))


## Years effect 
list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = ".rds")
MyVAR <- "yearsbetweensurveys"
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.TRUE.rds$"),full.names = T)
B <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.FALSE.rds$"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=PredVAR[c(1:2)]) ## be carreful with the order 

p.zt$data$species <- as.factor(p.zt$data$species)
pall <- p.zt+
  labs(title=c("a) Census interval effect (BIN models)"))+
  scale_color_manual("",values = Mycol20)+
  scale_shape_manual("",values = Myshape20)+
  scale_fill_manual("",values = Mycol20)+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  coord_cartesian(ylim=c(0,0.4))+
  
  p.bin+
  labs(title=c("b) Census interval effect (NB models)"))+
  scale_color_manual(name="",values = Mycol20)+
  scale_shape_manual(name="",values = Myshape20)+
  scale_fill_manual(name="",values = Mycol20)+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  coord_cartesian(ylim=c(0,150))+

  
  plot_layout(ncol=2)
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/",MyVAR,"Figure_S12.pdf"),
          plot = pall,base_width = 9, base_height = 7, dpi = 300 ,units = "in",nrow=1,ncol=2)


mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=B,y=PredVAR[c(1:2)]) ## be carreful with the order 
pall2 <- p.zt+
  scale_color_manual(name="",values = Mycol20)+
  scale_shape_manual(name="",values = Myshape20)+
  scale_fill_manual(name="",values = Mycol20)+
  labs(title=c("a) Census interval effect (BIN models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  p.bin+
  scale_color_manual(name="",values = Mycol20)+
  scale_shape_manual(name="",values = Myshape20)+
  scale_fill_manual(name="",values = Mycol20)+
  labs(title=c("b) Census interval effect (NB models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  plot_layout(ncol=2)
pall2
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/",MyVAR,"Effect_bothmodel.No.SE.png"),
          plot = pall2,base_width = 9, base_height = 7, dpi = 300 ,units = "in",nrow=1,ncol=2)


## treedensity 
list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = ".rds")
MyVAR <- "treeNbr"
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.TRUE.rds$"),full.names = T)
B <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.FALSE.rds$"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=PredVAR[c(1:2)]) ## be carreful with the order 
pall <- p.zt+
  scale_color_manual(name="",values = Mycol20)+
  scale_shape_manual(name="",values = Myshape20)+
  scale_fill_manual(name="",values = Mycol20)+
  labs(title=c("a) Tree density effect (BIN models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  #scale_y_continuous(expand = c(0.1,0)) + # for the loop (most of the figure)
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  theme(legend.position = c(0.01,0.99),legend.justification = c("left","top"))+guides(colour = guide_legend(ncol = 1))+
  #coord_cartesian(ylim=c(0,0.4))+
  p.bin+
  scale_color_manual(name="",values = Mycol20)+
  scale_shape_manual(name="",values = Myshape20)+
  scale_fill_manual(name="",values = Mycol20)+
  labs(title=c("b) Tree density effect (NB models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  theme(legend.position = c(0.01,0.99),legend.justification = c("left","top"))+guides(colour = guide_legend(ncol = 1,size=14))+
  coord_cartesian(ylim=c(0,350))+
  plot_layout(ncol=2)
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/",MyVAR,"S11.pdf"),
          plot = pall,base_width = 9, base_height = 7, dpi = 300 ,units = "in",nrow=1,ncol=2)



## dbh
list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = ".rds")
MyVAR <- "logdbh"
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.TRUE.rds$"),full.names = T)
B <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.FALSE.rds$"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=PredVAR[c(1:2)]) ## be carreful with the order 
pall <- p.zt+
  labs(title=c("a) Mean DBH effect (BIN models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  #scale_y_continuous(expand = c(0.1,0)) + # for the loop (most of the figure)
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  theme(legend.position = c(0.01,0.99),legend.justification = c("left","top"))+guides(colour = guide_legend(nrow = 3))+
  #coord_cartesian(ylim=c(0,0.4))+
  scale_color_manual(name="",values = Mycol20)+
  scale_shape_manual(name="",values = Myshape20)+
  scale_fill_manual(name="",values = Mycol20)+
  p.bin+
  scale_color_manual(name="",values = Mycol20)+
  scale_shape_manual(name="",values = Myshape20)+
  scale_fill_manual(name="",values = Mycol20)+
  labs(title=c("b) Mean DBH effect (NB models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  theme(legend.position = c(0.01,0.99),legend.justification = c("left","top"))+guides(colour = guide_legend(nrow = 2))+
  coord_cartesian(ylim=c(0,350))+
  plot_layout(ncol=2)
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/",MyVAR,"S10.pdf"),
          plot = pall,base_width = 9, base_height = 7, dpi = 300 ,units = "in",nrow=1,ncol=2)



### BAIJmean
list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = ".rds")
MyVAR <- "BAIj.plot.1.mean"
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.TRUE.rds$"),full.names = T)
B <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.FALSE.rds$"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=PredVAR[c(1:2)]) ## be carreful with the order 
pall <- p.zt+
  labs(title=c("a) Species mean growth rate (BIN models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  #theme(legend.position = c(0.01,0.99),legend.justification = c("left","top"))+guides(colour = guide_legend(ncol = 1))+
  #coord_cartesian(ylim=c(0,0.4))+
  p.bin+
  labs(title=c("b) Species mean growth rate (NB models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  #theme(legend.position = c(0.01,0.99),legend.justification = c("left","top"))+guides(colour = guide_legend(ncol = 1))+
  coord_cartesian(ylim=c(0,280))+
  plot_layout(ncol=2)
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/",MyVAR,"Effect_bothmodel.png"),
          plot = pall,base_width = 9, base_height = 7, dpi = 300 ,units = "in",nrow=1,ncol=2)



### Intra
list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = ".rds")
MyVAR <- "BAj"
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.TRUE.rds$"),full.names = T)
B <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.FALSE.rds$"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=PredVAR[c(1:2)]) ## be carreful with the order 
pall <- p.zt+
  labs(title=c("a) Intraspecific competition effect (BIN models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  #theme(legend.position = c(0.01,0.99),legend.justification = c("left","top"))+guides(colour = guide_legend(ncol = 1))+
  #coord_cartesian(ylim=c(0,0.4))+
  p.bin+
  labs(title=c("b) Intraspecific competition effect (NB models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  #theme(legend.position = c(0.01,0.99),legend.justification = c("left","top"))+guides(colour = guide_legend(ncol = 1))+
  coord_cartesian(ylim=c(0,200))+
  plot_layout(ncol=2)
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/",MyVAR,"Effect_bothmodel.png"),
          plot = pall,base_width = 9, base_height = 7, dpi = 300 ,units = "in",nrow=1,ncol=2)



## Mean spei
list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = ".rds")
MyVAR <- "mean_spei"
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.TRUE.rds$"),full.names = T)
B <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.FALSE.rds$"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=PredVAR[c(1:2)]) ## be carreful with the order 
pall <- p.zt+
  labs(title=c("a) Relative drought effect (BIN models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  #theme(legend.position = c(0.01,0.99),legend.justification = c("left","top"))+guides(colour = guide_legend(ncol = 1))+
  #coord_cartesian(ylim=c(0,0.4))+
  p.bin+
  labs(title=c("b) Relative drought effect (NB models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  #theme(legend.position = c(0.01,0.99),legend.justification = c("left","top"))+guides(colour = guide_legend(ncol = 1))+
  coord_cartesian(ylim=c(0,120))+
  plot_layout(ncol=2)
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/",MyVAR,"Effect_bothmodel.png"),
          plot = pall,base_width = 9, base_height = 7, dpi = 300 ,units = "in",nrow=1,ncol=2)




### Inter
list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = ".rds")
MyVAR <- "BA.O"
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.TRUE.rds$"),full.names = T)
B <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.FALSE.rds$"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=PredVAR[c(1:2)]) ## be carreful with the order 
pall <- p.zt+
  labs(title=c("a) Interspecific competition effect (BIN models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  #theme(legend.position = c(0.01,0.99),legend.justification = c("left","top"))+guides(colour = guide_legend(ncol = 1))+
  coord_cartesian(ylim=c(0,0.05))+
  p.bin+
  labs(title=c("b) Interspecific competition effect (NB models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  #theme(legend.position = c(0.01,0.99),legend.justification = c("left","top"))+guides(colour = guide_legend(ncol = 1))+
  coord_cartesian(ylim=c(0,250))+
  plot_layout(ncol=2)
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/",MyVAR,"Effect_bothmodel.png"),
          plot = pall,base_width = 9, base_height = 7, dpi = 300 ,units = "in",nrow=1,ncol=2)


## Check temp
list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = ".rds")
MyVAR <- "bio1"
A <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.TRUE.rds$"),full.names = T)
B <- list.files(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect"),pattern = paste0(MyVAR,".*.FALSE.rds$"),full.names = T)
# Load only the first half 
mapply(function(x,y){
  assign(paste0("p.",y),readRDS(x),envir = .GlobalEnv)
},x=A,y=PredVAR[c(1:2)]) ## be carreful with the order 
pall <- p.zt+
  labs(title=c("a) Annual mean temperature effect (BIN models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  #theme(legend.position = c(0.01,0.99),legend.justification = c("left","top"))+guides(colour = guide_legend(ncol = 1))+
  coord_cartesian(ylim=c(0,0.4))+
  p.bin+
  labs(title=c("b) Annual mean temperature effect (NB models)"))+
  theme(plot.title = element_text(size=16,hjust = 0))+
  theme(legend.text=element_text(size=10))+
  theme(legend.key.size = unit(1.2,"line"))+
  #theme(legend.position = c(0.01,0.99),legend.justification = c("left","top"))+guides(colour = guide_legend(ncol = 1))+
  coord_cartesian(ylim=c(0,150))+
  plot_layout(ncol=2)
pall
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/",MyVAR,"Effect_bothmodel.png"),
          plot = pall,base_width = 9, base_height = 7, dpi = 300 ,units = "in",nrow=1,ncol=2)








