rm(list = ls())
gc()
library(lattice)
library(tidyr)
library(usdm)

#### 1 #### Load our data 
Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
source(paste0(Dir,"Myscripts/Fundiv.project/function3.ModelBoot.R"))


i = 1
Allcode <- c("PINPINA","ABIALB","ACEPSE","ALNGLU","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINNIG","PINPIN","PINSYL","POPNIG","POPTRE","QUEILE","QUEPET",
             "QUEPYR","QUEROB","QUESUB")
Allmod <- c("Mbin13A.18","Mbin14A.19","Mbin13B.26","Mbin13A.22","M2bin13B.22","M2bin13A21","M2bin15B23","Mbin3C.26","M2bin7B.17","Mbin11B.19","Mbin3B.20","Mbin15B.21",
            "Mbin5B.17","Mbin3B.31","Mbin13A.27","M2bin1C.20","M2bin1B.23","Mbin13B.26","Mbin7A.26","Mbin1B.23")            
Allseuil <- c(0.8,0.7,0.55,0.6,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.7,0.8,0.7,0.7,0.8,0.7,0.7,0.8,0.7)
AllseuilC <- c(0.6,0.6,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)

for (i in 1:length(Allcode)){
  print(paste0(Allcode[i],Allmod[i]))
  Dir = (paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",Allcode[i])) # Directory 
  setwd(Dir)
  assign(paste0("dfplot"),readRDS(paste0("CLIMAP/dfplot2",Allcode[i],Allseuil[i],".rds"))) #Base de données plot
  load(paste0(Dir,"/",Allcode[i],"_summary.RData"),envir = globalenv())
  load(paste0(Dir,"/CLIMAP/Climap_",Allcode[i],"_",Allseuil[i],"_",AllseuilC[i],".RData"),envir = globalenv())
  ras <- rasterFromXYZ(test6[,c(1:2,26)],crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ") #New raster with our three categories 
  A=summary_all[summary_all$Clim_Var=="bio1",c(2:4)] # my plots
  M.Plotcat <- raster::extract(ras,dfplot[,2:3])
  Myplot = cbind(dfplot,M.Plotcat)
  z <- Myplot
  
  ##### 2 ### Extract the names of the variable for the VIF analysis 
  Dir <- paste0(Dir,"/CLIMAP/Models/binomial/",Allmod[i])
  setwd(Dir)
  assign(sub(".rda","",list.files(pattern = ".rda"), ignore.case = FALSE,fixed = T),get(load(file = list.files(pattern = ".rda")))) 
  rm("x")# Load the function 
  Dir <- c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/VIF.all")
  setwd(Dir)
  Explain <- Extraction(get(Allmod[i]))
  Explain[Explain=="Plotcat"] <- "M.Plotcat" # Replace the plotcat categorical variable by the continous one
  Resp <- c("sp.mortality.plot.count.yr")
  VIF <- vif(z[,Explain])
  VIFSTEP <- vifstep(z[,Explain],th=8)
  colnames(VIF)[1] <- paste0(Allcode[i]," Variables")
  #capture.output(print(VIF),file=paste0(Allcode[i],Resp,"_Vif.txt")) # Output as a latex wrapped in a txt file
  #capture.output(print(VIFSTEP),file=paste0(Allcode[i],Resp,"_Vifstep.txt")) # Output as a latex wrapped in a txt file
  write.csv(VIF, sep=",", file=paste0(Resp,"_All.Vif.csv",row.names = F),append = T) # Table for the paper with coefficient
  #write.table(c(Allcode[i],VIFSTEP), sep=",", file=paste0(Resp,"_All.Vifstep.csv",row.names = F),append = T) # Table for the paper with coefficient
  
  Resp <- c("sp.mort.bin")
  VIF <- vif(z[,Explain])
  VIFSTEP <- vifstep(z[,Explain],th=8)
  #capture.output(print(VIF),file=paste0(Allcode[i],Resp,"_Vif.txt")) # Output as a latex wrapped in a txt file
  #capture.output(print(VIFSTEP),file=paste0(Allcode[i],Resp,"_Vifstep.txt")) # Output as a latex wrapped in a txt file
  colnames(VIF)[1] <- paste0(Allcode[i]," Variables")
  write.table(VIF, sep=",", file=paste0(Resp,"_All.Vif.csv",row.names = F),append = T) # Table for the paper with coefficient
  #write.table(c(Allcode[i],VIFSTEP), sep=",", file=paste0(Resp,"_All.Vifstep.csv",row.names = F),append = T) # Table for the paper with coefficient
}

