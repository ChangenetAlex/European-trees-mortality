#########################################################################"
###########   Extract climatic data for any species       ###############"
###########       from a specified time period            ###############"
###########    And calculate the mean,max and min         ###############"
###########         At different intervalles              ###############"
#########################################################################"

SavingInfo = "ARGUMENTS
# 'Dir is a character. This is the path in which different folders and tables are stored
# 'path is a character. This is the path in which climatic data are stored. 
# 'CODE' character giving the species we want to extract. See AllspFinal for the full list of species. 
# 'yearINI' is the first year we want to extract
# 'yearFIN' is the last year we want to extract
# 'Intervalle' is the vector indicating how far in the past we want to summarize the climatic data. For instance '5' means we want the min the max and the mean for each climatic data on the last 5 years. By default i has the value 5,10 and 15
# 'files.wanted is a vector giving the climatic parameter we are intesrest in. The list of the 21 is called files. 
# 'save' logical. If TRUE, 2 tables will be stored (full path of the destination folder given and can be modified in the function).

# Exemple : ExtractClimate(CODE=species,yearINI=yearI,yearFIN=yearF,Intervalle=Inter,files.wanted=files.w,save=T)

"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))

#Required packages 
library(stringr)
library(parallel)
library(dplyr)
require(reshape2)


#Modif du 20/02/2018
Dir=c("~/SynologyDrive/FUNDIV - NFI - Europe/")
#load(paste0(Dir,"our-data/treefinal.RData"))
#load(paste0(Dir,"our-data/tfinal.biotic.June2018.RData"))
#treefinal <- tfinal.biotic
#Here is the code of the particular species we want to work with 
species = "QUERUB"
#Here are the years we want basically extract for all plots
yearI = 1960
yearF = 2014

# climatic variables provided by EuMedClim
bioclim=paste("bio",c(1,12,13,14,2,5,6),sep="") #
T.seas=paste0("tmean.",c("djf","jja","mam", "son"))#
P.seas=paste0("prec.",c("djf","jja","mam","son"))#
pet=paste0("pet.",c("max","mean","min")) #
ppet=paste0("ppet.",c("max","mean","min"))#
files=(c(bioclim, pet, ppet, P.seas, T.seas))
files.w=files[c(1:21)] #Select among the 21 parameters the ones that we want to extract

#This parameter define how many years we want from the surveys 
Inter=c(5,10,15,30)

#Fonction
ExtractClimate = function(Dir=c("~/SynologyDrive/FUNDIV - NFI - Europe/"),
                          path=paste0(Dir,"our-data/climate/RDATA"),
                          CODE,
                          yearINI = 1901,
                          yearFIN = 2014,
                          Intervalle = c(5,10,15),
                          files.wanted,
                          save=T){
  
  bioclim=paste("bio",c(1,12,13,14,2,5,6),sep="") #
  T.seas=paste0("tmean.",c("djf","jja","mam", "son"))#
  P.seas=paste0("prec.",c("djf","jja","mam","son"))#
  pet=paste0("pet.",c("max","mean","min")) #
  ppet=paste0("ppet.",c("max","mean","min"))#
  files=(c(bioclim, pet, ppet, P.seas, T.seas))
  
  if(is.null(files.wanted)) files.wanted=files else
    for(ki in files.wanted) if(all(ki!=files))
     stop(paste('the climatic variable',ki,'is not provided by EuMedClim'))
  
  if(yearINI < 1901)
    stop("Year minimum is 1901, please provide a year above that time you jerk !")
  
  if(yearFIN > 2014)
    stop("Year maximum is 2014, please provide a year below that time you jerk !")
  
  if(yearINI > yearFIN)
    stop("Year minimum cannot be supriori to the final year idiot ! Please provide something different or exchange the parameters. ")
  
  if(is.null(Dir))
   stop("the argument 'Dir' is empty. Precise the full path of tables folder")
  
  if(is.null(path))
   stop("the argument 'path' is empty. Precise the full path of climate data folder (and/or where to store downloaded files)")
  
  Allspfinal <- read.csv2(paste0(Dir,"/our-data/Allspfinal.csv"))  #These are all the referenced species
  
  if(species %in% unique(Allspfinal$code)==FALSE)
    stop("Common. We have 201 species in the list. Could you not find at least one which match the list ? Please try again")
  
  load(paste0(Dir,"our-data/tfinal.biotic.June2018.RData")) #This is the database we want from which we want to analyse the metrics
  treefinal <- tfinal.biotic
  
  #Here is extracted the species we want from the global individual database
  sp <- treefinal[treefinal$code==CODE&!is.na(treefinal$code),]
  
  #sp.df <- readRDS(paste0(Dir,"our-data/species/",CODE,"/Mydf_",CODE,".rds")) ## This line added on december 2020 to process invasive species
  sp.df <- readRDS(paste0(Dir,"our-data/species/",CODE,"/Mydf_FR_",CODE,".rds")) ## This line added on jan 2021 to process invasive species FR
  
  #Here are listed the files corresponding to the variable (on 21) we want to extact
  file_names<-list.files(path=path,pattern=paste0("_",files.wanted, collapse = "|"), full.names=T) 
  
  #Here these are loaded for ALL plots
  load_files<- sapply(file_names, function(x) get(load(x)), simplify = FALSE) 
  
  #Here are taken the years we wanted to extract from year INI to yearfin. (Both included)
  A=lapply(load_files, function(x) x[c(1,2,(-1898+yearINI):(-1898+yearFIN))])
  
  #And are precised what plotcode we want on the 154000 (basically the ones in which the species is present)
  B=lapply(A, function(x) x[c(rownames(x)%in%sp[!duplicated(sp[,c('plotcode')]),c(1)]),]) #seulement les plots d'intéret
  
  #These are bind together 
  rbind_load_files <- do.call("rbind",B) #ici déjà un tableau
  #str(rbind_load_files)
  
  #The plotocde are given. Here rownames = plotID  and the file_names = the number of parameters wanted. 
  rbind_load_files$plotcode<-rep(rownames(B[[1]]),length(file_names)) #modified 30/01/18
  
  #Here we name each climatic parameter for each plotcode and bind it with our list of data
  clim_files=rep(files.wanted,each=length(unique(sp$plotcode))) 
  all_climate<-cbind(clim_files, rbind_load_files)  
  #Here are extracted the date of survey for french and no-french data. (Need to take a year before). This varible is added on all clim
  treefinal[treefinal$country!="FR","year"] <- as.numeric(str_sub(treefinal[treefinal$country!="FR","surveydate2"],1,4))
  #08/03/2018 Added lines to take the median year between both s
  #treefinal[treefinal$country!="FR","year"] <-round((as.numeric(str_sub(treefinal[treefinal$country!="FR","surveydate2"],1,4))
                    #+as.numeric(str_sub(treefinal[treefinal$country!="FR","surveydate1"],1,4)))/2)
  
  all_climate$surveydate <- treefinal$year[match(all_climate$plotcode, treefinal$plotcode,incomparables = NA)]
  #A dataframe is created and for each climatic variable we calculate the mean, the max, and the min for each intervalle specified. 
  climate_var<-data.frame(matrix(NA,ncol=1+length(Intervalle)*3,nrow=length(clim_files)))
  climate_var[,1]<-clim_files
  for (i in 1:length(Intervalle)){
    climate_var[3*i-1]<-cbind(round(apply(all_climate, 1, function(x) mean(as.numeric(x[(length(x)-2-(yearFIN-as.numeric(last(x)))-Intervalle[i]):
                                                                                          (length(x)-2-(yearFIN-as.numeric(last(x))))]))),2))
    
    climate_var[3*i]<-cbind(round(apply(all_climate, 1, function(x) min(as.numeric(x[(length(x)-2-(yearFIN-as.numeric(last(x)))-Intervalle[i]):
                                                                                       (length(x)-2-(yearFIN-as.numeric(last(x))))]))),2))
    
    climate_var[3*i+1]<-cbind(round(apply(all_climate, 1, function(x) max(as.numeric(x[(length(x)-2-(yearFIN-as.numeric(last(x)))-Intervalle[i]):
                                                                                         (length(x)-2-(yearFIN-as.numeric(last(x))))]))),2))
    i=i+1
  }
  #Important to note here that for a survey in year 2001. We tke the average (for 5 years) in between 1996 to 2001 included. 
  #This is 6 years wether or not samples have been taken in january or december. 
  #If we want to take 1997 to 2001 we should had "+1" after "intervalle[i]".
  #If we want to take from 1996 to 2000 we should had "-1" after the second "last (x)"
  
  #Union of mean, max and min and rename colum according to the previous calculation
  summary_all=cbind(all_climate[,c(1:3,(ncol(all_climate)-1):ncol(all_climate))],climate_var[,c(2:ncol(climate_var))])
  namecol = c()
  for (j in 1:length(Intervalle)){
    namecol[3*j-2]<-paste0("climate_mean.",Intervalle[j])
    namecol[3*j-1]<-paste0("climate_min.",Intervalle[j])
    namecol[3*j]<-paste0("climate_max.",Intervalle[j])
    
  } 
  #Named are given automatically after the calculations
  names(summary_all) <- c("Clim_Var", "longitude", "latitude","plotcode","surveydate",namecol)
  
  #Merge climate data with inventory individual data
  Sp<-merge(sp, summary_all, by="plotcode")
  Sp.df<-merge(sp.df, summary_all, by="plotcode") # added on december 2020 to process invasive species 
  
  # Convert the DF with the format we want
  long1 <- melt(summary_all, id=c("plotcode", "latitude", "longitude","surveydate","Clim_Var")) # Melt the table first
  wide1 <- dcast(long1,plotcode+latitude+longitude+surveydate~Clim_Var+variable) # Rebuild the table 
  Sp2<-merge(sp, wide1, by="plotcode")
  Sp2.df<-merge(sp.df, wide1, by="plotcode") # added on december 2020 to process invasive species 
  
  # And save as RData if it is specified
  if(save == "TRUE") 
  { dir.create(path=paste0(Dir,"our-data/species/",CODE,"/"))
    #saveRDS(Sp, file=paste0(Dir,"our-data/species/",CODE,"/",CODE,"_all.rds"))
    #saveRDS(summary_all, file=paste0(Dir,"our-data/species/",CODE,"/",CODE,"_summary.rds"))
    saveRDS(Sp2.df, file=paste0(Dir,"our-data/species/",CODE,"/","Mydf_",CODE,"_allVariable.rds")) # added on december to process invasive species 
    saveRDS(Sp2, file=paste0(Dir,"our-data/species/",CODE,"/",CODE,"_allVariable.rds"))} #Here the database is saved in the same file
  }

### ARGUMENTS
# 'Dir is a character. This is the path in which different folders and tables are stored
# 'path is a character. This is the path in which climatic data are stored. 
# 'CODE' character giving the species we want to extract. See AllspFinal for the full list of species. 
# 'yearINI' is the first year we want to extract
# 'yearFIN' is the last year we want to extract
# 'Intervalle' is the vector indicating how far in the past we want to summarize the climatic data. For instance '5' means we want the min the max and the mean for each climatic data on the last 5 years. By default i has the value 5,10 and 15
# 'files.wanted is a vector giving the climatic parameter we are intesrest in. The list of the 21 is called "files". 
# 'save' logical. If TRUE, 2 tables will be stored (full path of the destination folder given and can be modified in the function).