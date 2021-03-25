## This script follow 36. 
## It plots the effect we want and wrap it for both model 
rm(list = ls())
gc()

##############
##         ##
##  Figure ##
##         ##
##############

# Here are the orginal variable possible to plot

#################################################
nameVAR <- c("sqrtBAIj.plot.1.mean",            #
             "sqrtBAIj.plot.1",                 #
             "logdbh.plot.mean",                #
             "sqrtdbh.plot.mean",               #
             "treeNbr",                         #
             "yearsbetweensurveys",             #
             "min_spei12",                      #
             "mean_spei12",                     #
             "bio1_climate_mean.30",            #
             "bio5_climate_mean.30",            #
             "tmean.djf_climate_mean.30",       #
             "bio14_climate_mean.30",           #
             "bio13_climate_mean.30",           #
             "ppet.mean_climate_mean.30",       #
             "bio12_climate_mean.30",           #
             "sqrtBA.ha.plot.1",                #
             "sqrtBA.O.plot.1",                 #
             "logBAj.plot.1")                   #
#################################################

# Here are the real name 

#################################################################
nameReal <- c("Species mean\ngrowth rate (square rooted)",      #
              "Species total\ngrowth rate (square rooted)",     #
              "Species mean DBH\n(log)",                     #
              "Species mean DBH\n(sqrt)",                    #
              "Tree density",                                #
              "Years between \nsurveys",                     #
              "Min relative \ndrought index",                #
              "Mean relative \ndrought index",               #
              "Annual mean temperature",                     #
              "Maximal temperature \nof the warmest month ", #
              "Winter mean temperature",                     #
              "Driest month precipitation",                  #
              "Wettest month precipitation",                 #
              "Annual water balance",                        #
              "Annual mean precipitation",                   #
              "Total \ncompetition (log)",                   #
              "Interspecific \ncompetition (square rooted)", #
              "Intraspecific \ncompetition (log)")           #
#################################################################

###############
PredVAR <- c( #
  "zt",       #
  "bin")      #
######################################################
Predname <- c(                                       #
  "Predicted mortality intensity (number of events)",#
  "Mortality occurrence probability")                #
######################################################

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

## First pick a variable and a model (or not)
VARIABLE <- nameVAR[1]
Mymod <- PredVAR[2] # maybe do both automatically 
FACET <- FALSE
SE <- FALSE

SimpleEffect <- function(VARIABLE = "sqrtBAIj.plot.1.mean",Mymod="bin",SE=F,save=F,FACET=F){
print(paste0(VARIABLE,".",Mymod))
  
## Load the database corresponding
assign("myDF",readRDS(paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/",VARIABLE,".",Mymod,".rds")))

## load the model significiance file and extract the species for which it is significant
if (Mymod==PredVAR[1]){tableresults <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/tableresultsNB.csv")
colnames(tableresults)[2] <- "variable"
}else {tableresults <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/tableresults.csv")}
## Then check that the database we have correponds to the species in which the effect is significant ! 
Signifspecies <- as.character(unique(tableresults[tableresults$variable%in%c(VARIABLE) & tableresults$signif=="*",
                                    "species"])) # Keep only significant values
Signifspecies <- sub(" *r","",x=Signifspecies,fixed=F)

myDF <-myDF[myDF$species%in%Signifspecies,] ### Keep only species for which the effect is significant (normally tghis is already the)

Myspecies <- unique(myDF$species)
MycolSignif <- Mycol20[which(Allcode%in%Myspecies)] # Here association of the signif species names with the associated colors (according to number)
MyshapeSignif <- Myshape20[which(Allcode%in%Myspecies)] # Here association of the signif species names with the associated colors (according to number)


## Now the plot 

try(rm(p1),silent = T)
p1<-ggplot(myDF,aes(x=get(VARIABLE), y=fit,color=species,group=species,shape=species))+
geom_line(size=1)+
geom_point(size=3,stroke=1)
if(FACET==T){p1 <- p1+facet_wrap(vars(species),scales = "free")+scale_color_manual(values = rep("black",length(Myspecies)))
}else{p1 <- p1}
if(SE==T){p1 <- p1+
  geom_ribbon(aes(ymin = predVar_0.05, ymax = predVar_0.95,fill=species),alpha=0.1,size=0.5,lty="blank") #to test
}else{p1 <- p1}
p1 <- p1 +
  labs(y=Predname[PredVAR==Mymod], x=nameReal[nameVAR==VARIABLE])+
  # title=paste0(Predname[PredVAR==Pred]," VS ",nameReal[nameVAR==VARIABLE])
  scale_color_manual(name="",values = MycolSignif)+
  scale_shape_manual(name="",values = MyshapeSignif)+
  scale_fill_manual(name="",values = MycolSignif)+
  #ylim(0,max(myDF[,Pred])*1.1)+
  #scale_y_continuous(expand = c(0.1,0)) + # for the loop (most of the figure)
  scale_y_continuous(expand = c(0.05,0)) + # after the loop
  scale_x_continuous(expand = c(0.03,0)) +
  theme_light(base_size = 15)+
  theme(text = element_text(face="bold"),legend.position = c(.01, .99),
        legend.justification = c("left", "top"),legend.direction ="horizontal",
        axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
        legend.text=element_text(size=12),
        legend.key.size = unit(1.5,"line"),
        legend.margin = margin(1,1,1,1),
        legend.box.margin = margin(1,1,1,1),
        legend.title = element_blank(),
        legend.background=element_rect(fill="white",colour="black",size=0.2),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        axis.line = element_line(colour="black"),
        plot.title = element_text(size=18,hjust = 0.5),
        plot.caption = element_text(face="bold.italic"))

print(p1)
if (save==T){
save_plot(filename = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/Plot_",VARIABLE,".",Mymod,"_",SE,".png"),
            plot = p1,base_width = 12, base_height = 7, dpi = 400 ,units = "in",1)
saveRDS(p1,file = paste0("~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/SimpleEffect/Plot_",VARIABLE,".",Mymod,"_",SE,".rds"))}
}


SimpleEffect(VARIABLE = nameVAR[5],Mymod=PredVAR[1],FACET=F,save=T,SE=F)



for (i in 1:length(nameVAR)){
  for (j in 1:2){
    SimpleEffect(VARIABLE = nameVAR[i],Mymod=PredVAR[j],FACET=F,save=T,SE=T)
  }
}  


### Figure cowplot ! 
# Ned to select a few








