# Coauthor script responses ! 

Dir=c("/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/")
load(paste0(Dir,"our-data/tfinal.biotic.Feb2019.RData")) #This is the database we want from which we want to analyse the metrics

load(paste0(Dir,"our-data/tfinal.biotic.Feb2019.RData")) #This is the database we want from which we want to analyse the metrics
tfinal <- readRDS(paste0(Dir,"our-data/tfinal.biotic.August2020.rds"))
A <- by(tfinal$country,tfinal$code,summary)
QUERUB <- tfinal[tfinal$code=="QUERUB",]
length(unique(QUERUB$plotcode))
summary(QUERUB$country)

ROBPSE <- tfinal[tfinal$code=="ROBPSE",]
length(unique(ROBPSE$plotcode))
summary(ROBPSE$country)

PRUSER <- tfinal[tfinal$code=="PRUSER",]
length(unique(PRUSER$plotcode))


AILALT <- tfinal[tfinal$code=="AILALT",]
length(unique(AILALT$plotcode))

PINSYL <- tfinal[tfinal$code=="PINSYL",]
length(unique(PINSYL$plotcode))


summary(as.factor(tfinal.biotic$yearsbetweensurveys))

## 29 years census intervalle 
colnames(tfinal.biotic)
test <- tfinal.biotic[tfinal.biotic$yearsbetweensurveys==2,c("yearsbetweensurveys","surveydate1","surveydate2","plotcode","longitude","latitude","treecode","country")]
test <- subset(tfinal.biotic,yearsbetweensurveys==2,select=c("yearsbetweensurveys","surveydate1","surveydate2","plotcode","longitude","latitude","treecode","country"))


## Wallony mortality 
test <- subset(tfinal.biotic,country=="WA",select=c("plotcode","longitude","latitude","sp.mortality.plot.rate"))
summary(as.factor(test$sp.mortality.plot.rate))
boxplot(test$sp.mortality.plot.rate)
boxplot(tfinal.biotic$sp.mortality.plot.rate~tfinal.biotic$country)


tfinal <- tfinal.biotic[tfinal.biotic$code%in%Allcode,]
tfinal$code <- as.character(tfinal$code)
Allcode <- c("LARDEC","QUEPUB","ABIALB","ACEPSE","BETPEN","CASSAT","FAGSYL","FRAEXC","PICABI","PINHAL","PINSYL","POPNIG","QUEILE","PINPINA","ALNGLU","PINNIG","PINPIN","POPTRE","QUEPET","QUEPYR","QUEROB","QUESUB")
A <- by(tfinal$yearsbetweensurveys,tfinal$code,summary)
?by
colnames(tfinal.biotic)
All <- do.call(cbind,A)
All2 <- All[c(1,4,6),]
All3 <- t(All2)
All3 <- round('