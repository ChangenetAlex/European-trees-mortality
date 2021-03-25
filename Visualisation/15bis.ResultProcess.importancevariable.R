library(reshape2)
rm(list = ls())
gc()
getwd()

tableresults <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/tableresults.csv")
summary(as.factor(tableresults$variable))
tableresults <- tableresults[-c(which(tableresults$variable=="(Intercept)")),] # On enlève l'intercept
tableresults <- tableresults[grep(tableresults$variable,pattern = ":",fixed = T,value=F,invert=F),]
#tableresults$estimates_normed <- as.vector(unlist(by(tableresults[,3],tableresults[,"species"],function(x) scale(x,center = F,scale = T))))
#tableresults$signif <- as.character(tableresults$signif)
tableresults <- tableresults[tableresults$signif=="*",] # Keep only significant values
tableresults$species <- paste0(tableresults$eco," ",tableresults$species)
tableresults <- tableresults[,c(1,2,3)] # remove conif and angio
tableresults$estimates <- round(tableresults$estimates,3)
#tableresults$eco <- as.character(tableresults$eco)
tableresults$variable <- as.character(tableresults$variable)
tableresults$variable <- sub(":Plotcat1"," : Marginality (TE)",tableresults$variable)
tableresults$variable <- sub("Plotcat1:","Marginality (TE) : ",tableresults$variable)
tableresults$variable <- sub(":Plotcat2"," : Marginality (LE)",tableresults$variable)
tableresults$variable <- sub("Plotcat2:","Marginality (LE) : ",tableresults$variable)
tableresults$variable <- sub("Plotcat2","Marginality (LE)",tableresults$variable)
tableresults$variable <- sub("Plotcat1","Marginality (TE)",tableresults$variable)
tableresults$variable <- sub("_climate_mean.30","",tableresults$variable)
tableresults$variable <- sub("_climate_mean.30","",tableresults$variable)
tableresults$variable <- sub("bio12","Precipitation",tableresults$variable) # 1
tableresults$variable <- sub("bio13","Precipitation",tableresults$variable)# 2
tableresults$variable <- sub("bio14","Precipitation",tableresults$variable)# 3
tableresults$variable <- sub("bio1","Temperature",tableresults$variable)# 1
tableresults$variable <- sub("bio5","Temperature",tableresults$variable)# 2
tableresults$variable <- sub("tmean.djf","Temperature",tableresults$variable)# 2
tableresults$variable <- sub("mean_spei12","Drought1",tableresults$variable)# 1
tableresults$variable <- sub("min_spei12","Drought2",tableresults$variable)# 2
tableresults$variable <- sub("ppet.mean","Drought3",tableresults$variable)# 3
tableresults$variable <- sub("sqrtBA.ha.plot.1","CompetitionINTER",tableresults$variable)# 1
tableresults$variable <- sub("sqrtBA.O.plot.1","CompetitionINTER",tableresults$variable)# 2
tableresults$variable <- sub("logBAj.plot.1","CompetitionINTRA",tableresults$variable)# 3
tableresults$variable <- sub("log","",tableresults$variable)
tableresults$variable <- sub("sqrt","",tableresults$variable)

tableresults$variable <- as.character(tableresults$variable)
datDF=split(tableresults[,],as.character(tableresults$species)) # split mon df
A=lapply(datDF, function(x) x[order(abs(x[,3]),decreasing=T),][1:14,])
All <- do.call(rbind,A)
All <- na.omit(All[,2])
All <- sub(" ","",All)# 3
All <- sub(" ","",All)# 3
All <- sub(" ","",All)# 3
All <- strsplit(All, ":",fixed = T)
All <- lapply(All, function(x) sort(x))
All <- lapply(All, function(x) paste0(x[1],":",x[2]))
All <- do.call(rbind,All)
sort(summary(as.factor(All))) # for our interpretation


A <- do.call(rbind,A)
A$Varank <- as.character(c(paste0("Estimate ",1:14)))
B <- unique(A$species,na.rm=T)
B <- B[-c(2)]
A$species <- rep(B,each = 14)
A <- recast(A, species ~ Varank + variable, measure.var = c("estimates", "variable"))
A
A <- A[,c("species",unlist(strsplit(paste0("Estimate ",1:14,"_variable:Estimate ",1:14,"_estimates"), ":",fixed = T)))]
write.table(A, sep=",", file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/CoefInter.binom.csv",row.names = F) # Table for the paper with coefficient
A


#### Idem ZTNB models 
tableresults <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/tableresultsNB.csv")
summary(as.factor(tableresults$variable))
tableresults <- tableresults[-c(which(tableresults$variable=="(Intercept)")),] # On enlève l'intercept
tableresults <- tableresults[grep(tableresults$variable,pattern = ":",fixed = T,value=F,invert=F),]
#tableresults$estimates_normed <- as.vector(unlist(by(tableresults[,3],tableresults[,"species"],function(x) scale(x,center = F,scale = T))))
#tableresults$signif <- as.character(tableresults$signif)
tableresults <- tableresults[tableresults$signif=="*",] # Keep only significant values
tableresults$species <- paste0(tableresults$eco," ",tableresults$species)
tableresults <- tableresults[,c(1,2,3)] # remove conif and angio
tableresults$estimates <- round(tableresults$estimates,3)
names(tableresults)[2] <- "variable" # supress the s
#tableresults$eco <- as.character(tableresults$eco)
tableresults$variable <- as.character(tableresults$variable)
tableresults$variable <- sub(":Plotcat1"," : Marginality (TE)",tableresults$variable)
tableresults$variable <- sub("Plotcat1:","Marginality (TE) : ",tableresults$variable)
tableresults$variable <- sub(":Plotcat2"," : Marginality (LE)",tableresults$variable)
tableresults$variable <- sub("Plotcat2:","Marginality (LE) : ",tableresults$variable)
tableresults$variable <- sub("Plotcat2","Marginality (LE)",tableresults$variable)
tableresults$variable <- sub("Plotcat1","Marginality (TE)",tableresults$variable)
tableresults$variable <- sub("_climate_mean.30","",tableresults$variable)
tableresults$variable <- sub("_climate_mean.30","",tableresults$variable)
tableresults$variable <- sub("bio12","Precipitation",tableresults$variable) # 1
tableresults$variable <- sub("bio13","Precipitation",tableresults$variable)# 2
tableresults$variable <- sub("bio14","Precipitation",tableresults$variable)# 3
tableresults$variable <- sub("bio1","Temperature",tableresults$variable)# 1
tableresults$variable <- sub("bio5","Temperature",tableresults$variable)# 2
tableresults$variable <- sub("tmean.djf","Temperature",tableresults$variable)# 2
tableresults$variable <- sub("mean_spei12","Drought1",tableresults$variable)# 1
tableresults$variable <- sub("min_spei12","Drought2",tableresults$variable)# 2
tableresults$variable <- sub("ppet.mean","Drought3",tableresults$variable)# 3
tableresults$variable <- sub("sqrtBA.ha.plot.1","CompetitionINTER",tableresults$variable)# 1
tableresults$variable <- sub("sqrtBA.O.plot.1","CompetitionINTER",tableresults$variable)# 2
tableresults$variable <- sub("logBAj.plot.1","CompetitionINTRA",tableresults$variable)# 3
tableresults$variable <- sub("log","",tableresults$variable)
tableresults$variable <- sub("sqrt","",tableresults$variable)

tableresults$variable <- as.character(tableresults$variable)
datDF=split(tableresults[,],as.character(tableresults$species)) # split mon df
A=lapply(datDF, function(x) x[order(abs(x[,3]),decreasing=T),][1:14,])
All <- do.call(rbind,A)
All <- na.omit(All[,2])
All <- sub(" ","",All)# 3
All <- sub(" ","",All)# 3
All <- sub(" ","",All)# 3
All <- strsplit(All, ":",fixed = T)
All <- lapply(All, function(x) sort(x))
All <- lapply(All, function(x) paste0(x[1],":",x[2]))
All <- do.call(rbind,All)
sort(summary(as.factor(All))) # for our interpretation

A <- do.call(rbind,A)
A$Varank <- as.character(c(paste0("Estimate ",1:14)))
B <- unique(A$species,na.rm=T)
B <- B[-c(2)]
A$species <- rep(B,each = 14)
A <- recast(A, species ~ Varank + variable, measure.var = c("estimates", "variable"))
A
A <- A[,c("species",unlist(strsplit(paste0("Estimate ",1:14,"_variable:Estimate ",1:14,"_estimates"), ":",fixed = T)))]
write.table(A, sep=",", file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/CoefInter.Negbin.csv",row.names = F) # Table for the paper with coefficient
A


tableresults <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/1_table_bin_gene.csv")
tableresults <- tableresults[,-c(6)]
tableresults[3:7] <- round(tableresults[,3:7],2)
write.table(tableresults, sep=",", file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/binomallmodels.csv",row.names = F)


tableresults <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/3_table_negbin_gene.csv")
tableresults <- tableresults[,-c(7)]
tableresults[,3:8] <- round(tableresults[,3:8],2)
write.table(tableresults, sep=",", file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/negbinmodels.csv",row.names = F)






