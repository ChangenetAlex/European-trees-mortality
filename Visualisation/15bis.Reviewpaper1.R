library(reshape2)
rm(list = ls())
gc()
getwd()

tableresults <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Results.all/tableresults.csv") # BIN
tableresults <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Results.all/tableresultsNB.csv") #ZT
colnames(tableresults)[2] <- "variable"
summary(as.factor(tableresults$variable))
tableresults <- tableresults[-c(which(tableresults$variable=="(Intercept)")),] # On enlève l'intercept
tableresults <- tableresults[grep(tableresults$variable,pattern = ":",fixed = T,value=F,invert=F),] ## keep onluy interactions 
tableresults <- tableresults[grep(tableresults$variable,pattern = "Plotcat0",invert=T),]

#tableresults$estimates_normed <- as.vector(unlist(by(tableresults[,3],tableresults[,"species"],function(x) scale(x,center = F,scale = T))))
#tableresults$signif <- as.character(tableresults$signif)


#tableresults <- tableresults[tableresults$signif=="*",] # Keep only significant values
tableresults <- tableresults[,-c(4,5)] # remove conif and angio
tableresults$estimates <- round(tableresults$estimates,2)
tableresults$t.value <- round(tableresults$t.value,2)
tableresults$SE <- round(tableresults$SE,2)

tableresults$species <- sub("r","",tableresults$species) # remove the r indication
tableresults$species <- sub(" ","",tableresults$species) # remove the " " indication
unique(tableresults$species)


tableresults$Values <- paste0(tableresults$estimates,"\n(",tableresults$SE,")\n",tableresults$t.value," ",tableresults$signif)
tableresults <- tableresults[,c(1,2,7)]
unique(tableresults$variable)
tableresults$variable <- as.character(tableresults$variable)


unique(tableresults$variable)


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
tableresults$variable <- sub("mean_spei12","meanSPEI",tableresults$variable)# 1
tableresults$variable <- sub("min_spei12","minSPEI",tableresults$variable)# 2
tableresults$variable <- sub("ppet.mean","Precipitation",tableresults$variable)# 3
tableresults$variable <- sub("sqrtBA.ha.plot.1","BA",tableresults$variable)# 1
tableresults$variable <- sub("sqrtBA.O.plot.1","BA.O",tableresults$variable)# 2
tableresults$variable <- sub("logBAj.plot.1","BAj",tableresults$variable)# 3
tableresults$variable <- sub("log","",tableresults$variable)
tableresults$variable <- sub("sqrt","",tableresults$variable)
tableresults$variable <- sub(" : ",":",tableresults$variable)
tableresults$variable <- sub(":"," X ",tableresults$variable)


tableresults$variable <- as.character(tableresults$variable)
unique(tableresults$variable)

tableresults[tableresults$variable=="Temperature X Marginality (TE)","variable"] <-  "Marginality (TE) X Temperature"
tableresults[tableresults$variable=="Precipitation X Marginality (TE)","variable"] <-  "Marginality (TE) X Precipitation"
tableresults[tableresults$variable=="minSPEI X Marginality (TE)","variable"] <-  "Marginality (TE) X minSPEI"
tableresults[tableresults$variable=="meanSPEI X Marginality (TE)","variable"] <-  "Marginality (TE) X meanSPEI"
tableresults[tableresults$variable=="BA X Marginality (TE)","variable"] <-  "Marginality (TE) X BA"
tableresults[tableresults$variable=="BAj X Marginality (TE)","variable"] <-  "Marginality (TE) X BAj"
tableresults[tableresults$variable=="BA.O X Marginality (TE)","variable"] <-  "Marginality (TE) X BA.O"

tableresults[tableresults$variable=="Temperature X Marginality (LE)","variable"] <-  "Marginality (LE) X Temperature"
tableresults[tableresults$variable=="Precipitation X Marginality (LE)","variable"] <-  "Marginality (LE) X Precipitation"
tableresults[tableresults$variable=="minSPEI X Marginality (LE)","variable"] <-  "Marginality (LE) X minSPEI"
tableresults[tableresults$variable=="meanSPEI X Marginality (LE)","variable"] <-  "Marginality (LE) X meanSPEI"
tableresults[tableresults$variable=="BA X Marginality (LE)","variable"] <-  "Marginality (LE) X BA"
tableresults[tableresults$variable=="BAj X Marginality (LE)","variable"] <-  "Marginality (LE) X BAj"
tableresults[tableresults$variable=="BA.O X Marginality (LE)","variable"] <-  "Marginality (LE) X BA.O"

tableresults[tableresults$variable=="Precipitation X Temperature","variable"] <-  "Temperature X Precipitation"
tableresults[tableresults$variable=="minSPEI X Temperature","variable"] <-  "Temperature X minSPEI"
tableresults[tableresults$variable=="meanSPEI X Temperature","variable"] <-  "Temperature X meanSPEI"
tableresults[tableresults$variable=="BA X Temperature","variable"] <-  "Temperature X BA"
tableresults[tableresults$variable=="BAj X Temperature","variable"] <-  "Temperature X BAj"
tableresults[tableresults$variable=="BA.O X Temperature","variable"] <-  "Temperature X BA.O"

tableresults[tableresults$variable=="minSPEI X Precipitation","variable"] <-  "Precipitation X minSPEI"
tableresults[tableresults$variable=="meanSPEI X Precipitation","variable"] <-  "Precipitation X meanSPEI"
tableresults[tableresults$variable=="BA X Precipitation","variable"] <-  "Precipitation X BA"
tableresults[tableresults$variable=="BAj X Precipitation","variable"] <-  "Precipitation X BAj"
tableresults[tableresults$variable=="BA.O X Precipitation","variable"] <-  "Precipitation X BA.O"

tableresults[tableresults$variable=="minSPEI X meanSPEI","variable"] <-  "meanSPEI X minSPEI"
tableresults[tableresults$variable=="BA X meanSPEI","variable"] <-  "meanSPEI X BA"
tableresults[tableresults$variable=="BAj X meanSPEI","variable"] <-  "meanSPEI X BAj"
tableresults[tableresults$variable=="BA.O X meanSPEI","variable"] <-  "meanSPEI X BA.O"

tableresults[tableresults$variable=="minSPEI X BA","variable"] <-  "BA X minSPEI"
tableresults[tableresults$variable=="BA.O X BA","variable"] <-  "BA X BA.O"
tableresults[tableresults$variable=="BAj X BA","variable"] <-  "BA X BAj"

tableresults[tableresults$variable=="minSPEI X BAj","variable"] <-  "BAj X minSPEI"
tableresults[tableresults$variable=="BA.O X BAj","variable"] <-  "BAj X BA.O"

tableresults[tableresults$variable=="minSPEI X BA.O","variable"] <-  "BA.O X minSPEI"

long1 <- melt(tableresults, id=c("species", "variable", "Values")) # Melt the table first
wide1 <- dcast(long1,species~variable) # Rebuild the table 
wide1[is.na(wide1)] <- "-"
#write.table(wide1, sep=",", file="~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/ZT.Table.Inter.Review.csv",row.names = F) # Table for the paper with coefficient
write.table(wide1, sep=",", file="~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/BIN2.Table.Inter.Review.csv",row.names = F) # Table for the paper with coefficient



### effet simples
library(reshape2)
rm(list = ls())
gc()
getwd()

tableresults <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Results.all/tableresults.csv") # BIN
tableresults <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Results.all/tableresultsNB.csv") #ZT
colnames(tableresults)[2] <- "variable"
summary(as.factor(tableresults$variable))
tableresults <- tableresults[-c(which(tableresults$variable=="(Intercept)")),] # On enlève l'intercept
tableresults <- tableresults[grep(tableresults$variable,pattern = ":",fixed = T,value=F,invert=T),] ## keep onluy simple
tableresults <- tableresults[grep(tableresults$variable,pattern = "Plotcat0",invert=T),]

#tableresults$estimates_normed <- as.vector(unlist(by(tableresults[,3],tableresults[,"species"],function(x) scale(x,center = F,scale = T))))
#tableresults$signif <- as.character(tableresults$signif)


#tableresults <- tableresults[tableresults$signif=="*",] # Keep only significant values
tableresults <- tableresults[,-c(4,5)] # remove conif and angio
tableresults$estimates <- round(tableresults$estimates,2)
tableresults$t.value <- round(tableresults$t.value,2)
tableresults$SE <- round(tableresults$SE,2)

tableresults$species <- sub("r","",tableresults$species) # remove the r indication
tableresults$species <- sub(" ","",tableresults$species) # remove the " " indication
unique(tableresults$species)


tableresults$Values <- paste0(tableresults$estimates,"\n(",tableresults$SE,")\n",tableresults$t.value," ",tableresults$signif)
tableresults <- tableresults[,c(1,2,7)]
unique(tableresults$variable)
tableresults$variable <- as.character(tableresults$variable)

unique(tableresults$variable)


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
tableresults$variable <- sub("mean_spei12","meanSPEI",tableresults$variable)# 1
tableresults$variable <- sub("min_spei12","minSPEI",tableresults$variable)# 2
tableresults$variable <- sub("ppet.mean","Precipitation",tableresults$variable)# 3
tableresults$variable <- sub("sqrtBA.ha.plot.1","BA",tableresults$variable)# 1
tableresults$variable <- sub("sqrtBA.O.plot.1","BA.O",tableresults$variable)# 2
tableresults$variable <- sub("logBAj.plot.1","BAj",tableresults$variable)# 3
tableresults$variable <- sub("log","",tableresults$variable)
tableresults$variable <- sub("sqrt","",tableresults$variable)
tableresults$variable <- sub(" : ",":",tableresults$variable)
tableresults$variable <- sub(":"," X ",tableresults$variable)


tableresults$variable <- as.character(tableresults$variable)
unique(tableresults$variable)

long1 <- melt(tableresults, id=c("species", "variable", "Values")) # Melt the table first
wide1 <- dcast(long1,species~variable) # Rebuild the table 
wide1[is.na(wide1)] <- "-"
write.table(wide1, sep=",", file="~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/BIN2.Table.Simple.Review.csv",row.names = F) # Table for the paper with coefficient
#write.table(wide1, sep=",", file="~/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/ZT.Table.Simple.Review.csv",row.names = F) # Table for the paper with coefficient





