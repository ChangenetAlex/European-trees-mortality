library(reshape2)
rm(list = ls())
gc()
getwd()

tableresults <- read.csv("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/species/Results.all/tableresults.csv")
tableresults <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/tableresultsNB.csv")
colnames(tableresults)[2] <- "variable"
summary(as.factor(tableresults$variable))

## Simple effects 
tableresults <- tableresults[tableresults$signif=="*",] # Keep only significant values
tableresults <- tableresults[,c(1,2,6,8)] # remove conif and angio
tableresults$t.value <- round(tableresults$t.value,3)
tableresults$SE <- round(tableresults$SE,3)

unique(tableresults$species)

tableresults$species <- sub("r","",tableresults$species) # remove the r indication
tableresults$species <- sub(" ","",tableresults$species) # remove the " " indication
unique(tableresults$species)

tableresults$variable <- as.character(tableresults$variable)

# Change names of interaction that are in fact the same:
tableresults$variable <- sub(" ","",tableresults$variable) # remove the " " indication
tableresults$variable <- sub(" ","",tableresults$variable) # remove the " " indication
tableresults$variable <- sub(" ","",tableresults$variable) # remove the " " indication
unique(tableresults$variable)

MargeVar <- c("mean_spei12","bio1_climate_mean.30","bio5_climate_mean.30","bio12_climate_mean.30","bio13_climate_mean.30","bio14_climate_mean.30",
                "bio1_climate_mean.30","tmean.djf_climate_mean.30","logBAj.plot.1","sqrtBA.ha.plot.1","logBAj.plot.1")
Marge <- c("Plotcat2","Plotcat1")
for (i in 1:length(Marge)){
  for (j in 1:length(MargeVar)){
    tableresults$variable <- sub(paste0(MargeVar[j],":",Marge[i]),paste0(Marge[i],":",MargeVar[j]),tableresults$variable)
    print(c(paste0(MargeVar[j],":",Marge[i]),paste0(Marge[i],":",MargeVar[j])))
  }
}

# Keep only interaction with margins + simple effect 
unique(tableresults$variable)
A <- NULL
for (i in 1:length(Marge)){
  for (j in 1:length(MargeVar)){
  A[length(MargeVar)*(i-1)+j] <- paste0(Marge[i],":",MargeVar[j])
  
    }
}

# Simple effect 
unique(tableresults$species)
Abis <- unique(tableresults$variable[grep(tableresults$variable,pattern = ":",fixed = T,value=F,invert=T)])

# All
Kept <- c(Abis,A)

# Keep the variable that we want (All or just simple effect or just interaction)
#tableresults <- tableresults[tableresults$variable%in%Kept,]
#tableresults <- tableresults[tableresults$variable%in%Abis,]
tableresults <- tableresults[tableresults$variable%in%A,]
unique(tableresults$species)
tableresults <- tableresults[,-c(4)]
unique(tableresults$species)

# Transform the name
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
tableresults$variable <- sub("mean_spei12","Drought mean",tableresults$variable)# 1
tableresults$variable <- sub("min_spei12","Drought min",tableresults$variable)# 2
tableresults$variable <- sub("ppet.mean","Waterbalance",tableresults$variable)# 3
tableresults$variable <- sub("sqrtBA.ha.plot.1","CompetitionTOTAL",tableresults$variable)# 1
tableresults$variable <- sub("sqrtBA.O.plot.1","CompetitionINTER",tableresults$variable)# 2
tableresults$variable <- sub("logBAj.plot.1","CompetitionINTRA",tableresults$variable)# 3
tableresults$variable <- sub("log","",tableresults$variable)
tableresults$variable <- sub("sqrt","",tableresults$variable)
tableresults$variable <- sub("BAIj.plot.1.mean","Species growth (plot)",tableresults$variable)# 3
tableresults$variable <- sub("BAIj.plot.1","Total species growth",tableresults$variable)# 3
tableresults$variable <- sub("dbh.plot.mean","dbh mean (plot)",tableresults$variable)# 3
tableresults$variable <- sub("treeNbr","density",tableresults$variable)# 3
tableresults$variable <- sub("yearsbetweensurveys","Census interval",tableresults$variable)# 3


long1 <- melt(tableresults, id=c("species", "variable", "t.value")) # Melt the table first
wide1 <- dcast(long1,species~variable) # Rebuild the table 
wide1[is.na(wide1)] <- "-"
colnames(wide1)
wide1 <- wide1[,-c(9,11)]
write.table(wide1, sep=",", file="/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/TableINTER.ZT.csv",row.names = F) # Table for the paper with coefficient
write.table(wide1, sep=",", file="/home/achangenet/SynologyDrive/FUNDIV - NFI - Europe/Redaction/Paper1/TableSIMPLE.Binom.csv",row.names = F) # Table for the paper with coefficient


