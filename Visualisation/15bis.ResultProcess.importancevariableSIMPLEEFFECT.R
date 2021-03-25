library(reshape2)

rm(list = ls())
gc()
getwd()

tableresults <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/tableresults.csv")
summary(as.factor(tableresults$variable))
#tableresults$estimates_normed <- as.vector(unlist(by(tableresults[,3],tableresults[,"species"],function(x) scale(x,center = F,scale = T))))
tableresults <- tableresults[grep(tableresults$variable,pattern = ":",fixed = T,value=F,invert=T),]
#tableresults$signif <- as.character(tableresults$signif)
tableresults <- tableresults[tableresults$signif=="*",] # Keep only significant values
tableresults$species <- paste0(tableresults$eco," ",tableresults$species)
tableresults <- tableresults[,c(1,2,3)] # remove conif and angio
tableresults$estimates <- round(tableresults$estimates,3)
#tableresults$eco <- as.character(tableresults$eco)
tableresults$variable <- as.character(tableresults$variable)
datDF=split(tableresults[,],as.character(tableresults$species)) # split mon df
A=lapply(datDF, function(x) x[order(abs(x[,3]),decreasing=T),][1:14,])
All <- do.call(rbind,A)
All <- na.omit(All[,2])
summary(as.factor(All)) # for our interpretation


A <- do.call(rbind,A)
A$Varank <- as.character(c(paste0("Estimate ",1:14)))
B <- unique(A$species,na.rm=T)
B <- B[-c(2)]
A$species <- rep(B,each = 14)
A <- recast(A, species ~ Varank + variable, measure.var = c("estimates", "variable"))
A
A <- A[,c("species",unlist(strsplit(paste0("Estimate ",1:14,"_variable:Estimate ",1:14,"_estimates"), ":",fixed = T)))]
write.table(A, sep=",", file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/CoefSimple.binom.csv",row.names = F) # Table for the paper with coefficient
A

#### Idem ZTNB models 
tableresults <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/tableresultsNB.csv")
summary(as.factor(tableresults$variable))
#tableresults$estimates_normed <- as.vector(unlist(by(tableresults[,3],tableresults[,"species"],function(x) scale(x,center = F,scale = T))))
tableresults <- tableresults[grep(tableresults$variable,pattern = ":",fixed = T,value=F,invert=T),]
#tableresults$signif <- as.character(tableresults$signif)
tableresults <- tableresults[tableresults$signif=="*",] # Keep only significant values
tableresults$species <- paste0(tableresults$eco," ",tableresults$species)
tableresults <- tableresults[,c(1,2,3)] # remove conif and angio
tableresults$estimates <- round(tableresults$estimates,3)
#tableresults$eco <- as.character(tableresults$eco)
tableresults$variable <- as.character(tableresults$variable)

# Lines to sub names 

datDF=split(tableresults[,],as.character(tableresults$species)) # split mon df
A=lapply(datDF, function(x) x[order(abs(x[,3]),decreasing=T),][1:14,])
All <- do.call(rbind,A)
All <- na.omit(All[,2])
All <- as.character(All)
summary(as.factor(All)) # for our interpretation

A <- do.call(rbind,A)
A$Varank <- as.character(c(paste0("Estimate ",1:14)))
B <- unique(A$species,na.rm=T)
B <- B[-c(2)]
A$species <- rep(B,each = 14)
A <- recast(A, species ~ Varank + variable, measure.var = c("estimates", "variable"))
A
A <- A[,c("species",unlist(strsplit(paste0("Estimate ",1:14,"_variable:Estimate ",1:14,"_estimates"), ":",fixed = T)))]
write.table(A, sep=",", file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/CoefSimple.negbin.csv",row.names = F)



tableresults <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/1_table_bin_gene.csv")
tableresults <- tableresults[,-c(6)]
tableresults[3:7] <- round(tableresults[,3:7],2)
write.table(tableresults, sep=",", file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/binomallmodels.csv",row.names = F)


tableresults <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/3_table_negbin_gene.csv")
tableresults <- tableresults[,-c(7)]
tableresults[,3:8] <- round(tableresults[,3:8],2)
write.table(tableresults, sep=",", file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/negbinmodels.csv",row.names = F)






