rm(list = ls())
gc()
#Script to do the calculation again for all the variable that interest us with the fundivData
Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
setwd(Dir)
load(paste0(Dir,"our-data/Fundiv_alltree_allmetric.2018.RData")) #Added 5/02/2018 #This is the clean database with all the new calculated metrics
load(paste0(Dir,"our-data/tfinal.biotic.June2018.RData"))

# Bind names for clusters
ftp.short[,"clusters"] <- paste0(ftp.short$country,ftp.short$cluster)
length(unique(ftp.short$cluster)) # Fewcluster names in common (13834 VS 14307 unique clusters in reality)
length(unique(ftp.short$clusters))
nrow(unique(ftp.short[c("cluster", "country")])) #14307

tfinal.biotic$clusters <- ftp.short$clusters[match(tfinal.biotic$treecode, ftp.short$treecode,incomparables = NA)]
save(tfinal.biotic,file="/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/tfinal.biotic.Feb2019.RData")

