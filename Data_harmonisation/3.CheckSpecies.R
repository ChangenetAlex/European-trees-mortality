library(dplyr)
Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")

#Edit 05/02/2018
# To run script after have made correction at step 1 or 2. Go directly to script 4. This one is just the creation of the database for species

#########################################################
##              Load the databases                     ##
#########################################################

#Those are the two lists on which we are going to work

binomial.fr <- read.delim(paste0(Dir,"our-data/species_list.txt"))
SpeciesList <- read.csv("~/Documents/FUNDIV - NFI - Europe/FunDivEUROPE_Inventory_data_7Oct/FunDivEUROPE_species.csv")

length(unique(binomial$espar))
length(unique(SpeciesList$id))

#Those are the others databases (french and Fundiv)

load(paste0(Dir,"our-data/Fundiv_alltree_allmetric.RData")) #This is the clean database with all the new calculated metrics
load(paste0(Dir,"our-data/Fundiv_alltree_allmetric.2018.RData")) #This is the clean database with all the new calculated metrics 2018 (BAIJ corrected)
load(paste0(Dir,"our-data/tree_data_IFN_france_with_competition-indices.bai_ha_forFUNDIVcomp.RData"))
ifn.fr$espar <- as.character(ifn.fr$espar)  
length(unique(ftp.short$speciesid)) #172
length(unique(ifn.fr$espar)) #170

#And this is the final  we are suposed to obtain 
load(paste0(Dir,"our-data/tree_data_harmonise_fr.fundiv.indexes.RData")) #Table final that we are trying to obain
length(unique(tree$speciesid)) #141

#########################################################
######    Solves problems in the IFN.fr databse:   ######
#########################################################


# 1 # Here we look the individuals without identification in the two databases
ifn.fr2 <- ifn.fr[is.na(ifn.fr$espar)|(!is.na(ifn.fr$espar)&!is.na(ifn.fr$binomial)),] #partie des données correctes ou sans code 
A <- ifn.fr[is.na(ifn.fr$espar)|is.na(ifn.fr$binomial),] # 11882 without name or species ID
A <-(ifn.fr[is.na(ifn.fr$espar),]) # among which 8168 individuals without names nor ID but counted for the mortality rate and Sp.mortl rate(in ifn.fr)
A <- ifn.fr[!is.na(ifn.fr$espar)&is.na(ifn.fr$binomial),] #3714 individuals that are probably mistakes for which just the binomial is missing (No Lot of them are species id = 6 instead of 06 and so on...
#Problem with 3 9 2 6 5 7 4 and 2.50E+06
ifn.fr[is.na(ifn.fr$espar)&!is.na(ifn.fr$binomial),] # 0 with a name but no SpID

# 2 # Now we look how many names there are in the ifn.fr and how much of them can be found in the specieslist
length(unique(ifn.fr$espar)) #170 unique codes (but mistakes)
#When we look we can see 
length(unique(ifn.fr$binomial)) #152 unique names

# 3 # Need to fix this by merging the part with the problem (3714 indiv)
ifn.fr3 <- ifn.fr[!is.na(ifn.fr$espar)&is.na(ifn.fr$binomial),] #3714 indiv for which names
#ifn.fr3 <- semi_join(ifn.fr3,binomial, by = "espar")
# 2 espèces étant juste pas marquées : 25E3 et 25E5 mais figurent pourtant dans le tableau 15 individuals
#ifn.fr3 <- anti_join(ifn.fr3,binomial, by = "espar") #3699 indiv
summary(ifn.fr3$espar)
ifn.fr3$espar <- as.character(ifn.fr3$espar)
ifn.fr3[ifn.fr3$espar=="2",6] <- "02"
ifn.fr3[ifn.fr3$espar=="3",6] <- "03"
ifn.fr3[ifn.fr3$espar=="5",6] <- "05"
ifn.fr3[ifn.fr3$espar=="9",6] <- "09"
ifn.fr3[ifn.fr3$espar=="6",6] <- "06"
ifn.fr3[ifn.fr3$espar=="7",6] <- "07"
ifn.fr3[ifn.fr3$espar=="4",6] <- "04"
ifn.fr3[ifn.fr3$espar=="29FI",6] = "29fi"

ifn.fr3$binomial <- binomial.fr$binomial[match(ifn.fr3$espar, binomial.fr$espar)]
#Cette commande permet de regarder ce qui match sur ue variable (ici espar) et de remplir une autre colonne correspondante (binomial)
A<- ifn.fr3[is.na(ifn.fr3$binomial),] #Reste 8 sur 3714 sans nom. Et les 8168 sans l'un ni l'autre. 

# 4 ### Merge de nouveau cette partie du tableau avec les autres
ifn.fr.renamed <- rbind(ifn.fr2,ifn.fr3)
length(unique(ifn.fr.renamed$espar)) #163 unique codes (but mistakes) (-7 codes). 4 erreurs a noter + NA a remove + quelques exemplaires
#When we look if all codes are found in the list 

ifn.fr.renamed$espar <- as.numeric(ifn.fr.renamed$espar)
A<- ifn.fr.renamed[is.na(ifn.fr.renamed$espar),]
test <- semi_join(ifn.fr.renamed,binomial.fr,by="espar")#772626 Ce que j'explique avec mon tableau
test <- anti_join(treefinal,Allspfinal,by="espar")#8176 Ce que je n'explique pas (pas de code)
test <-semi_join(ifn.fr.renamed,binomial.fr,by="binomial") #763055 avec un nom
test <-anti_join(ifn.fr.renamed,binomial.fr,by="binomial") #17747 (pas de nom)

test <-anti_join(binomial.fr,ifn.fr.renamed,,by="binomial") #51 espèces que l'on utilise pas
test <-semi_join(binomial.fr,ifn.fr.renamed,,by="binomial") #131 espèces que l'on utiliseras comme premier tableau
test <-anti_join(binomial.fr,ifn.fr.renamed,,by="espar") #24 code que l'on utilise pas
AllSpecies.fr <-semi_join(binomial.fr,ifn.fr.renamed,,by="espar") #158 codes que l'on utiliseras comme premier tableau !!!!!
#27 noms qui sont des noms en français a remplacer par le nom latin avant de reatribuer les noms d'espèces dans la table
# et cheker les 5 derniers codes (erreurs et NA)

length(unique(ifn.fr.renamed$binomial)) #155 unique names 
length(unique(ifn.fr.renamed$espar)) #163 unique code

#Need to still be working on this data base
save(ifn.fr.renamed,file="/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/ifn.fr.renamed.RData")
save(AllSpecies.fr,file="/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/AllSpecies.fr.RData")



#########################################################
######   Solves problems in the Fundiv database:   ######
#########################################################


SpeciesList <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/FunDivEUROPE_speciesbis.csv",stringsAsFactors = FALSE) #document corrigé
load(paste0(Dir,"our-data/Fundiv_alltree_allmetric.RData"))

B <-(ftp.short[is.na(ftp.short$speciesid),]) # 0 obs without names or spId in fundiv.data
colnames(SpeciesList)[1] <- "speciesid" #rename id in species id
length(unique(ftp.short$speciesid)) #172 uique codes
test2 <- semi_join(ftp.short,SpeciesList, by = "speciesid") # 100% 1351328 rows have a code
AllSpecies.Fundiv <- semi_join(SpeciesList,ftp.short, by = "speciesid") # The 172 names are found in the database. This is the second part !!!
#This is the second part of our global database


SpeciesList[which(SpeciesList$species==""),6] <- "sp." #replace empty space by sp.
binomial=paste0(SpeciesList$genus,"_",SpeciesList$species) #create new variable
SpeciesList[,"binomial"]=binomial


AllSpecies.Fundiv[which(AllSpecies.Fundiv$species==""),6] <- "sp." #replace empty space by sp.
fix(AllSpecies.Fundiv) #fix manually and replace white spaces by other conif or oaks and so on
binomial=paste0(AllSpecies.Fundiv$genus,"_",AllSpecies.Fundiv$species) #create new variable
AllSpecies.Fundiv[,"binomial"]=binomial

save(AllSpecies.Fundiv,file="/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/AllSpecies.Fundiv.RData")


#########################################################
######     List all species in a common table      ######
#########################################################

load(paste0(Dir,"our-data/AllSpecies.Fundiv.RData"))
load(paste0(Dir,"our-data/AllSpecies.fr.RData"))
AllSpecies.Fundiv = AllSpecies.Fundiv[order(AllSpecies.Fundiv$binomial),]
fix(AllSpecies.Fundiv)#ajout turbinata et laricio pour pinus et juniperus 
save(AllSpecies.Fundiv,file="/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/AllSpecies.Fundiv.RData") ###good table

AllSpecies.fr = AllSpecies.fr[order(AllSpecies.fr$binomial),]
write.csv(AllSpecies.fr,file=paste0(Dir,"our-data/AllSpecies.fr.csv"))
AllSpecies.Fundiv = AllSpecies.Fundiv[order(AllSpecies.Fundiv$binomial),]
write.csv(AllSpecies.Fundiv,file=paste0(Dir,"our-data/AllSpecies.Fundiv.csv"))
binomial.fr <- read.delim(paste0(Dir,"our-data/species_list.txt"))
binomial.fr = binomial.fr[order(binomial.fr$binomial),]
write.csv(binomial.fr,file=paste0(Dir,"our-data/binomial.fr.csv"))
SpeciesList = SpeciesList[order(SpeciesList$binomial),]
write.csv(SpeciesList,file=paste0(Dir,"our-data/SpeciesList.csv"))
#Sortir les 4 en tableau 
#Work on the Allspeciesfr => Ajout de multiple noms latins. 
#Juniperus oxycedrus et communis instead of 1 et 2 et bosser sur les 4 pinus
#renommé AllSpecies.fr.bis
AllSpecies.fr <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/AllSpecies.fr.bis.csv", header=TRUE) #good table

length(unique(AllSpecies.fr$espar))
length(unique(AllSpecies.fr$binomial))
#158 code VS 156 names
length(unique(AllSpecies.Fundiv$speciesid))
length(unique(AllSpecies.Fundiv$binomial))
#172 codes VS 172 names 
#probleme pinus nigra * 4 et fois 


#Ici test des noms en commun ou pas. 
#Technique la miuex : extraire les deux tableaux arranger les colonnes dans l'odre. Choper les noms fundiv pour les espèces fr commune
#mettre le tout ensemble et comparer avec la vraie liste fundiv. rajoutd'une colonne pour des précisions (avec 3 options a b ou c cf rmarkdown pour détails)
AllSpecies.fr <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/AllSpecies.fr.bis.csv", header=TRUE) #good table
AllSpecies.fr <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/AllSpecies.Fundiv.csv", header=TRUE) #good table

AllSpecies.fr$code <- NA
AllSpecies.fr$rank <- NA
AllSpecies.fr$speciesid <- NA
names(AllSpecies.fr)
Allsp.fr <- subset(AllSpecies.fr, select=c("speciesid","code","espar","nom","binomial","rank","origin"))
colnames(Allsp.fr)[4] <- "acceptedname"

names(AllSpecies.Fundiv)
AllSpecies.Fundiv$espar <- NA
AllSpecies.Fundiv$origin <- NA
Allsp.Fundiv <- subset(AllSpecies.Fundiv, select=c("speciesid","code","espar","acceptedname","binomial","rank","origin"))

Allsp<- rbind(Allsp.fr,Allsp.Fundiv) #table a trier avec toutes les espèces
Allsp <- Allsp[order(Allsp$binomial),]
write.csv(Allsp,file=paste0(Dir,"our-data/Allsp.csv")) #Table with all the species explaining all our data. This table must be checked if we want to analyse just the french data for instance. 
#But also the AllSpecies.Fundiv & AllSpecies.fr

#Modification by hand (see details) => Allspfinal.csv

test <- semi_join(Allspecies.fr,AllSpecies.Fundiv, by = "binomial") #97 en commun
test2 <- semi_join(AllSpecies.Fundiv,Allspecies.fr, by = "binomial") #93 en commun

test <- anti_join(Allspecies.fr,AllSpecies.Fundiv, by = "binomial") #97 en commun
test2 <- anti_join(AllSpecies.Fundiv,Allspecies.fr, by = "binomial") #93 en commun

test <- inner_join(AllSpecies.fr,SpeciesList, by = "binomial") #84 en commun
test <- left_join(AllSpecies.fr,SpeciesList, by = "binomial") #84 en commun


test <- anti_join(AllSpecies.fr,AllSpecies.Fundiv, by = "binomial") #93
test3 <- anti_join(AllSpecies.Fundiv,AllSpecies.fr, by = "binomial")


AllSpecies.fr <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/Allspecies.fr.csv", header=TRUE)


