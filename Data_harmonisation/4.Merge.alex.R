library(dplyr)
Dir <- c("/home/achangenet/Documents/FUNDIV - NFI - Europe/")
Dir <- c("/Users/alexandrechangenet/Dropbox/FUNDIV/")
load(paste0(Dir,"our-data/tree_data_harmonise_fr.fundiv.indexes.RData")) #tableau obtenu avant 1953568 lines
#load(paste0(Dir,"our-data/tree_data_IFN_france_with_competition-indices.bai_ha_forFUNDIVcomp.RData")) #IFN et metrics
load(paste0(Dir,"our-data/Fundiv_alltree_allmetric.2018.RData")) #Added 5/02/2018 #This is the clean database with all the new calculated metrics
load(paste0(Dir,"our-data/ifn.fr.metric.2018.RData")) #IFN with 3500 names corrected
Allspfinal <- read.csv2(paste0(Dir,"/our-data/Allspfinal.csv"), header=TRUE)




#Ajout d'un treecode et d'un country #modif 20/02/2018
treecode=paste0("FR",ifn.fr[,1],"_",ifn.fr[,5])
ifn.fr[,"treecode"]=treecode
ifn.fr$country <- "FR"
str(ifn.fr)
str(ftp.short)
#Ajout du BAI fr à Fundiv (ramener à l'année mais pas a l'hectare)

BAI=(ftp.short$ba2-ftp.short$ba1)/ftp.short$yearsbetweensurveys
ftp.short[,"BAI"]=BAI


kk <-subset(ftp.short, country=="FR") #Check no french data in the dunfiv database

#harmonize names of ifn.fr & fundiv.tree
ifn.fr.renamed <- ifn.fr

length(unique(ifn.fr.renamed$espar))
length(unique(ifn.fr.renamed$binomial))
#Correct problem de translation during the conversion
ifn.fr.renamed$espar <- as.character(ifn.fr.renamed$espar)
ifn.fr.renamed$binomial <- as.character(ifn.fr.renamed$binomial)
Allspfinal$espar <- as.character(Allspfinal$espar)
ifn.fr.renamed[!is.na(ifn.fr.renamed$espar)&ifn.fr.renamed$espar=="2,50E+06",7] <- "Salix_pentandra"
ifn.fr.renamed[!is.na(ifn.fr.renamed$binomial)&ifn.fr.renamed$binomial=="Salix_pentandra",6] <- "2500000"
ifn.fr.renamed[!is.na(ifn.fr.renamed$binomial)&ifn.fr.renamed$binomial=="Salix_trianda",6] <- "25000"

Allspfinal[!is.na(Allspfinal$espar)&Allspfinal$espar=="2",3] <- "02"
Allspfinal[!is.na(Allspfinal$espar)&Allspfinal$espar=="3",3] <- "03"
Allspfinal[!is.na(Allspfinal$espar)&Allspfinal$espar=="4",3] <- "04"
Allspfinal[!is.na(Allspfinal$espar)&Allspfinal$espar=="5",3] <- "05"
Allspfinal[!is.na(Allspfinal$espar)&Allspfinal$espar=="6",3] <- "06"
Allspfinal[!is.na(Allspfinal$espar)&Allspfinal$espar=="7",3] <- "07"
Allspfinal[!is.na(Allspfinal$espar)&Allspfinal$espar=="8",3] <- "08"
Allspfinal[!is.na(Allspfinal$espar)&Allspfinal$espar=="9",3] <- "09"

length(unique(ifn.fr.renamed$espar))  #157
length(unique(ifn.fr.renamed$binomial)) #151
length(unique(Allspfinal$espar)) # 159
length(unique(Allspfinal$binomial))
length(unique(Allspfinal$code))
length(unique(Allspfinal$speciesid))


#Here are the control 

test <- semi_join(ifn.fr.renamed,Allspfinal,by="espar")
test <- anti_join(ifn.fr.renamed,Allspfinal,by="espar")#All but 2 data are not explained by our data base species.

ftp.short$speciesid <- as.character(ftp.short$speciesid)
test <- semi_join(ftp.short,Allspfinal,by="speciesid")
test <- anti_join(ftp.short,Allspfinal,by="speciesid") #All are also explained for the fundiv data

colnames(ftp.short)
colnames(ifn.fr.renamed)
cor.test(ftp.short$ba_ha,ftp.short$BA.ha.plot)
cor.test(ftp.short$ba_ha,ftp.short$BA.ha.plot2)
cor.test(ftp.short$BA.ha.plot,ftp.short$BA.ha.plot2)

ftp.short$BAnei <- ftp.short$ba_ha - ftp.short$ba_ha2 #Remplacement de ba.haplot( calculated by me) by the real ba_ha (calcultaed buy sophie)

#########################################################
######   Put all the columns in the same order     ######
#########################################################

ftp.short$espar <- NA
ftp.short$binomial <- NA
ftp.short$dbh <- NA
ftp.short$ir5 <- NA
ftp.short$c13 <- NA
ftp.short$htot <- NA
ftp.short$age <- NA
ftp.short$year <- NA
ftp.short$origin <- NA

ifn.fr.renamed$speciesid <- NA
ifn.fr.renamed$dbh1 <- NA
ifn.fr.renamed$dbh2 <- NA
ifn.fr.renamed$height1 <- NA
ifn.fr.renamed$height2 <- NA
ifn.fr.renamed$ba1 <- NA
ifn.fr.renamed$ba_ha1 <- NA
ifn.fr.renamed$bachange_ha  <- NA
ifn.fr.renamed$recruitment.plot.rate <- NA
ifn.fr.renamed$sp.recruitment.plot.rate <- NA
ifn.fr.renamed$recruitment.plot.ba <- NA
ifn.fr.renamed$sp.recruitment.plot.ba <- NA
ifn.fr.renamed$weight1 <- NA
ifn.fr.renamed$yearsbetweensurveys <- NA
ifn.fr.renamed$surveydate1 <- NA
ifn.fr.renamed$surveydate2 <- NA
ifn.fr.renamed$dbh1_mod <- NA
ifn.fr.renamed$ba1_mod <- NA
ifn.fr.renamed$bachange_mod <- NA

##Subset target traits from fundiv
# 20/02/2018 Added two mortality ba and rates
names(ftp.short)
fundiv<- subset(ftp.short, select=c("plotcode", "longitude", "latitude", "treecode", "country", 
                                    "speciesid","espar","binomial","dbh","ir5","c13","htot","age", "dbh1", "dbh2", "height1", "height2",
                                    "ba1","ba_ha1","ba2","ba_ha2",
                                    "bachange_ha","BAI","bachange_ha_yr","speciesrichness",
                                    "ba_ha","BAnei","BAj.plot","BAneicon",
                                    "BAI.plot","BAInei","BAIj.plot","BAIneicon",
                                    "mortality.plot.rate","sp.mortality.plot.rate","mortality.plot.ba","sp.mortality.plot.ba",
                                    "recruitment.plot.rate","sp.recruitment.plot.rate","recruitment.plot.ba","sp.recruitment.plot.ba",    
                                    "treestatus_th",
                                    "weight1","weight2","year",
                                    "yearsbetweensurveys", "surveydate1", "surveydate2",
                                    "dbh1_mod","ba1_mod","bachange_mod","origin"))
##Subset target traits from ifn.fr
names(ifn.fr.renamed)
fr <- subset(ifn.fr.renamed, select =c( "idp", "lon", "lat","treecode","country", "speciesid",
                                     "espar", "binomial","dbh","ir5","c13","htot","age","dbh1", "dbh2", "height1", "height2",
                                     "ba1","ba_ha1","BA","BA.ha","bachange_ha","BAI","BAI.ha","richness.plot","BA.plot",
                                     "BAnei","BAj.plot","BAneicon",
                                     "BAI.plot","BAInei","BAIj.plot","BAIneicon",
                                     "mortality.plot.rate","sp.mortality.plot.rate","mortality.plot.ba","sp.mortality.plot.ba",
                                     "recruitment.plot.rate","sp.recruitment.plot.rate","recruitment.plot.ba","sp.recruitment.plot.ba", 
                                     "alive","weight1","w","year",
                                     "yearsbetweensurveys", "surveydate1", "surveydate2",
                                     "dbh1_mod","ba1_mod","bachange_mod","origin"))

#Variable C13 et htot et HG ori veget lib mortb sfcoeur etc etc  ?
#Need to add origin et C13 et htot et décaller les variable en dessous

colnames(fr)[25] <- "speciesrichness" #richness.plot
colnames(fr)[24] <- "bachange_ha_yr" #BAI.ha     #BAI individual
colnames(fr)[26] <- "BA.plot.ha" #BA.plot #sum des BA a l'hectare
colnames(fundiv)[26] <- "BA.plot.ha" #ba_ha #a checker# somme des BA a l'hectare pour le plot. Nom prete a confusion
colnames(fr)[1] <- "plotcode" 
colnames(fr)[20] <- "ba2" #BA
colnames(fr)[21] <- "ba_ha2"#BA.ha
colnames(fr)[2] <- "longitude" #lon
colnames(fr)[3] <- "latitude" #lat
colnames(fr)[42] <- "treestatus_th" #alive
colnames(fr)[44] <- "weight2" #w

#Here we have two tables with the 52 variables we are interested in. #52 and before 48

#########################################################
##  Add a new tree code for merging france and fundiv: ##
#########################################################

#fill in species id and origin in fr from the allspfinal file
fr$speciesid <- Allspfinal$speciesid[match(fr$espar, Allspfinal$espar, incomparables = NA)]
fr$origin <- Allspfinal$origin[match(fr$espar, Allspfinal$espar,incomparables = NA)]
fr$code <- Allspfinal$code[match(fr$espar, Allspfinal$espar,incomparables = NA)] #recup les codes
fr$binomial <- Allspfinal$binomial[match(fr$espar, Allspfinal$espar,incomparables = NA)] #a partir de ces codes corriger les noms latins

#fill in espar binomial and origin in fr from the allspfinal file
fundiv$espar <- Allspfinal$espar[match(fundiv$speciesid, Allspfinal$speciesid,incomparables = NA)]
fundiv$origin <- Allspfinal$origin[match(fundiv$speciesid, Allspfinal$speciesid,incomparables = NA)]
fundiv$binomial <- Allspfinal$binomial[match(fundiv$speciesid, Allspfinal$speciesid,incomparables = NA)]
fundiv$code <- Allspfinal$code[match(fundiv$speciesid, Allspfinal$speciesid,incomparables = NA)]

#Merge together both database : 2132130 obs and 50 variables. 2132130 - 142972 (>10 dbh) = 1989158
treefinal <- rbind(fundiv,fr)
test <- semi_join(treefinal,Allspfinal,by="speciesid") #8 indiv sans noms
test <- semi_join(treefinal,Allspfinal,by="binomial") #8 indiv sans noms

length(unique(Allspfinal$binomial))
length(unique(treefinal$binomial))#9 sp disparaisse


#check si on explique bien les données

test=anti_join(treefinal,Allspfinal,by="binomial") #Hop la    8170 pbmmmmm
test=anti_join(Allspfinal,treefinal,by="binomial") #3 indiv

test=anti_join(treefinal,Allspfinal,by="speciesid") #Hop la    8170 pbmmmmm
test=anti_join(Allspfinal,treefinal,by="speciesid") #3 indiv 

test=anti_join(treefinal,Allspfinal,by="code") #Hop la    8170 pbmmmmm
test=anti_join(Allspfinal,treefinal,by="code") #3 indiv


test <- semi_join(treefinal,Allspfinal,by="speciesid")#772626 Ce que j'explique avec mon tableau
test <- anti_join(treefinal,Allspfinal,by="binomial")#8170 Ce que je n'explique pas (pas de code)
length(unique(treefinal$speciesid))

A <-(treefinal[is.na(treefinal$speciesid)&is.na(treefinal$code),]) #8168 indiv sans info
A <-(treefinal[is.na(treefinal$speciesid)&is.na(treefinal$code)&is.na(treefinal$code),]) # 2 indiv mauvaise sp relevé
test=anti_join(Allspfinal,treefinal,by="code") #3 indiv plus présent dans treefinal car dbh < 10. 

# 12, 26, 98, 263, 266, 267, 314, 329, 335, 339, 349, 353 (apparaissent pas)
# 8 => Apparait car premier NA (8000 indiv) #checked 
# 181, 248 coup de chance bien attribué
# 234,236,455,347 les deux ont un espar donc okay pour la base fr. Mais lequel attribuer pour fundiv ?
treefinal[!is.na(treefinal$espar)&treefinal$espar=="49PC",]
treefinal[!is.na(treefinal$espar)&treefinal$espar=="49RT",]
treefinal[!is.na(treefinal$espar)&treefinal$espar=="25000",]

save(treefinal,file="/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/treefinal.RData")

#Trait to model: increment of basal area per species and per hectare
bachange_ha_yr = BAnei + mortality.plot + sp.mortality.plot + speciesrichness + AI + SWC






####check for 1 sps and 1 plot
kk1 <-subset(xxx, species= )
kk <-subset(kk1, plotcode == "10000_1")


