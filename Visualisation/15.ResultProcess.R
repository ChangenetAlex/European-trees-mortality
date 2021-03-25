rm(list = ls())
gc()
getwd()

tableresults <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/tableresults.csv")
summary(as.factor(tableresults$variable))
tableresults$estimates_normed <- as.vector(unlist(by(tableresults[,3],tableresults[,"species"],function(x) scale(x,center = F,scale = T))))
tableresults <- tableresults[-c(which(tableresults$variable=="(Intercept)")),] # On enlève l'intercept
tableresults$signif <- as.character(tableresults$signif)
tableresults$eco <- as.character(tableresults$eco)
tableresults$variable <- as.character(tableresults$variable)
tableresults[!is.na(tableresults$signif)&tableresults$signif=="*","signif"] <- T # Boleen
tableresults[!is.na(tableresults$signif)&tableresults$signif=="","signif"] <- F
tableresults <- tableresults[which(tableresults$variable%in%grep("Plotcat",tableresults$variable,value=T)),] #seulement interaction
tableresults$cat <- NA
tableresults$varOrder <- NA
tableresults$varCat <- NA
tableresults[which(tableresults$variable%in%grep("Plotcat1",tableresults$variable,value=T)),"cat"] <- "RearEdge" #
tableresults[which(tableresults$variable%in%grep("Plotcat2",tableresults$variable,value=T)),"cat"] <- "LeadingEdge"
tableresults$variable <- sub(":Plotcat1","",tableresults$variable)
tableresults$variable <- sub("Plotcat1:","",tableresults$variable)
tableresults$variable <- sub(":Plotcat2","",tableresults$variable)
tableresults$variable <- sub("Plotcat2:","",tableresults$variable)
tableresults$variable <- sub("Plotcat2","Plotcat",tableresults$variable)
tableresults$variable <- sub("Plotcat1","Plotcat",tableresults$variable)
tableresults[which(tableresults$variable%in%c("Plotcat")),"varCat"] <- "Localisation" 
tableresults[which(tableresults$variable%in%c("bio12_climate_mean.30","bio13_climate_mean.30","bio14_climate_mean.30")),
             "varCat"] <- "Precipitation" 
tableresults[which(tableresults$variable%in%c("bio1_climate_mean.30","bio5_climate_mean.30")),
             "varCat"] <- "Temperature" 
tableresults[which(tableresults$variable%in%c("mean_spei12","min_spei12")),
             "varCat"] <- "Drought" 
tableresults[which(tableresults$variable%in%c("sqrtBA.ha.plot.1","sqrtBA.O.plot.1","logBAj.plot.1")),
             "varCat"] <- "Competition"
n = 1 
for (j in c("RearEdge","LeadingEdge")){
  for (i in c("Plotcat","bio12_climate_mean.30","bio13_climate_mean.30","bio14_climate_mean.30","bio1_climate_mean.30","bio5_climate_mean.30",
              "mean_spei12","min_spei12","sqrtBA.ha.plot.1","sqrtBA.O.plot.1","logBAj.plot.1")){
    tableresults[tableresults$variable==i & tableresults$cat==j,"varOrder"] <- n
    n = n+1
  }
}




tableM <- tableresults[tableresults$eco=="M",]
tableT <- tableresults[tableresults$eco=="T",]
tableC <- tableresults[tableresults$type=="C",]
tableBL <- tableresults[tableresults$type=="BL",]


library(ggplot2)
#Box plots
#Mycol <- c("blue","darkblue","red","red2","green","darkgreen","gray10","gray40","brown","brown2")
Myshape = c(21,21,22,22,23,23,24,25,25)

# 1a
Myshape = c(21,21,22,22,23,23,24,25,25)
p <- ggplot(tableM, aes(x=variable,y=estimates_normed,shape=variable,group=cat,alpha=signif,size=signif,fill=cat)) + 
  facet_wrap(species~.) +
  #,scales = "free_x"
  #geom_point(position=position_dodge(1)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black") +
  geom_point(position = position_dodge2(width=0.7,preserve = "single")) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,8),labels=c("No significant","Significant"))
p
# 1b
p <- ggplot(tableM, aes(x=variable,y=estimates_normed,shape=variable,group=cat,alpha=signif,size=signif,fill=cat)) + 
  facet_wrap(species~.) +
  #,scales = "free_x"
  #geom_point(position=position_dodge(1)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black") +
  geom_point() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,8),labels=c("No significant","Significant"))
p
# 1c (double facet)
p <- ggplot(tableM, aes(x=variable,y=estimates_normed,shape=variable,group=cat,alpha=signif,size=signif,fill=cat)) + 
  facet_wrap(species~cat) +
  #,scales = "free_x"
  #geom_point(position=position_dodge(1)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black") +
  geom_point(position = position_dodge2(width=0.7,preserve = "single")) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,8),labels=c("No significant","Significant"))
p

# 1d
p <- ggplot(tableM, aes(x=varOrder,y=estimates_normed,shape=variable,group=cat,alpha=signif,size=signif,color=cat,fill=cat)) + 
  facet_wrap(species~.) +
  #,scales = "free_x"
  #geom_point(position=position_dodge(1)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black") +
  geom_point() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_color_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,8),labels=c("No significant","Significant"))
p


# 2d
p <- ggplot(tableT, aes(x=varOrder,y=estimates_normed,shape=variable,group=cat,alpha=signif,size=signif,color=cat,fill=cat)) + 
  facet_wrap(species~.) +
  #,scales = "free_x"
  #geom_point(position=position_dodge(1)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black") +
  geom_point() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_color_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,8),labels=c("No significant","Significant"))
p

# Change level order
levels(tableresults$species)
tableresults$species<-factor(tableresults$species,
                             levels = c("ABIALB", "BETPEN", "FAGSYL", "FRAEXC", "PICABI", "PINSYL", 
                                        "POPNIG", "POPTRE r", "QUEROB r","CASSAT", "PINPINAr", "PINHAL",
                                        "PINNIG r", "PINPIN r", "QUEILE", 
                                        "QUEPYR r", "QUESUB r"))
#dput(unique(as.character(tableM$species)))
#tableresults$species<-factor(tableresults$species,
# levels = rev(levels(tableresults$species)))


Myshape=c(4,21)
jitter <- position_jitter(width = 0.2, height = 0)
p <- ggplot(tableresults, aes(x=cat,y=species,shape=signif,group=variable,alpha=signif,color=cat,fill=cat,size=signif,label=estimates_normed)) + 
  facet_wrap(varCat~.) +
  theme(legend.position = c(0.8, 0.2)) +
  #,scales = "free_x"
  #geom_point(position = position_jitterdodge()) +
  geom_point(position=position_dodge(width=0.3))+
  geom_hline(yintercept=9.5) +
  #geom_hline(yintercept=0, linetype="dashed",color = "black") +
  #geom_point(position = jitter) +
  #theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_color_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,3),labels=c("No significant","Significant"))
p

install.packages("plotly")
library("plotly")

# Tout regrouper 
Myshape=c(4,21)
jitter <- position_jitter(width = 0.2, height = 0)
p <- ggplot(tableresults, aes(x=cat,y=species,shape=signif,group=variable,alpha=signif,color=cat,fill=cat,size=signif,label=estimates_normed)) + 
  #,scales = "free_x"
  #geom_point(position = position_jitterdodge()) +
  geom_point(position=position_dodge(width=0.3))+
  geom_hline(yintercept=9.5) +
  #geom_hline(yintercept=0, linetype="dashed",color = "black") +
  #geom_point(position = jitter) +
  #theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_color_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,3),labels=c("No significant","Significant"))
p



# Change level order
levels(tableresults$species)
tableresults$species<-factor(tableresults$species,
                             levels = c("ABIALB", "PICABI", "PINPINAr", "PINHAL", "PINNIG r", "PINPIN r", 
                                        "PINSYL","BETPEN", "CASSAT", "FAGSYL", "FRAEXC", "POPNIG", "POPTRE r", 
                                        "QUEILE", "QUEPYR r", "QUEROB r", "QUESUB r"))
dput(unique(as.character(tableBL$species)))
#tableresults$species<-factor(tableresults$species,
# levels = rev(levels(tableresults$species)))

Myshape=c(4,21)
jitter <- position_jitter(width = 0.2, height = 0)
p <- ggplot(tableresults, aes(x=cat,y=species,shape=signif,group=variable,alpha=signif,color=cat,fill=cat,size=signif,label=estimates_normed)) + 
  facet_wrap(varCat~.) +
  theme(legend.position = c(0.8, 0.2)) +
  #,scales = "free_x"
  #geom_point(position = position_jitterdodge()) +
  geom_point(position=position_dodge(width=0.3))+
  geom_hline(yintercept=8.5) +
  #geom_hline(yintercept=0, linetype="dashed",color = "black") +
  #geom_point(position = jitter) +
  #theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_color_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,3),labels=c("No significant","Significant"))
p


# Tout regrouper 
Myshape=c(4,21)
jitter <- position_jitter(width = 0.2, height = 0)
p <- ggplot(tableresults, aes(x=cat,y=species,shape=signif,group=variable,alpha=signif,color=cat,fill=cat,size=signif,label=estimates_normed)) + 
  #,scales = "free_x"
  #geom_point(position = position_jitterdodge()) +
  geom_point(position=position_dodge(width=0.3))+
  geom_hline(yintercept=8.5) +
  #geom_hline(yintercept=0, linetype="dashed",color = "black") +
  #geom_point(position = jitter) +
  #theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_color_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,3),labels=c("No significant","Significant"))
p





########################### NB ZT #####################################################


tableresults <- read.csv("~/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/tableresultsNB.csv")
tableresults$estimates_normed <- as.vector(unlist(by(tableresults[,3],tableresults[,"species"],function(x) scale(x,center = F,scale = T))))
tableresults <- tableresults[-c(which(tableresults$variable=="(Intercept)")),] # On enlève l'intercept
tableresults$signif <- as.character(tableresults$signif)
tableresults$eco <- as.character(tableresults$eco)
tableresults$variable <- as.character(tableresults$variable)
tableresults[!is.na(tableresults$signif)&tableresults$signif=="*","signif"] <- T # Boleen
tableresults[!is.na(tableresults$signif)&tableresults$signif=="","signif"] <- F
tableresults <- tableresults[which(tableresults$variable%in%grep("Plotcat",tableresults$variable,value=T)),] #seulement interaction
tableresults$cat <- NA
tableresults$varOrder <- NA
tableresults$varCat <- NA
tableresults[which(tableresults$variable%in%grep("Plotcat1",tableresults$variable,value=T)),"cat"] <- "RearEdge" #
tableresults[which(tableresults$variable%in%grep("Plotcat2",tableresults$variable,value=T)),"cat"] <- "LeadingEdge"
tableresults <- tableresults[-c(which(tableresults$variable%in%grep("Plotcat0",tableresults$variable,value=T))),] #seulement interaction avec leading and rear edge
tableresults$variable <- sub(":Plotcat1","",tableresults$variable)
tableresults$variable <- sub("Plotcat1:","",tableresults$variable)
tableresults$variable <- sub(":Plotcat2","",tableresults$variable)
tableresults$variable <- sub("Plotcat2:","",tableresults$variable)
tableresults$variable <- sub("Plotcat2","Plotcat",tableresults$variable)
tableresults$variable <- sub("Plotcat1","Plotcat",tableresults$variable)
tableresults[which(tableresults$variable%in%c("Plotcat")),"varCat"] <- "Localisation" 

tableresults[which(tableresults$variable%in%c("bio12_climate_mean.30","bio13_climate_mean.30","bio14_climate_mean.30")),
             "varCat"] <- "Precipitation" 
tableresults[which(tableresults$variable%in%c("bio1_climate_mean.30","bio5_climate_mean.30","tmean.djf_climate_mean.30")), # Temperature DJF en plus 
             "varCat"] <- "Temperature" 
tableresults[which(tableresults$variable%in%c("mean_spei12","min_spei12")),
             "varCat"] <- "Drought" 
tableresults[which(tableresults$variable%in%c("sqrtBA.ha.plot.1","sqrtBA.O.plot.1","logBAj.plot.1")),
             "varCat"] <- "Competition"


n = 1 
for (j in c("RearEdge","LeadingEdge")){
  for (i in c("Plotcat","bio12_climate_mean.30","bio13_climate_mean.30","bio14_climate_mean.30","bio1_climate_mean.30","bio5_climate_mean.30","tmean.djf_climate_mean.30",
              "mean_spei12","min_spei12","sqrtBA.ha.plot.1","sqrtBA.O.plot.1","logBAj.plot.1")){
    tableresults[tableresults$variable==i & tableresults$cat==j,"varOrder"] <- n
    n = n+1
  }
}

tableM <- tableresults[tableresults$eco=="M",]
tableT <- tableresults[tableresults$eco=="T",]
tableC <- tableresults[tableresults$type=="C",]
tableBL <- tableresults[tableresults$type=="BL",]


library(ggplot2)
#Box plots
#Mycol <- c("blue","darkblue","red","red2","green","darkgreen","gray10","gray40","brown","brown2")
Myshape = c(21,21,22,22,23,23,24,25,25)

# 1a
Myshape = c(21,21,22,22,23,23,24,25,25)
p <- ggplot(tableM, aes(x=variable,y=estimates_normed,shape=variable,group=cat,alpha=signif,size=signif,fill=cat)) + 
  facet_wrap(species~.) +
  #,scales = "free_x"
  #geom_point(position=position_dodge(1)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black") +
  geom_point(position = position_dodge2(width=0.7,preserve = "single")) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,8),labels=c("No significant","Significant"))
p
# 1b
p <- ggplot(tableM, aes(x=variable,y=estimates_normed,shape=variable,group=cat,alpha=signif,size=signif,fill=cat)) + 
  facet_wrap(species~.) +
  #,scales = "free_x"
  #geom_point(position=position_dodge(1)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black") +
  geom_point() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,8),labels=c("No significant","Significant"))
p
# 1c (double facet)
Myshape = c(21,21,22,22,23,23,24,25,25)
p <- ggplot(tableM, aes(x=variable,y=estimates_normed,shape=variable,group=cat,alpha=signif,size=signif,fill=cat)) + 
  facet_wrap(species~cat) +
  #,scales = "free_x"
  #geom_point(position=position_dodge(1)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black") +
  geom_point(position = position_dodge2(width=0.7,preserve = "single")) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,8),labels=c("No significant","Significant"))
p

p <- ggplot(tableT, aes(x=variable,y=estimates_normed,shape=variable,group=cat,alpha=signif,size=signif,fill=cat)) + 
  facet_wrap(species~cat) +
  #,scales = "free_x"
  #geom_point(position=position_dodge(1)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black") +
  geom_point(position = position_dodge2(width=0.7,preserve = "single")) +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,8),labels=c("No significant","Significant"))
p


# 1d
p <- ggplot(tableM, aes(x=varOrder,y=estimates_normed,shape=variable,group=cat,alpha=signif,size=signif,color=cat,fill=cat)) + 
  facet_wrap(species~.) +
  #,scales = "free_x"
  #geom_point(position=position_dodge(1)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black") +
  geom_point() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_color_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,8),labels=c("No significant","Significant"))
p

# 2d
p <- ggplot(tableT, aes(x=varOrder,y=estimates_normed,shape=variable,group=cat,alpha=signif,size=signif,color=cat,fill=cat)) + 
  facet_wrap(species~.) +
  #,scales = "free_x"
  #geom_point(position=position_dodge(1)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black") +
  geom_point() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_color_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,8),labels=c("No significant","Significant"))
p





install.packages("plotly")
library("plotly")

# Change the order of the levels to obtain the right order (med and then temp)
levels(tableresults$species)
dput(unique(as.character(tableM$species)))
dput(unique(as.character(tableT$species)))
tableresults$species<-factor(tableresults$species,
                             levels = c("ABIALB", "BETPEN", "FAGSYL", "FRAEXC", "PICABI", "PINSYL", 
                                        "POPTRE r","QUEPET r", "QUEROB r","CASSAT", "PINPINAr", "PINHAL",
                                        "PINNIG r", "PINPIN r", "QUEILE", 
                                        "QUEPYR r", "QUESUB r"))

tableresults$species<-factor(tableresults$species,
                             levels = rev(levels(tableresults$species))) # To reverse the order of the levels 




Myshape=c(4,21)
jitter <- position_jitter(width = 0.2, height = 0)
p <- ggplot(tableresults, aes(x=cat,y=species,shape=signif,group=variable,alpha=signif,color=cat,fill=cat,size=signif,label=estimates_normed)) + 
  facet_wrap(varCat~.) +
  #,scales = "free_x"
  #geom_point(position = position_jitterdodge()) +
  geom_point(position=position_dodge(width=0.3))+
  geom_hline(yintercept=9.5) +
  #geom_hline(yintercept=0, linetype="dashed",color = "black") +
  #geom_point(position = jitter) +
  #theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_color_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,3),labels=c("No significant","Significant"))
p


# Tout regrouper 
Myshape=c(4,21)
jitter <- position_jitter(width = 0.2, height = 0)
p <- ggplot(tableresults, aes(x=cat,y=species,shape=signif,group=variable,alpha=signif,color=cat,fill=cat,size=signif,label=estimates_normed)) + 
  #,scales = "free_x"
  #geom_point(position = position_jitterdodge()) +
  geom_point(position=position_dodge(width=0.3))+
  geom_hline(yintercept=9.5) +
  #geom_hline(yintercept=0, linetype="dashed",color = "black") +
  #geom_point(position = jitter) +
  #theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_color_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,3),labels=c("No significant","Significant"))
p



### Angio VS Gymno

# Change level order
levels(tableresults$species)
tableresults$species<-factor(tableresults$species,
                             levels = c("ABIALB", "PICABI", "PINPINAr", "PINHAL", "PINNIG r", "PINPIN r", 
                                        "PINSYL","BETPEN", "CASSAT", "FAGSYL", "FRAEXC","POPTRE r", 
                                        "QUEILE","QUEPET r","QUEPYR r", "QUEROB r", "QUESUB r"))
dput(unique(as.character(tableBL$species)))
dput(unique(as.character(tableC$species)))
#tableresults$species<-factor(tableresults$species,
# levels = rev(levels(tableresults$species)))

Myshape=c(4,21)
jitter <- position_jitter(width = 0.2, height = 0)
p <- ggplot(tableresults, aes(x=cat,y=species,shape=signif,group=variable,alpha=signif,color=cat,fill=cat,size=signif,label=estimates_normed)) + 
  facet_wrap(varCat~.) +
  theme(legend.position = c(0.8, 0.2)) +
  #,scales = "free_x"
  #geom_point(position = position_jitterdodge()) +
  geom_point(position=position_dodge(width=0.3))+
  geom_hline(yintercept=7.5) +
  #geom_hline(yintercept=0, linetype="dashed",color = "black") +
  #geom_point(position = jitter) +
  #theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_color_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,3),labels=c("No significant","Significant"))
p


# Tout regrouper 
Myshape=c(4,21)
jitter <- position_jitter(width = 0.2, height = 0)
p <- ggplot(tableresults, aes(x=cat,y=species,shape=signif,group=variable,alpha=signif,color=cat,fill=cat,size=signif,label=estimates_normed)) + 
  #,scales = "free_x"
  #geom_point(position = position_jitterdodge()) +
  geom_point(position=position_dodge(width=0.3))+
  geom_hline(yintercept=7.5) +
  #geom_hline(yintercept=0, linetype="dashed",color = "black") +
  #geom_point(position = jitter) +
  #theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  scale_color_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_fill_manual(values = c("blue","red"),labels=c("Leading Edge","Rear Edge")) +
  scale_shape_manual(values=Myshape) +
  scale_alpha_manual(values=c(0.3,1),labels=c("No significant","Significant")) +
  scale_size_manual(values=c(3,3),labels=c("No significant","Significant"))
p



