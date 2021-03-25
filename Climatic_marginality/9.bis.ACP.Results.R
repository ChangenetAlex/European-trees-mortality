library("xtable")
Allcode <- c("QUEILE","PINPINA","PINSYL","PICABI","FAGSYL","PINHAL","QUEROB","PINNIG","QUEPET","CASSAT","ABIALB","QUEPYR","FRAEXC","PINPIN",
        "QUESUB","BETPEN","POPTRE","POPNIG","ACEPSE","ALNGLU","LARDEC","QUEPUB")
ACP.all <- data.frame(matrix(nrow=21,ncol=3*22,NA))
names1 <- c("poids.tot","poids","poids.rank")
allnames <- expand.grid(names1,Allcode)
Tablenames <- paste0(allnames[,2],allnames[,1])
colnames(ACP.all) <- Tablenames

ACP.sum <- data.frame(matrix(nrow=22,ncol=6,NA))
colnames(ACP.sum) <- c("Species","Min Variable","Min %","Max Variable", "Max %","Mean %")

i = 1
for (code in Allcode){
acp10000 <- readRDS(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",code,"/PCA/acp10000.rds"))
ACP.all[,paste0(code,"poids.tot")] <- apply(acp10000$co,1,function(x) sum(abs(x)*c(acp10000$eig[1:2]/sum(acp10000$eig)*100)))
#poids2 <- poids[order(poids,decreasing=T)]
ACP.all[,paste0(code,"poids")] <- ACP.all[,paste0(code,"poids.tot")]/sum(ACP.all[,paste0(code,"poids.tot")])*100
ACP.all[,paste0(code,"poids.rank")] <- rank(-ACP.all[,paste0(code,"poids.tot")])

ACP.sum[i,1] <- code
ACP.sum[i,2] <- rownames(acp10000$co)[which(ACP.all[,paste0(code,"poids")]==min(ACP.all[,paste0(code,"poids")]))]
ACP.sum[i,3] <- round(min(ACP.all[,paste0(code,"poids")]),2)
ACP.sum[i,4] <- rownames(acp10000$co)[which(ACP.all[,paste0(code,"poids")]==max(ACP.all[,paste0(code,"poids")]))]
ACP.sum[i,5] <- round(max(ACP.all[,paste0(code,"poids")]),2)
ACP.sum[i,6] <- round(mean(ACP.all[,paste0(code,"poids")]),2)
i = i+1
}

ACP.all[,"Variable"] <- rownames(acp10000$co)
# Obtain the graph as a latex table 
#xtable(ACP.sum, caption = NULL, label = NULL, align = NULL, digits = 2)
write.table(ACP.sum, sep=",", file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/ACP.sum.csv",row.names = FALSE)



#### Table with variance explained by all axes for each species 
Allcode <- c("QUEILE","PINPINA","PINSYL","PICABI","FAGSYL","PINHAL","QUEROB","PINNIG","QUEPET","CASSAT","ABIALB","QUEPYR","FRAEXC","PINPIN",
             "QUESUB","BETPEN","POPTRE","POPNIG","ACEPSE","ALNGLU","LARDEC","QUEPUB")
ACP.all <- data.frame(matrix(nrow=21,ncol=22,NA))
#names1 <- c(" Eigenvalues"," Cumuleigenvalue")
#names1 <- c(" Cumuleigenvalue")
#allnames <- expand.grid(names1,Allcode)
#Tablenames <- paste0(allnames[,2],allnames[,1])
colnames(ACP.all) <- Allcode
i = 1
for (code in Allcode){
  acp <- readRDS(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",code,"/PCA/acp10000.rds"))
  #ACP.all[,paste0(code," Eigenvalues")] <- round(acp$eig/sum(acp$eig)*100,3) # Don't need it
  ACP.all[,paste0(code)] <- round(cumsum(acp$eig/sum(acp$eig)*100),3)
  i = i+1
}
ACP.all <- t(ACP.all)[,1:5]
colnames(ACP.all) <- c("Axe 1","Axe 2","Axe 3","Axe 4","Axe 5")
# Obtain the graph as a latex table 

#xtable(ACP.all, caption = NULL, label = NULL, align = NULL, digits = 2)
write.table(ACP.all, sep=",", file="/home/achangenet/Documents/FUNDIV - NFI - Europe/Redaction/Paper1/ACP.all.csv",row.names = TRUE)




