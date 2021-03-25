# Alex on the 01/06/2018
# Edited on the 25/07/2018
# Script to compute the confidence interval with spaMM for any interaction 
# Function extraction + function Bootstrap

SavingInfo = "ModelBoot( \n x = my model, \n N = the variable to vary as a number in extraction function, \n Ncat = the variable with three levels as a number in extraction function, \n LvL = number of step in the variable is by defaut 20 but can be something else, \n CAT = what part of the map as a factor beaing 0 = core, 1 = leading or 2 = rear edge, \n nBoot = number of bootstraps by defaut 100, \n saveboot = To save or not to save is a logical by defaut is TRUE, \n)
Require the dplyr, spaMM, parallel, ggplot2 and reshape2 packages.\n
!!! WARNING !!! Before running this function you need to be in the directory corresponding to the model you want to study. !!! WARNING !!!\n
Output: \n A interaction plot named after your model and your variables of interest \n
Example : 
Extraction(MnbZT80.60.Gbis4)
ModelBoot(MnbZT80.60.Gbis4,3,5,LvL=40,CAT,nBoot=20,saveboot=F)\n"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))

## Packages

library(spaMM)
library(dplyr)
library(plyr)
library(parallel)
library(ggplot2)
library(pbapply)
#library(plotly)
library(reshape2)
library(stringr)

# 
################### MY PARAMETERS ####################
LvL <- 20              # Number of levels           ##
N = 9                  # Variable to vary           ## # Can be the quadratic effects
Ncat = 7               # The three level variable   ##
CAT <- as.factor("0")  # Category by defaut is 0.   ##
nBoot = 10             # nbr of bootstraps          ##
saveboot = T           # Save or not the plot       ##
Yportion = 0.66        # Pourcentage a echantilloner##
nCoeur = 10            # Nbr of cores               ##
######################################################

### Extract my fixed effects names 
#x <- M2bin7B.17
Extraction <- function(x){
  A <- paste0(getCall(x)[2])
  A <- unlist(strsplit(A, "~",fixed = T))[2]
  A <- unlist(strsplit(A, "+",fixed = T))
  A <- grep(A,pattern = "|",fixed=T,value=T,invert=T)
  A <- unlist(strsplit(A, ":",fixed = T))
  #A <- grep(A,pattern = ":",fixed=T,value=T,invert=T)
  A <- grep(A,pattern = "e)",fixed=T,value=T,invert=T) #Remove the AC term which remains
  A <- grep(A,pattern = "^[ ]$",value=T,invert=T) #Remove empty characters
  A <- sub("I(","", A, ignore.case = FALSE,fixed = T)
  A <- sub("^2)","", A, ignore.case = FALSE,fixed = T)
  A <- sub("offset(log(","", A, ignore.case = FALSE,fixed = T)
  A <- sub("))","", A, ignore.case = FALSE,fixed = T)
  A <- sub("\n ","", A, ignore.case = FALSE,fixed = T) # New line
  A <- sub(" ","", A, ignore.case = FALSE,fixed = T)
  A <- sub(" ","", A, ignore.case = FALSE,fixed = T) #Twice
  A <- sub(" ","", A, ignore.case = FALSE,fixed = T) #Twice
  A <- sub(" ","", A, ignore.case = FALSE,fixed = T) #Twice
  A <- unique(A)
  A
}

## Bootstrap function 
ModelBoot <- function(x,N = N,Ncat = Ncat,LvL = LvL, CAT = CAT, nBoot = nBoot,Yportion=Yportion, saveboot = saveboot,nCoeur=nCoeur){
  A <- Extraction(x)
  Variation = Extraction(x)[N]
  VariationCAT = Extraction(x)[Ncat]
  Resp <- paste0(getCall(x)[2])
  Resp <- unlist(strsplit(Resp, " ~ ",fixed = T))[1]
  ### For each single one except one give average values
  data <- x$data[,colnames(x$data)%in%A[1:length(A)]] # Extract by colnames and in the right order
  data <- data[,c(A[1:length(A)])] # Reorder my data Need in all cases !!!
  # Creer table avec nouvelles valeurs 
  myDF = as.data.frame(matrix(ncol = length(A), nrow = LvL*3, NA)) # Df vierge
  colnames(myDF) <- colnames(data) # Bon noms de colonnes
  for (i in c(1:(length(A)))){
    if (!is.na(myDF[1,i])){next} # It is possible for Plotcatr to be a column before the interesting variable. Therefore, it is possible for a loop to fill two columns in the same iteration and then to erase the info we want (VariationCAT located after Plotcat). This line ensure that is a column is filled, it remains as it is. 
    if (is.numeric(data[,i])){myDF[,i] <- as.numeric(rep(mean(data[,i]),LvL))} # if numeric, put at the average from data
    else if (is.numeric(data[,VariationCAT])){ # If not, is "i" corresponds to Variation CAT ??? If VarCat is a numeric, it is not i ...
      myDF[,VariationCAT] <- c(as.numeric(rep(mean(data[,VariationCAT]),LvL)), # ... Then we need Lvl * the mean of the pop for this column whis is not i
                               rep(as.numeric(quantile(data[,VariationCAT],0.995)), length=LvL), # Lvl time the first quantile = Cold pop
                               mvar1_Low<-rep(as.numeric(quantile(data[,VariationCAT],0.005)), length=LvL)) # Lvl time the Warm pop = lowest quantile
      myDF[,i] <- CAT}  # Then i is has to be equal to CAT that we have chosen.
    #else myDF[,i] <- as.factor(c(rep("0",LvL),rep("1",LvL),rep("2",LvL)))}
    else myDF[,VariationCAT] <- as.factor(c(rep("0",LvL),rep("1",LvL),rep("2",LvL)))} # If i is not numeric neither does VariationCat, they are the same and need to varying
  
  myDF[,Variation] <- seq(min(data[,Variation]),max(data[,Variation]),length=LvL) ### The interest variable was at mean first but is replaced by a sequence of lenbgth LVL.
  
  #predict(Mod,newdata=numeric(0)) ## TO RUN IF THERE IS A BUG
  
  #colnames(myDF)[names(myDF)%in%(C[1:length(C)])] <- D  # If there is a quadratic term -> Need to rename the myDF before bootstrap !!
  pred.boot =rep(NA,nBoot)
  Bootpred <- mclapply(pred.boot,function(train){
    if (x$family$family=="binomial"){
      calibrate.data = sample(1:nrow(dfplot2[!is.na(dfplot2$sp.mort.bin), ]), Yportion*nrow(dfplot2[!is.na(dfplot2$sp.mort.bin), ])) # 66% of our data
      train = dfplot2[!is.na(dfplot2$sp.mort.bin), ][calibrate.data, ]
      #train <- rbind(sample_n(dfplot2[!is.na(dfplot2$sp.mort.bin) & dfplot2$Plotcat.80.60==0,],400),
                     #sample_n(dfplot2[!is.na(dfplot2$sp.mort.bin) & dfplot2$Plotcat.80.60==1,],400),
                     #sample_n(dfplot2[!is.na(dfplot2$sp.mort.bin) & dfplot2$Plotcat.80.60==2,],200))
    } else if (grepl('.ba', Resp)==F){ 
      calibrate.data = sample(1:nrow(dfplot2[dfplot2$sp.mortality.plot.count.yr >0, ]), Yportion * nrow(dfplot2[dfplot2$sp.mortality.plot.count.yr > 0, ])) # 66% of our data
      train = dfplot2[!is.na(dfplot2$BAI.O.plot.1) & dfplot2$sp.mortality.plot.count.yr > 0, ][calibrate.data, ]
      #train <- rbind(sample_n(dfplot2[dfplot2$sp.mortality.plot.count.yr>0 & !is.na(dfplot2$sp.mortality.plot.count.yr) & dfplot2$Plotcat.80.60==0,],300),
       #              sample_n(dfplot2[dfplot2$sp.mortality.plot.count.yr>0 & !is.na(dfplot2$sp.mortality.plot.count.yr) & dfplot2$Plotcat.80.60==1,],300),
        #             sample_n(dfplot2[dfplot2$sp.mortality.plot.count.yr>0 & !is.na(dfplot2$sp.mortality.plot.count.yr) & dfplot2$Plotcat.80.60==2,],30))
    } else { 
      calibrate.data = sample(1:nrow(dfplot2[dfplot2$sp.mortality.plot.count.yr >0, ]), Yportion * nrow(dfplot2[dfplot2$sp.mortality.plot.count.yr > 0, ])) # 66% of our data
      train = dfplot2[!is.na(dfplot2$BAI.O.plot.1) & dfplot2$sp.mortality.plot.count.yr > 0, ][calibrate.data, ]
      #train <- rbind(sample_n(dfplot2[dfplot2$sp.mortality.plot.ba>0 & !is.na(dfplot2$sp.mortality.plot.ba) & dfplot2$Plotcat.80.60==0,],300),
       #              sample_n(dfplot2[dfplot2$sp.mortality.plot.ba>0 & !is.na(dfplot2$sp.mortality.plot.ba) & dfplot2$Plotcat.80.60==1,],300),
        #             sample_n(dfplot2[dfplot2$sp.mortality.plot.ba>0 & !is.na(dfplot2$sp.mortality.plot.ba) & dfplot2$Plotcat.80.60==2,],30))
    }
    sub_binomial = update(x, data=train) # Fit my model on my sample 
    pred.boot=as.numeric(predict(sub_binomial, newdata=myDF,re.form=NA)) # Predict the values for each cases describe in MyDf
    },mc.cores=nCoeur,mc.silent=T)
  as.data.frame(pred.boot) # Does it is usefull for anything ? 
  
  Myerrors <- unlist(Bootpred[sapply(Bootpred, function(x) inherits(x, "try-error"))]) # Take the errors out 
  if (length(Myerrors)!=0) {print(Myerrors) # If any, print it on the screen 
    Bootpred <- Bootpred[sapply(Bootpred, function(x) !inherits(x, "try-error"))] # Keep the non-errors
    message(paste0(nBoot-length(Bootpred)," Bootstraps on ",nBoot," failed because of an error")) 
    }
  Bootpred <- as.data.frame(Bootpred) # Convert it as a dataframe 
  message(paste0("There is ",length(Bootpred[Bootpred>2000])," predicted values above 2000, the maximum being ",max(Bootpred),". These are going to be replaced by NA"))
  Above <- as.data.frame(table(unlist(sapply(as.list(Bootpred),function(x) which(x >2000)))))
  message(paste0(Above$Var1," * ",Above$Freq,sep="\n"),"Are the lines for which a huge values as been estimated and the number of times")
  Bootpred[Bootpred>2000] <- NA
  Bootpred0 <- as.data.frame(Bootpred[1:LvL,]) # Mean or Core
  Bootpred1 <- as.data.frame(Bootpred[(1+LvL):(2*LvL),]) # High values or rear edge
  Bootpred2 <- as.data.frame(Bootpred[(1+(LvL*2)):(LvL*3),]) # Low values or leading edge
  
  # compute mean and SD predicted values without a accounting for the NA 
  Means_Bootpred = rowMeans(Bootpred,na.rm=T) # mean of each row for the range
  SD_Bootpred = apply(Bootpred, 1, function(x) sd(x,na.rm=T)) # sd of each column for the range 
  yMAX <- max(Means_Bootpred)+max(SD_Bootpred) # Here is the maximum values that will be the limit for the plot
  Means_Bootpred0 = rowMeans(Bootpred0,na.rm=T) # mean of each row
  SD_Bootpred0 = apply(Bootpred0, 1, function(x) sd(x,na.rm=T)) # sd of each column
  Means_Bootpred1 = rowMeans(Bootpred1,na.rm=T) # mean of each row
  SD_Bootpred1 = apply(Bootpred1, 1, function(x) sd(x,na.rm=T)) # sd of each columnMeans_Bootpred = rowMeans(Bootpred) # mean of each row
  Means_Bootpred2 = rowMeans(Bootpred2,na.rm=T) # mean of each row
  SD_Bootpred2 = apply(Bootpred2, 1, function(x) sd(x,na.rm=T)) # sd of each column
  
  #yMAX <- max(c(Means_Bootpred0,Means_Bootpred1))+max(c(SD_Bootpred0,SD_Bootpred1)) # Here is the maximum values when removing one category (need to modify the number)
  
  # Retrieve original variation AND if log ou sqrt redefine the column to go 
  #if (grepl(Variation,pattern = "sqrt",fixed=T)==T){Variation <- sub("sqrt","", Variation, ignore.case = FALSE,fixed = T)
  #}else if (grepl(Variation,pattern = "log",fixed=T)==T) Variation <- sub("log","", Variation, ignore.case = FALSE,fixed = T)
  #original_Variation=myDF[1:LvL,Variation]*sd(dfplot[,Variation],na.rm=T)+mean(dfplot[,Variation],na.rm=T) #From myDF data to real values based on dfplot

  # Retrieve original variation AND if log ou sqrt redefine the column to go 
  
  # Check infinite and Nan values and remove the line. (Just infinite for the moment) and remove them after 
  l.NaN <- length(which(is.nan(dfplotbis[,Variation])))
  l.Inf <- length(which(is.infinite(dfplotbis[,Variation])))
  if (l.NaN!=0|l.Inf!=0){
    message(paste0("There is ",l.NaN," NaN values in the original database and ",l.Inf," Infinite values for the Variation variable : ",Variation,
                   " It will be removed from dfplotbis database to calculate mean, sd and to retrieve original values"))}
  
  ## Replace Inf values by NA
  
  if (l.Inf!=0){
    Z <- dfplotbis[,Variation]
    Z[is.infinite(Z)] <- NA
    dfplotbis[,Variation] <- Z}
  
  # If sqrt scaled
  if (grepl(Variation,pattern = "sqrt",fixed=T)==T){
  original_Variation=myDF[1:LvL,Variation]*sd(dfplotbis[,Variation],na.rm=T)+mean(dfplotbis[,Variation],na.rm=T) #Descale and denorm with tranformed but no scaled variables (dfplotbis)
  original_Variation=original_Variation*original_Variation # Detrasnformed
  Variation <- sub("sqrt","", Variation, ignore.case = FALSE,fixed = T)
  #If logscaled
  }else if (grepl(Variation,pattern = "log",fixed=T)==T){ #If 
  original_Variation=myDF[1:LvL,Variation]*sd(dfplotbis[,Variation],na.rm=T)+mean(dfplotbis[,Variation],na.rm=T) #Descale and denorm with tranformed but no scaled variables (dfplotbis)
  original_Variation <- exp(original_Variation)-1 # De transfo  
  Variation <- sub("log","", Variation, ignore.case = FALSE,fixed = T)
  #If not scaled
  }else original_Variation=myDF[1:LvL,Variation]*sd(dfplotbis[,Variation],na.rm=T)+mean(dfplotbis[,Variation],na.rm=T) #Descale and denorm with original variable (dfplotbis)
  
  # Create a longformat DF with the original variation parameter as well
  testDF <- data.frame(Means_Bootpred0,Means_Bootpred1,Means_Bootpred2,original_Variation)
  colnames(testDF)<-gsub("Means_","",colnames(testDF))
  test_data_long <- melt(testDF, id=c("original_Variation"), value.name=c("Means"), variable="Marginality")  # convert to long format
  testDF <- data.frame(SD_Bootpred0,SD_Bootpred1,SD_Bootpred2,original_Variation)
  colnames(testDF)<-gsub("SD_","",colnames(testDF))
  test_data_long_sd <- melt(testDF, id=c("original_Variation"), value.name="SD", variable="Marginality")  # convert to long format
  
  #Plot with ggplot2 package
  Esp <- getwd()
  Esp <- stringr::str_extract(Esp, "species/.{0,6}")
  Esp <- str_sub(Esp, 9)
  p1<-ggplot(data=merge(test_data_long, test_data_long_sd), # All
  #p1<-ggplot(data=merge(test_data_long[-c(31:60),], test_data_long_sd[-c(31:60),]), #If want to remove the red part (Need to change legend as well)
  #p1<-ggplot(data=merge(test_data_long[-c((1+(LvL*2)):(LvL*3)),], test_data_long_sd[-c((1+(LvL*2)):(LvL*3)),]), #If want to remove the last part (Need to change legend as well)
    aes(x=original_Variation, y=Means, col=Marginality,shape=Marginality,
    ymin=Means-SD, ymax=Means+SD))+
    geom_line()+
    geom_point(size=2)+
    geom_linerange()
  # Legend levels depend upon the used variable
  if (is.numeric(data[,VariationCAT]))
  {p1 <- p1 + scale_color_manual("",values = c("black", "blue", "red"),labels=(c("Average level","High level","Low level")))+
    scale_shape_manual("", values=c(1,2,3),labels=(c("Average level","High level","Low level")))
  }else p1 <- p1 + scale_color_manual("",values = c("black", "blue", "red"),labels=(c("Core","Leading Edge","Rear Edge")))+
    scale_shape_manual("", values=c(1,2,3),labels=(c("Core","Leading Edge","Rear Edge")))
  if (x$family$family=="binomial"){
    #setwd(dir = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/binomial/",deparse(substitute(x)),"/")) # All models
    setwd(dir = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/",deparse(substitute(x)),"/")) # One model
    p1 <- p1 + labs(title=paste0(VariationCAT," : ",Variation), y=paste0(CODE,": Mortality occurence probability"), x=Variation)
  } else {
    #setwd(dir = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/Negbin/",deparse(substitute(x)),"/")) # All models
    setwd(dir = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/",deparse(substitute(x)),"/")) # One model
    p1 <- p1 + labs(title=paste0(VariationCAT," : ",Variation), y=paste0(CODE,": Annual predicted mortality (in ‰)"), x=Variation)}
  p1 <- p1 + coord_cartesian(ylim=c(0,yMAX),expand = T)+
    #scale_y_continuous(breaks = waiver(),limits = c(0,yMAX)) + 
    theme_light(base_size = 15)+
    theme(text = element_text(face="bold"),legend.position = "none",          # legend.position = c(0.5,0.8),legend.direction ="horizontal"
          axis.text.x = element_text(size=17,color="black"),axis.text.y = element_text(size=17,color="black"),
          axis.title.x = element_text(size=17,color="black"),axis.title.y = element_text(size=17,color="black"),
          legend.background=element_rect(fill="white",colour="black",size=0.2),
          panel.border = element_rect(colour = "black", fill=NA, size=0.8),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=20,hjust = 0.5))
  print(p1)
  if (saveboot==T){ggsave(filename = paste0("GG_",CODE,"_",deparse(substitute(x)),"_",seuil,"_",Variation,"_",VariationCAT,".150.png"),plot = p1, width = 12.37, height = 7.04, dpi=150,units = "in")
    ggsave(filename = paste0("GG_",CODE,"_",deparse(substitute(x)),"_",seuil,"_",Variation,"_",VariationCAT,".300.png"),plot = p1, width = 12.37, height = 7.04, dpi=300,units = "in")
    ggsave(filename = paste0("GG_",CODE,"_",deparse(substitute(x)),"_",seuil,"_",Variation,"_",VariationCAT,".600.png"),plot = p1, width = 12.37, height = 7.04, dpi=600,units = "in")
    ggsave(filename = paste0("GG_",CODE,"_",deparse(substitute(x)),"_",seuil,"_",Variation,"_",VariationCAT,".800.png"),plot = p1, width = 12.37, height = 7.04, dpi=800,units = "in")
    ggsave(filename = paste0("GG_",CODE,"_",deparse(substitute(x)),"_",seuil,"_",Variation,"_",VariationCAT,".1000.png"),plot = p1, width = 12.37, height = 7.04, dpi=1000,units = "in")
    } # To modify dpi (300 ou 600 au lieu de 150)
  
  #Here is the final plot
  dev.off()
  gc()
}
