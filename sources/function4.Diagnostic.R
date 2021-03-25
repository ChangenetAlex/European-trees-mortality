# Alex on the 01/06/2018
# Script to compute the confidence interval with spaMM for any interaction 
# Function extraction + function Bootstrap

SavingInfo = " This function is to make the diagnostic and the cross validation for the mortality models 
Diagnostic(
x = Mymodel
Yportion = Percentage of the data on which the cross validation is done
AllInOne = Do you want all the figure to be on a single one or in different ones 
nAC = Number of plots we want to use for the autocorrelation spatial analyse
)
!!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! 

Before running this function you need to be in the directory corresponding to the model you want to study. 

!!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! !!! WARNING !!! 

Output:
Several plot and the a txt file with the values of the cross validation
"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))

## Packages 

## Packages
library(spaMM)
library(dplyr)
library(arm)
library(pROC) # For the Roc Curve 
library(lattice)

## The function 

Diagnostic <- function(x,Yportion,AllInOne=T){
  
  # Extract the response and the explicative variable 
  ExplVar <- paste0(getCall(x)[2])
  Response <- unlist(strsplit(ExplVar, " ~ ",fixed = T))[1]
  ExplVar <- unlist(strsplit(ExplVar, " ~ ",fixed = T))[2]
  ExplVar <- unlist(strsplit(ExplVar, "+",fixed = T))
  ExplVar <- grep(ExplVar,pattern = "|",fixed=T,value=T,invert=T)
  ExplVar <- grep(ExplVar,pattern = ":",fixed=T,value=T,invert=T)
  ExplVar <- grep(ExplVar,pattern = "e)",fixed=T,value=T,invert=T) #Remove the AC term which remains
  ExplVar <- grep(ExplVar,pattern = "^[ ]$",value=T,invert=T) #Remove empty characters
  ExplVar <- sub("I(","", ExplVar, ignore.case = FALSE,fixed = T)
  ExplVar <- sub("^2)","", ExplVar, ignore.case = FALSE,fixed = T)
  ExplVar <- sub("offset(log(","", ExplVar, ignore.case = FALSE,fixed = T)
  ExplVar <- sub("))","", ExplVar, ignore.case = FALSE,fixed = T)
  ExplVar <- sub("\n ","", ExplVar, ignore.case = FALSE,fixed = T) # New line
  ExplVar <- sub(" ","", ExplVar, ignore.case = FALSE,fixed = T)
  ExplVar <- sub(" ","", ExplVar, ignore.case = FALSE,fixed = T) #Twice
  ExplVar <- unique(ExplVar)
  
  
  #### AC Spatial on the plots on which mortality occurs
  
  #lctools
  if (x$family$family=="negbin"){
    sampled=sample(1:nrow(dfplot2[dfplot2$sp.mortality.plot.rate.yr>0,]),nrow(dfplot2[dfplot2$sp.mortality.plot.rate.yr>0,])) # All values above 0
    dfplotVARIO=dfplot2[sampled,]
    dfplotVARIO <- dfplotVARIO[!duplicated(dfplotVARIO[c("longitude","latitude")]),]
    Coords <- cbind(dfplotVARIO$longitude, dfplotVARIO$latitude) ## Coordonee
    bws <- c(2,3,5, 10, 15, 20,30,40,50,100) # Neighbors
    
    png(file=paste0(CODE,"_",deparse(substitute(x)),"VarioAC.LCtools.png"),width = 12.37, height = 7.04, res=600,units = "in")
    moran <- lctools::moransI.v(Coords, bws, dfplotVARIO$sp.mortality.plot.rate.yr, WType='Binary') #par defaut = binary et neighbors
    capture.output(lctools::moransI.v(Coords, bws, dfplotVARIO$sp.mortality.plot.rate.yr, WType='Binary'),file="MoranTest.LCTools")
    dev.off()
    bws <- c(10,20,30,40,42) # Kernels distance
    #moran <- lctools::moransI.v(Coords, bws, dfplotVARIO$sp.mortality.plot.rate.yr, WType='Binary',family="fixed")
    png(file=paste0(CODE,"_",deparse(substitute(x)),"ScatterAC.LCtools.png"),width = 12.37, height = 7.04, res=600,units = "in")
    l.moran <- lctools::l.moransI(Coords,6,dfplotVARIO$sp.mortality.plot.rate.yr) # scatterplot 
    dev.off()
    # Vario with Correlog
    #corD <- correlog(Coords, dfplotVARIO$sp.mortality.plot.rate.yr, method = "Moran")
    #jpeg(file=paste0(deparse(substitute(x)),"Vario.AC.Correlog.jpeg"),width = 12.37, height = 7.04, res=150,units = "in")
    #plot(corD)
    #dev.off()
    
  }else{ 
    sampled=sample(1:nrow(dfplot2[dfplot2$sp.mort.bin>0,]),nrow(dfplot2[dfplot2$sp.mort.bin>0,])) # All values above 0
    dfplotVARIO=dfplot2[sampled,]
    dfplotVARIO <- dfplotVARIO[!duplicated(dfplotVARIO[c("longitude","latitude")]),]
    Coords <- cbind(dfplotVARIO$longitude, dfplotVARIO$latitude) ## Coordonee
    bws <- c(2,3,5, 10, 15, 20,30,40,50,100) # Neighbors
    
    png(file=paste0(CODE,"_",deparse(substitute(x)),"VarioAC.LCtools.png"),width = 12.37, height = 7.04, res=600,units = "in")
    moran <- lctools::moransI.v(Coords, bws, dfplotVARIO$sp.mort.bin, WType='Binary') #par defaut = binary et neighbors
    capture.output(lctools::moransI.v(Coords, bws, dfplotVARIO$sp.mort.bin, WType='Binary'),file="MoranTest.LCTools")
    dev.off()
    bws <- c(10,20,30,40,42) # Kernels distance
    moran <- lctools::moransI.v(Coords, bws, dfplotVARIO$sp.mortality.plot.rate.yr, WType='Binary',family="fixed")
    png(file=paste0(CODE,"_",deparse(substitute(x)),"ScatterAC.LCtools.png"),width = 12.37, height = 7.04, res=600,units = "in")
    l.moran <- lctools::l.moransI(Coords,6,dfplotVARIO$sp.mort.bin) # scatterplot 
    dev.off()
    
    # Vario with Correlog (a enlever ou pas)
    corD <- correlog(Coords, dfplotVARIO$sp.mort.bin, method = "Moran")
    png(file=paste0(CODE,"_",deparse(substitute(x)),"Vario.AC.Correlog.png"),width = 12.37, height = 7.04, res=600,units = "in")
    plot(corD)
    dev.off()
  }
  
  
  # Add lines to do the same on the residuals 
  
  
  ### Cross validation ###
  
  if (x$family$family=="negbin"){
    #setwd(dir = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/Negbin/",deparse(substitute(x)),"/")) # All negbin model
    setwd(dir = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/",deparse(substitute(x)),"/")) # One model
    
    #### Define my data ##### That are all data I have not taken in my model
    calibrate.data=sample(1:nrow(dfplot2[dfplot2$sp.mortality.plot.rate.yr>0,]), Yportion*nrow(dfplot2[dfplot2$sp.mortality.plot.rate.yr>0,])) # 66% of our data 
    df1=dfplot2[!is.na(dfplot2$BAI.O.plot.1)&dfplot2$sp.mortality.plot.rate.yr>0,][calibrate.data,]
    df2=dfplot2[!is.na(dfplot2$BAI.O.plot.1)&dfplot2$sp.mortality.plot.rate.yr>0,][-calibrate.data,]
    Nbr <- nrow(df2)
    df2 <- df2[complete.cases(df2[,ExplVar]),]
    message(paste0(nrow(df2)," individuals with no NA in my explicative variables to do the predictions on ",Nbr))
    
    M1 = update(x, data=df1)
    M2 = predict (M1, newdata=df2,re.form=NA,type="response")
    
    ### Cor Test ###
    if (AllInOne == T){png(file=paste0(CODE,"_",deparse(substitute(x)),"Diagnostic.png"),width = 12.37, height = 7.04, res=600,units = "in");par(mfrow = c(2,2))
    } else {png(file=paste0(CODE,"_",deparse(substitute(x)),"_Response_VS_predicted.png"),width = 12.37, height = 7.04, res=600,units = "in")}
    
    D <- cor.test(df2$sp.mortality.plot.rate.yr, M2, method = ("pearson")) # To record and store somewhere
    capture.output(print(D), file=paste0(CODE,"_",deparse(substitute(x)),"_CrossValid.Table.txt")) # Output as a latex wrapped in a txt file
    plot(df2$sp.mortality.plot.count.yr, M2,xlab="Fitted values", ylab="Predicted", main="Predicted vs. fitted",
         cex.main=1.5, cex.lab =1.3) # To save 
    abline(0,1)
    if (AllInOne == F){dev.off()}
    
    ### Residuals ###
    ResM <- residuals(M1)
    FitM <- fitted(M1)
    if (AllInOne == F) png(file=paste0(CODE,"_",deparse(substitute(x)),"_Residuals_VS_fitted.png"),width = 12.37, height = 7.04, res=600,units = "in") # Save it 
    plot(ResM ~ FitM, xlab="Fitted values", ylab="Residuals", main="Residuals vs. fitted",
         cex.main=1.5, cex.lab =1.3)
    abline(h=0)
    if (AllInOne == F){dev.off()} # Save it 
    
    # Save this plot 
    if (AllInOne == F){png(file=paste0(CODE,"_",deparse(substitute(x)),"_QQplot.png"),width = 12.37, height = 7.04, res=600,units = "in")} # Save it 
    qqnorm(ResM, cex.main=1.5, cex.lab =1.3)
    qqline(ResM)
    if (AllInOne == F){dev.off()} # Save it 
    
    # Save this one 
    if (AllInOne == F){png(file=paste0(CODE,"_",deparse(substitute(x)),"_Observed_VS_fitted.png"),width = 12.37, height = 7.04, res=600,units = "in")} # Save it 
    plot(M1$data$sp.mortality.plot.count.yr ~ FitM, xlab="Fitted values", ylab="Observed", main="Observed vs. fitted",
         cex.main=1.5, cex.lab =1.3)
    abline(0,1)
    dev.off() # Save it 
    
  } else {
    
    ##################################
    ######    NegaBin models   #######
    ##################################
    
    #setwd(dir = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/binomial/",deparse(substitute(x)),"/")) # All model
    setwd(dir = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/Results.all/Species/",CODE,"/",deparse(substitute(x)),"/")) # One model
    
    calibrate.data=sample(1:nrow(dfplot2[!is.na(dfplot2$sp.mort.bin),]), Yportion*nrow(dfplot2[!is.na(dfplot2$sp.mort.bin),])) # 66% of our data 
    df1=dfplot2[!is.na(dfplot2$sp.mort.bin),][calibrate.data,]
    df2=dfplot2[!is.na(dfplot2$sp.mort.bin),][-calibrate.data,]
    
    # Added on the 13th of july
    png(file=paste0(CODE,"_",deparse(substitute(x)),"_AUroC.png"),width = 12.37, height = 7.04,units = "in",res=600) # Save it 
    myroc <- roc(response = x$data$sp.mort.bin, predictor=as.numeric(predict(x)),
                 density.controls, density.cases, percent=FALSE, na.rm=TRUE,
                 quiet = TRUE, smooth=FALSE, auc=TRUE, ci=TRUE, plot=TRUE,print.auc=T,auc.polygon=T,grid=T)
    myrocCI <- ci.se(myroc, boot.n=100,progress = "none") # Add the CI 
    plot(myrocCI) ## Plot the AUC for the ROC
    dev.off()
    
    M1 = update(x, data=df1)
    M2 = predict (M1, newdata=df2,re.form=NA,type="response")
    
    ### Cor Test ###
    df2$names <- rownames(df2)
    M2 <- as.data.frame(M2)
    M2$names <-rownames(M2)
    df2 <- semi_join(df2,M2,by="names")
    D <- cor.test(df2$sp.mort.bin, M2$V1, method = ("pearson")) # To record somewhere
    capture.output(print(D), file=paste0(CODE,"_",deparse(substitute(x)),"_CrossValid.Table.txt")) # Output as a latex wrapped in a txt file
    
    if (AllInOne == T){png(file=paste0(CODE,"_",deparse(substitute(x)),"Diagnostic.png"),width = 12.37, height = 7.04, res=600,units = "in");par(mfrow = c(2,2))
    } else {png(file=paste0(CODE,"_",deparse(substitute(x)),"_Response_VS_predicted.png"),width = 12.37, height = 7.04, res=600,units = "in")}
    plot(df2$sp.mort.bin, M2$V1,xlab="Fitted values", ylab="Predicted", main="Predicted vs. fitted",
         cex.main=1.5, cex.lab =1.3) # To save 
    if (AllInOne == F){dev.off()} # To save 
    
    ### Residuals ###
    ResM <- residuals(M1)
    FitM <- fitted(M1)
    PredM <- predict(M1)
    if (AllInOne == F) png(file=paste0(CODE,"_",deparse(substitute(x)),"_Residuals_VS_fitted.png"),width = 12.37, height = 7.04, res=600,units = "in") # Save it 
    par(mfrow=c(1,1))
    plot(ResM ~ FitM, xlab="Fitted values", ylab="Residuals", main="Residuals vs. fitted",
         cex.main=1.5, cex.lab =1.3)
    abline(h=0)
    if (AllInOne == F){dev.off()} # Save it 
    
    # Save this plot 
    if (AllInOne == F){png(file=paste0(CODE,"_",deparse(substitute(x)),"_QQplot.png"),width = 12.37, height = 7.04, res=600,units = "in")} # Save it 
    qqnorm(ResM, cex.main=1.5, cex.lab =1.3)
    qqline(ResM)
    if (AllInOne == F){dev.off()} # Save it 
    
    # Save this one 
    if (AllInOne == F){png(file=paste0(CODE,"_",deparse(substitute(x)),"_BinnedPlot.png"),width = 12.37, height = 7.04, res=600,units = "in")} # Save it 
    binnedplot(PredM,ResM)
    dev.off()
  }
  
  # Plot residuals against all my explicative variable
  
  c1 <- rep(ResM,times=length(ExplVar))
  M1$data$Plotcat <- as.numeric(as.character(M1$data$Plotcat))
  c3 <- rep(M1$data$Plotcat,times=length(ExplVar))
  c2 <- gather(M1$data[,c(ExplVar)]) # To remove the warning
  #c2 <- c2[-c(c2$key=="Plotcat"),]
  test <- cbind(c1,c2[,c(2,1)],c3)
  
  png(file=paste0(CODE,"_",deparse(substitute(x)),"Residuals_Scatterplots.png"),width = 12.37, height = 7.04, res=600,units = "in")
  par(mfrow=c(1,1))
  p <- xyplot(test[,1]~test[,2] | test[,3],grid = T,
              #group=test[,4],
              scales=list(alternating=TRUE,
                          x=list(relation="free"),
                          y=list(relation="same")),
              xlab="Explanatory variables",
              ylab="Response variable",
              panel=function(x,y){
                panel.grid(h=-1,v=2)
                panel.points(x,y,col=1,cex=0.2,pch=19)
                panel.loess(x,y,xol=1,lwd=2)
                #panel.groups=test[,4]
              })
  print(p)
  dev.off()
  png(file=paste0(CODE,"_",deparse(substitute(x)),"_VisualisationTEST.png"),width = 12.37, height = 7.04, res=600,units = "in") # Save it 
  plot(x,ask=F)
  dev.off()
  dev.off()
  #plot(x)
  #dev.print(file=paste0(deparse(substitute(x)),"_Visualisation2.jpeg"),device=jpeg,width=710) # Save it 
  #dev.off()
  #dev.print(file=paste0(deparse(substitute(x)),"_Visualisation.jpeg"),device=jpeg,width=710) # Save it 
  #par(ask=F)
  #par(mfrow = c(1,1))
}

