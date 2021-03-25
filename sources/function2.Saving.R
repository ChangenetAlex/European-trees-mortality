# Alex on the 23/05/2018
# Script to save the output of my models with package spaMM
# It also requires the function Rsquared

SavingInfo = "Saving 
A R function to save your model and various output  
Require the knitr & the spaMM packages. RandomFields might be needed if there is spatial autocorrelation

Input data: Your model (HLfit)

create and save four elements: 
(1) A directory named after your model 
(2) A csv file with the output results of the Rsquared.AC function 
(3) A csv file with the estimated coefs, standards deviations and significance levels 
(4) A txt file with these information as LaTex format \n (5) The pattern of the spatial autocorrelation is also plotted "
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))

library(knitr)
library(spaMM)
library(fields)
library(piecewiseSEM)
library(raster)
library(rworldmap)
library(rgdal)
library(rasterVis)
library(rgeos)
library(ggplot2)

Saving <- function(x){
  if (x$family$family=="binomial"){
    Dir <- paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/binomial")
  }else Dir <- paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models/Negbin")
  if (length(grep(substitute(x),pattern = "get",fixed=T,value=F,invert=F))!=0){ # Added to simplify the model selection process 
    z <- paste0(Mymod,num)
    dir.create(paste0(Dir,"/",z,"/"))
    save(x, file = paste0(Dir,"/",z,"/",z,".rda")) # save the model as an RDA file 
  }else {dir.create(paste0(Dir,"/",deparse(substitute(x)),"/"))
    save(x, file = paste0(Dir,"/",deparse(substitute(x)),"/",deparse(substitute(x)),".rda"))
    z <- deparse(substitute(x))} # save the model as an RDA file 
  
  # New section 
  r1 <- data.frame(x$data$longitude, x$data$latitude, predict(x))
  colnames(r1) <- c("X","Y","Z")
  if (x$family$family=="binomial") {r1$bin <- ifelse(r1$Z>0.20,"0",ifelse(r1$Z>0.1,"1","2")) ## Percentage of mortality as classes 
  }else r1$bin <- ifelse(r1$Z>40,"0",ifelse(r1$Z>20,"1","2"))
  par(mfrow=c(1,1))
  worldmap <- getMap(resolution = "high")
  europe <- worldmap[which(worldmap$REGION=="Europe"),]
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & as.character(worldmap$NAME) != "Russia"),] 
  europe <- worldmap[which(grepl("Europe", worldmap$GEO3)),] 
  europe <-spTransform(europe,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ") # Convert to the right system
  europe = crop(europe,c(-15, 45, 35, 70))
  p <- ggplot() + theme(panel.background = element_rect(fill="lightblue", colour="black", size=3, 
                                                        linetype=1, color="black"),legend.key.size = unit(1.5, "cm")) +
    theme(panel.grid.major = element_line(colour="blue", size=4, linetype=0,lineend="square", color="red"),
          panel.grid.minor = element_line(colour="blue", size=1, linetype=0,lineend="round", color="blue"))+
    theme(legend.background = element_rect(fill="black",size=0.5, linetype="solid", 
                                           colour ="black"))+
    geom_polygon(data = europe, aes(x=long, y = lat,group=group),fill="gray20",col="gray14") + 
    coord_fixed(1.3) + geom_point(data = r1, aes(x = X, y = Y, group=bin, col = bin, shape = bin, size=bin))
  if (x$family$family=="binomial") {p <- p + scale_colour_manual(values = c("red","white","blue"),name="Proportion %",labels = c(">20%", "10-20%","0-10%")) +
    scale_shape_manual(values = c(20,3,17),name="Proportion %",labels = c(">20%", "10-20%","0-10%")) +
    scale_size_manual(values= c(5,1,1),name="Proportion %",labels = c(">20%", "10-20%","0-10%"))
  }else {p <- p + scale_colour_manual(values = c("red","white","blue"),name="Number of events by plot",labels = c(">40", "20-40","0-20")) +
    scale_shape_manual(values = c(20,3,17),name="Number of events by plot",labels = c(">40", "20-40","0-20")) +
    scale_size_manual(values= c(1,1,1),name="Number of events by plot",labels = c(">40", "20-40","0-20"))}
  p <- p + guides(shape = guide_legend(override.aes = list(size = 5)))+
    labs(title=paste0('Predicted mortality'), y=paste0("Latitude"), x="Longitude", caption="Changenet et al. 2018")+
    theme(text = element_text(face="bold"),legend.direction ="vertical",
          axis.text.x = element_text(size=13,color="black"),axis.text.y = element_text(size=13,color="black"),
          legend.key = element_rect(fill = "gray20", colour = "white"),
          legend.background=element_rect(fill="white",colour="black",size=0.2),
          panel.border = element_rect(colour = "black", fill=NA, size=0.8),
          axis.line = element_line(colour="black"),
          plot.title = element_text(size=18,hjust = 0.5),
          plot.caption = element_text(face="bold.italic"))
  #ggsave(filename = paste0(Dir,"/",deparse(substitute(x)),"/",deparse(substitute(x)),"Pred.Mort_",deparse(substitute(x)),".png"),plot = p, width = 12.37, height = 7.04, dpi=150,units = "in")
  ggsave(filename = paste0(Dir,"/",z,"/",z,"Pred.Mort_",z,".png"),plot = p, width = 12.37, height = 7.04, dpi=150,units = "in")
  
  # End new section ##
  
  if (x$family$family=="binomial"){Namecat <- "Binomial"} else Namecat <- "Negbin" # Identidy the family
  y <- data.frame(matrix(unlist(rsquared.AC(x)), nrow=1,byrow=T,dimnames = list(c(NULL),c("family","link","method","Marginal","Conditional","Lik","AICm","Craw1","ICCadj1", "PCVran", "PCVObs","ICCraw2","ICCadj2","PCV1","PCV2"))))
  if (length(list.files(path=paste0(Dir,"/"),pattern=paste0("Models_",Namecat,CODE,"_",seuil,".csv")))==0){
    write.table(y,paste0(Dir,"/Models_",Namecat,CODE,"_",seuil,".csv"),append=T, quote = FALSE, sep=";",row.names = z,col.names=NA) # If file does not exist
  }else write.table(y, paste0(Dir,"/Models_",Namecat,CODE,"_",seuil,".csv"),append=T, quote = FALSE, sep=";",row.names = z,col.names=F)  # If it does exist                             
  # The rsquared output are stored at the end of each other in a huge Csv file
  
  A <- as.data.frame(summary(x)[2]) # All intercepts fixed
  
  if (length(grep("Matern",x[["predictor"]],fixed=T))>=1)    # If there is a spatial effect
  {B <- cbind("",t(t(as.numeric(unclass(x[["lambda"]])))),"") # Extract the lambda values in a table 
  colnames(B) <- colnames(A)
  rownames(B) <- sub(" .","", names(x[["lambda"]]), ignore.case = FALSE,fixed = T) # Name them correctly
  B <- rbind(B,B) # Duplicate this table twice 
  B[1:(nrow(B)/2),1] <- "Log Lambda" # First half of the first column are loglambda
  B[(1+(nrow(B)/2)):nrow(B),1] <- "Lambda" # The second half of the first column is lambda
  B[1:(nrow(B)/2),2] <- log(as.numeric(B[(1+(nrow(B)/2)):nrow(B),2])) # The second half of the second column is then equal to the exponential of the first half
  D <- cbind("Estimates",t(t(unlist(x$ranFix[[1]]))),"") # Extract the estimated parameters of the spatial effect
  rownames(D) <- sub("2.","",rownames(D), ignore.case = FALSE,fixed = T)
  B <- rbind(B,D) # Combine all together
  par(mfrow=c(1,1))
  d<- seq(0,40,length.out=80)
  y<- Matern(d, range=as.numeric(D[2,2]), smoothness=as.numeric(D[1,2])) # Plot my spatial effects
  matplot( d, y, type="l", lty=1, lwd=2,main=c(bquote(paste("Estimated parameters : ",rho," = ",.(round(as.numeric(D[2,2]),3))," and ",nu," = ",.(round(as.numeric(D[1,2]),3))))),ylab="Response covariance",xlab="Distance between plots")
  dev.print(file=paste0(Dir,"/",z,"/",z,"Autocor.jpeg"),device=jpeg,width=710) # Save it 
  
  } else { # If there is no spatial effect 
    
    B <- as.data.frame(summary(x)[3])[,2:4] # Extract the random effects 
    colnames(B) <- colnames(A)
    B <- rbind(B,B) # Duplicate this table twice
    B[,1] <- as.character(B[,1]) # Make the first column as a character
    B[1:(nrow(B)/2),1] <- "Log Lambda" # First half of the first column are loglambda
    B[(1+(nrow(B)/2)):nrow(B),1] <- "Lambda" # The second half of the first column is lambda
    B[(1+(nrow(B)/2)):nrow(B),2] <- exp(B[1:(nrow(B)/2),2]) # The first half of the second column is then equal to the log of the second half
    B[(1+(nrow(B)/2)):nrow(B),3] <- "" # The second half of the third column is empty
  }
  
  C <- rbind(A,B) # All together (random + fixed effects)
  colnames(C) <- c("Estimate","Cond. SE","t-value")
  C[,"Signif"] <- ""
  C[C$`t-value`!="" & (as.numeric(C$`t-value`)>=2 | as.numeric(C$`t-value`)<=-2),"Signif"] <- "*" # Add significiance 
  write.table(C,paste0(Dir,"/",z,"/",z,"_Table.csv"),sep=";",col.names = NA, row.names = T) # Variable outputs in a table.csv
  D <- kable(C,digits = 3,"latex",align=c("rrrc"))
  capture.output(print(D), file=paste0(Dir,"/",z,"/",z,"_Tex.Table.txt")) # Output as a latex wrapped in a txt file
}

