SavingInfo = " This function extact the SPEI indexes for the time period we want and for different indexes of calculated spei (from 1 to 48 months)

Example :  SPEI(nStart = 1981,nEnd=2015,ncname = 'spei12') 
To extract the spei indexes based on the last 12 months for my plots between 1981  and 2015
"
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak,SavingInfo,bannerBreak,"\n"))

#################################
#######                  ########
#######   Extract SPEI   ########
#######                  ########
##########################################
nStart = 1981   #    Starting year      ##
nEnd = 2015     # Final year (included) ##
ncname <- "spei12" # Whic Spei Indice   ##
##########################################
library("sp")              #
library("ncdf4")           #
library("chron")           #
library("RColorBrewer")    #
library("lattice")         #
library("raster")          #
library("rgdal")           #
library("maptools")        #
library("matrixStats")     #
library("stringr")         #
############################

SPEI <- function(nStart = nStart, nEnd=nEnd,ncname = ncname){
  load("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/tfinal.biotic.June2018.RData")
  df <- unique(tfinal.biotic[,c("plotcode","latitude","longitude","country","surveydate1","surveydate2","year")])
  Dir = paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/climate/SPEI/")
  setwd(Dir)
  

  
  ncfname <- paste0(Dir, ncname,".nc")
  dname <- "spei"  # note: tmp means temperature (not temporary)
  # open a netCDF file
  ncin <- nc_open(ncfname)
  
  # Obtain the longitude area we are interested in 
  lon.ti <- which(ncvar_get(ncin, "lon") == (round(min(df$longitude))-0.75))
  lon.tf <- which(ncvar_get(ncin, "lon") == (round(max(df$longitude))+0.75)) - which(ncvar_get(ncin, "lon") == (round(min(df$longitude))-0.75))
  lon <- ncvar_get(ncin, "lon", start=c(lon.ti), count=c(lon.tf))
  
  #Obtain the latitude area we are interested in 
  lat.ti <- which(ncvar_get(ncin, "lat") == (round(min(df$latitude))-0.75))
  lat.tf <- which(ncvar_get(ncin, "lat") == (round(max(df$latitude))+0.75)) - which(ncvar_get(ncin, "lat") == (round(min(df$latitude))-0.75))
  lat <- ncvar_get(ncin, "lat", start=c(lat.ti), count=c(lat.tf))
  
  # Obtain the intervalle of time we are interested in 
  time <- ncvar_get(ncin, "time", start=c((nStart-1901)*12+1), count=c(nEnd-nStart+1)*12)
  
  
  #####################################################
  ## Get the input variable (SPEI) and its attributes,# 
  ##        and verify the size of the array         ##
  #####################################################
  
  spei.array <- ncvar_get(ncin, dname, start=c(lon.ti, lat.ti,(nStart-1901)*12+1),count=c(lon.tf,lat.tf,(nEnd-nStart+1)*12), verbose = TRUE)
  dlname <- ncatt_get(ncin, dname,"long_name")
  dunits <- ncatt_get(ncin, dname, "units")
  fillvalue <- ncatt_get(ncin, dname, "_FillValue")
  
  #################################################
  ###CONVERSION OF THE NETCDF FILES TO DATA.FRAMES
  #################################################
  # Replace netCDF fillvalues with R NAs
  #values of a variable that are either missing or simply not available 
  #(i.e. ocean grid points in a terrestrial data set) are flagged 
  #using specific fill values (_FillValue) or missing values (missing_value)
  spei.array[spei.array == fillvalue$value] <- NA
  length(na.omit(as.vector(spei.array[, , 1])))
  
  #################################################################################
  ##Convert the whole array to a data frame, and calculate the annual mean
  #create a long vector tmp.vec.long using the as.vector() reshaping function, 
  #verify its length, which should be 4568760
  spei.vec.long <- as.vector(spei.array)
  
  #Then reshape that vector into a 4568760 by 12 matrix using the matrix() function, 
  ##verify its dimensions, which should be 4568760 by 12.
  
  nlon <- dim(lon)
  nlat <- dim(lat)
  nt <- dim(time)
  
  #spei ann has dimension of 10878 (no. pixels) by 420
  spei.ann <- matrix(spei.vec.long, nrow = nlon * nlat, ncol = nt)

  #there is the montly information for each pixel by 12 months by 35 years
  #head(na.omit(spei.ann))
  #Create data frame from the spei.ann matrix.
  lonlat <- as.matrix(expand.grid(lon, lat))
  spei.df02 <- data.frame(cbind(lonlat, spei.ann))
  mm<-c("speiJan", "speiFeb", "speiMar", "speiApr", "speiMay", 
        "speiJun", "speiJul", "speiAug", "speiSep", "speiOct", "speiNov", "speiDec")
  mmm<-rep(mm,nEnd-nStart+1) #nbr year
  names(spei.df02) <- c("lon", "lat", mmm)
  options(width = 110)
  #head(na.omit(spei.df02, 20))
  
  #Get annual mean
  #str(spei.df02)
  #dim(spei.df02)
  ss <- spei.df02[,3:(nt+2)]
  b <- matrix(nrow=nt/12, ncol=nlat*nlon) 
  b <-apply(b, 1,as.numeric)
  year_data <- data.frame(b)
  colnames(year_data) <- paste("yr", 1:(nt/12), sep = "")
  #str(year_data)
  
  # loop through the 35 years
  for(i in 1:(nt/12)){
    z <- i*12
    j <- z-11
    aa <- ss[,j:z]
    str(aa)
    year_data[i] <- rowMeans(aa)
    summary(year_data[,i])
  }
  long <- spei.df02[1]
  lat <- spei.df02[2]
  
  final_data <- data.frame(long, lat, year_data)
  save(final_data, file=paste0(ncname,"_",nStart,"_",nEnd,".RData"))
  
  
  # match the spei to the plots
  coords <- df[,c(3,2)]
  dfr <- rasterFromXYZ(final_data) 
  SPEIplot <- extract(dfr, coords)
  SPEIplot <- as.data.frame(SPEIplot)
  
  # select the years and bind to the plot data
  #SPEIplot <- SPEIplot[,c(1:35)]
  spei_plots_all <- cbind(df, SPEIplot)
  colnames(spei_plots_all)[8:(7+nt/12)] <- (nStart:nEnd)
  
  # calculate the min and mean spei of the years between the two surveys
  spei_plots_all$min_spei_survey_years <- NA
  spei_plots_all$mean_spei_survey_years <- NA
  # French plots missing start year
  spei_plots_all[spei_plots_all$country!="FR","start_year"] <- as.numeric(str_sub(spei_plots_all[spei_plots_all$country!="FR","surveydate1"],1,4))
  spei_plots_all[spei_plots_all$country!="FR","end_year"] <- as.numeric(str_sub(spei_plots_all[spei_plots_all$country!="FR","surveydate2"],1,4))
  spei_plots_all[spei_plots_all$country=='FR',"end_year"] <- spei_plots_all[spei_plots_all$country=='FR',"year"]
  spei_plots_all[spei_plots_all$country=='FR',"start_year"] <- spei_plots_all[spei_plots_all$country=='FR',"end_year"] - 5
  
  
  for(i in 1:nrow(spei_plots_all)){
    if(is.na(spei_plots_all$mean_spei_survey_years[i])){
      print(spei_plots_all$plotcode[i])
      yrs <- as.character(spei_plots_all$start_year[i]:spei_plots_all$end_year[i])
      spei_plots_all$mean_spei_survey_years[i] <- rowMeans(spei_plots_all[i, yrs], na.rm=TRUE)
      spei_plots_all$min_spei_survey_years[i] <- rowMins(as.matrix(spei_plots_all[i, yrs]), na.rm=TRUE)
    }
  }
  save(spei_plots_all, file=paste0("plots_",ncname,"_",nStart,"_",nEnd,".RData"))
}
  