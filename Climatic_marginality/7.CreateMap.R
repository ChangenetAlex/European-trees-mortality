#####################################################################"
######## Create climate maps with location of sampling points #######"
#####################################################################"

### Required packages
require(rgdal)
require(raster)
require(parallel)
require(fields)

### FUNCTION
clim.map<-function(clim.var,period,fun.map='mean',extent.map=c(-20,60,20,72),xy=NULL,dir.in,download=F,dir.out=NULL,n.core=1,plot.map=T) {
  
  # climatic variables provided by EuMedClim
  bioclim=paste0("bio",c(1,2,5,6,12,13,14))# WorldClim type parameters
  T.seas=paste0("tmean.",c("djf","mam","jja","son"))# seasonal mean temperature
  P.seas=paste0("prec.",c("djf","mam","jja","son"))# seasonal precipitation
  pet=paste0("pet.",c("mean","min","max"))# potential evapotranspiration
  ppet=paste0("ppet.",c("mean","min","max"))# water balance
  eumedclim.vars=c(bioclim,T.seas,P.seas,pet,ppet)# all
  
  # verify conditions
  if(is.null(clim.var)) stop("error in the 'clim.var' argument, select one EuMedClim climatic variable") else
    if(length(clim.var)>1) stop("error in the 'clim.var' argument, select one EuMedClim climatic variable") else
      if(all(clim.var!=eumedclim.vars))
        stop(paste('the climatic variable',clim.var,'is not provided by EuMedClim'))
  
  if(is.null(period))
    stop("the argument 'period' is empty. Select years within 1901-2014") else
      if(min(period)<1901 | max(period)>2014)
        stop("error in the 'period' argument. Select years within 1901-2014")
  
  if(length(extent.map)!=4 | (extent.map[1]<(-20) | extent.map[2]>60 | extent.map[3]<20 | extent.map[4]>72))
    stop("error in the 'extent.map' argument. Select extent within -20,60° longitude E and within 20,72° latitude N")
  
  if(length(period)>1 & all(fun.map!=c('mean','min','max')))
    stop("error in the 'fun.map' argument.")
  
  if(is.null(dir.in))
    stop("the argument 'dir.in' is empty. Precise the full path of climate data folder (or where to store downloaded files)") else
      if(substr(dir.in,nchar(dir.in),nchar(dir.in))!='/')
        dir.in=paste0(dir.in,'/')
  
  if(is.null(dir.out))
    warning("the argument 'dir.out' is empty. To save map, precise the full path of a destination folder") else
      if(substr(dir.out,nchar(dir.out),nchar(dir.out))!='/')
        dir.out=paste0(dir.out,'/')
  
  # function to select years
  sel.yr<-function(r,yrs)
    unlist(lapply(as.list(yrs),function(y) c(1:length(names(r)))[names(r)==paste("X",y,sep="")]))
  
  # tiles
  tiles=list(c(-20,0,20,45),c(0,20,20,45),c(20,40,20,45),c(40,60,20,45),
             c(-20,0,45,72),c(0,20,45,72),c(20,40,45,72),c(40,60,45,72))
  tiles=lapply(tiles,extent)
  # tiles to be selected
  tiles.in=c()
  for(t in 1:length(tiles))
    tiles.in[t] = any(extent.map[1:2]>=xmin(tiles[[t]])) & any(extent.map[1:2]<=xmax(tiles[[t]])) & any(extent.map[3:4]>=ymin(tiles[[t]])) & any(extent.map[3:4]<=ymax(tiles[[t]]))
  
  
  # compute tiled maps
  list.maps<-mclapply(as.list(c(1:8)[tiles.in]),function(t) {
    
    # tile t
    tile=tiles[[t]]
    
    # corresponding file
    file.nm=paste0(clim.var,"_1901-2014_1km_lon_",xmin(tile),"_",xmax(tile),"_lat_",ymin(tile),"_",ymax(tile),"_eumedclim.tif")
    
    # find folder
    dir.in.var=dir.in
    if(all(list.files(dir.in)!=file.nm))
      if(any(list.dirs(dir.in,full.names=F)==clim.var)) dir.in.var=paste0(dir.in,clim.var,'/')
    
    # download file
    if(all(list.files(dir.in.var)!=file.nm))
      
      if(download)
        
        download.file(url=paste0('http://gentree.data.inra.fr/climate/datasets/',clim.var,'/',file.nm),destfile=paste0(dir.in.var,file.nm),method='auto',cacheOK=F) else
          
          stop(paste("the following file",file.nm,"was not found in 'dir.in'. Verify the 'dir.in' argument or set download=T to automatically download climate files"))
    
    # load file
    stk=stack(paste0(dir.in.var,file.nm))
    names(stk)=paste0('X',1901:2014)
    
    # compute the mean, min or max values by pixel over the period
    FUN=get(fun.map)
    if(length(period)==1)
      map.t=stk[[sel.yr(stk,period)]] else
        map.t=calc(stk[[sel.yr(stk,period)]],FUN,na.rm=T)
    
    return(map.t)
    
  },mc.cores=n.core)
  
  # merge tiled maps
  for(i in 1:length(list.maps))
    if(i==1) Cmap=list.maps[[1]] else
      Cmap=merge(Cmap,list.maps[[i]])
  
  # mask on extent
  Cmap=crop(Cmap,extent.map)
  
  # convert units (to °C or mm)
  Cmap=0.1*round(Cmap)
  
  # plot map
  if(plot.map) {
    
    # map
    plot(Cmap,col=tim.colors(),bigplot=c(.1,.75,.1,.95),legend=F,axes=F,box=F)
    # add sampling points
    if(!is.null(xy)) if(ncol(xy)==2) points(xy,pch=3,cex=.5)
    # box
    axis(1,extent.map[1:2],labels=F,lwd.ticks=0,pos=extent.map[3])
    axis(2,extent.map[3:4],labels=F,lwd.ticks=0,pos=extent.map[1])
    axis(3,extent.map[1:2],labels=F,lwd.ticks=0,pos=extent.map[4])
    axis(4,extent.map[3:4],labels=F,lwd.ticks=0,pos=extent.map[2])
    # axes
    xmn.lab=seq(-20,60,5)[(seq(-20,60,5)-extent.map[1])>0][1]
    ymn.lab=seq(20,72,5)[(seq(20,72,5)-extent.map[3])>0][1]
    axis(1,seq(xmn.lab,extent.map[2],5),tck=-.02,lwd.ticks=1,lwd=0,cex.axis=.8,pos=extent.map[3])
    axis(2,seq(ymn.lab,extent.map[4],5),tck=-.02,lwd.ticks=1,lwd=0,cex.axis=.8,pos=extent.map[1])
    # legend
    if(any(substr(clim.var,1,5)==c("bio12","bio13","bio14","pet.m","prec.","ppet."))) unit.var="(mm)" else unit.var='(°C)'
    plot(Cmap,col=tim.colors(),legend.only=T,smallplot=c(.80,.82,.25,.75),legend.args=list(text=paste(clim.var,unit.var),side=3,adj=.2,font=2,line=1.5,cex=1))
    
  }
  
  # save map
  if(!is.null(dir.out)) {
    
    nm.map=paste0("map_",clim.var,"_",min(period),"-",max(period),ifelse(!is.null(xy),paste0("_with_",nrow(xy),"_sampling_points"),""),".tif")
    # write to the disk in GeoTiff format
    writeRaster(Cmap,filename=paste0(dir.out,nm.map),options="INTERLEAVE=BAND",overwrite=T)
  }
  
  return(Cmap)
  
}

### ARGUMENTS
# 'clim.var' character. Name of one EuMedClim climatic variable.
# 'fun.map' character. Name of the funtion to be applied among selected years among the 'mean' (by default), 'min' or 'max'.
# 'period' numeric. Vector of years to be selected within 1901-2014.
# 'extent.map' numeric. Vector of longitude and latitude boundaries (xmin,xmax,ymin,ymax) of the final map. By default the maximum geographical extent.
# 'xy' matrix or data.frame (optional). Coordinates of sampling points to plot on the climate map (two-columns coordinates of longitude and latitude in WGS84).
# 'dir.in' character. Full path of the folder where climate files are stored. If files are missing, set download=T to download them.
# 'download' logical. If TRUE, required climate files will be downloaded from the web and stored (full path of the destination folder given in 'dir.in'). FALSE by default.
# 'dir.out' character (optional). Full path of the folder where to save the climate map (GeoTiff format). Null (no saving) by default.
# 'n.core' numeric. Number of CPU cores to use (no parallelization by default). Useful when the map is generated from several tiles. Use 'detectCores()-1' to know the number of available cores.
# 'plot.map' logical.

detectCores()-1
### EXAMPLES
## Create map of 2000-2014 average annual mean temperature (bio1)
# Select a single climatic variable provided by EuMedClim (here annual mean temperature)
# Para to extract
bioclim=paste0("bio",c(1,2,5,6,12,13,14))# WorldClim type parameters
T.seas=paste0("tmean.",c("djf","mam","jja","son"))# seasonal mean temperature
P.seas=paste0("prec.",c("djf","mam","jja","son"))# seasonal precipitation
pet=paste0("pet.",c("mean","min","max"))# potential evapotranspiration
ppet=paste0("ppet.",c("mean","min","max"))# water balance
eumedclim.vars=c(bioclim,T.seas,P.seas,pet,ppet)# all
eumedclim.vars = eumedclim.vars[c(1:21)] #Here we can specify the files we want to extract

clim.var='bio1'
# Set the desired time period or years (within 1901-2014)
period=c(1901:2014)
# Set the function to compute mean, min or max values over the desired period
fun.map='mean'# by default
# Set the geographical extent of the map within [-20,60]°E [20,72]°N
extent.map=c(-20,60,20,72)# full extent
# Select sampling points (for plot)
#load("~/Documents/FUNDIV - NFI - Europe/our-data/treefinal.RData") #This is the database we want from which we want to analyse the metrics
#xy <- treefinal[treefinal$code==CODE&!is.na(treefinal$code),c("longitude","latitude")]
#xy=data.frame(longitude=c(-0.57,-3.70),latitude=c(44.82,40.41))
#rownames(xy)=c('bordeaux','madrid')
# path of the folder where climate data is stored (or will be stored if missing).
dir.in=c("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/climate/") # here your working directory
# path of the folder where you want to save map (NULL for no saving)
dir.out=paste0("~/SynologyDrive/FUNDIV - NFI - Europe/our-data/climate/") # here your working directory

# select number of available CPU cores
#n.core=detectCores() - 1

# run
#### On the 14th february 2020 it has been relaunched on all species with the period 1901-2014
for (i in c(1:length(eumedclim.vars))){
  try(clim.map(clim.var=eumedclim.vars[i],period=period,fun.map=fun.map,extent.map=extent.map,xy=NULL,dir.in=dir.in,download=F,dir.out=dir.out,n.core = 8,plot.map=F),silent=T) # Maybe not all the intervalle
}

clim.map(clim.var="bio5",period=period,fun.map=fun.map,extent.map=extent.map,xy=NULL,dir.in=dir.in,download=F,dir.out=dir.out,n.core = 8,plot.map = F)
bio1.map


