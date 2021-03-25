#################################################################################################"
###########  Get 1901-2014 climate data from geographic coordinates (WGS 84) ####################"
#################################################################################################"

### Required packages
require(rgdal)
require(raster)
require(parallel)


### FUNCTION
### Returns a dataframe with coordinates in row and yearly climatic values (1901-2014) in column. When more than one climatic parameter is selected, the function returns a list whose the length is the number of selected climatic parameters.
get.climate<-function(xy,clim.var,radius=NULL,dir.in,download=F,dir.out=NULL,n.core=20) {
  
  # climatic variables provided by EuMedClim
  bioclim=paste("bio",c(1,2,5,6,12,13,14),sep="")
  T.seas=paste0("tmean.",c("djf","mam","jja","son"))
  P.seas=paste0("prec.",c("djf","mam","jja","son"))
  pet=paste0("pet.",c("mean","min","max"))
  ppet=paste0("ppet.",c("mean","min","max"))
  eumedclim.vars=c(bioclim,T.seas,P.seas,pet,ppet)
  
  # verify conditions
  if(ncol(xy)!=2) stop('you should provide a two-column matrix of longitude and latitude data')
  if(!is.data.frame(xy)) xy=as.data.frame(xy)
  
  if(is.null(clim.var)) clim.var=eumedclim.vars else
    for(ki in clim.var) if(all(ki!=eumedclim.vars))
      stop(paste('the climatic variable',ki,'is not provided by EuMedClim'))
  
  if(is.null(dir.in))
    stop("the argument 'dir.in' is empty. Precise the full path of climate data folder (or where to store downloaded files)") else
      if(substr(dir.in,nchar(dir.in),nchar(dir.in))!='/')
        dir.in=paste0(dir.in,'/')
  
  if(is.null(dir.out))
    warning("the argument 'dir.out' is empty. To save map, precise the full path of a destination folder") else
      if(substr(dir.out,nchar(dir.out),nchar(dir.out))!='/')
        dir.out=paste0(dir.out,'/')
 
  # minimum radius
  if(!is.null(radius)) if(radius<5000) radius=5000
  
  # tiles
  tiles=list(c(-20,0,20,45),c(0,20,20,45),c(20,40,20,45),c(40,60,20,45),
             c(-20,0,45,72),c(0,20,45,72),c(20,40,45,72),c(40,60,45,72))
  tiles=lapply(tiles,extent)
  # tiles to be selected
  tiles.in=c()
  for(t in 1:length(tiles))
    tiles.in[t] = any(xy[,1]>=xmin(tiles[[t]]) & xy[,1]<xmax(tiles[[t]]) & xy[,2]>=ymin(tiles[[t]]) & xy[,2]<ymax(tiles[[t]]))

  # compute list of time series dataframes (one by parameter with coordinates in rows and yearly climate values in columns)
  list.ts=mclapply(as.list(clim.var),function(k) {
    
    clim.k=matrix(NA,nrow(xy),114,dimnames=list(rownames(xy),1901:2014))
    
    for(t in c(1:8)[tiles.in]) {
      
      # tile t
      tile=tiles[[t]]
      
      # corresponding coordinates
      in.tile <- xy[,1]>=xmin(tile) & xy[,1]<xmax(tile) & xy[,2]>=ymin(tile) & xy[,2]<ymax(tile)
      
      # get file
      file.nm=paste0(k,"_1901-2014_1km_lon_",xmin(tile),"_",xmax(tile),"_lat_",ymin(tile),"_",ymax(tile),"_eumedclim.tif")
      
      # find folder for variable k
      dir.in.k=dir.in
      if(all(list.files(dir.in)!=file.nm))
        if(any(list.dirs(dir.in,full.names=F)==k)) dir.in.k=paste0(dir.in,k,'/')
      
      # download file
      if(all(list.files(dir.in.k)!=file.nm))
        
        if(download)
          
          download.file(url=paste0('http://gentree.data.inra.fr/climate/datasets/',k,'/',file.nm),destfile=paste0(dir.in.k,file.nm),method='auto',cacheOK=F) else
            
            stop(paste("the following file",file.nm,"was not found in 'dir.in'. Verify the 'dir.in' argument or set download=T to automatically download climate files"))
      
      # load file
      stk=stack(paste0(dir.in.k,file.nm))
      names(stk)=1901:2014
      
      # extract yearly data at each point by applying a buffer (return average value)
      clim.k.t=extract(stk,xy[in.tile,],buffer=radius,fun=mean)
      clim.k[in.tile,]=clim.k.t
      
      # release memory
      rm(stk,clim.k.t)
      
    }
    
    # convert units (to °C or mm)
    clim.k=0.1*round(clim.k)
    
    # add coordinates
    clim.k=data.frame(xy,clim.k)
    
    # save
    if(!is.null(dir.out)) {
      
      nm.data=paste0(dir.out,"extraction_",k,"_",nrow(xy),"points",ifelse(!is.null(radius),paste0("_buffer",radius),""),"_eumedclim.RData")
      save(clim.k,file=nm.data)
      
    }
    
    return(clim.k)
    
  },mc.cores=n.core)
  
  # output
  names(list.ts)=clim.var
  return(list.ts)

}

### ARGUMENTS
# 'xy' matrix or data.frame. Coordinates of sampling points where to extract climate data (two-columns coordinates of longitude and latitude in WGS84).
# 'clim.var' character. Vector of selected EuMedClim climatic variables. All are selected by default.
# 'radius' numeric. The radius of a buffer (in meters) around each point from which to extract cell values. If the distance between the sampling point and the center of a cell is less than or equal to the buffer, the cell is included. If radius=NULL (by default), no buffer is applied.
# 'dir.in' character. Full path of the folder where climate files are stored. If files are missing, set download=T to download them.
# 'download' logical. If TRUE, required climate files will be downloaded from the web and stored (full path of the destination folder given in 'dir.in'). FALSE by default.
# 'dir.out' character (optional). Full path of the folder where to save climate data (RData format). If more than one climatic variable is selected, one file will be created for each. If NULL no saving (by default).
# 'n.core' numeric. Number of CPU cores to use (no parallelization by default). Useful when selecting several climatic variables. Use 'detectCores()-1' to know the number of available cores.


### EXAMPLES
## Get annual temperature (bio1) and precipitation (bio12) time series in Bordeaux and Madrid

# data.frame of sampling point coordinates
coord=data.frame(longitude=c(-0.57,-3.70),latitude=c(44.82,40.41))
rownames(coord)=c('bordeaux','madrid')

# select climatic variables
para=c('bio1','bio12')

# path of the folder where climate data is stored (or will be stored if missing).
dir.in=getwd() # here your working directory

# path of the folder where you want to save climate data (NULL for no saving)
dir.out=getwd() # here your working directory

# select number of available CPU cores
#n.core=detectCores() - 1

# run
data.emc=get.climate(xy=coord,clim.var=para,dir.in=dir.in,dir.out=dir.out,download=T)
head(data.emc)

# plot time series
par(mfrow=c(2,1))
par(mar=c(1.5,3,2,.5),mgp=c(2,.4,0))
for(k in para) {
  
  # load saved file (object named 'clim.k')
  #N=nrow(coord)
  #radius.val=ifelse(exists('radius'),ifelse(!is.null(radius),paste0("_buffer",ifelse(radius>5000,radius,5000)),""),"")
  #load(paste0(dir.out,"extraction_",k,"_",N,"points",radius.val,"_eumedclim.RData"))
  
  # get corresponding data.frame in the list
  clim.k=data.emc[[k]]
  
  # do not consider xy coordinates in computing the range
  ylim=range(clim.k[,-c(1:2)])

  # graph
  plot(c(1901,2014),ylim,t='n',xlab="",ylab=k)
  lines(1901:2014,clim.k['bordeaux',-c(1:2)],col=4)
  lines(1901:2014,clim.k['madrid',-c(1:2)],col=2)
  legend('topleft',bty='n',lty=1,col=c(4,2),legend=c('Bordeaux','Madrid'))
  
}
##################################
######### EXAMPLES ###############   #Added by Alex on the 15 December 2017
##################################

          #### How to use the function to extract the information for all my plotcodes ####

##### Obtain dataframe made of the plotcode, the latitude and the longitude and save it 

load("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/treefinal.RData") #the old version name old_with_duplicate now
treefinal$newplotcode = paste0(treefinal$country,treefinal$plotcode) #new line 30/01/2018
treefinal$newtreecode = paste0(treefinal$country,treefinal$treecode) #new line 30/01/2018
save(treefinal,file="/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/treefinal_Plot_newplot.RData") #this datatable is the only one having both the original plotcode and the new attributed one
treefinal$plotcode <- treefinal$newplotcode
treefinal$treecode <- treefinal$newtreecode
treefinal <- treefinal[,c(1:49)] 
save(treefinal,file="/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/treefinal.RData") #this datatable is the one we are going to use now. With the modified plot info 
# Problem with this. Several names in common # ll=split(treefinal[,c('plotcode','latitude','longitude')],treefinal[,'plotcode'])
ll=split(treefinal[,c('newplotcode','latitude','longitude')],treefinal[,'newplotcode']) #new line 30/01/2018
xy=do.call(rbind,lapply(ll,function(x) x[1,]))
save(xy,file="/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/plot.coord.RData")

# data.frame of sampling point coordinates
coord=data.frame(longitude=c(xy$longitude),latitude=c(xy$latitude))
rownames(coord)=c(xy$newplotcode)

# select climatic variables

bioclim=paste("bio",c(1,2,5,6,12,13,14),sep="")
T.seas=paste0("tmean.",c("djf","mam","jja","son"))
P.seas=paste0("prec.",c("djf","mam","jja","son"))
pet=paste0("pet.",c("mean","min","max"))
ppet=paste0("ppet.",c("mean","min","max"))
para=c(bioclim,T.seas,P.seas,pet,ppet)

# path of the folder where climate data is stored (or will be stored if missing).
dir.in=getwd() # here your working directory
setwd(dir = "/home/achangenet/Documents/FUNDIV\ -\ NFI\ -\ Europe/our-data/climate")
# path of the folder where you want to save climate data (NULL for no saving)
dir.out=getwd() # here your working directory

# select number of available CPU cores
#n.core=detectCores() - 1

# run
data.emc=get.climate(xy=coord,clim.var=para,dir.in=dir.in,dir.out=dir.out,download=T)
head(data.emc)

# plot time series
par(mfrow=c(2,1))
par(mar=c(1.5,3,2,.5),mgp=c(2,.4,0))
for(k in para) {
  
  # load saved file (object named 'clim.k')
  #N=nrow(coord)
  #radius.val=ifelse(exists('radius'),ifelse(!is.null(radius),paste0("_buffer",ifelse(radius>5000,radius,5000)),""),"")
  #load(paste0(dir.out,"extraction_",k,"_",N,"points",radius.val,"_eumedclim.RData"))
  
  # get corresponding data.frame in the list
  clim.k=data.emc[[k]]
  
  # do not consider xy coordinates in computing the range
  ylim=range(clim.k[,-c(1:2)])
  
  # graph
  plot(c(1901,2014),ylim,t='n',xlab="",ylab=k)
  lines(1901:2014,clim.k['bordeaux',-c(1:2)],col=4)
  lines(1901:2014,clim.k['madrid',-c(1:2)],col=2)
  legend('topleft',bty='n',lty=1,col=c(4,2),legend=c('Bordeaux','Madrid'))

