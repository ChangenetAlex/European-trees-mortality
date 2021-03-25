# Alex, 13/04/2018
# Autocor scripts based on several packages 

########################
##### Load my data #####
########################

rm(list = ls())
gc()
dev.off()
library(fields)
library(RandomFields)
library(lctools)
library(pgirmess)
CODE = "PINSYL" 
seuil = 0.875
Dir = c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/CLIMAP/Models")) # Directory 
setwd(Dir)
list.files(Dir,pattern = paste0(seuil,".rds"))
dfplot2 <- readRDS(paste0("dfplot2",CODE,seuil,".rds")) #Base de données plot
Dir = c(paste0("/home/achangenet/Documents/FUNDIV - NFI - Europe/our-data/species/",CODE,"/AutoCor")) # Directory 
setwd(Dir)
par(mfrow=c(1,2))

################################
##### Two useful functions #####
###### Euclidian distance  #####
################################

disteucl <- function(crds) {
  if(is.vector(crds)) crds <- as.matrix(crds)
  dist <- matrix(0, nrow(crds), nrow(crds))
  for(i in 1:ncol(crds)) {
    dist <- dist + outer(crds[, i], crds[, i], "-")^2
  }
  dist <- sqrt(dist)
  return(dist)}
# distance euclidienne entre un point p et un ensemble de points
distp <- function(crdsp, crds) {
  if(is.vector(crds)) crds <- as.matrix(crds)
  dist <- rep(0, nrow(crds))
  for(i in 1:ncol(crds)) {
    dist <- dist + (rep(crdsp[i],nrow(crds)) - crds[,i])^2
  }
  dist <- sqrt(dist)
  return(dist)}

#########################################
##### Visualise the mattern funtion #####
#########################################

# Effect of rho
par(mfrow=(c(1,2)))
d<- seq( 0,10,,200)
y<- Matern( d, range=0.5069, smoothness=16.6)
y1<- Matern( d, range=1.5, smoothness=2)
y2<- Matern( d, range=1.5, smoothness=5)
y3<- Matern( d, range=1.5, smoothness=10)
y4<- Matern( d, range=1.5, smoothness=100)
y<- cbind(y,y1,y2,y3,y4)
matplot( d, y, type="l", lty=1, lwd=2,main=c(expression(paste(rho,"= 1,5"))),ylab="Covariance",xlab="distance")
legend("topright", legend = c(expression(paste(nu," = 1")),expression(paste(nu," = 2"))
                              ,expression(paste(nu," = 5")),expression(paste(nu," = 10"))
                              ,expression(paste(nu," = 100"))), col = 1:5, lty = 1)

# Effect of mu
y<- Matern( d, range=0.1, smoothness=3)
y1<- Matern( d, range=0.5, smoothness=3)
y2<- Matern( d, range=1, smoothness=3)
y3<- Matern( d, range=1.5, smoothness=3)
y4<- Matern( d, range=3, smoothness=3)
y5<- Matern( d, range=10, smoothness=3)
y<- cbind(y,y1,y2,y3,y4,y5)
matplot( d, y, type="l", lty=1, lwd=2,main=c(expression(paste(nu,"= 3"))),ylab="Covariance",xlab="distance")
legend("topright", legend = c(expression(paste(rho," = 0.1")),expression(paste(rho," = 0.5"))
                              ,expression(paste(rho," = 1")),expression(paste(rho," = 1.5"))
                              ,expression(paste(rho," = 3")),expression(paste(rho," = 10"))), col = 1:6, lty = 1)

# Change the parameter values to understand effects
x <- 1:100
y <- 1:100
model <- RMmatern(nu=100,var=1,scale=10)
simu <- RFsimulate(model,x,y)
plot(simu,col=tim.colors(64))

########################
##### 1st approach #####
########################
# Plot my mortality values according to geographical distance
plot(dfplot2$latitude,dfplot2$longitude,typ="n")
symbols(dfplot2$latitude,dfplot2$longitude,circles=dfplot2$sp.mortality.plot.count/100,
        inches=F,add=T,bg=3)
# Same thing with a random sample of size= n 
n = 1000
sampled=sample(1:nrow(dfplot2), n)
dfplot3=dfplot2[sampled,]
plot(dfplot2$latitude,dfplot2$longitude,typ="n",xlim=c(min(dfplot2$latitude),max(dfplot2$latitude)),ylim=c(min(dfplot2$longitude),max(dfplot2$longitude)),
     main=paste0("Sample = ",n))
symbols(dfplot3$latitude,dfplot3$longitude,circles=dfplot3$sp.mortality.plot.count/100,
        inches=F,add=T,bg=3)
dev.print(file="MortSpatiale.jpeg",device=jpeg,width=1100)

# Variographie
# Distance and densities
par(mfrow=c(1,1))
distsite <- disteucl(dfplot3[,c(3,2)])
max(distsite)
z2 <- (disteucl(dfplot3$sp.mortality.plot.count)**2)/2.
plot(distsite,z2)
dev.print(file="Hist.Vario.jpeg",device=jpeg,width=1100)
# distclas definit la discretisation du variogramme experimental
# par des plages de distances constantes (classes de distances)
# variogram
distclas <- cut(distsite,breaks=c(-0.01,seq(1,46,2)))
distvario1 <- tapply(distsite,distclas,mean)
vario1 <- tapply(z2,distclas,mean,na.rm=T)
nbcouples1 <- tapply(z2,distclas,length)
vario1d <- cbind(vario1,distvario1,nbcouples1)
plot(vario1d[1:nrow(vario1d),2:1],typ="b") # This is my first vario 
dev.print(file="VarioManuel.jpeg",device=jpeg,width=1100)
# fonction du vario et cov asssociee
varexp <- function(h,c0,c1,scal)
{f <- c0 + c1 * (1-exp(-h/scal))
f[h==0] <- 0
f}
covexp <- function(h,c0,c1,scal)
{f <- c1 * (exp(-h/scal))
f[h==0] <- c0 + c1
f}
#variogramme exponentil
vario1dat <- as.data.frame(vario1d)
varioajust1 <- nls(vario1 ~ varexp(distvario1,c0,c1,scal),
                   data=vario1dat,
                   start=list(c0=1,c1=40,scal=0))
c0a <- coef(varioajust1)[1]
c1a <- coef(varioajust1)[2]
scala <- coef(varioajust1)[3]
# plot vario and ajust function 
plot(vario1d[1:nrow(vario1d),2:1],typ="b")
lines(seq(0.01,400,5),varexp(seq(0.01,400,5),c0a,c1a,scala),
      col=2,lwd=2)
############################################################
### FAIRE LA MEME AVEC FUNCTION MATERN (et CF spherique)####
############################################################


#####################
### Moran stats  ####
#####################

#Ape
distest <- as.matrix(dist(cbind(dfplot3$longitude, dfplot3$latitude)))
distestinv <- 1/distest
distestinv[is.infinite(distestinv)] <- 0
diag(distestinv) <- 0
Moran.I(dfplot3$mortality.plot.rate, distestinv) # Save these stats
capture.output(Moran.I(dfplot3$mortality.plot.rate, distestinv),file="MoranTest.Ape.txt")
A <- variog(coords = as.matrix(cbind(dfplot3$longitude, dfplot3$longitude)), data = dfplot3$mortality.plot.rate)
plot(A, type="l", lty=2)
dev.print(file="Vario.Ape.jpeg",device=jpeg,width=1100)

#lctools
n = 3000
sampled=sample(1:nrow(dfplot2), n)
dfplot3=dfplot2[sampled,]
Coords <- cbind(dfplot3$longitude, dfplot3$latitude) ## Coordonee
bws <- c(2,3,5, 10, 15, 20,30,40,50,100) # Neighbors
moran <- lctools::moransI.v(Coords, bws, dfplot3$mortality.plot.rate, WType='Binary') #par defaut = binary et neighbors
dev.print(file="Vario.LCtools.jpeg",device=jpeg,width=1100)
capture.output(lctools::moransI.v(Coords, bws, dfplot3$mortality.plot.rate, WType='Binary'),file="MoranTest.LCTools")
bws <- c(10,20,30,40,42) # Kernels distance
moran <- lctools::moransI.v(Coords, bws, dfplot3$mortality.plot.rate, WType='Binary',family="fixed")
l.moran<-lctools::l.moransI(Coords,6,dfplot3$mortality.plot.rate) # scatterplot 

# Vario with Correlog
corD <- correlog(Coords, dfplotVARIO$sp.mortality.plot.rate, method = "Moran")
plot(corD)
dev.print(file="Vario.Correlog.jpeg",device=jpeg,width=1100)
