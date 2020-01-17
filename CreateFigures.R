rm(list=ls())

library(mgcv)
library(mapproj)
library(RNetCDF)
library(ncdf.tools)
library(lattice)
library(lmvar)
library(spherstat) #devtools::install_github("baddstats/spherstat")
#library(GeoModels)
library(globe)
library(fields)
library(sphereplot)
library(latticeExtra) # for plotting 
library(maps)         # for ... maps 
library(RColorBrewer)


setwd("~/Desktop/R_codes/Covariance/Final_code_github") #Working directory

#######################################################################
#Loading data
open_dat <-  open.nc("./data/pr_wtr.eatm.2017.nc");
global_dat <- read.nc(open_dat); 
dat <- global_dat$pr_wtr;

###################################################################
##################### Preparing Data ##############################
###################################################################

lat = global_dat$lat;
lon = global_dat$lon;
for(j in 1:length(lon)){ if(lon[j]>180){lon[j] = lon[j]-360;} }
sq_lat = seq(1,length(lat),by=1);
sq_lon = seq(1,length(lon),by=1);

coord = expand.grid(lon[sq_lon],lat[sq_lat]);

dat0 = rep(0,length(coord[,1]))
for(j in 1:365){
  dat0 = dat0 + c( dat[sq_lon,sq_lat,  j] );
}
dat0 = dat0/365;

ind = which(dat0 != "NaN");
coord = coord[ind,];
ini.dat1 = dat0[ind];

###################################################################
##################### Filtering Data ##############################
###################################################################

west.lon.deg= -180;
east.lon.deg= 180; 
south.lat.deg=0;
north.lat.deg=70;

sq_lon = which(coord[,1]>= west.lon.deg & coord[,1]<= east.lon.deg);
sq_lat = which(coord[,2]>  south.lat.deg & coord[,2]< north.lat.deg);
set0 <- intersect(sq_lon,sq_lat)

coord_reg <- coord*(pi/90)
coord_reg <- coord_reg[set0,]
filter_dat1 <- ini.dat1[set0]
regresion <- lm( filter_dat1 ~ cos(coord_reg[,2]) + sin(coord_reg[,2]));summary(regresion)

pdf("./images/residuals1.pdf")
plot(fitted(regresion),residuals(regresion),xlab="",ylab="",pch = 16,ylim = c(-25,25),xaxt="none", yaxt="none")
mtext("Residuals", side=2, line=2.2, cex=2)
mtext("Fitted values", side=1, line=4, cex=2)
axis(1, seq(0,50, 5),las=1, font=1,cex.axis=1.25)
axis(2, seq(-20,20,4),las=1, font=1,cex.axis=1.25)
dev.off()

my_design <- cbind(cos(coord_reg[,2]),sin(coord_reg[,2]),sin(0.5*coord_reg[,1])*sin(0.5*coord_reg[,2]))

my_design_mu <- my_design[,-3]
my_design_sigma <- my_design
regresion2 <- lmvar(y =  filter_dat1, X_mu = my_design_mu, X_sigma = my_design_sigma);
summary(regresion2)
mus <- cbind(rep(1,length(coord_reg[,2])),my_design_mu)%*%regresion2$coefficients_mu
sigmas <- exp(cbind(rep(1,length(coord_reg[,2])),my_design_sigma)%*%c(regresion2$coefficients_sigma))
regresion2$coefficients_mu
regresion2$coefficients_sigma

pdf("./images/residuals2.pdf")
plot(mus,residuals(regresion2)/sigmas,xlab="",ylab="", ylim = c(-5,5),pch = 16,xaxt="none", yaxt="none")
mtext("Residuals", side=2, line=2.2, cex=2)
mtext("Fitted values", side=1, line=4, cex=2)
axis(1, seq(0,50, 5),las=1, font=1,cex.axis=1.25)
axis(2, seq(-4,4,1),las=1, font=1,cex.axis=1.25)
dev.off()
filter_dat2 <- residuals(regresion2)/sigmas

###########

source("Plot_functions.R")

coord_test <- coord_reg*(90/pi)
coord_test <- cbind(coord_test[,1]+179,coord_test[,2])
filter_test <- filter_dat1

#Static plots
www <- sphwin(type = "band", param=c((pi*(90-70)/180),pi/2))
pal = colorRampPalette(c("blue", "white", "red"))
eps = NULL;dimyx = NULL;main = "";
pdf("./images/data000to180.pdf")
eye0 <- c(pi/2,-pi/6)
eyeE <- c(pi/2,eye0[2] - pi)
Plot_function(values = filter_dat1,my_coords = coord_test,w=www, eye= eye0, top = c(0,0), eyeE = eyeE, topE = c(0,0),sigma0 = 0.015,col.image=pal(100), ribb = FALSE)
globedrawlat(lat = 0)
globedrawlat(lat = 70)
dev.off()

pdf("./images/data090to270.pdf")
eye0 <- eye0 + c(0,pi/2)
eyeE <- c(pi/2,eye0[2] - pi)
Plot_function(values = filter_dat1,my_coords = coord_test,w=www, eye= eye0, top = c(0,0), eyeE = eyeE, topE = c(0,0), sigma0 = 0.015,col.image=pal(100),ribb = FALSE)
globedrawlat(lat = 0)
globedrawlat(lat = 70)
dev.off()

pdf("./images/data180to360.pdf")
eye0 <- eye0 + c(0,pi/2)
eyeE <- c(pi/2,eye0[2] - pi)
Plot_function(values = filter_dat1,my_coords = coord_test,w=www, eye= eye0, top = c(0,0), eyeE = eyeE, topE = c(0,0), sigma0 = 0.015,col.image=pal(100), ribb = FALSE)
globedrawlat(lat = 0)
globedrawlat(lat = 70)
dev.off()
pdf("./images/data270to090.pdf")
eye0 <- eye0 + c(0,pi/2)
Plot_function(values = filter_dat1,my_coords = coord_test,w=www, eye= eye0, top = c(0,0), eyeE = eyeE, topE = c(0,0), sigma0 = 0.015,col.image=pal(100), ribb = FALSE)
globedrawlat(lat = 0)
globedrawlat(lat = 70)
dev.off()

#Animation 1
sseeqq <- seq(-pi,pi,0.015) #105
sseeqq <- sseeqq[-c(107:111,309:314)] #avoding numerical problems on -90 and 90
eyeE0 <- spherstat:::Convert.globe( cbind(pi/2,sseeqq+pi) )$lon
N0 <- length(sseeqq)
for(i in 1:N0){
  id <- formatC(i, width = 6, format = "d", flag = "0")
  filename <- paste("Animation1", id, ".jpeg", collapse = "", sep = "")
  jpeg(filename)
  eye0 <- c(pi/2,sseeqq[i])
  eyeE <- c(pi/2,sseeqq[i] - pi)
  Plot_function(values = filter_dat1,my_coords = coord_test,w=www, eye= eye0 , top = c(0,0), eyeE = eyeE, topE = c(0,0),sigma0 = 0.015,col.image=pal(100))
  globedrawlat(lat = 0)
  globedrawlat(lat = 70)
  degree <- format(round(eyeE0[i],digits = 2),nsmall = 2)
  title(paste("Longitude", degree, sep = " "))
  dev.off()
}
#install ImageMagick
system("convert -delay 02 *.jpeg Animation1.gif")
file.remove(list.files(pattern=".jpeg"))


pdf("./images/mymap.pdf",height = 3.00,width = 5.70)
world.map <- map('world', plot=FALSE, boundary=TRUE,
                 xlim=c(west.lon.deg,east.lon.deg), 
                 ylim=c(south.lat.deg,north.lat.deg)) 
world.df <- data.frame(lon=world.map$x, lat=world.map$y) 
l1 = levelplot(c(filter_dat1) ~ (90*coord_reg[,1]/pi)*(90*coord_reg[,2]/pi), col.regions=pal(100), xlab="",ylab="", margin = list(FUN = 'median'),aspect = "iso", colorkey = list(space = "bottom"))
update(l1) + xyplot(lat ~ lon, world.df, type='l', lty=1, lwd=1.5, col='black', xlab="",ylab="") 
dev.off()

pdf("./images/filter000to180.pdf")
eye0 <- c(1.5,-pi/6)
eyeE <- c(1.5,eye0[2] - pi)
Plot_function(values = filter_dat2,my_coords = coord_test,w=www, eye= eye0, top = c(0,0), eyeE = eyeE, topE = c(0,0), sigma0 = 0.015,col.image=pal(100), ribb = FALSE)
globedrawlat(lat = 0)
globedrawlat(lat = 70)
#Plot box 
line1x <- rep(180,l=200);line2x <- seq(180,120,l=200)
line1y <- seq(0,70,l=200);line2y <- rep(70,l=200)
line3x <- rep(120,0,l=200);line4x <- seq(120,180,l=200)
line3y <- seq(70,0,l=200);line4y <- rep(0,l=200)
globelines(loc=cbind(line1x,line1y),col="green")
globelines(loc=cbind(line2x,line2y),col="green")
globelines(loc=cbind(line3x,line3y),col="green")
globelines(loc=cbind(line4x,line4y),col="green")
dev.off()
pdf("./images/filter090to270.pdf")
eye0 <- eye0 + c(0,pi/2)
eyeE <- c(1.5,eye0[2] - pi)
Plot_function(values = filter_dat2,my_coords = coord_test,w=www, eye= eye0, top = c(0,0), eyeE = eyeE, topE = c(0,0), sigma0 = 0.015,col.image=pal(100), ribb = FALSE)
globedrawlat(lat = 0)
globedrawlat(lat = 70)
line1x <- rep(180,l=200);line2x <- seq(180,120,l=200)
line1y <- seq(0,70,l=200);line2y <- rep(70,l=200)
line3x <- rep(120,0,l=200);line4x <- seq(120,180,l=200)
line3y <- seq(70,0,l=200);line4y <- rep(0,l=200)
globelines(loc=cbind(line1x,line1y),col="green")
globelines(loc=cbind(line2x,line2y),col="green")
globelines(loc=cbind(line3x,line3y),col="green")
globelines(loc=cbind(line4x,line4y),col="green")
dev.off()
pdf("./images/filter180to360.pdf")
eye0 <- eye0 + c(0,pi/2)
eyeE <- c(1.5,eye0[2] - pi)
Plot_function(values = filter_dat2,my_coords = coord_test,w=www, eye= eye0, top = c(0,0), eyeE = eyeE, topE = c(0,0), sigma0 = 0.015,col.image=pal(100), ribb = FALSE)
globedrawlat(lat = 0)
globedrawlat(lat = 70)
dev.off()
pdf("./images/filter270to090.pdf")
eye0 <- eye0 + c(0,pi/2)
eyeE <- c(1.5,eye0[2] - pi)
Plot_function(values = filter_dat2,my_coords = coord_test,w=www, eye= eye0, top = c(0,0), eyeE = eyeE, topE = c(0,0), sigma0 = 0.015,col.image=pal(100), ribb = FALSE)
globedrawlat(lat = 0)
globedrawlat(lat = 70)
line1x <- rep(180,l=200);line2x <- seq(180,120,l=200)
line1y <- seq(0,70,l=200);line2y <- rep(70,l=200)
line3x <- rep(120,0,l=200);line4x <- seq(120,180,l=200)
line3y <- seq(70,0,l=200);line4y <- rep(0,l=200)
globelines(loc=cbind(line1x,line1y),col="green")
globelines(loc=cbind(line2x,line2y),col="green")
globelines(loc=cbind(line3x,line3y),col="green")
globelines(loc=cbind(line4x,line4y),col="green")
dev.off()

pdf("./images/mymap_filter.pdf",height = 3.00,width = 5.70)
world.map <- map('world', plot=FALSE, boundary=TRUE,
                 xlim=c(west.lon.deg,east.lon.deg), 
                 ylim=c(south.lat.deg,north.lat.deg)) 
world.df <- data.frame(lon=world.map$x, lat=world.map$y) 

l1 = levelplot(c(filter_dat2) ~ (90*coord_reg[,1]/pi)*(90*coord_reg[,2]/pi), col.regions=pal(100), xlab="",ylab="", margin = list(FUN = 'median'),aspect = "iso", colorkey = list(space = "bottom"), )
update(l1) + xyplot(lat ~ lon, world.df, type='l', lty=1, lwd=1.5, col='black', xlab="",ylab="") 
dev.off()

#Animation 2
sseeqq <- seq(-pi,pi,0.015) #105
sseeqq <- sseeqq[-c(107:111,309:314)] #avoding numerical problems on -90 and 90
eyeE0 <- spherstat:::Convert.globe( cbind(pi/2,sseeqq+pi) )$lon
N0 <- length(sseeqq)
for(i in 1:N0){
  id <- formatC(i, width = 6, format = "d", flag = "0")
  filename <- paste("Animation2", id, ".jpeg", collapse = "", sep = "")
  jpeg(filename)
  eye0 <- c(pi/2,sseeqq[i])
  eyeE <- c(pi/2,sseeqq[i] - pi)
  Plot_function(values = filter_dat2,my_coords = coord_test,w=www, eye= eye0 , top = c(0,0), eyeE = eyeE, topE = c(0,0),sigma0 = 0.015,col.image=pal(100))
  globedrawlat(lat = 0)
  globedrawlat(lat = 70)
  degree <- format(round(eyeE0[i],digits = 2),nsmall = 2)
  title(paste("Longitude", degree, sep = " "))
  #Plot box 
  line1x <- rep(180,l=200);line2x <- seq(180,120,l=200)
  line1y <- seq(0,70,l=200);line2y <- rep(70,l=200)
  line3x <- rep(120,0,l=200);line4x <- seq(120,180,l=200)
  line3y <- seq(70,0,l=200);line4y <- rep(0,l=200)
  globelines(loc=cbind(line1x,line1y),col="green")
  globelines(loc=cbind(line2x,line2y),col="green")
  globelines(loc=cbind(line3x,line3y),col="green")
  globelines(loc=cbind(line4x,line4y),col="green")
  dev.off()
}
#install ImageMagick
system("convert -delay 02 *.jpeg Animation2.gif")
file.remove(list.files(pattern=".jpeg"))
