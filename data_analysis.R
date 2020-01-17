rm(list=ls())

library(fields)
library(mgcv)
library(mapproj)
library(RNetCDF)
library(ncdf.tools)
library(lattice)

setwd("/Users/jeffpostdoc/Desktop/R_codes/Covariance/Upload_folder/")

system("rm ./c/matriz_cov.so")
system("rm ./c/matriz_cov.o")
system("rm ./c/mat_pred.so")
system("rm ./c/mat_pred.o")
system("R CMD SHLIB ./c/matriz_cov.c  -lgsl -larb -lflint")
system("R CMD SHLIB ./c/mat_pred.c  -lgsl -larb -lflint")

source("./r/Functions.R");
source("./r/MLE.R");
source("./r/matriz_cov.R");
source("./r/mat_pred.R");

dyn.load("./c/matriz_cov.so")
dyn.load("./c/mat_pred.so")
###############################################################
#############Loading and reparing data ########################
###############################################################
open_dat <-  open.nc("./data/pr_wtr.eatm.2017.nc");

global_dat <- read.nc(open_dat); 
dat <- global_dat$pr_wtr;

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

west.lon.deg= -180;
east.lon.deg= 180; 
south.lat.deg=0;
north.lat.deg=70;

sq_lon = which(coord[,1]>= west.lon.deg & coord[,1]<= east.lon.deg);
sq_lat = which(coord[,2]>  south.lat.deg & coord[,2]< north.lat.deg);
set0 <- intersect(sq_lon,sq_lat)

###################################################################
############### Filtering data ####################################
###################################################################

coord_reg <- coord*(pi/90)
coord_reg <- coord_reg[set0,]
filter_dat1 <- ini.dat1[set0]
regresion <- lm( filter_dat1 ~ cos(coord_reg[,2]) + sin(coord_reg[,2]));summary(regresion)

plot(fitted(regresion),residuals(regresion),xlab="Fitted values",ylab="Residuals",pch = 16,ylim = c(-25,25))
my_design <- cbind(cos(coord_reg[,2]),sin(coord_reg[,2]),sin(0.5*coord_reg[,1])*sin(0.5*coord_reg[,2]))
my_design_mu <- my_design[,-3]
my_design_sigma <- my_design

regresion2 <- lmvar::lmvar(y =  filter_dat1, X_mu = my_design_mu, X_sigma = my_design_sigma);

summary(regresion2)
mus <- cbind(rep(1,length(coord_reg[,2])),my_design_mu)%*%regresion2$coefficients_mu
sigmas <- exp(cbind(rep(1,length(coord_reg[,2])),my_design_sigma)%*%c(regresion2$coefficients_sigma))
regresion2$coefficients_mu
regresion2$coefficients_sigma

plot(mus,residuals(regresion2)/sigmas,xlab="Fitted values",ylab="Residuals", ylim = c(-5,5),pch = 16)
dev.off()
datos <- residuals(regresion2)/sigmas

###################################################################
############# Verifying anisotropy assumption #####################
###################################################################
#Memory consuming. Results were stored in the file "Data_results.Rdat"
#load("Data_results.Rdat")
source("Investigating_anisotropy.R")

###################################################################
###################################################################
#The training set 
set1 =  which( (90/pi)*coord_reg[,1] < 120 )

###################################################################
###################################################################
############# NOT RUN: it takes a lot of time (and memory). Results were stored in the file "Data_results.Rdat"
#Random choice
n.sim <- 512 
N <- rr1 <- 200;
niter <- 100;

scores <- list()
avg_ind <- matrix(0,nrow = n.sim,ncol=6)

#Setting server options

library("foreach")
library("parallel")
library("doParallel")
library("doRNG") #reproducible


cl <- makeCluster(64)
showConnections(all=TRUE)
registerDoParallel(cl)

set.seed(13) #Setting a seed after repeat several times

Iter_result_fix <- foreach(i = 1:n.sim)%dorng%{
  dyn.load("./c/matriz_cov.so")
  dyn.load("./c/mat_pred.so")
  ss1 = sample(set1, rr1)
  datos = dat1[ss1];
  coord2 = coord[set0,][ss1,]*(pi/180);
  lon = coord2[,1];
  lat = coord2[,2];
  ##################################################################
  ##################################################################
  #eleccion variable values:
  #0.- F-family (ARBlib); 1.- F-family (GSL); 2.- Chordal-Matern; 3.- Circular-Matern;
  inic0 <- c(1,0.5, 0.5);
  
  # ML estimation
  eleccion = 1;
  res11 <- optim(inic0, loglik , control=list(reltol=1e-14,maxit=30000,trace=F), hessian=TRUE);
  eleccion = 2;
  res22 <- optim(inic0, loglik , control=list(reltol=1e-14,maxit=30000,trace=F), hessian=TRUE);
  eleccion = 3;
  res33 <- optim(inic0, loglik , control=list(reltol=1e-14,maxit=30000,trace=F), hessian=TRUE);  
  
  #ML fix_estimate
  theta_vec <- rbind(res11$par,res22$par,res33$par)
  loglik_vec <- c(res11$value,res22$value,res33$value)
  ###################################################################
  ###################################################################
  
  predict_scores <- matrix(0,ncol=4,nrow=3)
  election_vec <- c(1,2,3)
  for(j in 1:3){
    mm1 = matriz_cov(eleccion = election_vec[j],N,theta_vec[j,],lon,lat);
    decompvarcov1 = chol(mm1);
    inv1 = chol2inv(decompvarcov1);
    predic1 = Prscores(datos,inv1,N);    #drop-one prediction
    predict_scores[j,] <- c(unlist(predic1))
  }
  
  mm1 = matriz_cov(eleccion = 1,N,theta_vec[1,],lon,lat);
  mm2 = matriz_cov(eleccion = 2,N,theta_vec[2,],lon,lat);
  mm3 = matriz_cov(eleccion = 3,N,theta_vec[3,],lon,lat);
  
  decompvarcov1 = chol(mm1 );
  decompvarcov2 = chol(mm2 );
  decompvarcov3 = chol(mm3 );
  
  inv1   = chol2inv(decompvarcov1 )
  inv2   = chol2inv(decompvarcov2 )
  inv3   = chol2inv(decompvarcov3 )
  
  RMSE <- MAES <- matrix(0,nrow=niter,ncol=3);    
  
  for(j in 1:niter){
    rr2 <- K <- 20;
    set2 = which(  (90/pi)*coord_reg[,1] >= 120) 
    
    ss2 = sample(set2, rr2)
    
    datos0 = dat1[ss2];
    coord0 = coord[set0,][ss2,]*(pi/180);
    lon0 = coord0[,1];
    lat0 = coord0[,2];
    
    mat1 =  mat_pred(1,N,K,theta_vec[1,],lon,lat,lon0,lat0)    
    mat2 =  mat_pred(2,N,K,theta_vec[2,],lon,lat,lon0,lat0)
    mat3 =  mat_pred(3,N,K,theta_vec[3,],lon,lat,lon0,lat0)
    
    predicciones1 = mat1%*%inv1%*%datos;
    predicciones2 = mat2%*%inv2%*%datos;
    predicciones3 = mat3%*%inv3%*%datos;
    
    aa1 = sum((predicciones1-datos0)^2)/K
    aa2 = sum((predicciones2-datos0)^2)/K
    aa3 = sum((predicciones3-datos0)^2)/K
    
    bb1 = sum(abs(predicciones1-datos0))/K
    bb2 = sum(abs(predicciones2-datos0))/K
    bb3 = sum(abs(predicciones3-datos0))/K
    
    RMSE[j,] =  c(sqrt(aa1),sqrt(aa2),sqrt(aa3))
    MAES[j,] = c(bb1,bb2,bb3)
  }
  
  avg_RMSE <- colMeans(RMSE)
  avg_MAES <- colMeans(MAES)
  
  result <- list(avg_MAES = avg_MAES , avg_RMSE = avg_RMSE,params_mat = theta_vec, score_mat = predict_scores, loglik = loglik_vec )
  return(result)
} 
closeAllConnections()
#stopCluster(cl)
#stopImplicitCluster()
#load("Data_results.Rdat")

#Parameter estimates
Ffamily.par <- do.call(rbind,lapply(Iter_result_fix,function(x) x$params_mat[1,]))
Matern.par <- do.call(rbind,lapply(Iter_result_fix,function(x) x$params_mat[2,]))
C_Matern.par <- do.call(rbind,lapply(Iter_result_fix,function(x) x$params_mat[3,]))

colMeans(Ffamily.par)
colMeans(Matern.par)
colMeans(C_Matern.par)

#log likelihood
loglik.all <- do.call(rbind,lapply(Iter_result_fix,function(x) x$loglik))
round(apply(loglik.all,2,mean),digits = 3)
round(apply(loglik.all,2,sd),digits = 3)
table(apply(loglik.all[,1:3],1,which.max))

#MAE and RMSE
avg_RMSE <- do.call(rbind,lapply(Iter_result_fix,function(x) x$avg_RMSE))
boxplot(avg_RMSE[,c(1,3,2)],ylab="", yaxt="none",names = c("F-Family","Circular Matérn","Chordal Matérn"),las = 1,cex.axis=1.25)
axis(2, seq(0,1, 0.1),las=1, font=1,cex.axis=1.25)
round(colMeans(avg_RMSE[,c(1,3,2)]),digits = 3)

avg_MAE <- do.call(rbind,lapply(Iter_result_fix,function(x) x$avg_MAES))
boxplot(avg_MAE[,c(1,3,2)],ylab="", yaxt="none",names = c("F-Family","Circular Matérn","Chordal Matérn"),las = 1,cex.axis=1.25)
axis(2, seq(0,1, 0.1),las=1, font=1,cex.axis=1.25)
round(colMeans(avg_MAE[,c(1,3,2)]),digits = 3)

########################
#### Model checking ####
########################
Ffam.par <- colMeans(Ffamily.par)[1:3]
Cmat.par <- colMeans(C_Matern.par)[1:3]
Mfam.par <- colMeans(Matern.par)[1:3]

#m0 = matriz_cov(eleccion = 0,100,Ffam.par,rep(0,100),seq(0,pi,l=100))[1,]; #Testing
m1 = matriz_cov(eleccion = 1,100,Ffam.par,rep(0,100),seq(0,pi,l=100))[1,];
m2 = matriz_cov(eleccion = 2,100,Cmat.par,rep(0,100),seq(0,pi,l=100))[1,];
m3 = matriz_cov(eleccion = 3,100,Mfam.par,rep(0,100),seq(0,pi,l=100))[1,];

#black = F family; red = Chordal Matern; green = Circular matern. 
#matplot(seq(0,pi,l=200),cbind(m1[1]-m1,m2[1]-m2,m3[1]-m3),type="l",col=c(1,2,3),lty = c(1,1,1),add = TRUE)
pdf("fitted_variogram.pdf")
plot(c(0,Variog.raw.new$hvec),c(0,Variog.raw.new$variog),type="l",lwd = 3, ylim = c(0,1.5), xlim = c(0,1), xlab = "", ylab = "", xaxt = "none", yaxt = "none")
matplot(seq(0,pi,l=100),cbind(m1[1]-m1,m2[1]-m2,m3[1]-m3),type="l",lty = c(1,1,1), xlim = c(0,1),add = TRUE,col=c("blue","green","red"),lwd = 3)
mtext(expression(gamma(theta)), side=2, line=2.2, cex=2)
mtext(expression(theta), side=1, line=4, cex=2)
axis(1, seq(0,1, l=5),las=1, font=1,cex.axis=1.25)
axis(2, seq(0,1.5,0.5),las=1, font=1,cex.axis=1.25)
dev.off()


