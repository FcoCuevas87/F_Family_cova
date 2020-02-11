#####################################
##### Global Rank Envelope Test #####
#####################################
load("Data_results.Rdat")
NNN <- nrow(coord0)
mm1 <- matriz_cov(eleccion = 1,NNN,Ffam.par,coord0[,1]*(pi/180),coord0[,2]*(pi/180))
mm2 <- matriz_cov(eleccion = 2,NNN,Cmat.par,coord0[,1]*(pi/180),coord0[,2]*(pi/180))
mm3 <- matriz_cov(eleccion = 3,NNN,Mfam.par,coord0[,1]*(pi/180),coord0[,2]*(pi/180))

chol01 <- chol(mm1)
chol02 <- chol(mm2)
chol03 <- chol(mm3)

Indexes1 <- findInterval(coord0[,2], c(35))  #Lat cut
Indexes2 <- findInterval(coord0[,1], c(-1.25)) #Long cut


Variog.iso <- Variog.lat.est <- Variog.lon.est <- list()
Variog.iso$ffam <- Variog.iso$cmat <- Variog.iso$mfam <- matrix(0,ncol=29,nrow = 12495)
Variog.lat.est$ffam <- Variog.lat.est$cmat <- Variog.lat.est$mfam <- list()
Variog.lon.est$ffam <- Variog.lon.est$cmat <- Variog.lon.est$mfam <- list()

Variog.lat.est$ffam[[1]] <- Variog.lat.est$cmat[[1]] <- Variog.lat.est$mfam[[1]] <- matrix(0,ncol=19,nrow = 12495)
Variog.lon.est$ffam[[1]] <- Variog.lon.est$cmat[[1]] <- Variog.lon.est$mfam[[1]] <- matrix(0,ncol=19,nrow = 12495)

Variog.lat.est$ffam[[2]] <- Variog.lat.est$cmat[[2]] <- Variog.lat.est$mfam[[2]] <- matrix(0,ncol=19,nrow = 12495)
Variog.lon.est$ffam[[2]] <- Variog.lon.est$cmat[[2]] <- Variog.lon.est$mfam[[2]] <- matrix(0,ncol=19,nrow = 12495)

RF <- list(X = dat1.tmp, longitude = coord0[,1], latitude = coord0[,2])
xxx <- Est.variog(RF,h.max = pi/3, l = 30)
plot(xxx$hvec,xxx$variog)

set.seed(123)

for(j in 1:12495){
  #Simulating GRF
  dat1.tmp <- chol01%*%c(rnorm(nrow(coord0)))
  dat2.tmp <- chol02%*%c(rnorm(nrow(coord0)))
  dat3.tmp <- chol03%*%c(rnorm(nrow(coord0)))
  
  #Creating lists
  RF1 <- list(X = dat1.tmp, longitude = coord0[,1], latitude = coord0[,2])
  RF2 <- list(X = dat2.tmp, longitude = coord0[,1], latitude = coord0[,2])
  RF3 <- list(X = dat3.tmp, longitude = coord0[,1], latitude = coord0[,2])
  
  #Computing isotropic variogram fixed
  Variog.iso$ffam[j,] <- Est.variog(RF1,h.max = pi/3, l = 30)$variog
  Variog.iso$cmat[j,] <- Est.variog(RF2,h.max = pi/3, l = 30)$variog
  Variog.iso$mfam[j,] <- Est.variog(RF3,h.max = pi/3, l = 30)$variog

  for(i in 0:1){
    RF1 <- list(X = (dat1.tmp[Indexes1==i]), longitude = coord0[Indexes1==i,1], latitude = coord0[Indexes1==i,2])
    RF2 <- list(X = (dat2.tmp[Indexes1==i]), longitude = coord0[Indexes1==i,1], latitude = coord0[Indexes1==i,2])
    RF3 <- list(X = (dat3.tmp[Indexes1==i]), longitude = coord0[Indexes1==i,1], latitude = coord0[Indexes1==i,2])
    Variog.lat.est$ffam[[i+1]][j,] <- Est.variog(RF1,h.max = pi/3, l = 20)$variog
    Variog.lat.est$cmat[[i+1]][j,] <- Est.variog(RF2,h.max = pi/3, l = 20)$variog
    Variog.lat.est$mfam[[i+1]][j,] <- Est.variog(RF3,h.max = pi/3, l = 20)$variog
  }
  
  for(i in 0:1){
    RF1 <- list(X = (dat1.tmp[Indexes2==i]), longitude = coord0[Indexes2==i,1], latitude = coord0[Indexes2==i,2])
    RF2 <- list(X = (dat2.tmp[Indexes2==i]), longitude = coord0[Indexes2==i,1], latitude = coord0[Indexes2==i,2])
    RF3 <- list(X = (dat3.tmp[Indexes2==i]), longitude = coord0[Indexes2==i,1], latitude = coord0[Indexes2==i,2])
    Variog.lon.est$ffam[[i+1]][j,] <- Est.variog(RF1,h.max = pi/3, l = 20)$variog
    Variog.lon.est$cmat[[i+1]][j,] <- Est.variog(RF2,h.max = pi/3, l = 20)$variog
    Variog.lon.est$mfam[[i+1]][j,] <- Est.variog(RF3,h.max = pi/3, l = 20)$variog
  }
  cat(j,"\n")
}
#base::save.image("Global_rank_env.Rdat")

#devtools::install_github('myllym/GET')
require("GET")

sim.to.env <- function(hvec = hvec, obs = obs, sims = sims){
  #set.curves <- create_curve_set( list(r = hvec[-1], obs = obs[-1] , sim_m = t(sims[,-1]) ))
  #envs.obs <- global_envelope_test(set.curves, type = "rank", alternative = "two.sided")
  envs.obs <- create_curve_set( list(r = hvec[-1], obs = obs[-1] , sim_m = t(sims[,-1]) ))
  return(envs.obs)
}

gre.curves.ffam <- gre.curves.cmat <- gre.curves.mfam <- list()
gre.curves.ffam[[1]] <- sim.to.env(hvec = Variog.raw.new$hvec, obs = Variog.raw.new$variog, sims = Variog.iso$ffam)
gre.curves.ffam[[2]] <- sim.to.env(hvec = Variog.fil.est$lat[[1]]$hvec, obs = Variog.fil.est$lat[[1]]$variog, sims = Variog.lat.est$ffam[[1]])
gre.curves.ffam[[3]] <- sim.to.env(hvec = Variog.fil.est$lat[[2]]$hvec, obs = Variog.fil.est$lat[[2]]$variog, sims = Variog.lat.est$ffam[[2]])
gre.curves.ffam[[4]] <- sim.to.env(hvec = Variog.fil.est$lon[[1]]$hvec, obs = Variog.fil.est$lon[[1]]$variog, sims = Variog.lon.est$ffam[[1]])
gre.curves.ffam[[5]] <- sim.to.env(hvec = Variog.fil.est$lon[[2]]$hvec, obs = Variog.fil.est$lon[[2]]$variog, sims = Variog.lon.est$ffam[[2]])

gre.curves.cmat[[1]] <- sim.to.env(hvec = Variog.raw.new$hvec, obs = Variog.raw.new$variog, sims = Variog.iso$cmat)
gre.curves.cmat[[2]] <- sim.to.env(hvec = Variog.fil.est$lat[[1]]$hvec, obs = Variog.fil.est$lat[[1]]$variog, sims = Variog.lat.est$cmat[[1]])
gre.curves.cmat[[3]] <- sim.to.env(hvec = Variog.fil.est$lat[[2]]$hvec, obs = Variog.fil.est$lat[[2]]$variog, sims = Variog.lat.est$cmat[[2]])
gre.curves.cmat[[4]] <- sim.to.env(hvec = Variog.fil.est$lon[[1]]$hvec, obs = Variog.fil.est$lon[[1]]$variog, sims = Variog.lon.est$cmat[[1]])
gre.curves.cmat[[5]] <- sim.to.env(hvec = Variog.fil.est$lon[[2]]$hvec, obs = Variog.fil.est$lon[[2]]$variog, sims = Variog.lon.est$cmat[[2]])

gre.curves.mfam[[1]] <- sim.to.env(hvec = Variog.raw.new$hvec, obs = Variog.raw.new$variog, sims = Variog.iso$mfam)
gre.curves.mfam[[2]] <- sim.to.env(hvec = Variog.fil.est$lat[[1]]$hvec, obs = Variog.fil.est$lat[[1]]$variog, sims = Variog.lat.est$mfam[[1]])
gre.curves.mfam[[3]] <- sim.to.env(hvec = Variog.fil.est$lat[[2]]$hvec, obs = Variog.fil.est$lat[[2]]$variog, sims = Variog.lat.est$mfam[[2]])
gre.curves.mfam[[4]] <- sim.to.env(hvec = Variog.fil.est$lon[[1]]$hvec, obs = Variog.fil.est$lon[[1]]$variog, sims = Variog.lon.est$mfam[[1]])
gre.curves.mfam[[5]] <- sim.to.env(hvec = Variog.fil.est$lon[[2]]$hvec, obs = Variog.fil.est$lon[[2]]$variog, sims = Variog.lon.est$mfam[[2]])

P_value.ffam <- global_envelope_test(curve_sets=list(gre.curves.ffam[[1]],  gre.curves.ffam[[2]], gre.curves.ffam[[3]], gre.curves.ffam[[4]], gre.curves.ffam[[5]]), type = "rank", alternative = "two.sided")
P_value.cmat <- global_envelope_test(curve_sets=list(gre.curves.cmat[[1]],  gre.curves.cmat[[2]], gre.curves.cmat[[3]], gre.curves.cmat[[4]], gre.curves.cmat[[5]]), type = "rank", alternative = "two.sided")
P_value.mfam <- global_envelope_test(curve_sets=list(gre.curves.mfam[[1]],  gre.curves.mfam[[2]], gre.curves.mfam[[3]], gre.curves.mfam[[4]], gre.curves.mfam[[5]]), type = "rank", alternative = "two.sided")


for(i in 1:5){
  pdf(paste("envelope",i,"_ffam.pdf",sep = ""))
  plot(global_envelope_test(gre.curves.ffam[[i]]),xlab = "", ylab = "",main = "",plot_style = "fv",xaxt="none", yaxt="none",legendmath=FALSE,lwd = 3,lty = 1)
  mtext(expression(gamma(theta)), side=2, line=2.2, cex=2)
  mtext(expression(theta), side=1, line=4, cex=2)
  axis(1,las=1, font=1,cex.axis=1.25)
  axis(2,,las=1, font=1,cex.axis=1.25)
  dev.off()
  
  pdf(paste("envelope",i,"_cmat.pdf",sep = ""))
  plot(global_envelope_test(gre.curves.cmat[[i]]),xlab = "", ylab = "",main = "",plot_style = "fv",xaxt="none", yaxt="none",legendmath=FALSE,lwd = 3,lty = 1)
  mtext(expression(gamma(theta)), side=2, line=2.2, cex=2)
  mtext(expression(theta), side=1, line=4, cex=2)
  axis(1,las=1, font=1,cex.axis=1.25)
  axis(2,,las=1, font=1,cex.axis=1.25)
  dev.off()
  
  pdf(paste("envelope",i,"_mfam.pdf",sep = ""))
  plot(global_envelope_test(gre.curves.mfam[[i]]),xlab = "", ylab = "",main = "",plot_style = "fv",xaxt="none", yaxt="none",legendmath=FALSE,lwd = 3,lty = 1)
  mtext(expression(gamma(theta)), side=2, line=2.2, cex=2)
  mtext(expression(theta), side=1, line=4, cex=2)
  axis(1,las=1, font=1,cex.axis=1.25)
  axis(2,,las=1, font=1,cex.axis=1.25)
  dev.off()
}
