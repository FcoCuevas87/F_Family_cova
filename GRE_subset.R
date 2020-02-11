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

table(Indexes1,Indexes2)
nsim.env <- 2499
max.dist <- pi/4
numbins <- 21

require("GeoModels")
dat0 <- dat1
vario <- GeoVariogram(data=dat0 ,coordx=coord0, maxdist=max.dist,distance="Geod" ,radius =1, numbins = numbins)
plot(vario$center,vario$variograms)
Variog.iso.latlon <- list()
Variogs.fit <- h.fit <- matrix(0, ncol = 4, nrow = numbins-1)
for(i in 1:2){
  for(j in 1:2){
    flag <- (Indexes1.1==(i-1))&(Indexes2.1==(j-1))
    sd4[j+(i-1)*2] <- sd(dat0[flag])
    RF <- list(X = dat0[flag], longitude = coord0[flag,1], latitude = coord0[flag,2])
    Variog.fil.est$latlon[[j+(i-1)*2]] <- GeoVariogram(data=dat0[flag] ,coordx=coord0[flag,], maxdist=max.dist,distance="Geod", radius =1, numbins = numbins)
    Variogs.fit[,j+(i-1)*2] <- Variog.fil.est$latlon[[j+(i-1)*2]]$variograms
    h.fit[,j+(i-1)*2] <- Variog.fil.est$latlon[[j+(i-1)*2]]$center
    Variog.iso.latlon$ffam[[j+(i-1)*2]] <- Variog.iso.latlon$cmat[[j+(i-1)*2]] <- Variog.iso.latlon$mfam[[j+(i-1)*2]] <- matrix(0,ncol=numbins-1,nrow = nsim.env)
  }
}


matplot(h.fit ,Variogs.fit,type="l")
set.seed(123)

for(k in 1:nsim.env){
  #Simulating GRF
  dat1.tmp <- chol01%*%c(rnorm(nrow(coord0)))
  dat2.tmp <- chol02%*%c(rnorm(nrow(coord0)))
  dat3.tmp <- chol03%*%c(rnorm(nrow(coord0)))
  
  #Creating lists
  
  #Computing isotropic variogram fixed

  for(i in 1:2){
    for(j in 1:2){
      flag <- (Indexes1==(i-1))&(Indexes2==(j-1))
      Variog.iso.latlon$ffam[[j+(i-1)*2]][k,] <- GeoVariogram(data=dat1.tmp[flag] ,coordx=coord0[flag,], maxdist=max.dist,distance="Geod" ,radius =1,numbins = numbins)$variograms
      Variog.iso.latlon$cmat[[j+(i-1)*2]][k,] <- GeoVariogram(data=dat2.tmp[flag] ,coordx=coord0[flag,], maxdist=max.dist,distance="Geod" ,radius =1,numbins = numbins)$variograms
      Variog.iso.latlon$mfam[[j+(i-1)*2]][k,] <- GeoVariogram(data=dat3.tmp[flag] ,coordx=coord0[flag,], maxdist=max.dist,distance="Geod" ,radius =1,numbins = numbins)$variograms
    }
  }
  
  cat(k,"\n")
}

#devtools::install_github('myllym/GET')
require("GET")


gre.curves.ffam <- gre.test.mfam <- gre.test.cmat <- list()
for(i in 1:4){
  gre.curves.ffam[[i]] <- create_curve_set( list(r = h.fit[,i], obs = Variogs.fit[,i] , sim_m = t(Variog.iso.latlon$ffam[[i]])) )
  gre.curves.mfam[[i]] <- create_curve_set( list(r = h.fit[,i], obs = Variogs.fit[,i] , sim_m = t(Variog.iso.latlon$mfam[[i]])) )
  gre.curves.cmat[[i]] <- create_curve_set( list(r = h.fit[,i], obs = Variogs.fit[,i] , sim_m = t(Variog.iso.latlon$cmat[[i]])) )

  gre.test.ffam <- global_envelope_test(gre.curves.ffam[[i]], type = "rank", alternative = "two.sided")
  gre.test.mfam <- global_envelope_test(gre.curves.mfam[[i]], type = "rank", alternative = "two.sided")
  gre.test.cmat <- global_envelope_test(gre.curves.cmat[[i]], type = "rank", alternative = "two.sided")
  
  print(attr(gre.test.ffam,"p_interval"))
  print(attr(gre.test.mfam,"p_interval"))
  print(attr(gre.test.cmat,"p_interval"))
  
  y.min = min(cbind(gre.test.ffam$lo, gre.test.ffam$lo, gre.test.ffam$lo))
  y.max = max(cbind(gre.test.ffam$hi, gre.test.ffam$hi, gre.test.ffam$hi))
  ylim = c(y.min,(1.1*y.max))
  
  pdf(paste("envelope",i,"_ffam.pdf",sep = ""))
  plot(gre.test.ffam,xlab = "", ylab = "",main = "",plot_style = "fv",xaxt="none", yaxt="none",legendmath=FALSE,lwd = 3,lty = 1, ylim = ylim)
  mtext(expression(gamma(theta)), side=2, line=2.2, cex=2)
  mtext(expression(theta), side=1, line=4, cex=2)
  axis(1,las=1, font=1,cex.axis=1.25)
  axis(2,las=1, font=1,cex.axis=1.25)
  dev.off()
  
  pdf(paste("envelope",i,"_mfam.pdf",sep = ""))
  plot(gre.test.mfam,xlab = "", ylab = "",main = "",plot_style = "fv",xaxt="none", yaxt="none",legendmath=FALSE,lwd = 3,lty = 1, ylim = ylim)
  mtext(expression(gamma(theta)), side=2, line=2.2, cex=2)
  mtext(expression(theta), side=1, line=4, cex=2)
  axis(1,las=1, font=1,cex.axis=1.25)
  axis(2,,las=1, font=1,cex.axis=1.25)
  dev.off()
  
  pdf(paste("envelope",i,"_cmat.pdf",sep = ""))
  plot(gre.test.cmat,xlab = "", ylab = "",main = "",plot_style = "fv",xaxt="none", yaxt="none",legendmath=FALSE,lwd = 3,lty = 1, ylim = ylim)
  mtext(expression(gamma(theta)), side=2, line=2.2, cex=2)
  mtext(expression(theta), side=1, line=4, cex=2)
  axis(1,las=1, font=1,cex.axis=1.25)
  axis(2,,las=1, font=1,cex.axis=1.25)
  dev.off()
}

#P_value.ffam <- global_envelope_test(curve_sets=list(gre.curves.ffam[[1]],  gre.curves.ffam[[2]], gre.curves.ffam[[3]], gre.curves.ffam[[4]]), type = "rank", alternative = "two.sided")
#P_value.cmat <- global_envelope_test(curve_sets=list(gre.curves.cmat[[1]],  gre.curves.cmat[[2]], gre.curves.cmat[[3]], gre.curves.cmat[[4]]), type = "rank", alternative = "two.sided")
#P_value.mfam <- global_envelope_test(curve_sets=list(gre.curves.mfam[[1]],  gre.curves.mfam[[2]], gre.curves.mfam[[3]], gre.curves.mfam[[4]]), type = "rank", alternative = "two.sided")

#base::save.image("Global_rank_env2.Rdat")

