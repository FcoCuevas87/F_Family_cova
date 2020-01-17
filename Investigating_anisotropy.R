dat0 <- datos
coord0 <- coord
RF <- list(X = dat0, longitude = coord0[,1], latitude = coord0[,2])
Variog.raw.new <- Est.variog(RF,h.max = pi/3, l = 30)

Indexes1 <- findInterval(coord0[,2], c(35))  #Lat cut
Indexes2 <- findInterval(coord0[,1], c(-1.25)) #Long cut
table(Indexes1)
table(Indexes2)


Variog.fil.est <- list()
Variog.fil.est$lat <- list()
sd2 <- c(0,0)
for(i in 0:1){
  sd2[i+1] <- sd(dat0[Indexes1==i])
  RF <- list(X = (dat0[Indexes1==i]), longitude = coord0[Indexes1==i,1], latitude = coord0[Indexes1==i,2])
  Variog.fil.est$lat[[i+1]] <- Est.variog(RF,h.max = pi/3, l = 20)
}

Variog.fil.est$lon <- list()
sd3 <- c(0,0)
for(i in 0:1){
  sd3[i+1] <- sd(dat0[Indexes2==i])
  RF <- list(X = (dat0[Indexes2==i]), longitude = coord0[Indexes2==i,1], latitude = coord0[Indexes2==i,2])
  Variog.fil.est$lon[[i+1]] <- Est.variog(RF,h.max = pi/3, l = 20)
}

Indexes1.1 <- findInterval(coord0[,2], c(35.5))  #Lat cut
Indexes2.1 <- findInterval(coord0[,1], c(-1.25)) #Long cut
table(Indexes1.1,Indexes2.1)

pdf("./images/iso_variogram.pdf")
plot(x = c(0,Variog.raw.new$hvec), y = c(0,Variog.raw.new$variog),type="l",xlab="",ylab="",ylim=c(0,1.5), lwd = 3,xaxt="none", yaxt="none")
mtext(expression(gamma(theta)), side=2, line=2.2, cex=2)
mtext(expression(theta), side=1, line=4, cex=2)
axis(1, seq(0,1, l=5),las=1, font=1,cex.axis=1.25)
axis(2, seq(0,1.5,0.5),las=1, font=1,cex.axis=1.25)
dev.off()

pdf("./images/variogram_latitude.pdf")
plot(x = c(0,Variog.fil.est$lat[[1]]$hvec), y = c(0,Variog.fil.est$lat[[1]]$variog),type="l",xlab="",ylab="",lty = 3,xlim=range(unlist(lapply(Variog.fil.est$lat,function(x) x$hvec))),ylim=c(0,1.5), lwd = 3, xaxt="none", yaxt="none")
mtext(expression(gamma(theta)), side=2, line=2.2, cex=2)
mtext(expression(theta), side=1, line=4, cex=2)
mtext(expression(theta), side=1, line=4, cex=2)
axis(1, seq(0,1, l=5),las=1, font=1,cex.axis=1.25)
axis(2, seq(0,1.5,0.5),las=1, font=1,cex.axis=1.25)
lines(x = c(0,Variog.fil.est$lat[[2]]$hvec), y = c(0,Variog.fil.est$lat[[2]]$variog),lty = 3, col="red", lwd = 3)
lines(c(0,Variog.raw.new$hvec),c(0,Variog.raw.new$variog), lwd = 3)
dev.off()

pdf("./images/variogram_longitude.pdf")
plot(x = c(0,Variog.fil.est$lon[[1]]$hvec), y = c(0,Variog.fil.est$lon[[1]]$variog),type="l",xlab="",ylab="",lty = 3,xlim=range(unlist(lapply(Variog.fil.est$lon,function(x) x$hvec))),ylim=c(0,1.5), lwd = 3, xaxt="none", yaxt="none")
mtext(expression(gamma(theta)), side=2, line=2.2, cex=2)
mtext(expression(theta), side=1, line=4, cex=2)
axis(1, seq(0,1, l=5),las=1, font=1,cex.axis=1.25)
axis(2, seq(0,1.5,0.5),las=1, font=1,cex.axis=1.25)
lines(x = c(0,Variog.fil.est$lon[[2]]$hvec), y = c(0,Variog.fil.est$lon[[2]]$variog),lty = 3, col="red", lwd = 3)
lines(c(0,Variog.raw.new$hvec),c(0,Variog.raw.new$variog), lwd = 3)
dev.off()
