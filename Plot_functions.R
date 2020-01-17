Plot_earth  <- function(eye,top){
  gdata <- globe::earth$coords
  runlen <- globe::earth$runlen
  spos <- spatialpos(gdata[, 1], gdata[, 2])
  mpos <- orthogproj(eye = eye, top = top, spos)
  x <- mpos[, 1]
  y <- mpos[, 2]
  ok <- (mpos[, 3] < 0)
  runlen <- runlen[runlen != 0]
  breaks <- cumsum(runlen)
  ok[breaks] <- FALSE
  s <- seq(x)[ok]
  segments(x[s], y[s], x[s + 1], y[s + 1])
  a <- seq(0, 2 * pi, length = 360)
  lines(cos(a), sin(a), lty = 2)
}

#i = 140;eye= c(pi/2,sseeqq[i])
#i = 30;eye= c(pi/2,sseeqq[i])
#values = filter_dat1; my_coords = coord_test; w=www; 
#eye= c(pi/2,sseeqq[i]); top = c(0,0); sigma0 = 0.015; col.image=pal(100)
Plot_function <- function(values,my_coords,w = sphwin(),eye,top,eyeE,topE,sigma0 = 0.2,col.image = NULL,ylim_o = NULL, ribb = FALSE){
  lon <- my_coords[,1]
  lat <- my_coords[,2]
  phi <- (lon/180) * pi
  theta <- ((90 - lat)/180) * pi
  Xgrid <- sp2(cbind(theta = theta, phi = phi))
  if (missing(w)) w <- x$X$win
  if (!is.null(w)) {
    inside <- in.W(Xgrid, win = w)
    Xgrid$X <- Xgrid$X[inside, , drop = FALSE]
    Xgrid$win <- w
  }
  #values <- predict(x, newdata = Xgrid)
  df <- cbind(spherstat:::Convert.globe(Xgrid$X), data.frame(values = values))
  if (!spherstat:::is.globe.point(eye)) eye <- spherstat:::Convert.globe(eye)
  if (!spherstat:::is.globe.point(top)) top <- spherstat:::Convert.globe(top)
  eye3 <- ensure3d(eye)
  top3 <- ensure3d(top)
  spos <- spatialpos(df[, 1], df[, 2])
  mpos <- orthogproj(eye3, top3, spos)
  xx <- mpos[, 1]
  yy <- mpos[, 2]
  ok <- (mpos[, 3] < 0)
  D <- disc(1)
  W <- disc(1, mask = TRUE, eps = eps, dimyx = dimyx)
  if (!is.null(w)) {
    wow <- spherstat:::sphwin2owin(w, eye = eye, top = top)
    W <- intersect.owin(W, wow)
  }
  X <- ppp(xx[ok], yy[ok], marks = df[ok, 3], window = W, check = FALSE)
  #print(X);str(X);mode(X)
  sigma <- sigma0
  Y <- Smooth(X, sigma = sigma, dimyx = dimyx,kernel="epanechnikov")
  Smooth.ppp(X = X, sigma = sigma, dimyx = dimyx, kernel = "epanechnikov")
  #print(Y);str(Y);mode(Y)
  do.call(plot.im, resolve.defaults(list(x = Y, col = col.image, zlim = range(values), ribbon = ribb), 
                                      list(main = main, box = FALSE), .MatchNull = FALSE,
                                      .StripNull = TRUE))
  if (!spherstat:::is.globe.point(eyeE)) eyeE <- spherstat:::Convert.globe(eyeE)
  if (!spherstat:::is.globe.point(topE)) topE <- spherstat:::Convert.globe(topE)
  Plot_earth(eyeE,topE)
}

plot.map<- function(database,center,...){
  Obj <- map(database,...,plot=F)
  coord <- cbind(Obj[[1]],Obj[[2]])
  
  # split up the coordinates
  id <- rle(!is.na(coord[,1]))
  id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=T)
  polygons <- apply(id,1,function(i){coord[i[1]:i[2],]})
  
  # split up polygons that differ too much
  polygons <- lapply(polygons,function(x){
    x[,1] <- x[,1] + center
    x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
    if(sum(diff(x[,1])>300,na.rm=T) >0){
      id <- x[,1] < 0
      x <- rbind(x[id,],c(NA,NA),x[!id,])
    }
    x
  })
  # reconstruct the object
  polygons <- do.call(rbind,polygons)
  Obj[[1]] <- polygons[,1]
  Obj[[2]] <- polygons[,2]
  
  map(Obj,...)
}

