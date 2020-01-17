mat_pred<-function(eleccion,N,K,param,lon,lat,lon0,lat0){
  storage.mode(param) <- "double"
  storage.mode(lon) <- "double"
  storage.mode(lat) <- "double"
  storage.mode(lon0) <- "double"
  storage.mode(lat0) <- "double"

  num<- N*K;
  o <- .C("mat_pred", eleccion=as.integer(eleccion),  N=as.integer(N), K=as.integer(K),  
                              param=param ,lon=lon, lat=lat, lon0=lon0, lat0=lat0,  
                              matrix = double(num))
  mat<-o$matrix
  M=t(array(dim=c(N,K),mat))
  return(M)
}

