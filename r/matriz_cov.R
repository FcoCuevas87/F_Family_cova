matriz_cov<-function(eleccion,N,param,lon,lat)
{
  storage.mode(param) <- "double"
  storage.mode(lon) <- "double"
  storage.mode(lat) <- "double"
  num<- N^2;
  o <- .C("matriz_cov", eleccion=as.integer(eleccion),  N=as.integer(N),  
                            param=param ,lon=lon,lat=lat, matrix = double(num))
  mat<-o$matrix
   
  M=t(array(dim=c(N,N),mat));
  return(M)
}
