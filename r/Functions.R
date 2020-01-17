Est.variog <- function(RF, h.max = (pi/8), l = 30){ #coords must be in radians
  h.seq  <- seq(0,h.max,l=l)
  hvec <- h.seq[-length(h.seq)] + 0.5*diff(h.seq)
  
  sph.ppobj <- spatstat.sphere::s2pp(list(longitude = RF$longitude, latitude = RF$latitude) , region = spatstat.sphere::s2(radius = 1))
  DIJ <- spatstat.sphere:::closepairs.s2pp(sph.ppobj,rmax = h.max)
  
  ddd2  <- (RF$X[DIJ$i] - RF$X[DIJ$j])^2
  
  ID <- base::findInterval(DIJ$d, h.seq, rightmost.closed = TRUE)
  
  N0 <- tabulate(ID)
  Variog <- c(unlist(lapply(split(ddd2,ID),sum))/N0)
  Variog <- Variog/2
  result <- list(variog = Variog, hvec = hvec, N0 = N0)
  return(result)
}
