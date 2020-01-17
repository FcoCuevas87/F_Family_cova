loglik<-function(param){  

    lik2 <- 1.0e8
    mmm = matriz_cov(eleccion,N,param,lon,lat);

    cholmmm <- try(chol(mmm),silent=TRUE)
    if(inherits(cholmmm, "try-error")) {return(lik2)}

    detvarcov <- sum(log(diag(cholmmm)))

    lik2 <- 0.5*( N*log(2*pi) + 2*detvarcov+sum(datos*backsolve(cholmmm,
                  forwardsolve(cholmmm,datos,transpose=TRUE,upper.tri=TRUE))))

    return(lik2)
}

loglik2_formula <- function(param,datos, lon, lat, design_mat = cbind(rep(1,length(lon)),lon,lat),eleccion = 0, nugget = FALSE){  
  N = length(datos)
  #sigma_for = rep(1,N) ~ sinlat + coslon
  if(!is.null(design_mat)){
    Sigma0 <- c(design_mat%*%param[1:ncol(design_mat)])
    Sigma0 <- sqrt(exp(Sigma0))
    Sigma0 <- diag(Sigma0,N)
    param <- param[-c(1:(ncol(design_mat)))]
  }
  
  if(nugget){
    tau <- param[1]*param[1]
    mm0 <- tau*diag(1,N)
    param <- param[-1]
  }
  
  
  lik2 <- 1.0e8
  mmm <- matriz_cov(eleccion,N,c(1,param),lon,lat);
  
  if(!is.null(design_mat)){
    mmm <- Sigma0%*%mmm%*%Sigma0
  }else if(nugget){
    mmm <- mmm + mm0
  }
  
  cholmmm <- try(chol(mmm),silent=TRUE)
  if(inherits(cholmmm, "try-error")) {return(lik2)}
  detvarcov <- sum(log(diag(cholmmm)))
  lik2 <- 0.5*( N*log(2*pi) + 2*detvarcov+sum(datos*backsolve(cholmmm,
                                                              forwardsolve(cholmmm,datos,transpose=TRUE,upper.tri=TRUE))))
  return(lik2)
}


loglik_fix<-function(param,fix = c(0,0)){  
  
  param <- c(fix,param)
  
  lik2 <- 1.0e8
  mmm = matriz_cov(eleccion,N,param,lon,lat);
  
  cholmmm <- try(chol(mmm),silent=TRUE)
  if(inherits(cholmmm, "try-error")) {return(lik2)}
  
  detvarcov <- sum(log(diag(cholmmm)))
  
  lik2 <- 0.5*( N*log(2*pi) + 2*detvarcov+sum(datos*backsolve(cholmmm,
                                                              forwardsolve(cholmmm,datos,transpose=TRUE,upper.tri=TRUE))))
  
  return(lik2)
}



Prscores<-function(data,inv,nsites)   {

vv=diag(inv)                                                                              
###########################
dime <- nsites;
###########################
D=diag(1/vv,dime,dime)
DD=sqrt(D)
temp=inv%*%data
z=D%*%temp
zz=DD%*%temp


MSE=(1/dime)*(t(z)%*%z)
MAE= (1/dime)*sum(abs(z))
LSCORE=(1/(2*dime))*(sum(log(2*pi/vv))+sum(zz^2))
CRPS=(1/dime)*(sum((1/vv)^0.5*zz*(2*pnorm(zz)-1))+2*sum((1/vv)^0.5*pnorm(zz))+sum((1/vv)^0.5)/sqrt(pi))
###########################
scores <- list(  MSE = MSE, MAE=MAE, LSCORE = LSCORE, CRPS = CRPS )
return(scores)
}







