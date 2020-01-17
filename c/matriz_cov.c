#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include "fun.h"
#include <Rinternals.h>
#include <R_ext/Applic.h>

#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>

#define LOW -1.0e15
#define MAXERR 1e-6

void matriz_cov(int *eleccion ,int *N, double *parametros,double *lon,double *lat, double *matrix) {
  int i,j,id=0;
  for(i=0;i<(*N);i++){
    for(j=0;j<(*N);j++){
      matrix[id] = modelo_de_covarianza(*eleccion,lon[i],lat[i],lon[j],lat[j],parametros);
      id=id+1;
    }      
  }
}









