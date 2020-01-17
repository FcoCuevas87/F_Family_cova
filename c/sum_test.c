#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <Rinternals.h>


#define LOW -1.0e15
#define MAXERR 1e-6
////////////////

void ele4(int *N, double *sequence, double *lon,double *lat, double *resultado) {
  int i,j,N0,id=0;
  double aux0, *ddd; 
  
  N0 = N[0]*(N[0]-1);
  
  ddd = (int *)malloc(sizeof(double)*N0);
  
  for(i=0;i<N0;i++){
    ddd[i]=0;
  }
  
  for(i=0;i<(*N-1);i++){
    for(j=0;(i+1)<(*N);j++){
      aux0 = pow(sin(0.5*(lat[i]-lat[j])),2) + cos(lat[i])*cos(lat[j])*pow(sin(0.5*(lon[i]-lon[j])),2);
      ddd[id] = 2*asin(sqrt(aux0));
      id=id+1;
    }      
  }

  for(i = 0;i < id;i++){
    for(j=0;j<1000;j++){ 
      resultado[i] = resultado[i] + sequence[j]*cos(j*ddd[i]); 
    }
  }

  
  
}