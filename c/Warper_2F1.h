/* F-Family of covariance functions using GSL or ARB */
/* Warper from ARB library to R was manually created*/
/* Author: Francisco Cuevas */
#include <R.h>
#include <Rmath.h>
#include <arb.h>
#include <arb_mat.h>
#include <acb_hypgeom.h>
#include <flint/profiler.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>
  
double warper_2f1(double alpha,double tau,double nu,double zz){
  slong prec = 4096;
  acb_t a, b, c, z, w;
  double re;
  /* initializing ARB values */
  acb_init(a); 
  acb_init(b);
  acb_init(c);
  acb_init(z);
  acb_init(w);
    
  /* Setting parameters of the model */
  acb_set_d_d(a,alpha,0);
  acb_set_d_d(b,tau,0);
  acb_set_d_d(c,nu,0);
  acb_set_d_d(z,zz,0); 
    
  /* Compute 2F1 function */
  acb_hypgeom_2f1(w, a, b, c, z, 0, prec);
    
  /* Conversion from ARB to C*/
  re = arf_get_d(arb_midref(acb_realref(w)), ARF_RND_NEAR);
  /* Result */
  return(re);
      
  /* Clean Values */
  acb_clear(a); 
  acb_clear(b);
  acb_clear(c);
  acb_clear(z);
  acb_clear(w);
  flint_cleanup();
}

/* Warper function */
double hypegeometric2f1_function_arb(double a,double b,double c,double x){    
  double result;
  result=warper_2f1(a,b,c,x);
  return(result);
}