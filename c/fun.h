#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <Rinternals.h>

#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>
#include "Warper_2F1.h"

#define LOW -1.0e15
#define MAXERR 1e-6
////////////////



double modelo_de_covarianza(int numero,double lon1,double lat1,double lon2,double lat2, double* par) ;


////////////////


double modelo_de_covarianza(int numero, double lon1,double lat1,double lon2,double lat2, double* par)  
{
    double aux, qq, suma, suma0 , eucl_dis, nu, covariance, factor1,factor2,factor3, alfa;
    int j;
    double s1,s2,aux1,aux2, constante;

    double  gc_dist, R_Earth = 1;
    double aux0 = pow(sin(0.5*(lat1-lat2)),2) + cos(lat1)*cos(lat2)*pow(sin(0.5*(lon1-lon2)),2);
    gc_dist = 2*R_Earth*asin(sqrt(aux0));
    eucl_dis = 2*sin(0.5*gc_dist);
    
    s1 = cos(lat1);
    s2 = cos(lat2);    

    if(numero==0){ // ARBlib
      aux1 = par[0];
      aux2 = par[0];
      
      alfa = 1/pow(sqrt(par[1]),2);
      nu =   pow(sqrt(par[2]),2);
      
      qq = cos(gc_dist);
      constante = (gammafn(alfa+0.5+nu)*gammafn(alfa+nu)) / (gammafn(2*alfa+0.5+nu)*gammafn(nu));
      if(lat1==lat2 & lon1==lon2){ 
        covariance = aux1*aux2 ; 
      }else{
        gsl_set_error_handler_off();
        covariance = aux1*aux2 * constante * warper_2f1(alfa,alfa+0.5,2*alfa+0.5+nu,qq);
      }
    }
    if(numero==1){
		    aux1 = par[0];
        aux2 = par[0];
                    
        alfa = 1/pow(sqrt(par[1]),2);
        nu =   pow(sqrt(par[2]),2);

        qq = cos(gc_dist);
        constante = (gammafn(alfa+0.5+nu)*gammafn(alfa+nu)) / (gammafn(2*alfa+0.5+nu)*gammafn(nu));
          if(lat1==lat2 & lon1==lon2){ 
            covariance = aux1*aux2 ; 
          }else{
				    gsl_set_error_handler_off();
            covariance = aux1*aux2 * constante * gsl_sf_hyperg_2F1(alfa,alfa+0.5,2*alfa+0.5+nu,qq);
          }
    }
    if(numero==2){
      //aux1  = pow(sqrt(par[0]),2)*sqrt(exp(par[1] + par[2]*s1)); 
      //aux2  = pow(sqrt(par[0]),2)*sqrt(exp(par[1] + par[2]*s2)); 		
      aux1 = par[0];
      aux2 = par[0];
      
		  alfa = pow(sqrt(par[1]),2);
		  nu = pow(sqrt(par[2]),2);
		       
      factor1  = pow(2.0, 1.0 - nu) / gammafn(nu);
      factor2  = pow( eucl_dis/alfa ,  nu );
      factor3  = bessel_k(eucl_dis/alfa ,nu,1);
           
      if(lat1==lat2 & lon1==lon2){ 
        covariance = aux1*aux2; 
        }else{
        covariance = aux1*aux2*factor1*factor2*factor3;
			 }
	}	
  if(numero==3){
    aux1 = par[0];
    aux2 = par[0];	      
      
    alfa = pow(sqrt(par[1]),2);
    nu = pow(sqrt(par[2]),2);	         
      
    suma0 = 0;
    suma = 0;
    for(j=0;j<1000;j++){ 
      aux = pow(pow(j,2) + pow(1/alfa,2),nu+0.5);
      suma0 = suma0 + 1/aux;
      suma = suma + cos(j*gc_dist)/aux; 
    }
    if(lat1==lat2 & lon1==lon2){
      covariance = aux1*aux2;
    }else{
      covariance = aux1*aux2*suma/suma0;
    }
  }
  if(numero==4){
    aux1 = par[0];
    aux2 = par[0];	      
    
    if(lat1==lat2 & lon1==lon2){
      covariance = aux1*aux2;
    }else{
      covariance = 0;
    }
  }
	return(covariance);
}






















