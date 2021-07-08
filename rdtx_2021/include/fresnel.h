#pragma once
#ifndef FRESNEL_H_GUARD
#define FRESNEL_H_GUARD

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <complex>
#include "global_constants.h"

const int MAXIT = 100;
const long double EPS =0.01;//numeric_limits<long double>::epsilon();
const long double FPMIN = 1e-8;//numeric_limits<long double>::min();
const long double sqrt_FPMIN=sqrt(FPMIN);
const long double BIG = std::numeric_limits<long double>::max()*EPS;
const long double PIBY2=(PI*0.5), PIBY6=(PI/6.0), XMIN=1.5;

const long double fn_Mielenz[12] = {
  0.318309844,
  9.34626e-08,
  9.34626E-08,
  0.000606222,
  0.325539361,
  0.325206461,
  -7.450551455,
  32.20380908,
  -78.8035274,
  118.5343352,
  -102.4339798,
  39.06207702};
const long double gn_Mielenz[12] = {
  0.0,
  0.101321519,
  -4.07292e-05,
  -0.152068115,
  -0.046292605,
  1.622793598,
  -5.199186089,
  7.477942354,
  -0.695291507,
  -15.10996796,
  22.28401942,
  -10.89968491};
inline void fresnel2(const long double x, std::complex <long double> &cs)
/*
  Computes the Fresnel integrals S(x) and C(x) for all real x. C(x) is
  returned as the real part of cs and S(x) is the imaginary part
  (Adapted from NIST publication J. Res. Natl. Inst. Stand. Technol. 105, 589 (2000)
  by Mielenz) IF |x|>3
*/
{
  long double invx=1.0/x,xsq=x*x;
  long double xmult;
  long double f_x=0.0,g_x=0.0; 
  for (int i=0;i<12;++i) {
    xmult=invx;
    for (int j=0;j<2*i;++j) xmult*=invx;
    f_x+=fn_Mielenz[i]*xmult;
    g_x+=gn_Mielenz[i]*xmult;
  }
  cs = std::complex<long double>(0.5+f_x*sin(PIBY2*xsq)-g_x*cos(PIBY2*xsq),
		       0.5-f_x*cos(PIBY2*xsq)-g_x*sin(PIBY2*xsq));
  return;
}

inline void fresnel(const long double x, std::complex <long double> &cs)
/*
  Computes the Fresnel integrals S(x) and C(x) for all real x. C(x) is
  returned as the real part of cs and S(x) is the imaginary part
  (Taken from Numerical Recipes in C++)
*/
{
  bool odd;
  int k,n;
  long double a,ax,fact,pix2,sign,sum,sumc,sums,term,test;
  std::complex<long double> b,cc,d,h,del;

  ax=sqrt(x*x);
  if (ax<sqrt_FPMIN) cs = ax;
  else if (ax <= XMIN) {
    sum=sums=0.0;
    sumc=ax;
    sign=1.0;
    fact=PIBY2*ax*ax;
    odd=true;
    term=ax;
    n=3;
    for (k=1; k<=MAXIT; ++k) {
      term *= fact/k;
      sum += sign*term/n;
      test=sqrt(sum*sum)*EPS;
      if (odd) {
	sign = -sign;
	sums=sum;
	sum=sumc;
      } else {
	sumc = sum;
	sum = sums;
      }
      if (term < test) break;
      odd = !odd;
      n += 2;
    }
    if (k > MAXIT) { std::cout << "\nSeries failed in Fresnel!\n"; exit(0);}
    cs = std::complex<long double>(sumc,sums);
  } else {
    pix2 = PI*ax*ax;
    b=std::complex<long double>(1.0,-pix2);
    cc=BIG;
    d=h=1.0L/b;
    n = -1;
    for (k=2; k<=MAXIT; ++k) {
      n += 2;
      a = -n*(n+1);
      b += 4.0;
      d = 1.0L/(a*d+b);
      cc = b+a/cc;
      del = cc*d;
      h *= del;
      if (fabs(real(del)-1.0) + fabs(imag(del)) <= EPS) break;
    }
    if (k > MAXIT) { std::cout << "\nContinued fraction failed in Fresnel!\n";
      std::cout<<"x = "<< x<<'\n';
      exit(0);}
    h *= std::complex<long double>(ax,-ax);
    cs = std::complex<long double>(0.5,0.5)
      *(1.0L-std::complex<long double>(cos(0.5*pix2),sin(0.5*pix2))*h);
  }
  if (x < 0.0) cs = -cs;
  // cs=complex<long double>(0.5+1.0/PI/x*sin(PIBY2*x*x),0.5-1.0/PI/x*cos(PIBY2*x*x));
  return;
}
inline void fresnel3(const long double x, std::complex <long double> &cs)
/*
  Computes the Fresnel integrals S(x) and C(x) for all real x. C(x) is
  returned as the real part of cs and S(x) is the imaginary part
  (Super fast but rough version calculated from lowest order expansion)
*/
{
  if (fabs(x)<1.0)
    cs=std::complex<long double>(x,PIBY6*x*x*x);
  else 
    cs=std::complex<long double>(0.5+1.0/(PI*x)*sin(PIBY2*x*x),0.5-1.0/(PI*x)*cos(PIBY2*x*x));
  return;
}

inline std::complex<long double> 
Spectral_Int_Exact( long double chi_1,  long double chi_2,  long double sign,  long double v_0, 
		    long double v_1, std::complex<long double> DFres,
		    long double sin_phi_p,  long double sin_phi_m, 
		    long double cos_phi_p,  long double cos_phi_m,
		    long double sC2,  long double sC3,  long double sC4,
		    long double sin_sC2,  long double cos_sC2)
{
  /*
    This function calculates the combination of fresnel functions that 
    yields exact spectral integral to quadratic function of tau
   */

   long double sC1=(2.0*chi_2*v_0-chi_1*v_1);
  
   long double Psi_p = sC1*cos_sC2*sC3;
   long double Psi_m = sC1*sin_sC2*sC3;
    
   long double Fresnels_real = (Psi_p*DFres.real()+Psi_m*DFres.imag());
   long double Fresnels_imag = (Psi_p*DFres.imag()-Psi_m*DFres.real());

   long double realpart=Fresnels_real+2.0*v_1*(sin_phi_p-sin_phi_m);
   long double imagpart=Fresnels_imag-2.0*v_1*(cos_phi_p-cos_phi_m);
   
  // Now combine all these terms
  return std::complex< long double>(realpart*sC4,imagpart*sign*sC4);
}

#endif
