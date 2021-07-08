#pragma once
#ifndef GLOBALCONST_H_GUARD
#define GLOBALCONST_H_GUARD

// GLOBAL CONSTANTS
//*****************************************************************
//Constants
 const long double PI=acos(-1.0);
 const long double oneover6=1.0/6.0;
 const long double e=1.602176e-19;
 const long double me=9.10938e-31;
 const long double mu0=4.0*PI*1e-7;
 const long double c=299792458;
 const long double hbar=1.05457172647e-34;
 const long double spectrum_const= e*e*c/PI/PI/PI/16.0*mu0/e; // in units of eV

//Messages and output
     const int Nmessages=10;
   const std::string version = "1.0";
 

//------------------------------------------------------
   const long double delta_step = 1.0e-6;
 // this is for calculating potentials by differential of EM tensor
// (NB. Can't be too small, or truncation errors arise - 1.0e-6 ok)

   const double model_threshold=0.001; // Threshold between Taylor series expansion and exact solution

   const long double damping_n_sigmas = 3.0; // number of std devs for endpoint filter

namespace twolaserpulses {

}
namespace onelaserpulse {

}
namespace otherlaserpulse {

}
namespace nolaserpulse {

}

#endif
