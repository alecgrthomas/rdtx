#pragma once
#ifndef CONST_H_GUARD
#define CONST_H_GUARD

#include "fourtens.h"
#include "input_params.h"

//------------------------------------------------------


// This is tau_0 in Rohrlich paper
   extern long double dtau0 ;//= taumax/(long double) Nmax;
   extern long double inv_dtau0 ;//= 1.0/dtau0;
extern    long double dtau0sq;// = dtau0*dtau0;
extern    long double rad_const;//mu0*e*e*omega_0/6.0/PI/c/me*radiation_force_on; // 0.62x10^-23 s * omega0
extern long double qed_const;//= hbar*omega_0/me/c/c; // ratio of compton wavelength to laser wavelength
//extern    long double spectrum_const;//= e*e*c/PI/PI/PI/16.0*mu0;
extern    long double shape_factor_constant;// = pow(omega_0*e*e*mu0/me/c,1.0L/3.0L);
extern    long double dtheta;//=(theta_scatter_max-theta_scatter_min)/(long double)(N_theta_bins-(1-!(N_theta_bins-1)))/180.0*PI;
extern    long double r_sheath;//=(r_bub+sheaththickness);

extern    long double inv_r_sheath_thick_square;//=1.0/(sheaththickness*sheaththickness);
extern    long double ratio_bub_sheath;// = r_bub/sheaththickness;
extern    long double inv_bub_sheath;// = sheaththickness/r_bub;
extern    long double one_plus_rbs;// = (1.0+ratio_bub_sheath);
extern    long double one_plus_rbs_2 ;//=(1.0+ratio_bub_sheath)*2.0;

extern long double r_sheath_half;//=r_sheath-0.5*sheaththickness;
extern long double r_bubsquared;//=r_bub*r_bub;
extern long double inv_r_bubsquared;//=1.0/r_bubsquared;
extern long double inv_r_bub ;//= 1.0/r_bub;
extern int spectrumerror;

extern    long double r_channelsquared;//=r_channel*r_channel;
extern    long double inv_bubble_aspect;//=1.0/(bubble_aspect);

extern    long double w0sq ;//= laserpulse1::w0*laserpulse1::w0;
extern    long double w0sq2 ;//= laserpulse2::w0*laserpulse2::w0;
  
extern long double invzRsq ;//= 1.0/laserpulse1::zR/laserpulse1::zR;
extern long double invzRsq2 ;//= 1.0/laserpulse2::zR/laserpulse2::zR;

extern long double inv_t0squared ;//= 1.0/laserpulse1::t0/laserpulse1::t0;
extern long double inv_t0squared2;// = 1.0/laserpulse2::t0/laserpulse2::t0;

extern long double a0_pulse1;// = laserpulse1::a0*sqrt((pow(cos(laserpulse1::phi_pol1)+cos(laserpulse1::phi_pol2),2.0L)
						    //+pow(sin(laserpulse1::phi_pol1)+sin(laserpulse1::phi_pol2),2.0L))*0.5);
extern long double a0_pulse2;// = laserpulse2::a0*sqrt((pow(cos(laserpulse2::phi_pol1)+cos(laserpulse2::phi_pol2),2.0L)
						    //+pow(sin(laserpulse2::phi_pol1)+sin(laserpulse2::phi_pol2),2.0L))*0.5);
extern int Ndumps;//=1000;
extern double yslice;
   extern int Nx_f;//=200,
extern int Nz_f; //Nz_f=400; //for outputting fields - no. points
extern long double xsize;//=100.0, 
extern long double zsize;//=150.0; // and spatial size
#endif
