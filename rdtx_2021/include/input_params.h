#pragma once
#ifndef INPARAMS_H_GUARD
#define INPARAMS_H_GUARD

#include "namespaces.h"
#include "global_constants.h"


/*
*********** RADAMPELTRAC *************

INPUT PARAMETERS

In this file, parameters for laser and electron bunch properties 
can be specified, in addition to the parameters of the 
'numerical spectrometer'.
 */
//*****************************************************************
//****** DUMP PARAMETERS ********
extern bool ifdumpnew ;// false; // switch to make new directory or overwrite

//*****************************************************************
//________________________________________________________________
// TIME PARAMETERS
extern int Nmax;// 10000; // Number of timesteps // Need 200000 for spectrum
extern long double taumax;//10.0; // maximum proper time (particles)
extern long double time_max ;// 200.0;// maximum lab time

//________________________________________________________________
// NORMALIZATION PARAMETERS 
extern long double omega_0;//2.36e15; //i.e. 800nm light
extern long double q_over_m;//-1.0; // -1.0 is electron, 1.0/931.5 proton etc.
extern long double rw0wp ;// 10.0; // ratio of w0 to wp

//________________________________________________________________
// RADIATION FORCE PARAMETERS 
extern long double radiation_force_on;//0.0; // 1 on, 0 off
//*****************************************************************



//*****************************************************************
//****** LASER PARAMETERS ********
// Should we use 2 pulses or not?

//using namespace twolaserpulses;

/*
  Monochromatic gaussian laser beam
 */
namespace laserpulse1
{
// relative frequency (to normalization)
  extern long double rww0 ;// 1.0; 

// group velocity of laser
  extern long double v_group ;// 1.0;//-1.5/(rw0wp*rw0wp*rww0*rww0); // +Etching rate - Lu PRSTAB
  //sqrt(1.0-1.0/(rw0wp*rw0wp*rww0*rww0));

// phase velocity of laser
  extern long double v_phase ;//1.0;///sqrt(1.0-1.0/(rw0wp*rw0wp*rww0*rww0));//1.0/v_group;

// field strength -> !!!!TIME AVERAGED!!!! (i.e. peak is sqrt(2)* for linear
  extern long double a0;//0.1;//3.5*0.0; 

//Propagation k fourvector - new - angle defined in zy plane from z axis
// IN radians
  extern long double theta_k0;//0.0; 

//________________________________________________________________

/* 
   DEFINITION OF POLARIZATION
   THese add phase and angle to the two waves, both 0 is linear pol,
   one of them ;//pi/2 -> circular polarization or elliptical
   The angles must be transverse (pi/2) for circular, parallel for linear
 */
  extern long double phase1;//0.0*PI*0.5,
  extern long double phase2;//0.0*PI*0.5;

//For phi, 0.0 is x, pi/2 is y
  extern long double phi_pol1 ;// 0.0*PI*0.5,
  extern long double phi_pol2 ;// 0.0*PI*0.5; 
//________________________________________________________________

/*
 Pulse length normalized to omega0
 (half of 1/e^2 width, FWHM ~ 1.1 t0)
*/
  extern long double t0;//65.0;//2.0*sqrt(5.0)*rw0wp; //65.0;

//Relative steepness of pulsefront
  extern long double pulse_front_steepness;//1.0;//3.0; 

/* Pulse waist (for intensity)

   w0;//0.0 means the waist is infinite (plane wave)

   Notes:
   remember, matched spot for bubble is 2sqrt(a)c/wp
*/
  extern long double w0;//20.0;//18.0;//5.0e10;//14.0*1000;//78.5;//2.0*sqrt(5.0)*rw0wp; 
// 
   extern long double delay;//0.5*laserpulse1::t0;//78.5;//2.0*sqrt(5.0)*rw0wp; 

//________________________________________________________________
// FOR FOCUSING PULSE
// Rayleigh range
  extern long double zR;//w0*w0*0.5; 
// position of focus
extern long double z_focus ;// 0.0;//1.5*laserpulse1::t0; 

}
//________________________________________________________________

//Switch on dumping of fields as a function of time
extern bool if_dump_fields;//true;
extern bool if_moving_frame ;// false;


// 4000
//*****************************************************************



//*****************************************************************
//****** ELECTROSTATIC WAVE PARAMETERS ********

//________________________________________________________________
// BUBBLE PROPERTIES
// Bubble traveling in +z direction

extern bool if_wakefield;//false;
extern bool if_bubble;//false;

//bubble radius
  extern long double r_bub;//2.0*sqrt(5.0)*rw0wp;//PI*rw0wp;//w0; 

// maximum potential, normalized to e/mc^2
  extern long double phi0_bub;//3.5*0.0;

// a number other than 1 makes it ellipsoidal in xz and xy planes
  extern long double bubble_aspect;//1.0; 

// phase velocity of bubble
  extern long double v_wake ;// laserpulse1::v_group;

// delay of bubble with respect to z;//0
  extern long double bub_delay;//25.0; 

/*
  Depletion modeling:
  These parameters model a bubble persisting for bub_width in t (not tau) 
  and then decreasing linearly in amplitude and increasing in width for 
  bub_ramp in t
 */
  extern long double bub_width;//4000.0, 
  extern long double bub_ramp;//1600.0; // ramp at end in lab time

//4000
//_______________________________________________________________
/*
  Modelling additional bubble features
  At the back and sides of the bubble is additional sheath formed by 
  electrons crossing. This results in a sheath potential which is 
  approximated by additional regions of phi with width and amplitude 
  specified below:
 */
// thickness of additional sheath.
  extern long double sheaththickness;//0.2*r_bub; 

// These next bits are slightly nonsensical

//additional phi at rear of bubble (due to electrons crossing)
  extern long double phi_bunch;//0.0;//0.0*2.0; 

//sheath thickness at rear of bubble
  extern long double ratio_sheath_to_phi_bunch;//1.0; 

//phi due to electron sheath at sides of bubble
  extern long double phi_sheath ;// 0.0*0.5;

//extra sheath thickness
  extern long double ratio_sheath_to_extrasheath;//2.0;  

//________________________________________________________________
//PARABOLIC CHANNEL PROPERTIES

extern bool if_channel;//false;
// maximum potential, normalized to e/mc^2
  extern long double phi0_chan;//0.0;//0.4*(r_bub0*r_bub0/rw0wp/rw0wp*0.5); 

// radius of channel
  extern long double r_channel;//sqrt(2.0);//
  extern long double r_bub0;//2.0*sqrt(3.5)*rw0wp;//sqrt(2.0); 

//*****************************************************************



//*****************************************************************
//****** ELECTRON BUNCH PARAMETERS ********
// Gaussian or tophat - check lower for namespace

// Number of particles in bunch
  extern int Npar ;// 1;

extern bool if_LWFA_dist;
// These are required for radiation emission calculations
// Use macroparticle approximation
  extern bool if_macroparticle ;// false;
// How many fundamental particles does a macro-particle represent  
extern long double sizeofbunch;//2.0;//1e9/(long double) Npar;
// what is the density of the particle bunch relative to the critical density?
  extern long double densityofbunch;//1.0;  

// width of bunch in x,y,z
  extern long double bunch_width[3] ;// {0.0,0.0,0.0}; 
//x,y,z of bunch peak/centre
  extern long double bunch_centre[3] ;// {0.0,0.0,2.0*laserpulse1::t0};

// Width in momentum space
extern long double bunch_momentum_spread[3] ;// {0.0,0.0,0.0}; 

// Peak/Center momentum
  extern long double bunch_momentum[3] ;// {0.0,0.0,-20.0}; 

// focal poisition of bunch
extern long double parzfocus;

// For large number of particles, may not want to dump data
  extern bool if_dump_par_dat ;// true; 

// Switch to be able to dump in lab time steps rather than particle
  extern bool dump_lab_time ;// true;

//*****************************************************************

//*****************************************************************
// IONIZATION
  extern bool ionization_on;//false;
 extern bool poincare_on;//false;

// Threshold PEAK a_0 for ionization
  extern long double ionization_threshold ;// 1.9;
//*****************************************************************


//*****************************************************************
//*********** UNIFORM MAGNETIC FIELD FOR DIAGNOSIS **********
// Uniform B field strength
extern bool if_bfield ;// false;

extern long double mag_Bfield ;// 0.0;
// B field direction unit vector, (x,y,z)
extern long double B_field_hat[3] ;// {0.0,0.0,1.0};
//*****************************************************************

//*****************************************************************
//***********  MAGNETIC WIGGLER FOR DIAGNOSIS **********
/*
  Wiggler period defined in z direction, but with infinite extent in
  x and y directions. Can either be sinusoidal ("not square") or
  square pulse shape. wiggle is in x-z plane - i.e. B=By. 
 */
extern bool if_wiggler;// = true;
extern bool wiggler_square;// = false;
// magnetic field strength
extern long double wiggler_Bfield;// = 800.0;
// wiggler k vector: wiggler period = 2pi/k
extern long double wiggler_k;// = 10.0;
// wiggler start and end positions 
extern long double wiggler_start;// = 0.01*PI;
extern long double wiggler_end;// = 0.21*PI;
extern long double wiggler_gap;// = 0.21*PI;


//*****************************************************************
//*********** NUMERICAL SPECTROMETER PARAMETERS **********
/*
  For calculating spectral power as a function of frequency for the 
  collection of particles
*/
// Switch spectrometer on/off (1/0)
  extern bool if_calc_spectrum;//false; 

// if particle radiation should be added coherently (1) or not (0)
  extern bool coherent_addition;//false; 

/*
  Spectral maximum/minimum frequency limits here
  Note these spectral limits are normalized to omega0
 */
  extern long double omega_max ;// 100000.0; //10.01;
  extern long double omega_min ;// 1.0; //0.01;
  extern int N_omega_bins ;// 500;

// Space omega bins linearly (0) or logarithmically (1)
  extern int if_omega_log ;// 1; 

/*
  Scattering angles of observation to be calculated
  Angles are defined in DEGREES, and theta is measured from z axis 
  rotating in z-y plane. Azimuth angle at x;//0
*/
extern long double theta_scatter_min;//180.0;//-0.573*20.0;
extern long double theta_scatter_max;//180.0;//+0.573*20.0;
  extern int N_theta_bins;//1;
 
// Dump orthoganal polarizations separately or combined
  extern bool separate_polarizations;//false; 

//________________________________________________________________
// This section is to do with the numerics of the Fourier transform

/*
 When switched on (1), the sinc functions arising due to the endpoints
 of the integral are subtracted, assuming that the velocity is 
    from -infinity to 0 and from taumax to infinity
 */
extern bool if_endpoint_integral;//false; 
extern int order_of_spectrometer; // determines method of calculation

/*
  For calculations where the endpoints are causing problems, a gaussian 
  shaped ramp can be applied to the endpoints to produce a window which
  damps artifacts from the sharp cutoff.
 */

//What damping do we want? 0 none, 1 end, 2 start and end
  extern bool if_endpoint_damping;//false; 

//Width of damping window ramp
  extern long double damping_width;//0.01*taumax;

// The ratio of the first ramp width to the second
  extern long double damping_endpoint_ratio;//1.0;

// A maximum time for calculation (0 means don't use)
  extern long double timemax;//0.0;
//*****************************************************************

//*****************************************************************
//****** LASER PARAMETERS ********
/*
  Monochromatic gaussian laser beam - LASER PULSE 2
*/

namespace laserpulse2 {
// relative frequency (to normalization)
  extern long double rww0 ;// 1.0; 

// group velocity of laser
  extern long double v_group ;// 1.0;//-1.5/(rw0wp*rw0wp*rww0*rww0); // +Etching rate - Lu PRSTAB
  //sqrt(1.0-1.0/(rw0wp*rw0wp*rww0*rww0));

// phase velocity of laser
  extern long double v_phase ;//1.0;///sqrt(1.0-1.0/(rw0wp*rw0wp*rww0*rww0));//1.0/v_group;

// field strength -> !!!!TIME AVERAGED!!!! (i.e. peak is sqrt(2)* for linear
  extern long double a0;//10.0;//3.5*0.0; 

//Propagation k fourvector - new - angle defined in zy plane from z axis
// IN radians
  extern long double theta_k0;//PI;
 
//________________________________________________________________

/* 
   DEFINITION OF POLARIZATION
   THese add phase and angle to the two waves, both 0 is linear pol,
   one of them ;//pi/2 -> circular polarization or elliptical
   The angles must be transverse (pi/2) for circular, parallel for linear
 */
  extern long double phase1;//0.0,
  extern long double phase2;//0.0*PI*0.5;

//For phi, 0.0 is x, pi/2 is y
  extern long double phi_pol1 ;//0.0*PI*0.5,
  extern long double phi_pol2 ;// 0.0*PI*0.5; 
//________________________________________________________________

/*
 Pulse length normalized to omega0
 (half of 1/e^2 width, FWHM ~ 1.1 t0)
*/
  extern long double t0;//65.0;//2.0*sqrt(5.0)*rw0wp; //65.0;

//Relative steepness of pulsefront
  extern long double pulse_front_steepness;//1.0;//3.0; 

/* Pulse waist (for intensity)

   w0;//0.0 means the waist is infinite (plane wave)

   Notes:
   remember, matched spot for bubble is 2sqrt(a)c/wp
*/
  extern long double w0;//5.0;//2.0*sqrt(5.0)*rw0wp; 
// 
  extern long double delay;//1.2*laserpulse1::t0;

//________________________________________________________________
// FOR FOCUSING PULSE
// Rayleigh range

  extern long double zR;//w0*w0*0.5; 
  extern long double z_focus ;// 0.0; 

}

//*****************************************************************

// Models (New defined here in rdt v1.0 agrt oct 2010)

extern bool if_one_pulses;
extern bool if_two_pulses;
extern bool if_with_focusing;
extern bool if_square_pulse;
extern bool if_trunc_pulse;
extern bool if_PLAD;
extern bool if_QPLAD;
extern bool if_uniform_bunch;

/*

using namespace onelaserpulse;
// laser model: with_focusing or no_focusing
using namespace with_focusing;

// laser model - temporal profile
using namespace asymmetric_gaussian_pulse;
//using namespace square_pulse;

// Radiation model: PLAD or LAD
using namespace LAD;

// Particle model: uniform_bunch or gaussian_bunch
using namespace gaussian_bunch;

*/


#endif
