#ifndef RUNPARAMS_C_
#define RUNPARAMS_C_

#include "fourtens.h"
#include "headers.h"





int set_all_parameters(void);

//*****************************************************************
// ABSOLUTE DEFAULTS

//****** DUMP PARAMETERS ********
 bool ifdumpnew = false; // switch to make new directory or overwrite


//*****************************************************************
//________________________________________________________________
// TIME PARAMETERS
 int Nmax = 3200; // Number of timesteps // Need 200000 for spectrum
 long double taumax=1.0; // maximum proper time (particles)
 long double time_max = 100.0;// maximum lab time

//________________________________________________________________
// NORMALIZATION PARAMETERS 
 long double omega_0=2.36e15; //i.e. 800nm light
 long double q_over_m=-1.0; // -1.0 is electron, 1.0/931.5 proton etc.
 long double rw0wp = 10.0; // ratio of w0 to wp

//________________________________________________________________
// RADIATION FORCE PARAMETERS 
 long double radiation_force_on=1.0; // 1 on, 0 off
//*****************************************************************



/*
	Standing wave for test cases
*/
bool if_standing_wave = false;
bool if_standing_circ = false;
long double standing_wave_k = 1.0;
long double standing_wave_a = 1.0;
long double standing_wave_phase = 0.0;
long double standing_wave_pol = 0.0; // polarization angle with respect to x axis
long double sinthetastanding = sin(standing_wave_pol);
long double costhetastanding = cos(standing_wave_pol);

	
/*
  Monochromatic gaussian laser beam
 */
namespace laserpulse1
{
// relative frequency (to normalization)
   long double rww0 = 1.0; 

// group velocity of laser
   long double v_group = 1.0;//-1.5/(rw0wp*rw0wp*rww0*rww0); // +Etching rate - Lu PRSTAB
  //sqrt(1.0-1.0/(rw0wp*rw0wp*rww0*rww0));

// phase velocity of laser
   long double v_phase =1.0;///sqrt(1.0-1.0/(rw0wp*rw0wp*rww0*rww0));//1.0/v_group;

// field strength -> !!!!TIME AVERAGED!!!! (i.e. peak is sqrt(2)* for linear
   long double a0=1.0;//3.5*0.0; 

//Propagation k fourvector - new - angle defined in zy plane from z axis
// IN radians
   long double theta_k0=0.0; 

//________________________________________________________________

/* 
   DEFINITION OF POLARIZATION
   THese add phase and angle to the two waves, both 0 is linear pol,
   one of them =pi/2 -> circular polarization or elliptical
   The angles must be transverse (pi/2) for circular, parallel for linear
 */
   long double phase1=0.0*PI*0.5,phase2=0.0*PI*0.5;

//For phi, 0.0 is x, pi/2 is y
   long double phi_pol1 = 0.0*PI*0.5,phi_pol2 = 0.0*PI*0.5; 
//________________________________________________________________

/*
 Pulse length normalized to omega0
 (half of 1/e^2 width, FWHM ~ 1.1 t0)
*/
   long double t0=50.0;//2.0*sqrt(5.0)*rw0wp; //65.0;

//Relative steepness of pulsefront
   long double pulse_front_steepness=1.0;//3.0; 

/* Pulse waist (for intensity)

   w0=0.0 means the waist is infinite (plane wave)

   Notes:
   remember, matched spot for bubble is 2sqrt(a)c/wp
*/
   long double w0=20.0;//18.0;//5.0e10;//14.0*1000;//78.5;//2.0*sqrt(5.0)*rw0wp; 
// 
    long double delay=1.0*laserpulse1::t0;//78.5;//2.0*sqrt(5.0)*rw0wp; 

//________________________________________________________________
// FOR FOCUSING PULSE
// Rayleigh range
   long double zR=w0*w0*0.5*rww0; 
// position of focus
 long double z_focus = 0.0;//1.5*laserpulse1::t0; 

}
//________________________________________________________________

//Switch on dumping of fields as a function of time
 bool if_dump_fields = false;
 bool if_moving_frame = false;
 int n_field_dumps = 20;

// 4000
//*****************************************************************



//*****************************************************************
//****** ELECTROSTATIC WAVE PARAMETERS ********

//________________________________________________________________
// BUBBLE PROPERTIES
// Bubble traveling in +z direction

 bool if_wakefield=false;
 bool if_bubble=false;
//bubble radius
   long double r_bub=2.0*sqrt(4.0)*rw0wp;//PI*rw0wp;//w0; 

// maximum potential, normalized to e/mc^2
   long double phi0_bub=2.0; // a0/2  

// a number other than 1 makes it ellipsoidal in xz and xy planes
   long double bubble_aspect=1.0; 

// phase velocity of bubble
   long double v_wake = 1.0-1.5/(rw0wp*rw0wp*laserpulse1::rww0*laserpulse1::rww0);//laserpulse1::v_group;

// delay of bubble with respect to z=0
   long double bub_delay=0.0; 

/*
  Depletion modeling:
  These parameters model a bubble persisting for bub_width in t (not tau) 
  and then decreasing linearly in amplitude and increasing in width for 
  bub_ramp in t
 */
   long double bub_width=1000.0, bub_ramp=100.0; // ramp at end in lab time

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
   long double sheaththickness=0.1*r_bub; 

// These next bits are slightly nonsensical

//additional phi at rear of bubble (due to electrons crossing)
   long double phi_bunch=1.0;//0.0*2.0; 

//sheath thickness at rear of bubble
   long double ratio_sheath_to_phi_bunch=1.0; 

//phi due to electron sheath at sides of bubble
   long double phi_sheath = 0.0*0.5;

//extra sheath thickness
   long double ratio_sheath_to_extrasheath=1.0;  

//________________________________________________________________
//PARABOLIC CHANNEL PROPERTIES

 bool if_channel=false;
// maximum potential, normalized to e/mc^2
   long double phi0_chan=1.0;//0.4*(r_bub0*r_bub0/rw0wp/rw0wp*0.5); 

// radius of channel
   long double r_channel=sqrt(2.0);//
   long double r_bub0=2.0*sqrt(1)*rw0wp;//sqrt(2.0); 

//*****************************************************************



//*****************************************************************
//****** ELECTRON BUNCH PARAMETERS ********
// Gaussian or tophat - check lower for namespace

// Number of particles in bunch
   int Npar = 1;
  bool if_LWFA_dist = false;
// These are required for radiation emission calculations
// Use macroparticle approximation
   bool if_macroparticle = false;
// How many fundamental particles does a macro-particle represent  
 long double sizeofbunch=1.0;//1e9/(long double) Npar;
// what is the density of the particle bunch relative to the critical density?
   long double densityofbunch=1.0;  

// width of bunch in x,y,z
   long double bunch_width[3] = {0.0,0.0,0.0}; 
//x,y,z of bunch peak/centre
   long double bunch_centre[3] = {0.0,0.0,1.0*laserpulse1::t0};

// Width in momentum space
 long double bunch_momentum_spread[3] = {0.0,0.0,0.0}; 

// Peak/Center momentum
   long double bunch_momentum[3] = {0.0,0.0,-100.0}; 

long double parzfocus = 0.0L;

// For large number of particles, may not want to dump data
   bool if_dump_par_dat = true; 

// Switch to be able to dump in lab time steps rather than particle
   bool dump_lab_time = true;

//*****************************************************************

//*****************************************************************
// IONIZATION
   bool ionization_on=false;
   bool poincare_on=false;
// Threshold PEAK a_0 for ionization
   long double ionization_threshold = 1.0;
//*****************************************************************


//*****************************************************************
//*********** UNIFORM MAGNETIC FIELD FOR DIAGNOSIS **********
// Uniform B field strength
 bool if_bfield = false;

 long double mag_Bfield = 0.0;
// B field direction unit vector, (x,y,z)
 long double B_field_hat[3] = {0.0,0.0,1.0};
//*****************************************************************


//*****************************************************************
//***********  MAGNETIC WIGGLER FOR DIAGNOSIS **********
/*
  Wiggler period defined in z direction, but with infinite extent in
  x and y directions. Can either be sinusoidal ("not square") or
  square pulse shape. wiggle is in x-z plane - i.e. B=By. 
 */
bool if_wiggler = false;
bool wiggler_square = false;
// magnetic field strength
long double wiggler_Bfield = 1.0;
// wiggler k vector: wiggler period = 2pi/k
long double wiggler_k = 1.0;
// wiggler start and end positions 
 long double wiggler_start = 0.1*PI;
 long double wiggler_end = 10.0*PI; // this is wiggler length now!
  long double wiggler_gap = 0.1*PI;

//*****************************************************************
//*********** NUMERICAL SPECTROMETER PARAMETERS **********
/*
  For calculating spectral power as a function of frequency for the 
  collection of particles
*/
// Switch spectrometer on/off (1/0)
   bool if_calc_spectrum=false; 

// if particle radiation should be added coherently (1) or not (0)
   bool coherent_addition=false; 

/*
  Spectral maximum/minimum frequency limits here
  Note these spectral limits are normalized to omega0
 */
   long double omega_max = 100000.0; //10.01;
   long double omega_min = 1000.0; //0.01;
   int N_omega_bins = 500;

// Space omega bins linearly (0) or logarithmically (1)
   int if_omega_log = 1; 

/*
  Scattering angles of observation to be calculated
  Angles are defined in DEGREES, and theta is measured from z axis 
  rotating in z-y plane. Azimuth angle at x=0
*/
 long double theta_scatter_min=PI;//-0.573*20.0;
 long double theta_scatter_max=PI;//+0.573*20.0;
   int N_theta_bins=1;
 
// Dump orthoganal polarizations separately or combined
   bool separate_polarizations=false; 

//________________________________________________________________
// This section is to do with the numerics of the Fourier transform

/*
 When switched on (1), the sinc functions arising due to the endpoints
 of the integral are subtracted, assuming that the velocity is 
   ant from -infinity to 0 and from taumax to infinity
 */
   bool if_endpoint_integral=false; 

 int order_of_spectrometer=2; // determines method of calculation
// 2 -> second order

/*
  For calculations where the endpoints are causing problems, a gaussian 
  shaped ramp can be applied to the endpoints to produce a window which
  damps artifacts from the sharp cutoff.
 */

//What damping do we want? 0 none, 1 end, 2 start and end
   bool if_endpoint_damping=false; 

//Width of damping window ramp
   long double damping_width=0.01*taumax;

// The width of the damping window ramp on the lower end
   long double damping_endpoint_ratio=0.01*taumax;

// A maximum time for calculation (0 means don't use)
   long double timemax=0.0;
//*****************************************************************

//*****************************************************************
//****** LASER PARAMETERS ********
/*
  Monochromatic gaussian laser beam - LASER PULSE 2
*/

namespace laserpulse2 {
// relative frequency (to normalization)
   long double rww0 = 1.0; 

// group velocity of laser
   long double v_group = 1.0;//-1.5/(rw0wp*rw0wp*rww0*rww0); // +Etching rate - Lu PRSTAB
  //sqrt(1.0-1.0/(rw0wp*rw0wp*rww0*rww0));

// phase velocity of laser
   long double v_phase =1.0;///sqrt(1.0-1.0/(rw0wp*rw0wp*rww0*rww0));//1.0/v_group;

// field strength -> !!!!TIME AVERAGED!!!! (i.e. peak is sqrt(2)* for linear
   long double a0=1.0;//3.5*0.0; 

//Propagation k fourvector - new - angle defined in zy plane from z axis
// IN radians
   long double theta_k0=PI;
 
//________________________________________________________________

/* 
   DEFINITION OF POLARIZATION
   THese add phase and angle to the two waves, both 0 is linear pol,
   one of them =pi/2 -> circular polarization or elliptical
   The angles must be transverse (pi/2) for circular, parallel for linear
 */
   long double phase1=0.0,phase2=0.0*PI*0.5;

//For phi, 0.0 is x, pi/2 is y
   long double phi_pol1 =0.0*PI*0.5,phi_pol2 = 0.0*PI*0.5; 
//________________________________________________________________

/*
 Pulse length normalized to omega0
 (half of 1/e^2 width, FWHM ~ 1.1 t0)
*/
   long double t0=50.0;//2.0*sqrt(5.0)*rw0wp; //65.0;

//Relative steepness of pulsefront
   long double pulse_front_steepness=1.0;//3.0; 

/* Pulse waist (for intensity)

   w0=0.0 means the waist is infinite (plane wave)

   Notes:
   remember, matched spot for bubble is 2sqrt(a)c/wp
*/
   long double w0=20.0;//2.0*sqrt(5.0)*rw0wp; 
// 
   long double delay=1.0*laserpulse2::t0;

//________________________________________________________________
// FOR FOCUSING PULSE
// Rayleigh range

   long double zR=w0*w0*0.5*rww0; 
   long double z_focus = 0.0; 

}



/*------------------------------------------------------------------------------

		OTHER THINGS THAT ARE NOW NO LONGER CONSTANT!
		
-----------------------------------------------------------------------------*/
//------------------------------------------------------------
// ****** OTHER    CONSTANTS **********
//------------------------------------------------------------
//these just make calculation quicker
// This is tau_0 in Rohrlich paper
long double dtau0 = taumax/(long double) Nmax;
long double inv_dtau0 = 1.0L/dtau0;
long double dtau0sq = dtau0*dtau0;
long double rad_const=mu0*e*e*omega_0/6.0L/PI/c/me*radiation_force_on; // 0.62x10^-23 s * omega0
long double qed_const = hbar*omega_0/me/c/c; // ratio of compton wavelength to laser wavelength
//long double spectrum_const= e*e*c/PI/PI/PI/16.0*mu0;
long double shape_factor_constant = pow(omega_0*e*e*mu0/me/c,1.0L/3.0L);
long double dtheta=(theta_scatter_max-theta_scatter_min)/(long double)(N_theta_bins-(1-!(N_theta_bins-1)));
long double r_sheath=(r_bub+sheaththickness);

long double inv_r_sheath_thick_square=1.0/(sheaththickness*sheaththickness);
long double ratio_bub_sheath = r_bub/sheaththickness;
long double inv_bub_sheath = sheaththickness/r_bub;
long double one_plus_rbs = (1.0+ratio_bub_sheath);
long double one_plus_rbs_2 =(1.0+ratio_bub_sheath)*2.0;

long double r_sheath_half=r_sheath-0.5*sheaththickness;
long double r_bubsquared=r_bub*r_bub;
long double inv_r_bubsquared=1.0/r_bubsquared;
long double inv_r_bub = 1.0/r_bub;


long double r_channelsquared=r_channel*r_channel;
long double inv_bubble_aspect=1.0/(bubble_aspect);

long double w0sq = laserpulse1::w0*laserpulse1::w0;
long double w0sq2 = laserpulse2::w0*laserpulse2::w0;
  
long double invzRsq = 1.0/laserpulse1::zR/laserpulse1::zR;
long double invzRsq2= 1.0/laserpulse2::zR/laserpulse2::zR;

long double inv_t0squared = 1.0/laserpulse1::t0/laserpulse1::t0;
long double inv_t0squared2 = 1.0/laserpulse2::t0/laserpulse2::t0;

// old and incorrect
//long double a0_pulse1 = laserpulse1::a0*sqrt((pow(cos(laserpulse1::phi_pol1)+cos(laserpulse1::phi_pol2),2.0L)
//						   +pow(sin(laserpulse1::phi_pol1)+sin(laserpulse1::phi_pol2),2.0L))*0.5);
//long double a0_pulse2 = laserpulse2::a0*sqrt((pow(cos(laserpulse2::phi_pol1)+cos(laserpulse2::phi_pol2),2.0L)
//						    +pow(sin(laserpulse2::phi_pol1)+sin(laserpulse2::phi_pol2),2.0L))*0.5);
long double a0_pulse1 = laserpulse1::a0/sqrt(1.0L + (sin(laserpulse1::phi_pol1)*sin(laserpulse1::phi_pol2)
					      +cos(laserpulse1::phi_pol1)*cos(laserpulse1::phi_pol2))
					      *cos(laserpulse1::phase1 - laserpulse1::phase2));
long double a0_pulse2 = laserpulse2::a0/sqrt(1.0L + (sin(laserpulse2::phi_pol1)*sin(laserpulse2::phi_pol2)
					      +cos(laserpulse2::phi_pol1)*cos(laserpulse2::phi_pol2))
					      *cos(laserpulse2::phase1 - laserpulse2::phase2));
int spectrumerror=0;
 int Ndumps=1000;
double yslice=0.0;
 int Nx_f=200,Nz_f=400; //for outputting fields - no. points
 long double xsize=100.0, zsize=150.0; // and spatial size

// DERIVED

long double particleisfree=1.0*!ionization_on;

long double theta_scatter=theta_scatter_min;

 fourvec *Kspec = NULL;

long double sinthetascatter = sin(theta_scatter);
long double costhetascatter = cos(theta_scatter);
long double excur_y=costhetascatter*costhetascatter;
long double excur_z = sinthetascatter*sinthetascatter;
long double domega = 1.0*(omega_max-omega_min)/N_omega_bins;
 
namespace laserpulse1 {
   fourvec * K0 = NULL;
   fourvec *dir0 = NULL;
   fourvec *rproj = NULL;
   fourvec *pol_1 = NULL;
   fourvec *pol_2 = NULL;
   fourvec *Kphase =  NULL;
   fourvec *Kgroup = NULL;
  long double focus_theta_1=0.5;
  
}
namespace laserpulse2 {
   fourvec *K0 =  NULL;
   fourvec *dir0 =  NULL;
   fourvec *rproj = NULL;
   fourvec *pol_1 = NULL;
   fourvec *pol_2 = NULL;
	fourvec *Kphase = NULL;
   fourvec *Kgroup = NULL;
  long double focus_theta_2=0.5;
}

// ARRAYS... Initialize as NULL and the set dynamically when running

// Parallel and perpendicular polarizations of radiation
long double **SpecInt_para=NULL;//[N_omega_bins][N_theta_bins];
long double **SpecInt_perp=NULL;//[N_omega_bins][N_theta_bins];

// These are a function of frequency and have Real and Imaginary parts
// Also x and y components, these are all combined to give final result

long double **Spe_In_Im_x=NULL;//[N_omega_bins][N_theta_bins];
long double **Spe_In_Re_x=NULL;//[N_omega_bins][N_theta_bins];
long double **Spe_In_Im_y=NULL;//[N_omega_bins][N_theta_bins];
long double **Spe_In_Re_y=NULL;//[N_omega_bins][N_theta_bins];
long double **Spe_In_Im_z=NULL;//[N_omega_bins][N_theta_bins];
long double **Spe_In_Re_z=NULL;//[N_omega_bins][N_theta_bins];
// To damp high freq errors.... numerical differencing
long double *sincwdt2=NULL;//[N_omega_bins];
long double *coswdt2=NULL;//[N_omega_bins];

// To represent finite bunching of electrons -> shape factor
long double *particle_shape_k=NULL;//[N_omega_bins];

long double *omega=NULL;//[N_omega_bins];
 
fourvec *int_x=NULL;//,[Nmax+1]
fourvec *int_v=NULL;//,[Nmax+1]
fourvec *int_S=NULL;//,[Nmax+1]
fourvec *int_x_m=NULL;//,[Nmax+1]
fourvec *int_v_m=NULL;//;[Nmax+1]
fourvec *int_S_m=NULL;//;[Nmax+1]

// spectral functions
long double *damping=NULL;//[Nmax+1];
long double *omega_r1=NULL;//[N_omega_bins];
long double *omega_r2=NULL;//[N_omega_bins];

// From main
long double **par_x=NULL;//[4][Npar];
long double **par_v=NULL;//[4][Npar];
long double **par_S=NULL;//[4][Npar];

int (*initialize_particles)(long double **par_x,long double **par_v,long double **par_S)=NULL;
long double (*a_laser_t_env)(fourvec &x,long double phase_in_exponent,long double inv_t0squared)=NULL;
fourvec (*Lambda)(fourvec& x,fourvec& v,fourvec& vdot,fourvec& vdotdot)=NULL;
int (*a_laser)(fourvec x,fourvec &ans)=NULL;
fourvec (*a_wakefield)(fourvec &x)=NULL;
long double (*ramp)(fourvec x)=NULL;
long double (*beam_waist_sq)(fourvec &x)=NULL;
long double (*inv_beam_Radius)(fourvec &x)=NULL;
long double (*Gouy_phase)(fourvec &x)=NULL;
long double (*beam_waist_sq2)(fourvec &x)=NULL;
long double (*inv_beam_Radius2)(fourvec &x)=NULL;
long double (*Gouy_phase2)(fourvec &x)=NULL;

// model switching
bool if_two_pulses=false;
bool if_one_pulses=false;
bool if_with_focusing=true;
bool if_square_pulse=false;
bool if_trunc_pulse=true;
bool if_PLAD=false;
bool if_QPLAD=false;
bool if_uniform_bunch=false;

bool if_import_data=false;

// PARAMETER LIST
namespace paramlist
{
const int param_length = 128;
std::string parameter_list[param_length] = {"ifdumpnew","Nmax","taumax","time_max","omega_0","Nx_f","Nz_f","xsize","zsize","rw0wp","radiation_force_on","if_one_pulses","laserpulse1::rww0","laserpulse1::v_group","laserpulse1::v_phase","laserpulse1::a0","laserpulse1::theta_k0","laserpulse1::phase1","laserpulse1::phase2","laserpulse1::phi_pol1","laserpulse1::phi_pol2","laserpulse1::t0","laserpulse1::pulse_front_steepness","laserpulse1::w0","laserpulse1::delay","laserpulse1::z_focus","if_dump_fields","if_moving_frame","n_field_dumps","if_wakefield","if_bubble","r_bub","phi0_bub","bubble_aspect","v_wake","bub_delay","bub_width","bub_ramp","sheaththickness","phi_bunch","ratio_sheath_to_phi_bunch","phi_sheath","ratio_sheath_to_extrasheath","if_channel","phi0_chan","r_channel","r_bub0","Npar","if_macroparticle","if_LWFA_dist","sizeofbunch","densityofbunch","bunch_width[0]","bunch_width[1]","bunch_width[2]","bunch_centre[0]","bunch_centre[1]","bunch_centre[2]","bunch_momentum_spread[0]","bunch_momentum_spread[1]","bunch_momentum_spread[2]","bunch_momentum[0]","bunch_momentum[1]","bunch_momentum[2]","Ndumps","q_over_m","if_dump_par_dat","dump_lab_time","ionization_on","ionization_threshold","if_bfield","mag_Bfield","B_field_hat[0]","B_field_hat[1]","B_field_hat[2]","if_wiggler","wiggler_square","wiggler_Bfield","wiggler_k","wiggler_start","wiggler_end","wiggler_gap","if_calc_spectrum","coherent_addition","omega_max","omega_min","N_omega_bins","if_omega_log","theta_scatter_min","theta_scatter_max","N_theta_bins","separate_polarizations","if_endpoint_integral","order_of_spectrometer","if_endpoint_damping","damping_width","damping_endpoint_ratio","timemax","if_two_pulses","laserpulse2::rww0","laserpulse2::v_group","laserpulse2::v_phase","laserpulse2::a0","laserpulse2::theta_k0","laserpulse2::phase1","laserpulse2::phase2","laserpulse2::phi_pol1","laserpulse2::phi_pol2","laserpulse2::t0","laserpulse2::pulse_front_steepness","laserpulse2::w0","laserpulse2::delay","laserpulse2::z_focus","if_with_focusing","if_square_pulse","if_trunc_pulse","if_PLAD","if_QPLAD","if_uniform_bunch","yslice","parzfocus","if_standing_wave","if_standing_circ","standing_wave_k","standing_wave_a","standing_wave_phase","standing_wave_pol"};
long double parameters[param_length] = {ifdumpnew,  Nmax,   taumax,   time_max ,   omega_0,   Nx_f,   Nz_f,   xsize,   zsize,   rw0wp,   radiation_force_on,     if_one_pulses,     laserpulse1::rww0,     laserpulse1::v_group,      laserpulse1::v_phase,      laserpulse1::a0,      laserpulse1::theta_k0,   laserpulse1::phase1,  laserpulse1::phase2,   laserpulse1::phi_pol1,  laserpulse1::phi_pol2,     laserpulse1::t0,     laserpulse1::pulse_front_steepness,     laserpulse1::w0,      laserpulse1::delay,   laserpulse1::z_focus ,    if_dump_fields ,    if_moving_frame ,    n_field_dumps ,   if_wakefield ,  if_bubble ,     r_bub,     phi0_bub,     bubble_aspect,     v_wake ,     bub_delay,   bub_width,   bub_ramp,     sheaththickness,   phi_bunch,   ratio_sheath_to_phi_bunch,  phi_sheath,     ratio_sheath_to_extrasheath,   if_channel,     phi0_chan,     r_channel,     r_bub0,     Npar ,     if_macroparticle,        if_LWFA_dist,   sizeofbunch,     densityofbunch,     bunch_width[0],     bunch_width[1],     bunch_width[2],     bunch_centre[0] ,     bunch_centre[1] ,     bunch_centre[2] ,   bunch_momentum_spread[0]  ,   bunch_momentum_spread[1]  ,   bunch_momentum_spread[2]  ,     bunch_momentum[0]  ,     bunch_momentum[1]  ,     bunch_momentum[2]  ,     Ndumps  ,   q_over_m,     if_dump_par_dat ,     dump_lab_time ,     ionization_on ,     ionization_threshold ,   if_bfield ,   mag_Bfield ,   B_field_hat[0] ,   B_field_hat[1] ,   B_field_hat[2] ,  if_wiggler ,  wiggler_square,  wiggler_Bfield,  wiggler_k ,   wiggler_start ,   wiggler_end ,    wiggler_gap ,   if_calc_spectrum,     coherent_addition,     omega_max ,     omega_min ,     N_omega_bins ,     if_omega_log ,   theta_scatter_min,   theta_scatter_max,     N_theta_bins,     separate_polarizations,     if_endpoint_integral,   order_of_spectrometer,     if_endpoint_damping,     damping_width,     damping_endpoint_ratio,   timemax,     if_two_pulses,   laserpulse2::rww0,     laserpulse2::v_group,      laserpulse2::v_phase,      laserpulse2::a0,      laserpulse2::theta_k0,   laserpulse2::phase1,  laserpulse2::phase2,   laserpulse2::phi_pol1,  laserpulse2::phi_pol2,     laserpulse2::t0,     laserpulse2::pulse_front_steepness,     laserpulse2::w0,      laserpulse2::delay,   laserpulse2::z_focus ,  if_with_focusing,  if_square_pulse,  if_trunc_pulse,  if_PLAD,  if_QPLAD,  if_uniform_bunch,  yslice,  parzfocus,if_standing_wave,if_standing_circ,standing_wave_k,standing_wave_a,standing_wave_phase,standing_wave_pol};
}

int findwhichslot(std::string stringtemp) {
for (int ii=0;ii<paramlist::param_length;++ii) 
	{
		if (stringtemp==paramlist::parameter_list[ii]) {
			return ii;
			}
	}
	return -1;
}
int load_generic(std::string FileName) {
std::ifstream defaultfilein(FileName.c_str());
std::string stringtemp;
std::string stringinputtemp;
long double doubletemp;
int whichslot=-1;
int whichline=0;
std::string line;
while (std::getline(defaultfilein, line))
{
++whichline;
//while (defaultfilein >> stringtemp >> doubletemp)
//{
	std::istringstream iss(line);
	std::istringstream iss2(line);
	std::istringstream iss3(line);
	stringtemp="";
   // ignore comments!
   	iss3 >> stringtemp; 
   //	cout<<stringtemp<<'\n';
	if (stringtemp[0] == '#'||stringtemp.empty()) {
	whichslot=-1;
	} else {
	
	if (!(iss >> stringtemp >> doubletemp)) { 
	if (!(iss2 >> stringtemp >> stringinputtemp)) { 
	std::cout<<"Error in input file on line "<<whichline<<'\n';
	exit(0); 
	} else {
		if (stringinputtemp=="true" || stringinputtemp=="TRUE") {
			doubletemp = 1;
		} else {
			doubletemp=0;
		}
	}
	}
   whichslot = findwhichslot(stringtemp);
   }
	
   if (whichslot>0){
   paramlist::parameters[whichslot] = doubletemp;
   }
		
}

/*for (int ii=0;ii<paramlist::param_length;++ii){
cout<<paramlist::parameter_list[ii]<<"\t"<<paramlist::parameters[ii]<<"\n";
}
*/
//exit(0);
set_all_parameters();

 defaultfilein.close();
return 0;
}

int save_generic(std::string FileName) {
std::ofstream defaultfileout(FileName.c_str());
for (int ii=0;ii<paramlist::param_length;++ii){
defaultfileout<< paramlist::parameter_list[ii]<<"\t"<<paramlist::parameters[ii]<<"\n";
}
 defaultfileout.close();
return 0;
}

int save_defaults() {
	save_generic("defaults.dat");
	return 0;
}
int load_defaults() {
std::string unix_com = "test -f defaults.dat";
int result = system(unix_com.c_str());

 if (!result){
 load_generic("defaults.dat");
 } else {
 save_defaults();
 }
 return 0;
 }
 
 
int set_all_parameters(void)
{
ifdumpnew=paramlist::parameters[0];
Nmax=paramlist::parameters[1];
taumax=paramlist::parameters[2];
time_max=paramlist::parameters[3];
omega_0=paramlist::parameters[4];
Nx_f=paramlist::parameters[5];
Nz_f=paramlist::parameters[6];
xsize=paramlist::parameters[7];
zsize=paramlist::parameters[8];
rw0wp=paramlist::parameters[9];
radiation_force_on=paramlist::parameters[10];
if_one_pulses=paramlist::parameters[11];
laserpulse1::rww0=paramlist::parameters[12];
laserpulse1::v_group=paramlist::parameters[13];
laserpulse1::v_phase=paramlist::parameters[14];
laserpulse1::a0=paramlist::parameters[15];
laserpulse1::theta_k0=paramlist::parameters[16];
laserpulse1::phase1=paramlist::parameters[17];
laserpulse1::phase2=paramlist::parameters[18];
laserpulse1::phi_pol1=paramlist::parameters[19];
laserpulse1::phi_pol2=paramlist::parameters[20];
laserpulse1::t0=paramlist::parameters[21];
laserpulse1::pulse_front_steepness=paramlist::parameters[22];
laserpulse1::w0=paramlist::parameters[23];
laserpulse1::delay=paramlist::parameters[24];
laserpulse1::z_focus=paramlist::parameters[25];
if_dump_fields=paramlist::parameters[26];
if_moving_frame=paramlist::parameters[27];
n_field_dumps=paramlist::parameters[28];
if_wakefield=paramlist::parameters[29];
if_bubble=paramlist::parameters[30];
r_bub=paramlist::parameters[31];
phi0_bub=paramlist::parameters[32];
bubble_aspect=paramlist::parameters[33];
v_wake=paramlist::parameters[34];
bub_delay=paramlist::parameters[35];
bub_width=paramlist::parameters[36];
bub_ramp=paramlist::parameters[37];
sheaththickness=paramlist::parameters[38];
phi_bunch=paramlist::parameters[39];
ratio_sheath_to_phi_bunch=paramlist::parameters[40];
phi_sheath=paramlist::parameters[41];
ratio_sheath_to_extrasheath=paramlist::parameters[42];
if_channel=paramlist::parameters[43];
phi0_chan=paramlist::parameters[44];
r_channel=paramlist::parameters[45];
r_bub0=paramlist::parameters[46];
Npar=paramlist::parameters[47];
if_macroparticle=paramlist::parameters[48];
if_LWFA_dist=paramlist::parameters[49];
sizeofbunch=paramlist::parameters[50];
densityofbunch=paramlist::parameters[51];
bunch_width[0]=paramlist::parameters[52];
bunch_width[1]=paramlist::parameters[53];
bunch_width[2]=paramlist::parameters[54];
bunch_centre[0]=paramlist::parameters[55];
bunch_centre[1]=paramlist::parameters[56];
bunch_centre[2]=paramlist::parameters[57];
bunch_momentum_spread[0]=paramlist::parameters[58];
bunch_momentum_spread[1]=paramlist::parameters[59];
bunch_momentum_spread[2]=paramlist::parameters[60];
bunch_momentum[0]=paramlist::parameters[61];
bunch_momentum[1]=paramlist::parameters[62];
bunch_momentum[2]=paramlist::parameters[63];
Ndumps=paramlist::parameters[64];
q_over_m=paramlist::parameters[65];
if_dump_par_dat=paramlist::parameters[66];
dump_lab_time=paramlist::parameters[67];
ionization_on=paramlist::parameters[68];
ionization_threshold=paramlist::parameters[69];
if_bfield=paramlist::parameters[70];
mag_Bfield=paramlist::parameters[71];
B_field_hat[0]=paramlist::parameters[72];
B_field_hat[1]=paramlist::parameters[73];
B_field_hat[2]=paramlist::parameters[74];
if_wiggler=paramlist::parameters[75];
wiggler_square=paramlist::parameters[76];
wiggler_Bfield=paramlist::parameters[77];
wiggler_k=paramlist::parameters[78];
wiggler_start=paramlist::parameters[79];
wiggler_end=paramlist::parameters[80];
wiggler_gap=paramlist::parameters[81];
if_calc_spectrum=paramlist::parameters[82];
coherent_addition=paramlist::parameters[83];
omega_max=paramlist::parameters[84];
omega_min=paramlist::parameters[85];
N_omega_bins=paramlist::parameters[86];
if_omega_log=paramlist::parameters[87];
theta_scatter_min=paramlist::parameters[88];
theta_scatter_max=paramlist::parameters[89];
N_theta_bins=paramlist::parameters[90];
separate_polarizations=paramlist::parameters[91];
if_endpoint_integral=paramlist::parameters[92];
order_of_spectrometer=paramlist::parameters[93];
if_endpoint_damping=paramlist::parameters[94];
damping_width=paramlist::parameters[95];
damping_endpoint_ratio=paramlist::parameters[96];
timemax=paramlist::parameters[97];
if_two_pulses=paramlist::parameters[98];
laserpulse2::rww0=paramlist::parameters[99];
laserpulse2::v_group=paramlist::parameters[100];
laserpulse2::v_phase=paramlist::parameters[101];
laserpulse2::a0=paramlist::parameters[102];
laserpulse2::theta_k0=paramlist::parameters[103];
laserpulse2::phase1=paramlist::parameters[104];
laserpulse2::phase2=paramlist::parameters[105];
laserpulse2::phi_pol1=paramlist::parameters[106];
laserpulse2::phi_pol2=paramlist::parameters[107];
laserpulse2::t0=paramlist::parameters[108];
laserpulse2::pulse_front_steepness=paramlist::parameters[109];
laserpulse2::w0=paramlist::parameters[110];
laserpulse2::delay=paramlist::parameters[111];
laserpulse2::z_focus=paramlist::parameters[112];
if_with_focusing=paramlist::parameters[113];
if_square_pulse=paramlist::parameters[114];
if_trunc_pulse=paramlist::parameters[115];
if_PLAD=paramlist::parameters[116];
if_QPLAD=paramlist::parameters[117];
if_uniform_bunch=paramlist::parameters[118];
yslice=paramlist::parameters[119];
parzfocus=paramlist::parameters[120];
if_standing_wave = paramlist::parameters[121];
if_standing_circ = paramlist::parameters[122];
standing_wave_k = paramlist::parameters[123];
standing_wave_a = paramlist::parameters[124];
standing_wave_phase = paramlist::parameters[125];
standing_wave_pol = paramlist::parameters[126];
return 0;	
}
 
 
#endif
