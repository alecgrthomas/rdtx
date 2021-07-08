

/*
  ********* RDTX *****************
  This is an electron tracking code **including radiation reaction** 
  forces and high order differencing, used for calculating radiation
  spectra in laser-particle interactions

  The finite difference equations use 4th order Runga Kutta
  
  Variables:
  v is four velocity i.e. (E,p)
  x is (ct,x)
  proper time tau = int (1/gamma )dt 

  Note this is in normalized units -> to laser frequency and c/w0 space
  
  Radiation spectrum calculation through fourier like integral
  Integral done numerically by taking quadratic interpolation and using
  exact analytic integral between timesteps

 The radiation can be observed in the z-y plane, with angle from z axis 
 specified - rotate the polarization of the laser to get azimuthal info

 Input parameters specified in input_params.h
 
 New version includes Spin dynamics through the BMT equation
 
 
*/

#pragma once

#ifndef HEADER_H
#define HEADER_H

#ifndef __cplusplus
#error A C++ compiler is required!
#endif 

#include <iostream>
#include <cmath>
#include <limits>
#include <complex>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "fourvec.h"
#include "fourtens.h"
#include "global_constants.h"
#include "constants.h"

//using namespace std;

//**********************************************************************
// Function declarations
//**********************************************************************
int rdtx_run(fourvec *x0,fourvec *x1,fourvec *x2,
             fourvec *v0,fourvec *v1);

long double string_eval(std::string expr) throw(int);

extern int n_field_dumps;

//int set_initial(RunPanel *OutputPanel);
int set_initial();
int calculatedamping();
int dump_fields(std::ofstream &outfile_fields,std::ofstream &outfile_potentials,long double time);
int defineshapefactor();

int ResetBins();
int ResetSpec();
long double par_traj(fourvec &v0,fourvec &x0,fourvec &S0, int par_label,std::ofstream &outfile_summary);
int SpectralEndpoint(int theta, long double tau, fourvec v,fourvec x,long double sgn);
int SpectralCalculation(int theta, int parnum);

int SpectralCalculationImport(int theta, int parnum,   fourvec *&x0_in,
                              fourvec *&x1_in,fourvec *&x2_in,fourvec *&v0_in,
                              fourvec *&v1_in,double dt);
int Calc_Spec(int theta);
int Write_Spec(int parnum,int theta);
void operator<< (std::ofstream& str_out, fourvec vec);
int write_data(std::ofstream &outfile, fourvec x, fourvec v,fourvec S, double Pave, fourvec a_las);
int write_parhead(std::ofstream &outfile,int par_label,int Ndumps_mod,double v_wake);

namespace onelaserpulse {
  int a_laser(fourvec x,fourvec &ans);
}
namespace twolaserpulses {
  int a_laser(fourvec x,fourvec &ans);
}
int a_field(fourvec x,fourvec &ans);

 
//**********************************************************************

// External variable definitions
extern  fourvec *int_x, *int_v, *int_x_m, *int_v_m;
extern  long double particleisfree;

extern  long double domega;
//**********************************************************************
// Global constants
//**********************************************************************
extern  long double theta_scatter;
extern  long double sinthetascatter;
extern  long double costhetascatter;
extern  long double excur_y;
extern  long double excur_z;

// Parallel and perpendicular polarizations of radiation
extern long double ** SpecInt_para;
extern long double ** SpecInt_perp;

// These are a function of frequency and have Real and Imaginary parts
// Also x and y components, these are all combined to give final result

extern  long double ** Spe_In_Im_x;
extern  long double ** Spe_In_Re_x;
extern  long double ** Spe_In_Im_y;
extern  long double ** Spe_In_Re_y;
extern  long double ** Spe_In_Im_z;
extern  long double ** Spe_In_Re_z;
// To damp high freq errors.... numerical differencing
extern  long double * sincwdt2;
extern  long double * coswdt2;

// To represent finite bunching of electrons -> shape factor
extern  long double * particle_shape_k;

extern  long double * omega;
extern  long double phasew;

extern std::string GetstrTime(long double time);


extern  long double totalErad;

extern std::string Data_directory;
extern std::string Static_Run_directory;
extern std::string Run_directory;
extern std::string Field_directory;
extern std::string filename_field;
extern std::string filename_potential;
extern std::string filename;
extern std::string filename_summary;
extern std::string filename_spectrum;
extern std::string filename_spec;

extern fourtens metrictens; // The metric tensor

extern long double max_excursion[3], min_excursion[3];

extern fourvec * int_x, * int_v, * int_S, * int_x_m, * int_v_m, * int_S_m;
extern long double particleisfree;

extern long double totalErad;
extern long double phasew;
extern long double particleisfree;

extern fourvec * Kspec;
extern long double theta_scatter;
extern long double sinthetascatter;
extern long double costhetascatter;
extern long double excur_y;
extern long double excur_z;
extern long double domega;

namespace laserpulse1 {
  extern fourvec *K0;
  extern fourvec *dir0;
  extern fourvec *rproj;
  extern fourvec *pol_1; //polarized in x
  extern fourvec *pol_2; //polarized in y
  extern fourvec *Kphase; // laser phase k four vector
  extern fourvec *Kgroup; // laser group k four vector
  extern long double focus_theta_1; // divergence of laser in far field
}
namespace laserpulse2 {
  extern fourvec *K0;
  extern fourvec *dir0;
  extern fourvec *rproj;
  extern fourvec *pol_1; //polarized in x
  extern fourvec *pol_2; //polarized in y
  extern fourvec *Kphase; // laser phase k four vector
  extern fourvec *Kgroup; // laser group k four vector
  extern long double focus_theta_2; // divergence of laser in far field
}
namespace with_focusing
{
  long double beam_waist_sq(fourvec &x);
  long double inv_beam_Radius(fourvec &x);
  long double Gouy_phase(fourvec &x);
  long double beam_waist_sq2(fourvec &x);
  long double inv_beam_Radius2(fourvec &x);
  long double Gouy_phase2(fourvec &x);
}
namespace no_focusing
{
  long double beam_waist_sq(fourvec &x);
  long double inv_beam_Radius(fourvec &x);
  long double Gouy_phase(fourvec &x);
  long double beam_waist_sq2(fourvec &x);
  long double inv_beam_Radius2(fourvec &x);
  long double Gouy_phase2(fourvec &x);
}

extern bool if_standing_wave;
extern bool if_standing_circ;
extern long double standing_wave_k;
extern long double standing_wave_a;
extern long double standing_wave_phase;
extern long double standing_wave_pol; // polarization angle with respect to x axis
extern long double costhetastanding;
extern long double sinthetastanding;


namespace onelaserpulse {
  int a_laser(fourvec x,fourvec &ans);
}
namespace twolaserpulses {
  int a_laser(fourvec x,fourvec &ans);
}
namespace otherlaserpulse {
  int a_laser(fourvec x,fourvec &ans);
}
namespace nolaserpulse {
  int a_laser(fourvec x,fourvec &ans);
}
fourvec a_channel(fourvec &x);
namespace bubblewake {
fourvec a_wakefield(fourvec &x);
}
namespace linearwake {
fourvec a_wakefield(fourvec &x);
}
namespace withramp{
long double ramp(fourvec x);
}
namespace withnoramp{
long double ramp(fourvec x);
}
fourvec a_Bfield(fourvec &x);
long double a_laser_envelope(fourvec(&x));
long double a_laser_envelope2(fourvec(&x));

extern fourtens metrictens; // The metric tensor
double Kirk_g(fourvec);
fourtens F(fourvec& x);

// Spin functions
fourvec BMT(long double &tau,fourvec& x,fourvec& v,fourvec& S);


namespace norad // no radiation force
{
  fourvec Lambda(fourvec& x,fourvec& v,fourvec& vdot,fourvec& vdotdot);
}
namespace LAD // Lorentz-Abraham-Dirac force
{
  fourvec Lambda(fourvec& x,fourvec& v,fourvec& vdot,fourvec& vdotdot);
}
namespace PLAD // Physical Lorentz-Abraham-Dirac force - Rohlich PRE 08
{
  fourvec Lambda(fourvec& x,fourvec& v,fourvec& vdot,fourvec& vdotdot);
}
namespace QPLAD // Quantum Corrected Physical Lorentz-Abraham-Dirac force - Rohlich PRE 08
{
  fourvec Lambda(fourvec& x,fourvec& v,fourvec& vdot,fourvec& vdotdot);
}



namespace uniform_bunch
{
  int initialize_particles(long double **par_x,long double **par_v,long double **par_S);
}
namespace gaussian_bunch
{
  int initialize_particles(long double **par_x,long double **par_v,long double **par_S);
}

namespace asymmetric_gaussian_pulse
{
  long double a_laser_t_env(fourvec &x,long double phase_in_exponent,long double inv_t0squared);
}
namespace square_pulse
{
  long double a_laser_t_env(fourvec &x,long double phase_in_exponent,long double inv_t0squared);
}

/*From spectral functions*/
extern long double *damping;//[Nmax+1];
extern long double *omega_r1;//[N_omega_bins];
extern long double *omega_r2;//[N_omega_bins];

// From main
extern  long double **par_x;//[4][Npar];
extern  long double **par_v;//[4][Npar];
extern  long double **par_S;//[4][Npar];

// Function pointers for namespaces switches
extern int (*initialize_particles)(long double **par_x,long double **par_v,long double **par_S);
extern long double (*a_laser_t_env)(fourvec &x,long double phase_in_exponent,long double inv_t0squared);
extern fourvec (*Lambda)(fourvec& x,fourvec& v,fourvec& vdot,fourvec& vdotdot);
extern int (*a_laser)(fourvec x,fourvec &ans);
extern long double (*ramp)(fourvec x);
extern fourvec (*a_wakefield)(fourvec &x);
extern long double (*beam_waist_sq)(fourvec &x);
extern long double (*inv_beam_Radius)(fourvec &x);
extern long double (*Gouy_phase)(fourvec &x);
extern long double (*beam_waist_sq2)(fourvec &x);
extern long double (*inv_beam_Radius2)(fourvec &x);
extern long double (*Gouy_phase2)(fourvec &x);

extern bool if_import_data;
#include "input_params.h"

#endif
