#ifndef RDTX_C_
#define RDTX_C_

#include "headers.h"
 #include <sys/stat.h>
int save_defaults();
int load_defaults();
int save_generic(std::string FileName);
int load_generic(std::string FileName);

int main(int argnum, char **argstr_in)
{
std::string fileName = "defaults.dat";
 std::cout<<"\n---- RDTX2 v"<<version<<" ----\n";
if (argnum>1){
	fileName =  argstr_in[1];
	}
 struct stat buffer;   
bool ifexist = (stat (fileName.c_str(), &buffer) == 0); 	

if (ifexist) {
  std::cout<<"Loading parameters from file '"<<fileName<<"'\n\n";
	} else {
	if (argnum>1){
		std::cout<< "File '"<< fileName<< "' does not exist\n\n";
		exit(0);
		} else {
		std::cout<<"Setting default parameters to file '"<<fileName<<"'\n\n";
		save_generic(fileName);
		exit(0);
		}
	}

	


	// First set up simulation!
 	load_generic(fileName);
 
	// Set derived params ------------------------------------------------------
	dtau0 = taumax/(long double) Nmax;
	inv_dtau0 = 1.0L/dtau0;
	dtau0sq = dtau0*dtau0;
	rad_const=mu0*e*e*omega_0/6.0L/PI/c/me*radiation_force_on; // 0.62x10^-23 s * omega0
	//spectrum_const= e*e*c/PI/PI/PI/16.0*mu0;
        spectrumerror=0;
	shape_factor_constant = pow(omega_0*e*e*mu0/me/c,1.0L/3.0L);
	dtheta=(theta_scatter_max-theta_scatter_min)/(long double)(N_theta_bins-(1-!(N_theta_bins-1)));
	r_sheath=(r_bub+sheaththickness);
	inv_r_sheath_thick_square=1.0L/(sheaththickness*sheaththickness);
	ratio_bub_sheath = r_bub/sheaththickness;
	inv_bub_sheath = sheaththickness/r_bub;
	one_plus_rbs = (1.0L+ratio_bub_sheath);
	one_plus_rbs_2 =(1.0L+ratio_bub_sheath)*2.0L;
	r_sheath_half=r_sheath-0.5L*sheaththickness;
	r_bubsquared=r_bub*r_bub;
	inv_r_bubsquared=1.0L/r_bubsquared;
	inv_r_bub = 1.0L/r_bub;
	r_channelsquared=r_channel*r_channel;
	inv_bubble_aspect=1.0L/(bubble_aspect);
	w0sq = laserpulse1::w0*laserpulse1::w0;
	w0sq2 = laserpulse2::w0*laserpulse2::w0;
	laserpulse1::zR=laserpulse1::w0*laserpulse1::w0*0.5*laserpulse1::rww0;
	laserpulse2::zR=laserpulse2::w0*laserpulse2::w0*0.5*laserpulse2::rww0;
	invzRsq = 1.0L/laserpulse1::zR/laserpulse1::zR;
	invzRsq2= 1.0L/laserpulse2::zR/laserpulse2::zR;
	inv_t0squared = 1.0L/laserpulse1::t0/laserpulse1::t0;
	inv_t0squared2 = 1.0L/laserpulse2::t0/laserpulse2::t0;
  	a0_pulse1 = laserpulse1::a0/sqrt(1.0L + (sin(laserpulse1::phi_pol1)*sin(laserpulse1::phi_pol2)
					      +cos(laserpulse1::phi_pol1)*cos(laserpulse1::phi_pol2))
					      *cos(laserpulse1::phase1 - laserpulse1::phase2));
	a0_pulse2 = laserpulse2::a0/sqrt(1.0L + (sin(laserpulse2::phi_pol1)*sin(laserpulse2::phi_pol2)
					      +cos(laserpulse2::phi_pol1)*cos(laserpulse2::phi_pol2))
					      *cos(laserpulse2::phase1 - laserpulse2::phase2));

	particleisfree=1.0L*!ionization_on;
	theta_scatter=theta_scatter_min;
 	Kspec = new fourvec(1.0L,0.0L,sin(theta_scatter),cos(theta_scatter));
	sinthetascatter = sin(theta_scatter);
	costhetascatter = cos(theta_scatter);
	excur_y=costhetascatter*costhetascatter;
	excur_z = sinthetascatter*sinthetascatter;
	domega = 1.0L*(omega_max-omega_min)/N_omega_bins;
	
	sinthetastanding = sin(standing_wave_pol);
	costhetastanding = cos(standing_wave_pol);

	
	laserpulse1::K0 = new fourvec(1.0L,0.0L,sin(laserpulse1::theta_k0),cos(laserpulse1::theta_k0));
    laserpulse1::dir0 = new fourvec(0.0L,0.0L,-sin(laserpulse1::theta_k0),-cos(laserpulse1::theta_k0)); // this defines the propagatino direction
    laserpulse1::rproj = new fourvec(0.0L,0.0L,-cos(laserpulse1::theta_k0),sin(laserpulse1::theta_k0)); // this defines the radial direction
    laserpulse1::pol_1 = new fourvec(0.0L,cos(laserpulse1::phi_pol1),sin(laserpulse1::phi_pol1),0.0L); //polarized in x
    laserpulse1::pol_2 = new fourvec(0.0L,cos(laserpulse1::phi_pol2),sin(laserpulse1::phi_pol2),0.0L); //polarized in y
    laserpulse1::Kphase = new fourvec((*laserpulse1::K0)*laserpulse1::rww0); // laser phase  k four vector
    laserpulse1::Kgroup = new fourvec((*laserpulse1::K0)); // laser group  k four vector
	
    laserpulse1::focus_theta_1=laserpulse1::w0/laserpulse1::zR;
    
   laserpulse2::K0 = new fourvec(1.0L,0.0,sin(laserpulse2::theta_k0),cos(laserpulse2::theta_k0));
   laserpulse2::dir0 = new fourvec(0.0L,0.0L,-sin(laserpulse2::theta_k0),-cos(laserpulse2::theta_k0)); // this defines the propagatino direction
   laserpulse2::rproj = new fourvec(0.0L,0.0L,-cos(laserpulse2::theta_k0),sin(laserpulse2::theta_k0)); // this defines the radial direction
   laserpulse2::pol_1 = new fourvec(0.0L,cos(laserpulse2::phi_pol1),sin(laserpulse2::phi_pol1),0.0L); //polarized in x
   laserpulse2::pol_2 = new fourvec(0.0L,cos(laserpulse2::phi_pol2),sin(laserpulse2::phi_pol2),0.0L); //polarized in y
	laserpulse2::Kphase = new fourvec((*laserpulse2::K0)*laserpulse2::rww0); // laser phase  k four vector
   laserpulse2::Kgroup = new fourvec((*laserpulse2::K0)); // laser group  k four vector
	
	laserpulse2::focus_theta_2=laserpulse2::w0/laserpulse2::zR;
	
	// Initialize Empty Arrays -------------------------------------------------
	SpecInt_para= new long double *[N_omega_bins];//[N_theta_bins];
	SpecInt_perp=new long double *[N_omega_bins];//[N_theta_bins];
	Spe_In_Im_x=new long double *[N_omega_bins];//[N_theta_bins];
	Spe_In_Re_x=new long double *[N_omega_bins];//[N_theta_bins];
	Spe_In_Im_y=new long double *[N_omega_bins];//[N_theta_bins];
	Spe_In_Re_y=new long double *[N_omega_bins];//[N_theta_bins];
	Spe_In_Im_z=new long double *[N_omega_bins];//[N_theta_bins];
	Spe_In_Re_z=new long double *[N_omega_bins];//[N_theta_bins];
	
	for (int i=0;i<N_omega_bins;++i){
		SpecInt_para[i]= new long double [N_theta_bins];
		SpecInt_perp[i]=new long double [N_theta_bins];
		Spe_In_Im_x[i]=new long double [N_theta_bins];
		Spe_In_Re_x[i]=new long double [N_theta_bins];
		Spe_In_Im_y[i]=new long double [N_theta_bins];
		Spe_In_Re_y[i]=new long double [N_theta_bins];
		Spe_In_Im_z[i]=new long double [N_theta_bins];
		Spe_In_Re_z[i]=new long double [N_theta_bins];
	}

	sincwdt2=new long double [N_omega_bins];
	coswdt2=new long double [N_omega_bins];
	particle_shape_k=new long double [N_omega_bins];
	omega=new long double [N_omega_bins];

	int_x=new fourvec [Nmax+1];
	int_v=new fourvec [Nmax+1];
	int_S=new fourvec [Nmax+1];
	int_x_m=new fourvec [Nmax+1];
	int_v_m=new fourvec [Nmax+1];
	int_S_m=new fourvec [Nmax+1];

	damping=new long double [Nmax+1];
	omega_r1=new long double [N_omega_bins];
	omega_r2=new long double [N_omega_bins];
	
	par_x=new long double *[4];//[Npar];
	par_v=new long double *[4];//[Npar];
	par_S=new long double *[4];//[Npar];
	
	for (int i=0;i<4;++i){
		par_x[i]=new long double [Npar];
		par_v[i]=new long double [Npar];
		par_S[i]=new long double [Npar];
	}
	
	// Set up the models
	if (if_two_pulses&&if_one_pulses){
		a_laser=&twolaserpulses::a_laser;}
	else{
		if (if_two_pulses) {
			a_laser=&otherlaserpulse::a_laser;
			} else {
				if (if_one_pulses) {
				a_laser=&onelaserpulse::a_laser;
				} else {
				a_laser=&nolaserpulse::a_laser;
				}
			}
	}
	if (if_with_focusing){
		beam_waist_sq=&with_focusing::beam_waist_sq;
		inv_beam_Radius=&with_focusing::inv_beam_Radius;
		Gouy_phase=&with_focusing::Gouy_phase;
		beam_waist_sq2=&with_focusing::beam_waist_sq2;
		inv_beam_Radius2=&with_focusing::inv_beam_Radius2;
		Gouy_phase2=&with_focusing::Gouy_phase2;}
	else{
		beam_waist_sq=&no_focusing::beam_waist_sq;
		inv_beam_Radius=&no_focusing::inv_beam_Radius;
		Gouy_phase=&no_focusing::Gouy_phase;
		beam_waist_sq2=&no_focusing::beam_waist_sq2;
		inv_beam_Radius2=&no_focusing::inv_beam_Radius2;
		Gouy_phase2=&no_focusing::Gouy_phase2;
	}
	if(if_square_pulse){
		a_laser_t_env=&square_pulse::a_laser_t_env;}
	else{
		a_laser_t_env=&asymmetric_gaussian_pulse::a_laser_t_env;
	}

/*
	if(if_PLAD){
		Lambda=&PLAD::Lambda;}
	else{
		if(if_QPLAD){
			Lambda=&QPLAD::Lambda;}
		else{
			Lambda=&LAD::Lambda;
		}
	}
	
 	New version AGRT Aug 2017		
	*/
	Lambda=&norad::Lambda;
	if(radiation_force_on) {
		Lambda=&LAD::Lambda;
	if(if_PLAD){
		Lambda=&PLAD::Lambda;}
	if(if_QPLAD){
		Lambda=&QPLAD::Lambda;}
	}
	
	
	if(if_uniform_bunch){
		initialize_particles=&uniform_bunch::initialize_particles;}
	else{
		initialize_particles=&gaussian_bunch::initialize_particles;
	}
	
	if(if_bubble) {
		a_wakefield = &bubblewake::a_wakefield;	
	} else {
		a_wakefield = &linearwake::a_wakefield;
	}
	if (bub_width>0.0L) {
		ramp = &withramp::ramp;
	} else {
		ramp = &withnoramp::ramp;	
	}
	// RUN THE SIM -------------------------------------------------------------
		
        fourvec *x0=NULL;
        fourvec *x1=NULL;
        fourvec *x2=NULL;
        fourvec *v0=NULL;
        fourvec *v1=NULL;
       
    // Run program
	rdtx_run(x0, x1,x2,v0,v1);

	// Free up some memory
	delete [] sincwdt2;
	delete [] coswdt2;
	delete [] particle_shape_k;
	delete [] omega;
	delete [] int_x;
	delete [] int_v;
	delete [] int_S;
	delete [] int_x_m;
	delete [] int_v_m;
	delete [] int_S_m;
	delete [] damping;
	delete [] omega_r1;
	delete [] omega_r2;
	
	for (int i=0;i<N_omega_bins;++i){
		delete [] SpecInt_para[i];
		delete [] SpecInt_perp[i];
		delete [] Spe_In_Im_x[i];
		delete [] Spe_In_Re_x[i];
		delete [] Spe_In_Im_y[i];
		delete [] Spe_In_Re_y[i];
		delete [] Spe_In_Im_z[i];
		delete [] Spe_In_Re_z[i];
	}
	delete [] SpecInt_para;
	delete [] SpecInt_perp;
	delete [] Spe_In_Im_x;
	delete [] Spe_In_Re_x;
	delete [] Spe_In_Im_y;
	delete [] Spe_In_Re_y;
	delete [] Spe_In_Im_z;
	delete [] Spe_In_Re_z;
	
	for (int i=0;i<4;++i){
		delete [] par_x[i];
		delete [] par_v[i];
	}
	delete [] par_x;
	delete [] par_v;
	
	delete Kspec;
	
	
   delete  laserpulse1::K0 ;
   delete  laserpulse1::dir0;
   delete  laserpulse1::rproj;
   delete  laserpulse1::pol_1;
   delete  laserpulse1::pol_2 ;
   delete  laserpulse1::Kphase;
   delete  laserpulse1::Kgroup;


   delete laserpulse2::K0 ;
   delete laserpulse2::dir0;
   delete laserpulse2::rproj;
   delete laserpulse2::pol_1;
   delete laserpulse2::pol_2;
	delete laserpulse2::Kphase;
   delete laserpulse2::Kgroup;
	
	
	SpecInt_para=NULL;//[N_omega_bins][N_theta_bins];
	SpecInt_perp=NULL;//[N_omega_bins][N_theta_bins];
	Spe_In_Im_x=NULL;//[N_omega_bins][N_theta_bins];
	Spe_In_Re_x=NULL;//[N_omega_bins][N_theta_bins];
	Spe_In_Im_y=NULL;//[N_omega_bins][N_theta_bins];
	Spe_In_Re_y=NULL;//[N_omega_bins][N_theta_bins];
	Spe_In_Im_z=NULL;//[N_omega_bins][N_theta_bins];
	Spe_In_Re_z=NULL;//[N_omega_bins][N_theta_bins];

	sincwdt2=NULL;//[N_omega_bins];
	coswdt2=NULL;//[N_omega_bins];

	particle_shape_k=NULL;//[N_omega_bins];

	omega=NULL;//[N_omega_bins];
 
	int_x=NULL;//,[Nmax+1]
	int_v=NULL;//,[Nmax+1]
	int_S=NULL;//,[Nmax+1]
	int_x_m=NULL;//,[Nmax+1]
	int_v_m=NULL;//;[Nmax+1]
	int_S_m=NULL;//;[Nmax+1]
	
	damping=NULL;//[Nmax+1];
	omega_r1=NULL;//[N_omega_bins];
	omega_r2=NULL;//[N_omega_bins];

	par_x=NULL;//[4][Npar];
	par_v=NULL;//[4][Npar];
	
	Kspec=NULL;
	

	laserpulse1::K0 =NULL;
	laserpulse1::dir0 = NULL;
	laserpulse1::rproj =  NULL;
    laserpulse1::pol_1 = NULL;
    laserpulse1::pol_2 = NULL;
    laserpulse1::Kphase = NULL;
    laserpulse1::Kgroup = NULL;

   laserpulse2::K0 = NULL;
   laserpulse2::dir0 = NULL;
   laserpulse2::rproj = NULL;
   laserpulse2::pol_1 = NULL;
   laserpulse2::pol_2 = NULL;
	laserpulse2::Kphase = NULL;
   laserpulse2::Kgroup = NULL;

	return 0;
}

#endif



