#ifndef SPECFUNCS_C_
#define SPECFUNCS_C_
#include <iostream>
#include <fstream>

#include "fresnel.h"
#include "constants.h"

#include "headers.h"

/*
  This file contains the spectral integration funcitons for RDT
 */

//------------------------------------------------------
// First a couple of utility functions....
// sinc function
long double sinc(long double x)
{
    if (x) return sin(x)/x;
    else return 1.0L;
}
/*
  In this function we define the shape factor for particles to determine coherency 
  properties. Assume a gaussian bunch with density so many times ambient density
  
  Shape factor S(omega/c)^2 is fourier transform of gaussian bunch structure, and so is
  also gaussian. But S^2 needs to tend towards N for incoherent radiation. So the combination
used is 
                  N^2*S^2 + (1-S^2)*N

 */
int defineshapefactor()
{
  // The shape factor constant is the normalization of n_c with respect to (omega_0/c)^3
  long double shape_width_square 
    = 0.5L*pow(sizeofbunch/densityofbunch*shape_factor_constant,2.0L/3.0L);
  long double shape_bunch_square=1.0L;
  if (if_macroparticle)
    {
      if (coherent_addition)
	for (int i=0;i<N_omega_bins;++i)
	  {
	    shape_bunch_square=exp(-omega[i]*omega[i]*0.5*shape_width_square);
	    particle_shape_k[i]=shape_bunch_square*sizeofbunch*sizeofbunch+(1.0-shape_bunch_square)*sizeofbunch;
	  }
      else for (int i=0;i<N_omega_bins;++i) particle_shape_k[i]=sizeofbunch;
    }
  else for (int i=0;i<N_omega_bins;++i) particle_shape_k[i]=1.0;
  return 0;
}

// This cleans all the bins used for integration
int ResetBins()
{
  for (int i=0;i<N_omega_bins;++i)
    { 
      for (int j=0;j<N_theta_bins;++j)
	{
	  Spe_In_Re_x[i][j] =0.0; 
	  Spe_In_Im_x[i][j] =0.0; 
	  Spe_In_Re_y[i][j] =0.0; 
	  Spe_In_Im_y[i][j] =0.0; 
	  Spe_In_Re_z[i][j] =0.0; 
	  Spe_In_Im_z[i][j] =0.0;  
	  
	}
    }
  return 0;
}
int ResetSpec()
{
  for (int i=0;i<N_omega_bins;++i)
    { 
            for (int j=0;j<N_theta_bins;++j)
	{
	  SpecInt_para[i][j]=0.0;
	  SpecInt_perp[i][j]=0.0;
	}
    }
  return 0;
}

/*
    To damp the significsnce of the endpoints
   */

int calculatedamping()
{
  // i0 starts at DW damping widths away so that damping -> 10^-7 at endpoint
  const long double DW=damping_n_sigmas;
  int i0=Nmax-DW*damping_width/dtau0;
  
  if (i0<0) {
	i0=0;}
  for (int i=0;i<Nmax;++i) {
	damping[i]=1.0L;}

  if (if_endpoint_damping)
    {
	if (damping_width>0.0L) {
	      for (int i=i0;i<Nmax;++i) {
		damping[i]=exp(-((double)i-i0)*((double)i-i0)
		       /damping_width/damping_width*dtau0*dtau0);
  		}
	} // end if damping_width

      i0=DW*damping_endpoint_ratio/dtau0;
      if (i0>Nmax-1) {
	i0=Nmax-1;}

	if (damping_endpoint_ratio>0.0L) {
      		for (int i=0;i<i0;++i) {
			damping[i]*=exp(-((double)i-i0)*((double)i-i0)
		      	/damping_endpoint_ratio/damping_endpoint_ratio*dtau0*dtau0);
		}
	}// end if
    }

  std::ofstream damping_info;
  damping_info.open((Run_directory+"/dampinginfo.dat").c_str());
  int dumpcounter=0;
  const int dumpnow=(Nmax)/Ndumps;
  damping_info<<Ndumps<<"\n\n";
  for (int i=0;i<=Nmax;++i) 
    {
      ++dumpcounter;
      if (dumpcounter>dumpnow)
	{	  
	  dumpcounter=1;
	  damping_info<<dtau0*i<<'\t'<<damping[i]<<'\n';
	}
    }
  damping_info.close();
  return 0;
}

int SpectralCalculation(int theta, int parnum)
{


int messagecounter=0; // to limit no. messages
  const int messagenow=Nmax/Nmessages;
  /*
Follows spectral integrals in corresponding paper but over 2dtau
   */
  const long double inv_dtau=1.0/dtau0;
  const long double inv_2dtau_sq = 2.0/dtau0/dtau0;
  long  double cos_phasew=0.0,sin_phasew=0.0;
  fourvec v0; 
 // First order approximation
  long double r1=0.0;
  fourvec v1;
  // Second order of approximation
  long double r2=0.0;
  fourvec v2;

  long  double v0x,v0y,v0z,v1x,v1y,v1z,sC2,sC3,sC4;
  long  double phase=(*Kspec)*int_x[0];
  long  double sign,sin_sC2,cos_sC2; 
  long  double I0,I1,I2;
  int iiplus=0;
  std::complex <long  double> DFres;
  std::complex<long  double> Fres_p, Fres_m;
  // The spectral integral has real and imaginary parts
  std::complex<long  double> Int_x,Int_y,Int_z;
  long  double Theta_p, Phi_p, sin_phi_p, cos_phi_p;
  long  double Theta_m, Phi_m, sin_phi_m, cos_phi_m;
  long  double start_time_step=clock(),end_time_step=clock();
  int omega_changeover=N_omega_bins; // change from Taylor expansion to Fresnel integral

   
   // CHECK RATIO OF R1 TO R2!!!!
  long double r3=0.0L, r2mean=0.0L,r3mean=0.0L,r1whenr3max=0.0L,r2max=0.0L,r1max=0.0L;

  for (int ii=1;ii<Nmax-1;++ii)
    {
      iiplus=ii+1;
      // work out 3rd order differential
      r1=(*Kspec)*(int_x[iiplus]-int_x[ii])*inv_dtau;
      r2=(*Kspec)*(int_x[iiplus]+int_x[ii]-int_x_m[ii]*2.0L)*inv_2dtau_sq;
      r3=(*Kspec)*(int_x[iiplus+1]-int_x_m[iiplus]*2.L+int_x_m[ii]*2.L-int_x[ii])*inv_2dtau_sq*inv_dtau*0.5L;
      if (fabs(r1)>r1max) {r1max=fabs(r1);}    
      if (fabs(r2)>r2max) {r2max=fabs(r2);}
      if (fabs(r3)>r3mean) {
	r3mean=fabs(r3);
	r2mean=fabs(r2);
	r1whenr3max=fabs(r1);
	}
    }
  // outout.close();

/*outstrstream<<"\n\nMaximum size of {r1*omegamax,r2*omegamax,r3*omegamax} when r3 max = {"<<r2mean*dtau0*dtau0*dtau0/3.0<<','<<r3mean*dtau0*dtau0*dtau0*dtau0*0.25<<','<<r3mean*dtau0*dtau0*dtau0*dtau0*0.25*omega_max<<"}\n\n";	
   *control<<wxString(outstrstream.str().c_str(),wxConvUTF8);
   control->Update();	
  std::cout<<"\n\nMaximum size of {r1*omegamax,r2*omegamax,r3*omegamax} when r3 max = {"<<r2mean*dtau0*dtau0*dtau0/3.0<<','<<r3mean*dtau0*dtau0*dtau0*dtau0*0.25<<','<<r3mean*dtau0*dtau0*dtau0*dtau0*0.25*omega_max<<"}\n\n";
*/

  std::cout<<"\n\nMaximum size of {r1*omegamax,r2*omegamax,r3*omegamax} when r3 max = {"<<r1whenr3max*dtau0*omega_max<<','<<r2mean*dtau0*dtau0*omega_max<<','<<r3mean*dtau0*dtau0*dtau0*omega_max<<"}\n\n";	
std::cout<<"\n\nMaximum size of {r1*omegamax,r2*omegamax} = {"<<r1max*dtau0*omega_max<<','<<r2max*dtau0*dtau0*omega_max<<"}\n\n";	


// Throw error messages
if (spectrumerror==0) {
if(order_of_spectrometer==2) {
if (r3mean*dtau0*dtau0*dtau0*omega_max>1.0L) {
	std::cout<<"Error in spectrum: 3rd order coefficient x omegamax too big (>1)";
}
}
if(order_of_spectrometer==1) {
	if (r2max*dtau0*dtau0*omega_max>1.0L) {
		std::cout<<"Error in 1st order spectrum: 2rd order coefficient x omegamax too big (>1)";
	}
}
if(order_of_spectrometer==0) {
	if (r1max*dtau0*omega_max>1.0L) {
		std::cout<<"Error in 0th order spectrum: 1rd order coefficient x omegamax too big (>1).";
	}
}
spectrumerror=1;
} // end error messages

  // Find limit for changeover
if(order_of_spectrometer==2) {
  for (int ii=0;ii<Nmax;++ii)
    {
      iiplus=ii+1;
      r2=(*Kspec)*(int_x[iiplus]+int_x[ii]-int_x_m[ii]*2.0)*inv_2dtau_sq;
      for (int i=0;i<N_omega_bins;++i) 
	if (fabs(omega[i]*r2)>model_threshold*inv_2dtau_sq) 
	  {omega_changeover=i; goto checkdone;} //only to exit loop
    }
}
checkdone:

 // omega_changeover=N_omega_bins; 
  for (int ii=0;ii<Nmax;++ii) {

    if (timemax) if (int_x[ii].gett()>timemax) break;
    // These sort out temproal boundary conditions
    iiplus=ii+1;

    // --------------------------------------------------------------------
    //Don't subtract phase!!! or radiation will be artificially COHERENT!!!
    // That is to say don't comment out the bit at the end of the below line
    // phase is basically r0 in the terminolgy I use below
    phase=(*Kspec)*(int_x_m[ii]);//-(*Kspec)*(int_x_m[0]);
    v0=int_v_m[ii];
        
     /*
    Use quadratic INTERPOLATION to find coefficients!
       (Using interpolation)
     */
    v1=(int_v[iiplus]-int_v[ii])*inv_dtau;
    r1=(*Kspec)*(int_x[iiplus]-int_x[ii])*inv_dtau;
    r2=(*Kspec)*(int_x[iiplus]+int_x[ii]-int_x_m[ii]*2.0)*inv_2dtau_sq;   
   
  // --------------------------------------------------------------------
 

    v0x=v0.getx(); v0y=v0.gety(); v0z=v0.getz();
    v1x=v1.getx(); v1y=v1.gety(); v1z=v1.getz();

    sign = 2.0*(0.5-(r2<0.0));   
    for (int i=0;i<N_omega_bins;++i) {
      omega_r1[i]=omega[i]*r1*sign;
      omega_r2[i]=fabs(omega[i]*r2);
    }
    for (int i=0;i<N_omega_bins;++i) {	  
      sincwdt2[i]=sinc(omega_r1[i]*dtau0*0.5);
      coswdt2[i]=cos(omega_r1[i]*dtau0*0.5);
    }

     /*
      Taylor expansion of quadratic expansion
    */

    for (int i=0;i<omega_changeover;++i) {
      phasew = omega[i]*phase;
      cos_phasew=cos(phasew)*damping[ii]*dtau0;
      sin_phasew=sin(phasew)*damping[ii]*dtau0;      
      // These integrals for 0th 1st and 2nd order terms
      I0=sincwdt2[i];
      if (omega_r1[i])
	I1=(coswdt2[i]-sincwdt2[i])/omega_r1[i];
      else I1=0.0;
      if (omega_r1[i])
	I2=2.0*I1/omega_r1[i]+dtau0*dtau0*0.25*I0;
      else I2=dtau0*dtau0*0.25*I0;
  
      if (order_of_spectrometer==1) {
        I2=I1=0.0; 
      }
      if (order_of_spectrometer==0) {
        I0=1.0;
      }
      //Remember: dtau = dt/gamma and vx=gamma*betax
      Spe_In_Re_x[i][theta] +=v0x*I0*cos_phasew
	-(-v1x*I1+sign*omega_r2[i]*I2*v0x)*sin_phasew; 
      //------------------
      Spe_In_Im_x[i][theta] +=v0x*I0*sin_phasew
	+(-v1x*I1+sign*omega_r2[i]*I2*v0x)*cos_phasew; 
      //------------------
      Spe_In_Re_y[i][theta] +=v0y*I0*cos_phasew
	-(-v1y*I1+sign*omega_r2[i]*I2*v0y)*sin_phasew; 
      //------------------
      Spe_In_Im_y[i][theta] +=v0y*I0*sin_phasew
	+(-v1y*I1+sign*omega_r2[i]*I2*v0y)*cos_phasew;
      //------------------
      Spe_In_Re_z[i][theta] +=v0z*I0*cos_phasew
	-(-v1z*I1+sign*omega_r2[i]*I2*v0z)*sin_phasew;
      //------------------
      Spe_In_Im_z[i][theta] +=v0z*I0*sin_phasew
	+(-v1z*I1+sign*omega_r2[i]*I2*v0z)*cos_phasew;
    }
    /*
      Full integral of quadratic expansion
    */
  
    for (int i=omega_changeover;i<N_omega_bins;++i) {
      if (fabs(omega[i]*r2)>model_threshold*inv_2dtau_sq){
      phasew = omega[i]*phase;
      cos_phasew=cos(phasew)*damping[ii]; // NB not * by dtau0
      sin_phasew=sin(phasew)*damping[ii];
      
      // Get Fresnel integrals
      Theta_p=(omega_r1[i]+omega_r2[i]*dtau0)/sqrt(2.0*PI*omega_r2[i]);
      Theta_m=(omega_r1[i]-omega_r2[i]*dtau0)/sqrt(2.0*PI*omega_r2[i]);
      Phi_p = dtau0sq*0.25*omega_r2[i]+dtau0*0.5*omega_r1[i];
      Phi_m = dtau0sq*0.25*omega_r2[i]-dtau0*0.5*omega_r1[i];
      
      fresnel(Theta_p,Fres_p);
      fresnel(Theta_m,Fres_m);
      DFres = Fres_p-Fres_m; // now Fres_p is Delta Fres
      
      sin_phi_p=sin(Phi_p);sin_phi_m=sin(Phi_m);
      cos_phi_p=cos(Phi_p);cos_phi_m=cos(Phi_m);
     
      sC2 = omega_r1[i]*omega_r1[i]/(4.0*omega_r2[i]);
      sin_sC2=sin(sC2);
      cos_sC2=cos(sC2);
      sC3 = sqrt(2.0*PI/omega_r2[i]);
      sC4=0.25/omega_r2[i];
      
      Int_x = 
	Spectral_Int_Exact(omega_r1[i],omega_r2[i],sign,v0x,v1x,DFres,
			   sin_phi_p,sin_phi_m,cos_phi_p,cos_phi_m,sC2,sC3,sC4,sin_sC2,cos_sC2);
      Int_y =
	Spectral_Int_Exact(omega_r1[i],omega_r2[i],sign,v0y,v1y,DFres,
			   sin_phi_p,sin_phi_m,cos_phi_p,cos_phi_m,sC2,sC3,sC4,sin_sC2,cos_sC2);
      Int_z =
	Spectral_Int_Exact(omega_r1[i],omega_r2[i],sign,v0z,v1z,DFres,
			   sin_phi_p,sin_phi_m,cos_phi_p,cos_phi_m,sC2,sC3,sC4,sin_sC2,cos_sC2);
      
      //Remember: dtau = dt/gamma and vx=gamma*betax
      Spe_In_Re_x[i][theta] += Int_x.real()*cos_phasew-Int_x.imag()*sin_phasew; 
      //------------------
      Spe_In_Im_x[i][theta] += Int_x.real()*sin_phasew+Int_x.imag()*cos_phasew; 
      //------------------
      Spe_In_Re_y[i][theta] += Int_y.real()*cos_phasew-Int_y.imag()*sin_phasew; 
      //------------------
      Spe_In_Im_y[i][theta] += Int_y.real()*sin_phasew+Int_y.imag()*cos_phasew;
      //------------------
      Spe_In_Re_z[i][theta] += Int_z.real()*cos_phasew-Int_z.imag()*sin_phasew;
      //------------------
      Spe_In_Im_z[i][theta] += Int_z.real()*sin_phasew+Int_z.imag()*cos_phasew;
      } else{
	phasew = omega[i]*phase;
	cos_phasew=cos(phasew)*damping[ii]*dtau0;
	sin_phasew=sin(phasew)*damping[ii]*dtau0;   
      // These integrals for 0th 1st and 2nd order terms
      I0=sincwdt2[i];
      if (omega_r1[i])
	I1=(coswdt2[i]-sincwdt2[i])/omega_r1[i];
      else I1=0.0;
      if (omega_r1[i])
	I2=2.0*I1/omega_r1[i]+dtau0*dtau0*0.25*I0;
      else I2=dtau0*dtau0*0.25*I0;

      //Remember: dtau = dt/gamma and vx=gamma*betax
      Spe_In_Re_x[i][theta] +=v0x*I0*cos_phasew
	-(-v1x*I1+sign*omega_r2[i]*I2*v0x)*sin_phasew; 
      //------------------
      Spe_In_Im_x[i][theta] +=v0x*I0*sin_phasew
	+(-v1x*I1+sign*omega_r2[i]*I2*v0x)*cos_phasew; 
      //------------------
      Spe_In_Re_y[i][theta] +=v0y*I0*cos_phasew
	-(-v1y*I1+sign*omega_r2[i]*I2*v0y)*sin_phasew; 
      //------------------
      Spe_In_Im_y[i][theta] +=v0y*I0*sin_phasew
	+(-v1y*I1+sign*omega_r2[i]*I2*v0y)*cos_phasew;
      //------------------
      Spe_In_Re_z[i][theta] +=v0z*I0*cos_phasew
	-(-v1z*I1+sign*omega_r2[i]*I2*v0z)*sin_phasew;
      //------------------
      Spe_In_Im_z[i][theta] +=v0z*I0*sin_phasew
	+(-v1z*I1+sign*omega_r2[i]*I2*v0z)*cos_phasew;
      }
      } 

  ++messagecounter;
  if (messagecounter>messagenow) {
 
    messagecounter=1;
    end_time_step=clock();
    std::cout<<"Current timestep "<<ii<<" of "<<Nmax<<", theta = "
	  <<theta_scatter*1000.0<<" mrad of "<<(theta_scatter_max-theta_scatter_min)*1000.0<<" mrad, particle "<<parnum+1<<" of "<<Npar;
      std::cout<<"\ntime (this step/total est.) = ("
	  <<GetstrTime(end_time_step-start_time_step)<<"/"
	  <<GetstrTime((end_time_step-start_time_step)*Nmessages) <<")\n\n"; 
     

	

      start_time_step=end_time_step;
  }
  }
  // outout.close();
  return 0;
}

int SpectralEndpoint(int theta, long double tau, fourvec v,fourvec x,long double sgn)
{
  /*
    In taking spectral integral, have integrated to T_0 - actually need
    to integrate to infinity - hence need to correct: Assume constant 
    velocity after T_0, add term:
   */
  if (if_endpoint_integral)
    {
     long  double cosphase=0.0L,sinphase=0.0L,vecphase=(*Kspec)*(v*tau),phasezero=(*Kspec)*x;
      for (int i=0;i<N_omega_bins;++i)
	{
	  phasew = omega[i]*vecphase;

	  if (phasew) {	  
	    sinphase=sin(phasew+phasezero*omega[i])*tau/phasew;
	    cosphase=cos(phasew+phasezero*omega[i])*tau/phasew;
	  } else {
	    sinphase=cosphase=0.0L;
	  }

	  //Remember: dtau = dt/gamma and vx=gamma*betax
	  
	  Spe_In_Re_x[i][theta] +=(v.getx())*sinphase*sgn; 
	  Spe_In_Im_x[i][theta] -=(v.getx())*cosphase*sgn; 
	  Spe_In_Re_y[i][theta] +=(v.gety())*sinphase*sgn;
	  Spe_In_Im_y[i][theta] -=(v.gety())*cosphase*sgn; 
	  Spe_In_Re_z[i][theta] +=(v.getz())*sinphase*sgn;
	  Spe_In_Im_z[i][theta] -=(v.getz())*cosphase*sgn; 
	  
    	}
    }
  return 0;
}

int Calc_Spec(int theta)
{
  costhetascatter = cos(theta_scatter);
  sinthetascatter = sin(theta_scatter);
  // NEW first calculate spectrum components
  fourvec Kspec0(1.0L,0.0L,sinthetascatter,costhetascatter);
  *Kspec=Kspec0;

  // Now put it all together and dump
  long double integration_temp=0.0; // to make the integration quicker
 
  for (int i=0;i<N_omega_bins;++i)
    {
      // Calculate perpendicular spectral intensity
      // line below changed to += AGRT 2010
      SpecInt_perp[i][theta]+=Spe_In_Im_x[i][theta]*Spe_In_Im_x[i][theta];
      SpecInt_perp[i][theta]+=Spe_In_Re_x[i][theta]*Spe_In_Re_x[i][theta];

      // Calculate parallel spectral intensity
      SpecInt_para[i][theta]+=Spe_In_Im_y[i][theta]*Spe_In_Im_y[i][theta]
	*costhetascatter*costhetascatter;
      SpecInt_para[i][theta]+=Spe_In_Re_y[i][theta]*Spe_In_Re_y[i][theta]
	*costhetascatter*costhetascatter;
      SpecInt_para[i][theta]+=Spe_In_Im_z[i][theta]*Spe_In_Im_z[i][theta]
	*sinthetascatter*sinthetascatter;
      SpecInt_para[i][theta]+=Spe_In_Re_z[i][theta]*Spe_In_Re_z[i][theta]
	*sinthetascatter*sinthetascatter;
      integration_temp=-2.0*(Spe_In_Re_z[i][theta]*Spe_In_Re_y[i][theta]
			     +Spe_In_Im_z[i][theta]*Spe_In_Im_y[i][theta])
	                     *sinthetascatter*costhetascatter;
      SpecInt_para[i][theta]+=integration_temp;   
    }

 return 0;
}
int Write_Spec(int parnum,int theta)
{
  // int to string
  std::ostringstream buffer;
  buffer << parnum;

  if (!separate_polarizations) {
    std::ofstream outfile_spectrum;
    std::string file_string = filename_spec+buffer.str()+".dat";
    if (parnum<0) file_string = filename_spectrum;
    if (N_theta_bins>1) 
      {
      buffer.str("");
      buffer<<theta;
      file_string=Run_directory+"/spectrum_theta_"+buffer.str()+".dat";
      }
    outfile_spectrum.open(file_string.c_str());
    outfile_spectrum<<"Final spectrum radiated at "<<theta_scatter<<" rad to +z direction\n";
     outfile_spectrum<<"E\n"<<N_omega_bins<<'\n';
    outfile_spectrum<<"omega\tI_spec\n";
    for (int i=0;i<N_omega_bins;++i)
      outfile_spectrum<<omega[i]<<'\t'<<SpecInt_perp[i][theta]+SpecInt_para[i][theta]<<'\n';
    outfile_spectrum.close();
  } else {
  buffer.str("");
      buffer<<theta;
      std::string file_string_perp=Run_directory+"/spectrum_perp_theta_"+buffer.str()+".dat";
      std::string file_string_para=Run_directory+"/spectrum_para_theta_"+buffer.str()+".dat";
      std::ofstream outfile_spectrum_perp, outfile_spectrum_para;
      outfile_spectrum_perp.open(file_string_perp.c_str());
      outfile_spectrum_para.open(file_string_para.c_str());
      {
	outfile_spectrum_perp<<"Final spectrum radiated at "<<theta_scatter<<" rad to +z direction\n";
		outfile_spectrum_perp<<"E\n"<<N_omega_bins<<'\n';
	outfile_spectrum_perp<<"omega\tI_spec\n";
      }
 {
	outfile_spectrum_para<<"Final spectrum radiated at "<<theta_scatter<<" rad to +z direction\n";
	outfile_spectrum_para<<"E\n"<<N_omega_bins<<'\n';
	outfile_spectrum_para<<"omega\tI_spec\n";
      }
 for (int i=0;i<N_omega_bins;++i){
      outfile_spectrum_para<<omega[i]<<'\t'<<SpecInt_para[i][theta]<<'\n';
      outfile_spectrum_perp<<omega[i]<<'\t'<<SpecInt_perp[i][theta]<<'\n';}
    outfile_spectrum_para.close();
    outfile_spectrum_perp.close();
  }
  return 0;
}

#endif
