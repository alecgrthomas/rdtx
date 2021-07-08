#ifndef IMPORT_C_
#define IMPORT_C_

#include "headers.h"
#include "fresnel.h"

void ferr(const char * file)
{
  std::cout<< "ERROR: "<<file<<" not found"<<std::endl;
}
long double sinc(long double x);

int import_data(fourvec *&x0,fourvec *&x1,fourvec *&x2,fourvec *&v0,fourvec *&v1,long double &dt,std::string &name) 
{
  const int Nfields=16; // number of data columns
  
  char nullchar[20];
  // long double dt=1.0L;
  int ii;
  long double tempdouble=0.0L;
  std::ifstream infile;
  infile.open(name.c_str());
  if (!infile) ferr(name.c_str());
  
  infile>>nullchar; // "Nt"
  infile>>Nmax;

  infile>>nullchar; // "dt"
  infile>>dt;

  x0=new fourvec[Nmax];
  x1=new fourvec[Nmax];
  x2=new fourvec[Nmax];
  v0=new fourvec[Nmax];
  v1=new fourvec[Nmax];

  for (ii=0;ii<Nfields; ++ii)
    {
      infile>>nullchar; // "dumps file titles"
    }
   std::string temp;
  for (ii=0;ii<Nmax; ++ii) // cycle over rows
    {
      std::cout<<"Timestep "<<ii<<'\n'; 

      
      // one time component
      infile>>tempdouble; 
      x0[ii].sett(tempdouble);
      x1[ii].sett(dt);
      // no t_2
      
      infile>>tempdouble;
      x0[ii].setx(tempdouble);
      infile>>tempdouble; 
      x1[ii].setx(tempdouble);
      infile>>tempdouble; 
      x2[ii].setx(tempdouble);
      
      infile>>tempdouble; 
      x0[ii].sety(tempdouble);
      infile>>tempdouble; 
      x1[ii].sety(tempdouble);
      infile>>tempdouble; 
      x2[ii].sety(tempdouble);
      
      infile>>tempdouble; 
      x0[ii].setz(tempdouble);
      infile>>tempdouble; 
      x1[ii].setz(tempdouble);
      infile>>tempdouble; 
      x2[ii].setz(tempdouble);
      
      infile>>tempdouble; 
      v0[ii].setx(tempdouble);
      infile>>tempdouble; 
      v1[ii].setx(tempdouble);
      
      infile>>tempdouble; 
      v0[ii].sety(tempdouble);
      infile>>tempdouble; 
      v1[ii].sety(tempdouble);
    
      infile>>tempdouble; 
      v0[ii].setz(tempdouble);
      infile>>tempdouble; 
      v1[ii].setz(tempdouble);
    }
  
  infile.close();
  
  return 0;
}


int SpectralCalculationImport(int theta, int parnum,   fourvec *&x0_in,
                              fourvec *&x1_in,fourvec *&x2_in,fourvec *&v0_in,
                              fourvec *&v1_in,double dt)
{
 
int messagecounter=0; // to limit no. messages
  const int messagenow=Nmax/Nmessages;
  /*
Follows spectral integrals in corresponding paper but over 2dtau
   */
   long double inv_dtau=1.0/dt;
   long double inv_2dtau_sq = 2.0/dt/dt;
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
     // if (fabs(r2)>r2mean) {r2mean=fabs(r2);}
      if (fabs(r3)>r3mean) {
	r3mean=fabs(r3);
	r2mean=fabs(r2);
	r1whenr3max=fabs(r1);
	}
    }
  // outout.close();
  std::cout<<"\n\nMaximum size of {r1*omegamax,r2*omegamax,r3*omegamax} when r3 max = {"<<r1whenr3max*dtau0*omega_max<<','<<r2mean*dtau0*dtau0*omega_max<<','<<r3mean*dtau0*dtau0*dtau0*omega_max<<"}\n\n";	
std::cout<<"\n\nMaximum size of {r1*omegamax,r2*omegamax} = {"<<r1max*dtau0*omega_max<<','<<r2max*dtau0*dtau0*omega_max<<"}\n\n";	

/*
// Throw error messages
if(order_of_spectrometer==2) {
if (r3mean*dtau0*dtau0*dtau0*omega_max>1.0L) {
	wxMessageDialog *dial = new wxMessageDialog(NULL,wxT("Error in spectrum: 3rd order coefficient x omegamax too big (>1). Continue?"), wxT("Spectrum Error"), wxYES_NO | wxICON_ERROR);
	if(dial->ShowModal()==wxID_NO) return 1;
	dial->Destroy();
}
}
if(order_of_spectrometer==1) {
	if (r2max*dtau0*dtau0*omega_max>1.0L) {
		wxMessageDialog *dial = new wxMessageDialog(NULL,wxT("Error in 1st order spectrum: 2rd order coefficient x omegamax too big (>1). Continue?"), wxT("Spectrum Error"), wxYES_NO | wxICON_ERROR);
		if(dial->ShowModal()==wxID_NO) return 1;
		dial->Destroy();
	}
}
if(order_of_spectrometer==0) {
	if (r1max*dtau0*omega_max>1.0L) {
		wxMessageDialog *dial = new wxMessageDialog(NULL,wxT("Error in 0th order spectrum: 1rd order coefficient x omegamax too big (>1). Continue?"), wxT("Spectrum Error"), wxYES_NO | wxICON_ERROR);
		if(dial->ShowModal()==wxID_NO) return 1;
		dial->Destroy();
	}
}
*/
  // Find limit for changeover
if(order_of_spectrometer==2) {
  for (int ii=0;ii<Nmax;++ii)
    {
      iiplus=ii+1;
      r2=(*Kspec)*(x2_in[ii]);
      for (int i=0;i<N_omega_bins;++i) 
	if (fabs(omega[i]*r2)>model_threshold*inv_2dtau_sq) 
	  {omega_changeover=i; goto checkdone;}
    }
 }
 checkdone:

  //omega_changeover=N_omega_bins; 
  for (int ii=0;ii<Nmax;++ii) {

    if (timemax) if (x0_in[ii].gett()>timemax) break;
    // These sort out temproal boundary conditions
    iiplus=ii+1;

    // --------------------------------------------------------------------
    //Don't subtract phase!!! or radiation will be artificially COHERENT!!!
    // That is to say don't comment out the bit at the end of the below line
    // phase is basically r0 in the terminolgy I use below
    phase=(*Kspec)*(x0_in[ii]);//-(*Kspec)*(int_x_m[0]);
    v0=v0_in[ii];
        
     /*
    Use quadratic INTERPOLATION to find coefficients!
       (Using interpolation)
     */
    v1=v1_in[ii];
    r1=(*Kspec)*x1_in[ii];
    r2=(*Kspec)*x2_in[ii];   
   
  // --------------------------------------------------------------------
 

    v0x=v0.getx(); v0y=v0.gety(); v0z=v0.getz();
    v1x=v1.getx(); v1y=v1.gety(); v1z=v1.getz();

    sign = 2.0*(0.5-(r2<0.0));   
    for (int i=0;i<N_omega_bins;++i) {
      omega_r1[i]=omega[i]*r1*sign;
      omega_r2[i]=fabs(omega[i]*r2);
    }
    for (int i=0;i<N_omega_bins;++i) {	  
      sincwdt2[i]=sinc(omega_r1[i]*dt*0.5);
      coswdt2[i]=cos(omega_r1[i]*dt*0.5);
    }

     /*
      Taylor expansion of quadratic expansion
    */

    for (int i=0;i<omega_changeover;++i) {
      phasew = omega[i]*phase;
      cos_phasew=cos(phasew)*damping[ii]*dt;
      sin_phasew=sin(phasew)*damping[ii]*dt;      
      // These integrals for 0th 1st and 2nd order terms
      I0=sincwdt2[i];
      if (omega_r1[i])
	I1=(coswdt2[i]-sincwdt2[i])/omega_r1[i];
      else I1=0.0;
      if (omega_r1[i])
	I2=2.0*I1/omega_r1[i]+dt*dt*0.25*I0;
      else I2=dt*dt*0.25*I0;
  
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
      Theta_p=(omega_r1[i]+omega_r2[i]*dt)/sqrt(2.0*PI*omega_r2[i]);
      Theta_m=(omega_r1[i]-omega_r2[i]*dt)/sqrt(2.0*PI*omega_r2[i]);
      Phi_p = dt*dt*0.25*omega_r2[i]+dt*0.5*omega_r1[i];
      Phi_m = dt*dt*0.25*omega_r2[i]-dt*0.5*omega_r1[i];
      
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
	cos_phasew=cos(phasew)*damping[ii]*dt;
	sin_phasew=sin(phasew)*damping[ii]*dt;   
      // These integrals for 0th 1st and 2nd order terms
      I0=sincwdt2[i];
      if (omega_r1[i])
	I1=(coswdt2[i]-sincwdt2[i])/omega_r1[i];
      else I1=0.0;
      if (omega_r1[i])
	I2=2.0*I1/omega_r1[i]+dt*dt*0.25*I0;
      else I2=dt*dt*0.25*I0;

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
	  <<theta_scatter<<", particle "<<parnum<<" of "<<Npar;
      std::cout<<"\ntime (this step/total est.) = ("
	  <<GetstrTime(end_time_step-start_time_step)<<"/"
	  <<GetstrTime((end_time_step-start_time_step)*Nmessages) <<")\n\n"; 
     

  
      start_time_step=end_time_step;
  }
  }
  // outout.close();
  return 0;
}

#endif

