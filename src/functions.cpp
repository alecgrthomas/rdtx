
#ifndef FUNCS_C_
#define FUNCS_C_
#include "headers.h"
/*#include <iostream>
#include <fstream>
#include <sstream>
#include "fourvec.h"
#include "fourtens.h"
#include "constants.h"
*/
#include "radiation_force.h"




//------------------------------------------------------




// Right hand side of equation dx/dtau = f1(tau,x,v,vdot,vdotdot)
fourvec fx(long double &tau,fourvec& x,fourvec& v,fourvec& vdot,fourvec& vdotdot)
{
  fourvec ans;
  ans=v;
  return ans;
}
// Right hand side of equation dv/dtau = f2(tau,x,v,vdot,vdotdot)
fourvec fv(long double &tau,fourvec& x,fourvec& v,fourvec& vdot,fourvec& vdotdot)
{
  fourvec ans(0.0L,0.0L,0.0L,0.0L);
  ans=F(x)*v; // Lorentz force
  ans=ans*(-q_over_m); //charge to mass ratio
  // Now add radiation reaction force
  ans=ans+Lambda(x,v,vdot,vdotdot);
 
  return ans;
}
// Poincare invariant (E^2-B^2)/2 - not implemented yet!
bool Poincare(fourvec& x) {

	bool ans = false;
 	int zero=0;
 
	fourtens FS = F(x);
	long double FPoincare=0.0L;
	if (poincare_on) {
		FPoincare = (FS*FS)*0.5L; // not /4, so that condition is thresh ~ E
	} else {
		for (int ii=1;ii<4;++ii){
			FPoincare +=  FS(zero,ii)*FS(zero,ii);	//E components
		}
	}
	if (sqrt(FPoincare)>ionization_threshold) ans = true;
	return ans;
}

long double par_traj(fourvec &v0,fourvec &x0,fourvec &S0,int par_label,std::ofstream &outfile_summary)
{
 
  fourvec S=S0; // four Spin vector
  fourvec v=v0; //four velocity equiv to energy/momentum
  fourvec x=x0; // four vector - (ct,x)
  fourvec a_las;

  fourvec vdot(0.0,0.0,0.0,0.0); //four acceleration
  fourvec vdotdot(0.0,0.0,0.0,0.0); // differential of acceleration

  long double tau=0.0; // tau is proper time  

  
  fourvec kx[4]; //For Runga Kutta step - in x
  fourvec kv[4]; //For Runga Kutta step - in velocity
  fourvec kS[4]; //For Runga Kutta step - in spin
  
  fourvec kvm(0.0,0.0,0.0,0.0); // midpoint velocity
  fourvec kxm(0.0,0.0,0.0,0.0);
  fourvec kSm(0.0,0.0,0.0,0.0);
  
  
  //fourvec kvdot[4]; //For Runga Kutta step - in x
  //fourvec kvdotdot[4]; //For Runga Kutta step - in velocity
  
  fourvec xtemp=x; //temporary stores
  fourvec Stemp=S;
  fourvec vtemp=v;
  fourvec vdottemp=vdot;
  fourvec vdotdottemp=vdotdot;
  
  const long double dtau=dtau0;
  long double tautemp = tau;
  for (int ii=0;ii<4;++ii)
    kx[ii]=fourvec(0.0,0.0,0.0,0.0);
  for (int ii=0;ii<4;++ii)
    kv[ii]=fourvec(0.0,0.0,0.0,0.0);
   for (int ii=0;ii<4;++ii)
    kS[ii]=fourvec(0.0,0.0,0.0,0.0);
 
 //main loop
  //*********************************************
  int dumpcounter=0; // to limit no. dumps
  const int dumpnow=Nmax/(Ndumps-1);

 
  long double next_time_step=0.0;
  const long double time_step_size = time_max/(Ndumps-1);
  int messagecounter=0; // to limit no. messages
  const int messagenow=Nmax/Nmessages;
  std::ofstream outfile;

  // int to string
  std::ostringstream buffer;
  buffer << par_label;
  std::string file_string = filename+buffer.str()+".dat";
  
  int Ndumps_mod=Ndumps;
  if (Ndumps>Nmax) Ndumps_mod=Nmax;

  // Set up data write 
  a_laser(x,a_las);
  if (if_dump_par_dat){  
    outfile.open(file_string.c_str());
    write_parhead(outfile,par_label,Ndumps_mod,v_wake);
    write_data(outfile,x, v,S, 0.0,a_las);
  }
  
 


  long double Erad=0.0,Pmax=0.0;

 
  std::ofstream tempout;

  tempout.open((Run_directory+"/a_temp_out.dat").c_str());
  long double start_time_step=clock(),end_time_step=clock();
  //***********************************************************
  // MAIN LOOP
  //***********************************************************
  for (int ii=0;ii<=Nmax;++ii)
    {

      //***********************************************************
      // START OF CALCULATION
      //***********************************************************
      // Store velocity and position
      int_v[ii]=v;
      int_x[ii]=x;
      int_S[ii]=S;
      
      //dtau=dtau0*v.gett();
      // First work out Runga Kutta coefficients
      //-----------------------------------------
      kx[0] = fx(tau,x,v,vdot,vdotdot);
      kv[0] = fv(tau,x,v,vdot,vdotdot);
      kS[0] = BMT(tau,x,v,S);

      xtemp=x+kx[0]*dtau*0.5;
      vtemp=v+kv[0]*dtau*0.5;
      Stemp=S+kS[0]*dtau*0.5;
      
      vdottemp = vdot+vdotdot*0.5;
      
      // don't have good approx for vdotdot - use kv0 (dv/dt)      
      vdotdottemp = (kv[0]-vdot)*(1.0/dtau);

      tautemp=tau+dtau*0.5;
      kx[1] = fx(tautemp,xtemp,vtemp,vdottemp,vdotdottemp);
      kv[1] = fv(tautemp,xtemp,vtemp,vdottemp,vdotdottemp);
	  kS[1] = BMT(tautemp,xtemp,vtemp,Stemp);
	  
      xtemp=x+kx[1]*dtau*0.5;
      vtemp=v+kv[1]*dtau*0.5;
      Stemp=S+kS[1]*dtau*0.5;
            
      vdottemp = vdot+vdotdot*0.5;
      // don't have good approx for vdotdot - use kv1 (dv/dt)      
      vdotdottemp = (kv[1]-vdot)*(1.0/dtau);
      kx[2] = fx(tautemp,xtemp,vtemp,vdottemp,vdotdottemp);
      kv[2] = fv(tautemp,xtemp,vtemp,vdottemp,vdotdottemp);
  	  kS[2] = BMT(tautemp,xtemp,vtemp,Stemp);
	
      xtemp=x+kx[2]*dtau;
      vtemp=v+kv[2]*dtau;
      Stemp=S+kS[2]*dtau;
            
      vdottemp = vdot+vdotdot*0.5;
      // don't have good approx for vdotdot - use kv2 (dv/dt)      
      vdotdottemp = (kv[2]-vdot)*(1.0/dtau);
       tautemp=tau+dtau;
      kx[3] = fx(tautemp,xtemp,vtemp,vdottemp,vdotdottemp);
      kv[3] = fv(tautemp,xtemp,vtemp,vdottemp,vdotdottemp); 
 	  kS[3] = BMT(tautemp,xtemp,vtemp,Stemp);
	
      vdottemp = vdot+vdotdot*0.5;
      // don't have good approx for vdotdot - use kv2 (dv/dt)      
      vdotdottemp = (kv[3]-vdot)*(1.0/dtau);

             //-----------------------------------------


      // intermediate! NEW!
      //------------------------------------------
      //update interpolation variables for use in spectral integration
      kvm = (kv[0]+kv[1])*0.5*dtau*0.5;//1st 0.5 for averaging, 2nd for 1/2 step
      kxm = (kx[0]+kx[1])*0.5*dtau*0.5;      
	  kSm = (kS[0]+kS[1])*0.5*dtau*0.5;      

      int_v_m[ii]=int_v[ii]+kvm;
      int_x_m[ii]=int_x[ii]+kxm;
      int_S_m[ii]=int_S[ii]+kSm;
      //------------------------------------------
      
      vdot=(kv[0]+(kv[1]+kv[2])*2.0+kv[3])*oneover6;
      
      // Now update variables x,v and S
   // Check "ionization" state
	if (Poincare(x)) {
		particleisfree = 1.0L;}
		
	if(particleisfree) {
      x=x+(kx[0]+(kx[1]+kx[2])*2.0+kx[3])*oneover6*dtau;
      v=v+vdot*dtau;
       S=S+(kS[0]+(kS[1]+kS[2])*2.0+kS[3])*oneover6*dtau;
	} else {
	x.sett(v.gett()*ii*dtau);	
	}

      //now update time
      tau+=dtau;


      // Now update radiated energy
      Erad-=Lambda(x,v,vdot,vdotdot).gett();//(v*rad_const*vdot.square()*dtau).gett();

      // Now for dumping data
      ++dumpcounter;
      
      if (!dump_lab_time) {
	if (dumpcounter>dumpnow) {
          a_laser(x,a_las);	  
	  dumpcounter=1;
	  write_data(outfile,x, v,S, Erad/x.gett(),a_las);
	  
	  
	  }
      }
      else if (x.gett()>next_time_step) {
	if (x.gett()<=time_max) {
          a_laser(x,a_las);
	  write_data(outfile,x, v,S, Erad/x.gett(),a_las);
	  next_time_step+= time_step_size;
	} 
      }
      
      if (Erad/x.gett()>Pmax&&x.gett()) Pmax=Erad*dtau0/x.gett();
      ++messagecounter;
      if (messagecounter>messagenow)
	{
	  messagecounter=1;
	  end_time_step=clock();
	  std::cout<<"Current timestep "<<ii<<" of "<<Nmax<<"(tau="<<tau<<", current particle "<<par_label<<")"<<"\ntime (this step/total est.) = ("<<GetstrTime(end_time_step-start_time_step)<<"/"<<GetstrTime((end_time_step-start_time_step)*Nmessages) <<"), P_rad_ave="<<Erad/x.gett()<<"\n\n"; 
	  start_time_step=end_time_step;
  
	}
    } // END OF MAIN LOOP
  tempout.close();
  if (dump_lab_time) { outfile<<"1.2345";}
  //-------------------------------------------------------
  /*
    In case you think this is wrong Alec, you already checked for synchrotron radiation that the radiation force 
    term exactly compensates for the radiated power! (April 09)
   */
  std::cout<<"\nTotal radiated energy for particle "<< par_label<< " = "<<Erad*dtau0*511.0<<" keV\n";
  std::cout<<"Final energy for particle "<< par_label<< " = "<<v.gett()*0.511<<" MeV\n\n";
  totalErad+=Erad*dtau0;
  long double Plas=0.5*a0_pulse1*a0_pulse1/.64/8.6e-10/8.6e-10*PI*(laserpulse1::w0*laserpulse1::w0/4.0*c*c/omega_0/omega_0);
  long double nelec=1e19*PI*(laserpulse1::w0*laserpulse1::w0/4*c*c/omega_0/omega_0)*laserpulse1::t0/omega_0*c*1e6;
  std::cout <<"P_rad/P_las (laser 1) for t0*w0^2*1e19 electrons = "<<Pmax*nelec*5.11e6*e*omega_0/Plas<<"\n\n--------------------------------------------\n";
 
    if (if_dump_par_dat){
    outfile_summary<<par_label<<'\t';
    outfile_summary<<x0; outfile_summary<<'\t';
    outfile_summary<<v0; outfile_summary<<'\t';
    outfile_summary<<x; outfile_summary<<'\t';
    outfile_summary<<v; outfile_summary<<'\t';
    outfile_summary<<Erad*511.<<'\t'<<Pmax*nelec*5.11e6*e*omega_0/Plas<<'\n';
    outfile.close();
  }
  // if (!par_label&&!if_dump_fields) dump_fields(1,tau);
  return x.gett();
}


int set_initial()
{
 
  // Create metric tensor
  metrictens.Minkowski();

  // Set K vectors to group/phase velocity
  laserpulse1::Kphase->sett(laserpulse1::v_phase*laserpulse1::Kphase->gett());
  laserpulse1::Kgroup->sett(laserpulse1::v_group*laserpulse1::Kgroup->gett());
  laserpulse2::Kphase->sett(laserpulse2::v_phase*laserpulse2::Kphase->gett());
  laserpulse2::Kgroup->sett(laserpulse2::v_group*laserpulse2::Kgroup->gett());

  srand ( time(NULL) ); 
  //Initialize spectral integrals
  //----------------------------------------
  ResetBins();  
   
  if (!if_omega_log)
    {
      omega[0]=omega_min;
      for (int i=1;i<N_omega_bins;++i)
	omega[i]=omega[i-1]+domega;
    }
  else
    {
      if (omega_min<=0.0) 
	{
	  std::cerr<<"ERROR: Logarithic omega scale requested but omega_min = "<<omega_min<<'\n';
 	  exit(0);
	}
      /*
	For a logaritmic scale, domega appears as:
	omega = omega_min *exp( n domega );
	This is equally spaced on a logarithmic plot
       */
      domega =  1.0*log(omega_max/omega_min)/(N_omega_bins-1);
      omega[0]=omega_min;
      for (int i=1;i<N_omega_bins;++i)
	omega[i]=omega[0]*exp(i*domega);
    }
  for (int i=0;i<N_omega_bins;++i)
    {	  
      sincwdt2[i]=1.0;
      coswdt2[i]=1.0;
    }
  defineshapefactor();
  return 0;
}
std::string GetstrTime(long double time)
{
  std::ostringstream outtime;
  int days,hours,mins;
  long double secs=time/CLOCKS_PER_SEC;
  days=int (secs)/86400;
  secs-=86400.0*days;
  hours=int (secs)/3600;
  secs-=3600.0*hours;
  mins=int (secs)/60;
  secs-=60.0*mins;
if (days==1) outtime<<days<<" day, ";
  if (hours==1) outtime<<hours<<" hour, ";
  if (mins==1) outtime<<mins<<" min, ";
  if (days>1) outtime<<days<<" days, ";
  if (hours>1) outtime<<hours<<" hours, ";
  if (mins>1) outtime<<mins<<" mins, ";
  outtime<<secs<<" s";
 return outtime.str();
}

#endif
