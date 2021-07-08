
#ifndef FIELDPAR_H_GUARD_
#define FIELDPAR_H_GUARD_
#include "fourvec.h"
#include "constants.h"
#include "headers.h"

//*******************************************************
// NEW! june 2011 - fields updated to include higher order corrections.
fourvec a_laser_phase(fourvec(&x))
{
   long double rsquared = ((*laserpulse1::rproj)*x)*((*laserpulse1::rproj)*x);
   long double phase_total = (*laserpulse1::Kphase)*x-rsquared*0.5L*inv_beam_Radius(x)*laserpulse1::rww0;
  phase_total+=Gouy_phase(x);
  //circular polarization

  fourvec ans = ((*laserpulse1::pol_1)*sin(phase_total+laserpulse1::phase1)
		 +(*laserpulse1::pol_2)*sin(phase_total+laserpulse1::phase2));
 // new bit - second order correction, add scalar potential
  ans.sett((cos(phase_total+laserpulse1::phase1-Gouy_phase(x))
            +cos(phase_total+laserpulse1::phase2-Gouy_phase(x))));
   return ans;
}

long double Psi_function(fourvec(&x))
{
  long double rsquared = x.getx()*x.getx()+((*laserpulse1::rproj)*x) * ((*laserpulse1::rproj)*x);

  long double zpos=(*laserpulse1::Kgroup)*x;
  //Temporal envelope
  long double phase_in_exponent=zpos-laserpulse1::delay;

  long double ans = 1.0L;
  if (laserpulse1::t0)
    {
      ans = a_laser_t_env(x,phase_in_exponent,inv_t0squared);
    }

  //Transverse profile
  if (laserpulse1::w0>0.0L) 
    ans = ans*laserpulse1::w0/sqrt(beam_waist_sq(x))*exp(-rsquared/beam_waist_sq(x));
  return ans;
}
long double a_laser_envelope(fourvec(&x))
{
  long double ans = Psi_function(x);
  long double zpos=(*laserpulse1::Kgroup)*x;
  long double rsquared = x.getx()*x.getx()+((*laserpulse1::rproj)*x) * ((*laserpulse1::rproj)*x);
    
  // New bit AGRT June 2011 - 2nd order corrections to fields
 if (laserpulse1::w0>0.0L) {
  long double f_factor = laserpulse1::w0/sqrt(beam_waist_sq(x))*sin(Gouy_phase(x));
  long double correctionfactor=1.0L;
  correctionfactor+=0.5L* laserpulse1::focus_theta_1*laserpulse1::focus_theta_1*(1.0L-0.5L*(1.0L-zpos*zpos/laserpulse1::zR/laserpulse1::zR)/(1.0L+zpos*zpos/laserpulse1::zR/laserpulse1::zR)*rsquared/laserpulse1::w0/laserpulse1::w0*f_factor);
  ans=ans*correctionfactor;
 }
  return ans;
}
long double phi_laser_envelope(fourvec(&x))
{

  long double ans = a_laser_envelope(x);
 long double zpos=(*laserpulse1::Kgroup)*x;
 //long double rsquared = x.getx()*x.getx()+((*laserpulse1::rproj)*x) * ((*laserpulse1::rproj)*x);
    
  ans -= laserpulse1::focus_theta_1*laserpulse1::focus_theta_1*(1.0L-zpos*zpos/laserpulse1::zR/laserpulse1::zR)/(1.0L+zpos*zpos/laserpulse1::zR/laserpulse1::zR)* Psi_function(x);
  //
  ans *= laserpulse1::focus_theta_1*((*laserpulse1::pol_1)*x)/sqrt(beam_waist_sq(x));
  return ans;
}

fourvec a_laser_phase2(fourvec(&x))
{
   long double rsquared = ((*laserpulse2::rproj)*x)*((*laserpulse2::rproj)*x);
   long double phase_total = (*laserpulse2::Kphase)*x-rsquared*0.5L*inv_beam_Radius2(x)*laserpulse2::rww0;
  phase_total+=Gouy_phase2(x);
  //circular polarization

  fourvec ans = ((*laserpulse2::pol_1)*sin(phase_total+laserpulse2::phase1)
		 +(*laserpulse2::pol_2)*sin(phase_total+laserpulse2::phase2));
 // new bit - second order correction, add scalar potential
  ans.sett((cos(phase_total+laserpulse2::phase1-Gouy_phase2(x))
            +cos(phase_total+laserpulse2::phase2-Gouy_phase2(x))));
   return ans;
}

long double Psi_function2(fourvec(&x))
{
  long double rsquared = x.getx()*x.getx()+((*laserpulse2::rproj)*x) * ((*laserpulse2::rproj)*x);

  long double zpos=(*laserpulse2::Kgroup)*x;
  //Temporal envelope
  long double phase_in_exponent=zpos-laserpulse2::delay;

  long double ans = 1.0L;
  if (laserpulse2::t0)
    {
      ans = a_laser_t_env(x,phase_in_exponent,inv_t0squared);
    }

  //Transverse profile
  if (laserpulse2::w0>0.0L) 
    ans = ans*laserpulse2::w0/sqrt(beam_waist_sq2(x))*exp(-rsquared/beam_waist_sq2(x));
  return ans;
}
long double a_laser_envelope2(fourvec(&x))
{
  long double ans = Psi_function2(x);
  long double zpos=(*laserpulse2::Kgroup)*x;
  long double rsquared = x.getx()*x.getx()+((*laserpulse2::rproj)*x) * ((*laserpulse2::rproj)*x);
    
  // New bit AGRT June 2011 - 2nd order corrections to fields
 if (laserpulse2::w0>0.0L) {
  long double f_factor = laserpulse2::w0/sqrt(beam_waist_sq2(x))*sin(Gouy_phase2(x));
  long double correctionfactor=1.0L;
  correctionfactor+=0.5L* laserpulse2::focus_theta_2*laserpulse2::focus_theta_2*(1.0L-0.5L*(1.0L-zpos*zpos/laserpulse2::zR/laserpulse2::zR)/(1.0L+zpos*zpos/laserpulse2::zR/laserpulse2::zR)*rsquared/laserpulse2::w0/laserpulse2::w0*f_factor);
  ans=ans*correctionfactor;
}
  return ans;
}
long double phi_laser_envelope2(fourvec(&x))
{

  long double ans = a_laser_envelope2(x);
 long double zpos=(*laserpulse2::Kgroup)*x;
 //long double rsquared = x.getx()*x.getx()+((*laserpulse2::rproj)*x) * ((*laserpulse2::rproj)*x);
    
  ans -= laserpulse2::focus_theta_2*laserpulse2::focus_theta_2*(1.0L-zpos*zpos/laserpulse2::zR/laserpulse2::zR)/(1.0L+zpos*zpos/laserpulse2::zR/laserpulse2::zR)* Psi_function2(x);
  //
  ans *= laserpulse2::focus_theta_2*((*laserpulse2::pol_1)*x)/sqrt(beam_waist_sq2(x));
  return ans;
}

namespace nolaserpulse {
  int a_laser(fourvec x,fourvec &ans)
  {

  	ans.reset();
    return 0;
  }
}
namespace onelaserpulse {
  int a_laser(fourvec x,fourvec &ans)
  {

      ans=a_laser_phase(x);
      long double a0temp=(a_laser_envelope(x)*a0_pulse1);
      ans.setx(ans.getx()*a0temp);
      ans.sety(ans.gety()*a0temp);
      ans.setz(ans.getz()*a0temp);
      ans.sett(ans.gett()*phi_laser_envelope(x)*a0_pulse1);

    return 0;
  }
}
namespace otherlaserpulse {
  int a_laser(fourvec x,fourvec &ans)
  {
      ans=a_laser_phase2(x);
      long double a0temp=(a_laser_envelope2(x)*a0_pulse2);
      ans.setx(ans.getx()*a0temp);
      ans.sety(ans.gety()*a0temp);
      ans.setz(ans.getz()*a0temp);
      ans.sett(ans.gett()*phi_laser_envelope2(x)*a0_pulse2);
    return 0;
  }
}
namespace twolaserpulses {
  int a_laser(fourvec x,fourvec &ans)
  {
    onelaserpulse::a_laser(x,ans);
    fourvec ans2;
    otherlaserpulse::a_laser(x,ans2);
    ans = ans+ans2;
    return 0;
  }
}

//*******************************************************
// Plasma structures
// Channels or bubble or wakefield
namespace withramp {
long double ramp(fourvec x)
{
	long double ans = 0.5L-0.5L*tanh(2.0L*(x.getz()-bub_width)/bub_ramp);
	return ans;
}
}
namespace withnoramp {
long double ramp(fourvec x)
{
	return 1.0L;
}
}

fourvec a_channel(fourvec &x)
{
  // Add ion channel
  fourvec ans2(0.0L,0.0L,0.0L,0.0L);
  long double rsquared=(x.getx()*x.getx()+x.gety()*x.gety());
  // long double rsquared=(x.getx()*x.getx()+x.getz()*x.getz());
  //long double multiplier = 0.0;
  if (rsquared<r_channelsquared)// multiplier=1.0;
    ans2.sett(phi0_chan*rsquared/r_channelsquared);
  return ans2;
}
namespace bubblewake {
fourvec a_wakefield(fourvec &x)
{
  /*
    Bubble consists of a focusing electric field and azimuthal B field
  */

	// ramp value
long double density = ramp(x);

  // Electromagentic potential
// phi is proportional to density
  fourvec ans1(1.0L,0.0L,0.0L,v_wake*inv_bubble_aspect); //z term is magnetic field
  fourvec ans2;
  /*
    B=E*vphase in these units
// radius gets smaller by sqrt(density)
   */
  long double rsquared=(x.getx()*x.getx()+x.gety()*x.gety())*density;
  long double z_pos = (x.getz()-v_wake*x.gett()+r_bub)*inv_bubble_aspect+bub_delay;
	z_pos = z_pos*((z_pos<0.0L)*sqrt(density)+(z_pos>=0.0L)*1.0L);
  long double zsquared = z_pos*z_pos;
  long double r_pos=sqrt(rsquared);
  long double ebunch_profile =((z_pos+r_sheath_half)*(z_pos+r_sheath_half))*(inv_r_sheath_thick_square*4.0L) //*ratio_sheath_to_phi_bunch
    +(4.0L*rsquared+zsquared*0.0L)*inv_r_bubsquared; 
  // if (ebunch_profile<0.0) ebunch_profile=0.0;
  // This next piece adds the effect of an electron bunch at the rear of
  long double ans3 =exp(-sqrt(ebunch_profile))*phi_bunch; //exp(ebunch_profile)*phi_bunch;
  
/*
  // This adds an additional electron sheath
  ebunch_profile =(r_pos-r_sheath_half)*(r_pos-r_sheath_half)*inv_r_sheath_thick_square*ratio_sheath_to_extrasheath;
  if (z_pos<0.0) ebunch_profile+=(zsquared)*inv_r_bubsquared;
  else ebunch_profile+=(zsquared)*inv_r_bubsquared*4.0;

  long double ans4 = exp(-sqrt(ebunch_profile))*phi_sheath; //exp(ebunch_profile)*phi_sheath;

*/


  //now for z part of bubble
  rsquared+=zsquared;
  r_pos=sqrt(rsquared);
  if (rsquared<(r_bubsquared)) 
    ans2=ans1*(rsquared*inv_r_bubsquared);
   else 
     ans2=ans1*(one_plus_rbs_2*inv_r_bub*r_pos-one_plus_rbs-ratio_bub_sheath*inv_r_bubsquared*rsquared);

  if (r_pos>r_sheath) ans2=ans1*(1.0L+inv_bub_sheath);
  
  //ans2.sett(ans2.gett()+ans3+ans4);

   ans2=ans2-ans1*(1.0L+inv_bub_sheath);
  ans2=ans2*phi0_bub+ans1*ans3;

  /* long double multiplier = ans2.gett()*(1.0+0.5*a0*a0*(1.0-a_laser_phase(x).square3())*a_laser_envelope(x)*a_laser_envelope(x));
  if (a0_0)
    ans2.sett(multiplier);*/
  return ans2;
}
}
namespace linearwake {
// Linear wakefield
fourvec a_wakefield(fourvec &x)
{	
     long double density = ramp(x);
     // radius squared normalized to wake radius
     long double rsquared=(x.getx()*x.getx()+x.gety()*x.gety())*inv_r_bubsquared*density;
     // this is k_p (z-v_p t)
     long double kpxi = (x.getz()-v_wake*x.gett()+bub_delay)*(v_wake/rw0wp*sqrt(density));
     kpxi = kpxi*(kpxi<0.0L); // this is so wake has a front!
     long double myphi = phi0_bub*exp(-2.0L*rsquared)*sin(kpxi);
     fourvec ans(myphi,0.0L,0.0L,v_wake*myphi);
     return ans;
}
}

// Bfield, uniform

fourvec a_Bfield(fourvec &x)
{
  fourvec ans(0.0L,0.0L,0.0L,0.0L);
  ans.setx(x.getz()*B_field_hat[1]-x.gety()*B_field_hat[2]);
  ans.sety(x.getx()*B_field_hat[2]-x.getz()*B_field_hat[0]);
  ans.setz(x.gety()*B_field_hat[0]-x.getx()*B_field_hat[1]);
  ans=ans*(0.5L*mag_Bfield);
  return ans;
}

fourvec a_wiggler(fourvec &x)
{
  /*// New wiggler field!
  // Parameterized by magnet spacing (normalized to c/omega) and
  // magnet strength (normalized B) and
  // wiggler length (normalized to c/omega) -> start to end
  fourvec ans(0.0L,0.0L,0.0L,0.0L);
  long double ztemp=x.getz();
  long double periodic=0.0L;
  long double delta = 0.0L;
  if (ztemp>wiggler_start && ztemp<wiggler_end)
    {
      periodic = sin((ztemp-wiggler_start)*wiggler_k);
      if (wiggler_square) 
		{
		  periodic = asin(periodic);
		  // add gap between wiggler b fields
	  		if (wiggler_gap) {
	   		 delta = (0.5L*PI-wiggler_gap*wiggler_k);
	   		 if (periodic>delta){ periodic = delta;}
	   		 if (periodic<-delta){ periodic = -delta;}
	    
	 		 }
		}
      else
		{
		  periodic = (periodic)/wiggler_k;
		}
		
      ans.setx(-periodic*wiggler_Bfield);
    }
  */
  
  // New new wiggler - using B field given in paper by 
  // Sagan in PAC conference proceedings 2003 (particle accelerator conf).
  fourvec ans(0.0,0.0,0.0,0.0);

  long double xtemp = x.getx();
  long double ytemp = x.gety();
  long double zdis = x.getz()-wiggler_start;

// wiggler_k is now actually the period 
  const long double localkz=(2.0L*PI/wiggler_k);

// calculated from wiggler_gap - which is magnetic width!
  const long double localkx=(PI/wiggler_gap);
 

// localky chosen to satisfy Maxwell's equations
  const long double localky=sqrt(localkx*localkx+localkz*localkz); 

// wiggler_k is now actually the period not k
  const long double kysquare = localky*localky;
  long double par = cosh(localky*ytemp);

  // wiggle_Bfield is C
  //equations used: Bx = -C*kx/ky*sin(kx*x)*sinh(ky*y)*cos(kz*z);
  //By = C*cos(kx*x)*cosh(ky*y)*cos(kz*z); Bz = -kz/ky*cos(kx*x)*sinh(ky*y)*sin(kz*z);

long double wiggler_xsize = 0.5L-0.5L*tanh(2.0*(fabs(xtemp)/(0.5L*wiggler_gap) - 1.0L));
long double wiggler_zsize = 0.25L*(1.0L+tanh(2.0*localkz*(zdis-0.5L*wiggler_k)))*(1.0L-tanh(2.0*localkz*(zdis-wiggler_end+0.5L*wiggler_k)));
  //if (zdis>=(0.0L-0.25L*wiggler_k) && zdis<=(wiggler_end+0.25L*wiggler_k))
    //{
     long double Az = -(localkx/kysquare)*sin(localkx*xtemp)*par*cos(localkz*zdis);
      ans.setz(Az*wiggler_Bfield*wiggler_xsize*wiggler_zsize);
    //}
  //if (zdis>=0.0L && zdis<=(wiggler_end))
    //{
      long double Ax = (localkz/kysquare)*cos(localkx*xtemp)*par*sin(localkz*zdis);
      ans.setx(Ax*wiggler_Bfield*wiggler_xsize*wiggler_zsize);
    //}

  return ans; 
}

// Plane standing wave for test problems
fourvec a_standing(fourvec &x)
{
  fourvec ans(0.0L,0.0L,0.0L,0.0L);
  long double stand_phase=sin(standing_wave_k*(x.getz())+standing_wave_phase);
  //stand_phase=stand_phase*sin(standing_wave_k*(x.gett()));
  if (if_standing_circ) {
   ans.setx(standing_wave_a*stand_phase*sin(standing_wave_k*(x.gett())+standing_wave_pol));
   ans.sety(standing_wave_a*stand_phase*cos(standing_wave_k*(x.gett())+standing_wave_pol));
  } else {
   ans.setx(costhetastanding*standing_wave_a*stand_phase*sin(standing_wave_k*(x.gett())));
   ans.sety(sinthetastanding*standing_wave_a*stand_phase*sin(standing_wave_k*(x.gett())));
  }
  return ans;
}



int a_field(fourvec x,fourvec &ans)
{
  a_laser(x,ans);
  if (if_standing_wave) ans = ans+a_standing(x);
  if (if_channel) ans=ans+a_channel(x);
  if (if_wakefield) ans=ans+a_wakefield(x);
  if (if_bfield) ans=ans+a_Bfield(x);
  if (if_wiggler) ans=ans+a_wiggler(x);
  return 0;
}

#endif
