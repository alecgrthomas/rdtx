#pragma once
#ifndef RADF_C_GUARD
#define RADF_C_GUARD
#include "fourvec.h"
#include "fourtens.h"
#include "constants.h"


// Radiation force
namespace norad
{
fourvec Lambda(fourvec& x,fourvec& v,fourvec& vdot,fourvec& vdotdot)
  {
    fourvec ans(0.0,0.0,0.0,0.0);

    return ans;
  }
}
namespace LAD // Lorentz-Abraham-Dirac force
{
  fourvec Lambda(fourvec& x,fourvec& v,fourvec& vdot,fourvec& vdotdot)
  {
    fourvec ans=(vdotdot+v*(vdot*vdot))*rad_const*(q_over_m*q_over_m);

    return ans;
  }
}
namespace PLAD // Physical Lorentz-Abraham-Dirac force - Rohlich PRE 08
{
  fourvec Lambda(fourvec& x,fourvec& v,fourvec& vdot,fourvec& vdotdot)
  {
// fourtens P_Rohlich(v,v); // P tensor in Rohlich PRE 08
  //  P_Rohlich=P_Rohlich+metrictens;
//maintain forward differencing
    fourvec xminus=x-v*delta_step; // delta step same space step dr as in field solver
    fourtens F_Rohlichminus;    
   // F_Rohlichplus=F(x); // was xplus
    F_Rohlichminus=F(xminus);

    // Lorentz force
    fourvec lforce=F(x)*v*(q_over_m); //charge to mass ratio

    // This is Schott term - small
    fourvec F_Rohlich=(lforce-F_Rohlichminus*v*q_over_m)*(0.5/delta_step) + F(x)*vdot*q_over_m; 
    // alternative formulation
   // fourvec ans=(P_Rohlich*(F_Rohlich))*(-rad_const*(q_over_m*q_over_m));

// second term is simple damping
   fourvec ans=(v*(lforce*lforce)-F_Rohlich)*(rad_const);
    return ans;
  }
}

// LAD model  with Kirk/Bell/Ridgers correction g
/* 
AGRT - although the function is called QPLAD, because the 
LAD form is more efficient, it is used
*/
namespace QPLAD // Quantum corrected Lorentz-Abraham-Dirac force 
{

 fourvec Lambda(fourvec& x,fourvec& v,fourvec& vdot,fourvec& vdotdot)
  {
    fourvec ans=(vdotdot+v*(vdot*vdot))*rad_const*(q_over_m*q_over_m);
    ans = ans*Kirk_g(vdot);
    return ans;
  }
 /*
OLD FORM - QPLAD:
 fourvec Lambda(fourvec& x,fourvec& v,fourvec& vdot,fourvec& vdotdot)
  {

//maintain forward differencing
    fourvec xminus=x-v*delta_step; // delta step same space step dr as in field solver
    fourtens F_Rohlichminus;    

    F_Rohlichminus=F(xminus);

    // Lorentz force
    fourvec lforce=F(x)*v*(q_over_m); //charge to mass ratio

    // This is Schott term - small
    fourvec F_Rohlich=(lforce-F_Rohlichminus*v*q_over_m)*(0.5/delta_step) + F(x)*vdot*q_over_m; 

// second term is simple damping
   fourvec ans=(v*(lforce*lforce)-F_Rohlich)*(rad_const*Kirk_g(lforce));


    return ans;
  }*/
}

/*namespace simple // Just damping
{
  fourvec Lambda(fourvec& x,fourvec& v,fourvec& vdot,fourvec& vdotdot)
  {
  fourvec lforce(0.0,0.0,0.0,0.0);
  lforce=F(x)*v; // Lorentz force
  lforce=lforce*(-q_over_m); //charge to mass ratio
    fourvec ans=(v*(lforce*lforce))*rad_const;
    return ans;
  }
}*/

/*namespace simple_w_qed_correct // Just damping with Kirk/Bell/Ridgers correction
{
  fourvec Lambda(fourvec& x,fourvec& v,fourvec& vdot,fourvec& vdotdot)
  {
  fourvec lforce(0.0,0.0,0.0,0.0);
  lforce=F(x)*v; // Lorentz force
  lforce=lforce*(-q_over_m); //charge to mass ratio
    fourvec ans=(v*(lforce*lforce))*rad_const*Kirk_g(lforce);
    return ans;
  }
}*/

// This is Kirk PPCF 2008 g function as fitted by CP Ridgers
double Kirk_g(fourvec lforce){
double chi = sqrt(fabs(lforce*lforce))*qed_const; // chi parameter
double ans = pow(3.7*chi*chi*chi+31.0*chi*chi+12.0*chi+1.0,-4.0/9.0);
return ans;
}
// This is chi parameter
double chi_param(fourvec &x,fourvec &v){
  fourvec lforce(0.0,0.0,0.0,0.0);
  lforce=F(x)*v; // Lorentz force
  lforce=lforce*(q_over_m); //charge to mass ratio
double chi = sqrt(fabs(lforce*lforce))*qed_const; // chi parameter
return chi;
}
#endif
