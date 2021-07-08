
#ifndef SPIN_H_GUARD_
#define SPIN_H_GUARD_
#include "constants.h"
#include "headers.h"

// Currently accepted NIST g-factor values
const long double g_e = -2.00231930436182;

// S is spin four vector (beta.S, S)

// BMT equation dS/dtau = f(tau,x,v)
fourvec BMT(long double &tau,fourvec& x,fourvec& v,fourvec& S)
{
  fourvec term1(0.0L,0.0L,0.0L,0.0L);
  term1=F(x)*S; 
  term1=term1*(-q_over_m*g_e*0.5L); //charge to mass ratio * g factor
  
  fourvec term2(0.0L,0.0L,0.0L,0.0L);
  
  long double term2scalar = S*(F(x)*v); // should be a scalar!
  
  term2 = v;
  term2=term2*(-q_over_m*(g_e*0.5L-1.0L)*term2scalar); //charge to mass ratio * g factor
  
  fourvec myans(0.0L,0.0L,0.0L,0.0L);
  myans = term1+term2;
  
  return myans;
}


int initialize_spins(long double **par_v,long double **par_S)
{

long double rest_s[4];
for (int i=0;i<Npar;++i)
   {
   		// Define spin vector in rest frame
   		rest_s[0] = 0.0L;
   		rest_s[1] = sqrt(2.0L)*0.5L;//0.0L;
   		rest_s[2] = 0.0L;
   		rest_s[3] = 0.5L;
   		

   		// fix 0 component
   		par_S[0][i] = 0.0L;
		for (int j=0;j<4;++j)
		{
			par_S[0][i] = par_S[0][i] + par_v[j][i]*rest_s[j];
		}	 
		
		for (int j=1;j<4;++j)
		{
			par_S[j][i] = rest_s[j]+par_S[0][i]*par_v[j][i]/(par_v[0][i]+1.0);
		}		
   				
   }
return 0;
}


#endif
