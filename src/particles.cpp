
#ifndef FIELDPAR_H_GUARD_
#define FIELDPAR_H_GUARD_
#include "fourvec.h"
#include "constants.h"
#include "headers.h"

//*******************************************************

// Gaussian random distribution function
long double gaussian()
{
	long double x, y;
	do
	{
		x = 6.0*rand()/RAND_MAX-3.0;  // +/-3 sigma range
		y = exp(-x*x/2.0); // modified compared to earlier versions - now like in Thomas PRX 2012
	} while (1.0*rand()/RAND_MAX > y);
	return x;
}

// uniform distribution - in sphere, gives x,y,z
fourvec uniformdist()
{
	fourvec ans;
	// first find r,theta,phi
	// spherical coords so uniform in r^3
	long double r = pow((1.0/RAND_MAX)*rand(),1.0/3.0);// 0->1 range
	long double costheta = ((2.0/RAND_MAX)*rand()-1.0);// acos (0 to pi);
	long double phi = (2.0*PI/RAND_MAX)*rand();// 0 to 2pi

	ans.setx(r*sqrt(1.0-costheta*costheta)*cos(phi));
	ans.sety(r*sqrt(1.0-costheta*costheta)*sin(phi));
	ans.setz(r*costheta);
return ans;
}

// lwfa distribution as in Esirkepov PRL 2006 random distribution function
// I have modified to include finite width deltaE because otherwise 
// distribution is singular
long double lwfadist(long double deltaE)
{
	long double x, y;
	if (deltaE<=0.0){
		deltaE =1e-9;	
	}
	do
	{
		x = (1.0/RAND_MAX)*rand();  // 0->1 range
		// the maximum of this function is 1
		if (x<(1.0-deltaE)) {
			y = sqrt(1.0-x)*(1.0-sqrt(1.0-deltaE/(1.0-x)))/sqrt(deltaE);
		} else { // x 1-dE -> 1
			y = sqrt((1.0-x)/deltaE);
		}
	} while (1.0*rand()/RAND_MAX > y);
	return x;
}
// this function adds focusing effect to particles -
// particle momenta modified as if they have eminated from single point
// z focus away. This function adds spatially dependent momentum
// make sure function not called for parzfocus=0!

int focus_particles(long double **par_x,long double **par_v){
// only for more than 1 particle!
 if (Npar>1&&parzfocus)
    for (int i=0;i<Npar;++i)
      {
	// 1 is x and 2 is y - minus so focusing for positive parzfocus
	par_v[1][i]-=par_v[3][i]*par_x[1][i]/parzfocus;
	par_v[2][i]-=par_v[3][i]*par_x[2][i]/parzfocus;
	
	/*
	  Now fix gamma
	  Gamma^2 = 1+p^2
	*/
	par_v[0][i]=1.0;
	for (int j=1;j<4;++j)
	  par_v[0][i]+=par_v[j][i]*par_v[j][i];
	par_v[0][i]=sqrt(par_v[0][i]);
      }
return 0;
}

int initialize_spins(long double **par_v,long double **par_S);


namespace gaussian_bunch
{
int initialize_particles(long double **par_x,long double **par_v,long double **par_S)
{
  // Distribute particles in a gaussian profile randomly
  if (Npar>1)
    for (int i=0;i<Npar;++i)
      {
	//time =0.0 for all
	par_x[0][i]=0.0;
	
	/* How do we want x,y, z distributed? */
	for (int j=1;j<4;++j)
	  par_x[j][i]=bunch_centre[j-1];
	// Add on a gaussian spread
	for (int j=1;j<4;++j)
	  if (bunch_width[j-1]>0.0) 
	    par_x[j][i]+=bunch_width[j-1]*gaussian();
	
	/* How do we want momenta distributed? */
	for (int j=1;j<3;++j)
	  par_v[j][i]=bunch_momentum[j-1]; 
	// Add on a gaussian spread
	for (int j=1;j<3;++j)
	  if (bunch_momentum_spread[j-1]>0.0)
	    par_v[j][i]+=bunch_momentum_spread[j-1]*gaussian();
	
	if (!if_LWFA_dist) {
	par_v[3][i]=bunch_momentum[2]; 
	if (bunch_momentum_spread[2]>0.0) {
	    par_v[3][i]+=bunch_momentum_spread[2]*gaussian();
	    }
	} else {
		par_v[3][i]=bunch_momentum[2]*lwfadist(fabs(bunch_momentum_spread[2]/bunch_momentum[2]));
	}
	/*
	  Now fix gamma
	  Gamma^2 = 1+p^2
	*/
	par_v[0][i]=1.0;
	for (int j=1;j<4;++j)
	  par_v[0][i]+=par_v[j][i]*par_v[j][i];
	par_v[0][i]=sqrt(par_v[0][i]);
      }
  else
    {
      par_x[0][0] = 0.0;
      /* How do we want x,y, z distributed? */
      for (int i=1;i<4;++i)
	par_x[i][0]=bunch_centre[i-1];
      /* How do we want momenta distributed? */
      for (int j=1;j<4;++j)
	par_v[j][0]=bunch_momentum[j-1]; 
      // Gamma^2 = 1+p^2
      par_v[0][0]=1.0;
      for (int j=1;j<4;++j)
	par_v[0][0]+=par_v[j][0]*par_v[j][0];
      par_v[0][0]=sqrt(par_v[0][0]);
    }
    
    initialize_spins(par_v,par_S);
  return 0;
}
}
namespace uniform_bunch
{
  int initialize_particles(long double **par_x,long double **par_v,long double **par_S)
{
  fourvec deltavector;
  // Distribute particles in a uniform profile randomly
  if (Npar>1)
    for (int i=0;i<Npar;++i)
      {
	//time =0.0 for all
	par_x[0][i]=0.0;
	
	/* How do we want x,y, z distributed? */
	for (int j=1;j<4;++j)
	  par_x[j][i]=bunch_centre[j-1];
	
	// Add on a uniform spread
	
	// calculate dx,dy,dz
	deltavector = uniformdist();

	for (int j=1;j<4;++j)
	  if (bunch_width[j-1]>0.0) 
	    par_x[j][i]+=bunch_width[j-1]*deltavector.get(j);
	
	/* How do we want momenta distributed? */
	for (int j=1;j<4;++j)
	  par_v[j][i]=bunch_momentum[j-1]; 
	  
	// Add on a uniform spread
	// calculate dx,dy,dz
	deltavector = uniformdist();
	for (int j=1;j<4;++j)
	  if (bunch_momentum_spread[j-1]>0.0)
	    par_v[j][i]+=bunch_momentum_spread[j-1]*deltavector.get(j);
	
	if (!if_LWFA_dist) {
	par_v[3][i]=bunch_momentum[2]; 
	if (bunch_momentum_spread[2]>0.0) {
	    par_v[3][i]+=bunch_momentum_spread[2]*deltavector.get(3);
	    }
	} else {
		par_v[3][i]=bunch_momentum[2]*lwfadist(fabs(bunch_momentum_spread[2]/bunch_momentum[2]));
	}
	/*
	  Now fix gamma
	  Gamma^2 = 1+p^2
	*/
	par_v[0][i]=1.0;
	for (int j=1;j<4;++j)
	  par_v[0][i]+=par_v[j][i]*par_v[j][i];
	par_v[0][i]=sqrt(par_v[0][i]);
      }
  else
    {
      par_x[0][0] = 0.0;
      /* How do we want x,y, z distributed? */
      for (int i=1;i<4;++i)
	par_x[i][0]=bunch_centre[i-1];
      /* How do we want momenta distributed? */
      for (int j=1;j<4;++j)
	par_v[j][0]=bunch_momentum[j-1]; 
      // Gamma^2 = 1+p^2
      par_v[0][0]=1.0;
      for (int j=1;j<4;++j)
	par_v[0][0]+=par_v[j][0]*par_v[j][0];
      par_v[0][0]=sqrt(par_v[0][0]);
    }
        initialize_spins(par_v,par_S);
  return 0;
}
}

namespace line_bunch
{
  int initialize_particles(long double **par_x,long double **par_v,long double **par_S)
{
  // Distribute particles in a square profile sequentially
  int centre=(Npar)/2; // subtractor
  long double multiplier=-1.0/(long double)centre;

  if (Npar>1)
    for (int i=0;i<Npar;++i)
      {
	//time =0.0 for all
	par_x[0][i]=0.0;
	
	/* How do we want x,y, z distributed? */
	for (int j=1;j<4;++j)
	  par_x[j][i]=bunch_centre[j-1];
	// Add on a gaussian spread
	for (int j=1;j<4;++j)
	  if (bunch_width[j-1]>0.0) 
	    par_x[j][i]+=bunch_width[j-1]*multiplier*(i-centre);
	
	/* How do we want momenta distributed? */
	for (int j=1;j<4;++j)
	  par_v[j][i]=bunch_momentum[j-1]; 
	// Add on a gaussian spread
	for (int j=1;j<4;++j)
	  if (bunch_momentum_spread[j-1]>0.0)
	    par_v[j][i]+=bunch_momentum_spread[j-1]*multiplier*(i-centre);
	
	/*
	  Now fix gamma
	  Gamma^2 = 1+p^2
	*/
	par_v[0][i]=1.0;
	for (int j=1;j<4;++j)
	  par_v[0][i]+=par_v[j][i]*par_v[j][i];
	par_v[0][i]=sqrt(par_v[0][i]);
      }
  else
    {
      par_x[0][0] = 0.0;
      /* How do we want x,y, z distributed? */
      for (int i=1;i<4;++i)
	par_x[i][0]=bunch_centre[i-1];
      /* How do we want momenta distributed? */
      for (int j=1;j<4;++j)
	par_v[j][0]=bunch_momentum[j-1]; 
      // Gamma^2 = 1+p^2
      par_v[0][0]=1.0;
      for (int j=1;j<4;++j)
	par_v[0][0]+=par_v[j][0]*par_v[j][0];
      par_v[0][0]=sqrt(par_v[0][0]);
    }
        initialize_spins(par_v,par_S);
  return 0;
}
}

#endif
