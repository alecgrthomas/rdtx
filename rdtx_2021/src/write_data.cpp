/*

  data write functions for rdt
  agrt oct 2010
  
 */
#include "headers.h"
double chi_param(fourvec &x,fourvec &v);
void operator<< (std::ofstream& str_out, fourvec vec)
{
  str_out<<vec.gett()<<'\t'<<vec.getx()<<'\t'<<vec.gety()<<'\t'<<vec.getz();
}

int write_data(std::ofstream &outfile, fourvec x, fourvec v,fourvec S, double Pave, fourvec a_las)
{
if (!if_moving_frame) {
  outfile<<x;
} else {
  outfile<<x.speedoflight();
}
  outfile<<"\t";
  outfile<<v;
  outfile<<'\t'<<Pave;
  outfile<<'\t';
  outfile<<a_las;
  outfile<<'\t';
  outfile<<chi_param(x,v);
  outfile<<'\t';  
  outfile<<S;
  outfile<<'\n';
  return 0;
}
int write_parhead(std::ofstream &outfile,int par_label,int Ndumps_mod,double v_wake)
{
  outfile<<"Particle data\n0\nPar_num_";
  outfile<<par_label<<"\nxy\n";
  outfile<<Ndumps_mod<<'\n'<<v_wake<<"\n0\n";
  outfile<<"t\tx\ty\tz\tgamma\tpx\tpy\tpz\tPave\tphi\tax\tay\taz\tchi\tS0\tSx\tSy\tSz\n";
  return 0;
}

fourtens F(fourvec& x);

// get electric field component from EM tensor - use 1,2,3
double Efield(fourvec &x,int &component) {
	fourtens FS = F(x);
	int zero=0;	
	return (FS(zero,component));
}
// get magnetic field component from EM tensor - use 1,2,3
double Bfield(fourvec &x,int &component) {
	fourtens FS = F(x);
	// permute indices
	int index2=component+1,index1=component+2;	
	if (index1>3) index1-=3;
	if (index2>3) index2-=3;
	// should be 3,2 or 1,3 or 2,1
	return -1.0*(FS(index1,index2));
}
int dump_fields(std::ofstream &outfile_fields,std::ofstream &outfile_potentials,long double time)
{
  outfile_potentials<<"Potential data, slice through y="<<yslice<<" plane\n0\n";
  outfile_potentials<<"|A|_and_phi\nzx\n";
  outfile_potentials<<Nx_f<<'\n'<<Nz_f<<"\n"<<xsize<<"\n"<<zsize<<"\n";
  outfile_potentials<<time<<"\n";
  fourvec xtemp(time,0.0,yslice,0.0);
  fourvec a_field_temp;
  long double xstart=-xsize*0.5, zstart=-zsize*0.5-+time*(int) if_moving_frame;
  const long double dxstart=fabs(xsize)/Nx_f, dzstart=fabs(zsize)/Nz_f; 
 

 for (int i=0;i<Nx_f;++i)
    {
      xtemp.setx(xstart);
      
      for (int xx=0;xx<4;++xx)
	{
	  zstart=-zsize*0.5+time*(int) if_moving_frame;

	 for (int j=0;j<Nz_f;++j)
	    {	      
	      xtemp.setz(zstart);
	      a_field(xtemp,a_field_temp);
	      
	      outfile_potentials<<a_field_temp.get(xx)<<'\t';

	      zstart+=dzstart;
	    }
	}
      outfile_potentials<<'\n';    
      xstart+=dxstart;
    }
 
 outfile_fields<<"Field data, slice through y="<<yslice<<" plane\n0\n";
  outfile_fields<<"E_and_B\nzx\n";
  outfile_fields<<Nx_f<<'\n'<<Nz_f<<"\n"<<xsize<<"\n"<<zsize<<"\n";
  outfile_fields<<time<<"\n";
//reset
xstart=-xsize*0.5;
 zstart=-zsize*0.5-+time*(int) if_moving_frame;

 for (int i=0;i<Nx_f;++i)
    {
      xtemp.setx(xstart);
      
      for (int xx=1;xx<4;++xx) // 3 components each - of E and B
	{
	  zstart=-zsize*0.5+time*(int) if_moving_frame;
	 for (int j=0;j<Nz_f;++j)
	    {	      
	      xtemp.setz(zstart);
	      outfile_fields<<Efield(xtemp,xx)<<'\t';
	      zstart+=dzstart;
	    }
	zstart=-zsize*0.5+time*(int) if_moving_frame;
	 for (int j=0;j<Nz_f;++j)
	    {	      
	      xtemp.setz(zstart);
	      outfile_fields<<Bfield(xtemp,xx)<<'\t';
	      zstart+=dzstart;
	    }
	}
      outfile_fields<<'\n';    
      xstart+=dxstart;
    }
  return 0;
}
