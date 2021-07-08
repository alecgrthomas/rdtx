#pragma once
#ifndef FOURVEC_H_GUARD
#define FOURVEC_H_GUARD
#include <iostream>
#include <fstream>
#include <math.h>

class fourvec
{
 private:
 long double coord[4];  //private data member
  
 public:
  //default constructor
  fourvec() {}
  
  //cartesian constructor
  fourvec(long double t,long double x, long double y, long double z)
    {
      coord[1] = x;
      coord[2] = y;
      coord[3] = z;
      coord[0] = t;
    }

  /*//method to write data to file (suitable for graphs)
    void graphwrite (ofstream & outfile)
    {
    outfile << xcoord << "\t" 
    << ycoord << "\t"
    << zcoord << endl;
    }
  */
  void reset()
  {
      for (int ii=0;ii<4;++ii)
      set(ii,0.0L);
  }
  
    void view ()
    {
     for (int ii=0;ii<4;++ii)
       std::cout<<coord[ii]<<',';
     std::cout<<'\n';
    }
  void copy(fourvec& vec)
  {
    for (int ii=0;ii<4;++ii)
      set(ii,vec(ii));
  }
  //access functions
  //access xcoord
  long double getx()
  {
    return coord[1];
  }
  
  //access ycoord
  long double gety()
  {
    return coord[2];
  }	
  
		//access zcoord
  long double getz()
  {
    return coord[3];
  }
  //access tcoord
  long double gett()
  {
    return coord[0];
  }
  
  //access function for square of 4vector - Minkowski space
  long double square()
  {
    return coord[0]*coord[0]-(coord[1]*coord[1] +
			      coord[2]*coord[2] +coord[3]*coord[3]);
  }

    //access function for square of 3vector
  long double square3()
  {
    return (coord[1]*coord[1]+coord[2]*coord[2] +coord[3]*coord[3]);
  }

  //access function for magnitude of 4vector
  long double mag()
  {
    return sqrt(square());
  }
    //access function for magnitude of 3vector
  long double mag3()
  {
    
    return sqrt(square3());
  }

  // This method turns dp^mu/dt to four force
  void fourforce(fourvec& v)
  {
    sett(gett()*v.mag3());
    setx(getx()*v.gett());
    sety(gety()*v.gett());
    setz(getz()*v.gett());
  } 
  void Minkowskii()
  {
    setx(-getx());
    sety(-gety());
    setz(-getz());
  }
  //methods to allow programs to rest vectors
  void setx(long double a)
  {
    coord[1] = a;
  }
  
  void sety(long double b)
  {
    coord[2] = b;
  }
  
  void setz(long double c)
  {
    coord[3] = c;
  }
  void sett(long double c)
  {
    coord[0] = c;
  }

  void hat(int ii) // makes into a directionvector
  {
    for (int i=0;i<4;++i)
      coord[i]=0.0;
    coord[ii]=1.0;
  }
  long double operator() (int index)
  {
    return coord[index];
  }
  long double get(int index)
  {
    return coord[index];
  }
  void set(int index,long double value)
  {
    coord[index]=value;
  }

    //pointwise multiply product

  fourvec pointwise(fourvec v_old)
  {
    fourvec v_new;
    for (int ii=0;ii<4;++ii)
      v_new.set(ii,v_old(ii)*get(ii));
    return v_new;
  }

// frame moving at speed of light in z
fourvec speedoflight()
{
	fourvec ans(coord[0],coord[1],coord[2],coord[3]-coord[0]);
	return ans;
}
  //operator overloadings (to allow maths with vectors)
 
  
  //vector addition
  fourvec operator+(fourvec v_old)
  {
    fourvec v_sum;
    for (int ii=0;ii<4;++ii)
      v_sum.set(ii,v_old(ii)+get(ii));
    return v_sum;
  }
	
  //vector subtraction
  fourvec operator-(fourvec v_old)
  {
    fourvec v_sum;
    for (int ii=0;ii<4;++ii)
      v_sum.set(ii,get(ii)-v_old(ii));
    return v_sum;
  }

  //multiplication of vector by scalar
  fourvec operator*(long double scalar)
    {
      fourvec v_sum;
      for (int ii=0;ii<4;++ii)
	v_sum.set(ii,get(ii)*scalar);
      return v_sum;
    }
  //multiplication of vector by scalar
  fourvec operator/(long double scalar)
    {
      fourvec v_sum;
      for (int ii=0;ii<4;++ii)
	v_sum.set(ii,get(ii)/scalar);
      return v_sum;
    }
  //method to perform vector multiplication (cross product)
  /* vector operator ^(vector b)
  {
    vector y (  ycoord*b.getz() - zcoord*b.gety(),
		zcoord*b.getx() - xcoord*b.getz(),
		xcoord*b.gety() - ycoord*b.getx() );
    return y;
    }*/
  //inner product
  long double operator*(fourvec vec)
  {
    long double v_product=get(0)*vec(0);
    for (int ii=1;ii<4;++ii)
      v_product-=get(ii)*vec(ii);
    return v_product;
  }

  //method to set one vector equal to another

 void operator =(fourvec b)
    {
      for (int ii=0;ii<4;++ii)
	coord[ii]=b(ii);
    }

};	

#endif
