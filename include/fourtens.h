#pragma once
#ifndef FOURTENS_H_GUARD
#define FOURTENS_H_GUARD
#include "fourvec.h"

class fourtens
{
 private:
  long double coord[4][4];  //private data member
  
 public:
  //default constructor
  fourtens() {}
  fourtens(fourvec& vec1,fourvec& vec2) // Dyadic constructor
  {
    for (int ii=0;ii<4;++ii)
      for (int jj=0;jj<4;++jj)
	coord[ii][jj]=vec1(ii)*vec2(jj);
  }
  long double operator() (int& i,int& j)
  {
    return coord[i][j];
  }
  void dyadic(fourvec& vec1,fourvec& vec2) // Dyadic
  {
    for (int ii=0;ii<4;++ii)
      for (int jj=0;jj<4;++jj)
	coord[ii][jj]=vec1(ii)*vec2(jj);

  }

  fourvec operator*(fourvec vec)
    {
      fourvec ans(0.0,0.0,0.0,0.0);
      fourvec vec_mod=vec;
      vec_mod.Minkowskii();
      for (int ii=0;ii<4;++ii)
	for (int jj=0;jj<4;++jj)
	  ans.set(ii,ans(ii)+coord[ii][jj]*vec_mod(jj));
      return ans;
    }
  void set(int& ii,int& jj,long double value)
  {
    coord[ii][jj]=value;
  }
  long double get(int& ii,int& jj)
  {
    return coord[ii][jj];
  }
  fourtens operator+(fourtens& tens)
  {
    fourtens ans;
      for (int ii=0;ii<4;++ii)
	for (int jj=0;jj<4;++jj)
	  ans.set(ii,jj,get(ii,jj)+tens.get(ii,jj));
      return ans;
  }  
fourtens operator-(fourtens& tens)
  {
    fourtens ans;
      for (int ii=0;ii<4;++ii)
	for (int jj=0;jj<4;++jj)
	  ans.set(ii,jj,get(ii,jj)-tens.get(ii,jj));
      return ans;
  }

long double operator*(fourtens& tens) 
  {
    long double ans=0.0L;
    int ii,jj;
    jj = 0;
      for (ii=0;ii<4;++ii)
	ans=ans+get(ii,jj)*tens.get(ii,jj);
    ii=0;
      for (jj=1;jj<4;++jj)
	ans=ans+get(ii,jj)*tens.get(ii,jj);

      for (ii=1;ii<4;++ii)
	for (jj=1;jj<4;++jj)
	  ans=ans-get(ii,jj)*tens.get(ii,jj);
      return ans;
  }

 void operator=(fourtens tens)
  {
      for (int ii=0;ii<4;++ii)
	for (int jj=0;jj<4;++jj)
	  set(ii,jj,tens.get(ii,jj));

  }  
 void Null()
 {
   for (int ii=0;ii<4;++ii)
     for (int jj=0;jj<4;++jj)
	coord[ii][jj]=0.0;
 }
  void Identity()
  {
    for (int ii=0;ii<4;++ii)
      for (int jj=0;jj<4;++jj)
	coord[ii][jj]=0.0;
    for (int ii=0;ii<4;++ii)
      coord[ii][ii]=1.0;
  }
  void Minkowski()
  {
    for (int ii=0;ii<4;++ii)
      for (int jj=0;jj<4;++jj)
	coord[ii][jj]=0.0;
    for (int ii=1;ii<4;++ii)
      coord[ii][ii]=-1.0;
    coord[0][0]=1.0;
  }
  void view()
  {
    for (int ii=0;ii<4;++ii){
      for (int jj=0;jj<4;++jj)
	std::cout<<coord[ii][jj]<<',';
      std::cout<<'\n';
    }
    std::cout<<'\n';
  }
};
  
#endif
