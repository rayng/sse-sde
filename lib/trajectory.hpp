#ifndef traj_HPP
#define traj_HPP


#include <fstream>
#include <string>
#include <iomanip>
#include <complex>
#include <math.h>
#include <iostream>
#include <RandomLib/Random.hpp>
#include <RandomLib/NormalDistribution.hpp>
#include "globalparams.hpp"
#include "timebin.hpp"
#include "bin.hpp"

using std::cout;


extern RandomLib::Random r;
extern RandomLib::NormalDistribution<double> normdist;


class trajectory
{
 public:
  // member data
  complex<double> z[N], zp[N];
  complex<double> y[N], yp[N];
  complex<double> zetafcn[N];
  int SS[N], SSp[N];
  timebin tbin;
  
  // member functions
  trajectory(): tbin(nT) {};   // Not sure why this works.
  
  void evolve();
  void evolve_y();
  void calcobservables();
  
  void init();
  void init_y();
  void write_variables(ofstream& fp, int t, string var);   // var is either z or y 
  void genzfromy(bool write, ofstream& fp);
  
  
 private:  
  // private member function
  void calc_zetafcn(complex<double>* y, complex<double>* yp);
  
};


#endif
