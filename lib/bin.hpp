#ifndef Bin_HPP
#define Bin_HPP



#include <fstream>
#include <iomanip>
#include "globalparams.hpp"


using namespace std;

class bin
{

public:
  long double obs[nobs];
  long double physobs[nobs];   // hard coded size of array. 
  bin();
  void operator= ( bin &rhs);
  void operator+= ( bin &rhs);
  void dump(long double O, int i );
  void average(int nn);
  void parityaverage(int nn);
  void zero();
  void generate_physical_obs(double& hh);
  void write_to_file(ofstream& fp, const double  (&epsilonn)[Nb]);
  
};

#endif

