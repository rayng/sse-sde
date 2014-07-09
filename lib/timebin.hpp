#ifndef timebin_HPP
#define timebin_HPP


#include <fstream>
#include <iomanip>
#include <complex>
#include "globalparams.hpp"

using namespace std;


class timebin
{

public:
  
  // Constructor
  timebin();
  timebin(int);
  
  // Destructor
  virtual ~timebin() {};
  

  // member data
  complex<double> *Sx, *Sy, *Sz;
  int size;      //variable size of bin. Same size as number of time steps.
  
  //member functions
  void zero();
  void dumpSx(complex<double> val, int pos);
  void dumpSy(complex<double> val, int pos);
  void dumpSz(complex<double> val, int pos);
  void average(int num);
  
};


#endif
