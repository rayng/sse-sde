#include "timebin.hpp"


timebin::timebin()
{
  size=1000;
  Sx = new complex<double> [size];
  Sy = new complex<double> [size];
  Sz = new complex<double> [size];
  zero();
}


timebin::timebin(int n )
{
  size=n;
  // allocate arrays
  Sx = new complex<double> [size];
  Sy = new complex<double> [size];
  Sz = new complex<double> [size];
  zero();
}



void timebin::zero()
{
  for(int i=0; i<size; i++)
    {
      Sx[i]=0.;
      Sy[i]=0.;
      Sz[i]=0.;
    }
}

void timebin::average(int num)
{
  for(int i=0; i<size; i++)
    {
      Sx[i]/=num;
      Sy[i]/=num;
      Sz[i]/=num;
    }
}
void timebin::dumpSx(complex<double> val, int pos)
{
  Sx[pos] = val;
}

void timebin::dumpSy(complex<double> val, int pos)
{
  Sy[pos] = val;
}

void timebin::dumpSz(complex<double> val, int pos)
{
  Sz[pos] = val;
}
