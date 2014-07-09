//#define nobs 14
#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include <RandomLib/Random.hpp>
#include <RandomLib/NormalDistribution.hpp>
#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>
#include <string>
#include <sstream>
#include <assert.h>
//#include <Eigen/Dense>
//#include <boost/filesystem.hpp> 
#include <sys/stat.h>
#include <sys/types.h>
#include <iomanip>

using namespace std;

#include "globalparams.hpp"
#include "bin.hpp"
#include "timebin.hpp"
#include "trajectory.hpp"
#include "lattice.hpp"


#include <RandomLib/Random.hpp>
#include <RandomLib/NormalDistribution.hpp>

//#include "bin.hpp"
using namespace std;

RandomLib::Random r;
RandomLib::NormalDistribution<double> normdist;
int psamples, msamples;
double beta, hb=2., epsilonn[Nb]; 
ofstream fp;
double tlast=0., dt=0.001, dto=0.001;

int main()
{
bin A;
trajectory worm;

timebin C;
  for(int i=0;i<10; i++)
    //cout << C.Sx[i] << endl;
    cout << worm.tbin.Sx[i] << endl;
  
  
  fp.open("output", ios::out);
  
  worm.init_y();
  worm.write_variables(fp, 0, "y");
  for(int t=1; t<100; t++)
    {
      tlast+=dt;
      worm.evolve_y();
      worm.write_variables(fp, t, "y");
    }
  
  //A.write_to_file(fp);
  
  fp.close();
  
  
  return 0;
}
