//#define nobs 14
#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <Eigen/Dense>
//#include <boost/filesystem.hpp> 

#include <RandomLib/Random.hpp>
#include <RandomLib/NormalDistribution.hpp>

#include <sys/stat.h>
#include <sys/types.h>
#include <iomanip>

#include "../lib/globalparams.hpp"
#include "../lib/bin.hpp"
#include "../lib/timebin.hpp"
#include "../lib/trajectory.hpp"
#include "../lib/lattice.hpp"
#include "../lib/optimise.h"


using namespace std;
using namespace optimise;
//using namespace Eigen;

RandomLib::Random r;
RandomLib::NormalDistribution<double> normdist;

int psamples;
int  msamples;
double beta=2.*Lx;
double  hb;
double  epsilonn[Nb]; 
double tlast=0.;
double dt=0.001;
double dto=0.001;
ofstream fp1, fp2, ferr;
int GLOBAL_RANK;
int GLOBAL_REAL;
int batch;
int h_iter;
int n_per_interval;
string dir="../";

void SETUP_FILE_STREAMS();
string convertInt (int t);
string convertDouble (double t);
int convertString(const string &phrase);
double convertStringD(const string &phrase);


int main(int argc, char* argv[])
{
  int i, num;
  lattice alpha;  
  int Mo=20;
  int nHo=0;
  int Nlo=0; 
  bin obs;
  double t;
  
  h_iter=convertString(argv[1]);     
  batch = convertString(argv[2]);
  hb = hbi+h_iter*dh;                  // field values
  n_per_interval=1;
  
  SETUP_FILE_STREAMS();
  r.Reseed();                // Every batch must get a new seed.   
  
  alpha.init_state="hightemp";    // random initial state
  alpha.set_M(Mo);
  alpha.set_nH(nHo);
  
  
  alpha.buildlattice();
  alpha.generate_probtable_simplex();
  alpha.check_probtable(ferr);
  alpha.init_config("hightemp");
  //alpha.print_probtable(fp);
  
  alpha.Nl_sum=0;
  alpha.n_sum=0;
  
  for(i=0; i<4*Ne; i++)
    {
      alpha.diagonal_update();
      alpha.adjust_truncate();
      alpha.loop_update_eqm(i,ferr);
    }
  alpha.set_Nl( (alpha.Nl_sum/(4*Ne)) );
  
  
  for(num=0; num<samples; num++)
      {
	psamples=0;       	
	msamples=0;            // Note that when taking parity into consideration 
	                       // Nm is such that Nm = msamples + psamples
	for( i=0; i < Nm; i++) // Most of the sampling is done here. 
	  {
	    alpha.diagonal_update();
	    alpha.loop_update(i,ferr);
#if(measureobs)	    
	    alpha.measure();
#endif
	  } 
	
	
#if(measureobs)
	assert(Nm==psamples+msamples);
	alpha.XPS.parityaverage(Nm);
	obs+=alpha.XPS;  // dump a bin
	alpha.XPS.zero();
#endif

#if(measureobs)
  obs.average(samples);
  obs.write_to_file(fp1, alpha.epsilonn);
  fp1 << endl;
  fp1.flush();
#endif
  alpha.free_arrays();
  

  fp1.close();
  fp2.close();
  ferr.close();
  
  return 0;
}
  
  
  
void SETUP_FILE_STREAMS()
{
  string sys = convertInt(Lx)+"x"+convertInt(Ly);
  // Name of file:
  string infile = dir+"data/"+sys+"/DIST/TIM-"+sys+"h"+convertInt(h_iter)+"p"+convertInt(batch)+"beta"+convertDouble(beta);
  string errfile= dir+"data/"+sys+"/ERR/ERR_SDE-TIM-h"+convertInt(h_iter)+"p"+convertInt(batch)+"beta"+convertDouble(beta);
  string sdefile= dir+"data/"+sys+"/SDE/SDE-TIM-h"+convertInt(h_iter)+"p"+convertInt(batch)+"beta"+convertDouble(beta);
  
  int ret=-1;
  struct stat buff;
  ret=stat(infile.c_str(), &buff );
  // This will test if the file exists..
  // If it does not then create it and start fresh. 
  
#if(!append)
  fp1.open(infile.c_str(), ios::out);
  ferr.open ( errfile.c_str(), ios:: out  );  
#else
  if(ret!=0) {
    fp1.open(infile.c_str(), ios::out);
    ferr.open ( errfile.c_str(), ios:: out  );
  }
  // If the file exists then we will append to it. 
  else if (ret==0) {
    fp1.open(infile.c_str(), ios::app);
    ferr.open ( errfile.c_str(), ios:: out  );
  }
#endif
  
  fp2.open(sdefile.c_str(), ios::out);
  if(fp2==NULL)
    cout << "Cannot open SDE file" <<endl;
  
#if saveconfig
  fp2.open  ( (dir+"CONFIG/init"+convertInt(Lx)+"x"+convertInt(Ly)+"h"+convertInt(h_iter)+"p"+convertInt(batch)+"beta"+convertDouble(beta)).c_str(), ios:: out );
  assert(fp2!=NULL);
#endif
  
  assert(ferr!=NULL);
  assert(fp1!=NULL);
  
}



string convertInt (int t)
{
  stringstream ss;
  ss << t;
  return ss.str();
}

// =======================================================================================================


string convertDouble (double t)
{
  stringstream ss;
  ss << t;
  return ss.str();
}

// =======================================================================================================


int convertString(const string &phrase)
{
  stringstream ss(phrase);
  int result;
  return ss >> result ? result : 0;
}

// =======================================================================================================


double convertStringD(const string &phrase)
{
  stringstream ss(phrase);
  double result;
  return ss >> result ? result : 0;
}

// =======================================================================================================


