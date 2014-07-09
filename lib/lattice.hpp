#ifndef lattice_HPP
#define lattice_HPP

#include <fstream>
#include <string>
#include <iomanip>
#include <complex>
#include <math.h>
#include <iostream>
#include <assert.h>


#include <Eigen/Dense>
#include <RandomLib/Random.hpp>
#include <RandomLib/NormalDistribution.hpp>


#include "globalparams.hpp"
#include "bin.hpp"
#include "trajectory.hpp"

#include "optimise.h"

using std::cout;
using std::endl;
using std::vector;
using namespace optimise;
using namespace Eigen;


extern RandomLib::Random r;
extern RandomLib::NormalDistribution<double> normdist;
extern int GLOBAL_REAL;
extern int GLOBAL_RANK;


class lattice
{
  
public:
  //member data
  long int M;                              // imaginary time direction
  long int nH;                             // number of operators
  long int Nl;                             // number of loops
  string init_state;                  // hightemp, FMup, FMdown
  int **bsite;
  int spin[N];
  
  
private:
  int first[N], last[N];     // lattice configuration
  int *sm;                            // initialised in DIAGONAL UPDATE
  int *linklist;                      // modifed in    
  int *vtx;
  int *vb;
  
public:
  
  bin XPS;                            // I just called this bin object: XPS
  
  long int Nl_sum=0;
  long int n_sum=0;
  double epsilonn[Nb];
  double h[N];
  int phi[N];
  double W[Nb][nV]; 
  
private:
  
  // internal loop variables
  double eps_bool[nV];
  double prob[Nb][4][4][nV];
  double prob_keynumbers[Nb][4][4][nV][2];   
  int legspin[nV][4];
  int optype[nV];
  int vlist[nV][4];
  int vtx_id[3][3][3][3]; // we don't use element [1].
  int vtx_flip[nV][4][4]; // one half of the array is redundant really
  
public:  
  // member functions
  void update_spin_config();
  void diagonal_update();
  void adjust_truncate();
  void construct_linklist();
  void extend_arrays();
  void loop_update(int nsamples, ofstream& ferr);
  void loop_update_eqm(int nsamples, ofstream& ferr);
  void measure();
  void init_config(string state);
  void write_to_file(ofstream& fp);
  void free_arrays();
  void buildlattice();
  void generate_probtable_simplex();

  // Diagonostics
  void check_probtable(ofstream& ferr);
  void print_probtable(ofstream& fp);
  
  // stochastic trajectories
  trajectory traj;
  void init_traj_z();
  void init_traj_y();
  
  // default constructor
  lattice() {   
    XPS.zero();
  }
  long int set_M( long int M_){
    return M=M_;
  }
  long int set_nH(long int n_){
    return nH=n_;
  }
  
  long int set_Nl(long int Nl_){
    return Nl=Nl_;
  }

  string set_init_state(string init){
    return init_state = init;
  }
  
};

#endif
