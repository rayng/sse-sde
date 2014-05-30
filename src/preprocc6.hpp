//#include "mpi.h"
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
#include <Eigen/Dense>
//#include <boost/filesystem.hpp>
#include <sys/stat.h>
#include <sys/types.h>
#include <iomanip>



// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// ----------------------------------------------------------
//
// THIS SPECIFIES THE PHYSICAL MODEL
// AND THE SIMULATION PARAMETERS.
//
// ----------------------------------------------------------
// ----------------------------------------------------------


#define loadfield 0
#define FS 0
#define bounceA 1  // This is the old way of calculating bounce weights
#define saveconfig 0
#define WUsweeps 5000
#define Nbins 50
#define MCsweeps 100
#define append 0


#define system "uniform"
#define Lx 8  /***/
#define Ly 8  /***/
#define Lz 1  
#define d 2
#define N  (Lx*Ly*Lz)
#define Nb (d*N)
#define nV 8
#define Pi 3.145159265
#define Delta 0.
#define nobs 14

// These are the key parameters
#define no_real 50000        /***/
#define nsteps 40             /***/
#define Np 1                /***/
#define hbi (1. /(2.*d))
#define hbf (3. /(2.*d))    
//#define dh  ((hbf- hbi) / (double) nsteps )
#define dh  (0.25 / (2.*d) )

// ------------------------
// beta doubling parameters
// ------------------------
#define nmax 10
#define Ne 256
#define Nm (2*Ne)
#define samples 1000

// ------------------------
// Simplex parameters
// ------------------------
#define  MMAX  1000
#define  NMAX  1000
#define  REAL  double
#define  EPSMIN 0.25
#define  zerotol 1e-13

typedef REAL MAT[MMAX][NMAX];

MAT  A;
int  IPOSV[MMAX], IZROV[NMAX];
int  i,j,ICASE;
int  NN = nV*16+1;         // Number of variables.
int  MM;                   // Total number of constraints
int  M1 =  0 ;             // Number of '<=' constraints
int  M2 =  1 ;             // Number of '>=' constraints
int  M3 =  nV*6 + nV*4;    // Number of '=' constraints
REAL R;




// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// -------------------------------------------------------------------------------------------
//
// GLOBAL PARAMETERS/ ARRAYS
//
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

using namespace std;
using namespace Eigen;

RandomLib::Random r;
//double beta, T, hb; 
double hb, beta;

// ------------------------
// Variable sized arrays
// ------------------------
int **bsite;    // Actualy there is no reason to malloc this bsite array

// ------------------------
// Static Arrays
// ------------------------
double W[Nb][nV], eps_bool[nV], epsilonn[Nb];
double h[N];
// -----------------------------------
// Global parameters for SSE algorithm
// -----------------------------------

long int n_sum=0, Nl_sum;  // number of loops
double prob[Nb][4][4][nV], prob_keynumbers[Nb][4][4][nV][2];   
int legspin[nV][4];
int optype[nV];
int phi[N];
int vlist[nV][4];
int vtx_id[3][3][3][3]; // we don't use element [1].
int vtx_flip[nV][4][4]; // one half of the array is redundant really

// ------------------------
// Output streams
// ------------------------
ofstream fp1,fp2, fp3, feq, ferr, fprobout;
ifstream frd, feps, flog, fprobin;


// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// -------------------------------------------------------------------------------------------
//
// Define OBJECTS/CLASSES.
//
//
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

#include "../lib/stats.h"
#include "../lib/myclasses.h"



// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// -------------------------------------------------------------------------------------------
//
// SSE SUBROUTINES OR FUNCTIONS
//
//
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

// ------------------------
// QMC SSE subroutines
// ------------------------
void BUILD_LATTICE();
void INIT_CONFIG(lattice& a);
void DIAGONAL_UPDATE(lattice& a);
void LOOP_UPDATE(lattice & a);
void LOOP_UPDATE_EQM(lattice & a, int nsamples);
void CONSTRUCT_LINKLIST(lattice& a);
void ADJUST_TRUNCATE(lattice& a);
void GENERATE_PROBTABLE_SIMPLEX();
void SETUP_FILE_STREAMS();
void MEASURE(lattice& a);
void EXTEND_ARRAYS(lattice&a);
void CHECK_PROBTABLE();
void WRITE_HEADER(ofstream& fp);
void FREE_ARRAYS(lattice& a);

// ------------------------
// SIMPLEX METHOD SUBROUT-
// INES
// ------------------------
void simplx(MAT a,int m,int n,int m1,int m2,int m3,int *icase,int *izrov, int *iposv);
void simp1(MAT,int,int *,int,int,int *,REAL *);
void simp2(MAT,int,int,int *,int,int *,int,REAL *);
void simp3(MAT,int,int,int,int);


// ------------------------
// PRINTING SUBROUTINES
// ------------------------
void PRINT_BONDLIST(ofstream& fp);
void PRINT_SPINLIST(lattice& a, ofstream& fp);
void PRINT_OPSTRING(ofstream& fp);
void PRINT_VTXLIST(ofstream& fp);
void PRINT_PROBTABLE(ofstream& fp);
void PRINT_WEIGHTS(ofstream& fp);
void PRINT_ARRAYTABLE(ofstream& fp);
void PRINT_LINKLIST(ofstream& fp);
void PRINT_STAG(ofstream& fp);
void PRINT_ARRAYTABLE2(ofstream& fp);
void PRINT_WEIGHTS2(ofstream& fp);
void PRINT_PROBTABLE2(ofstream& fp);
void PRINT_VLIST(ofstream &fp);
  
// ------------------------
// SSE/ MATH FUNCTIONS
// ------------------------
double CALC_MAG(int* s);
bool is_even(int i);
string convertInt (int t);
string convertDouble (double t);
int convertString(const string &phrase);
double convertStringD(const string &phrase);
string EQ_FILENAME, NAME, FILENAME, ERROR, FOLDER;




