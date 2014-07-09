#ifndef global_HPP
#define global_HPP

// MC toggles
#define loadfield 0
#define FS 1
#define bounceA 1  // This is the old way of calculating bounce weights
#define saveconfig 0
#define append 0
#define parityinclusion 1
#define measureobs 1   // Don't need to measure things here. 

//SDE toggles
#define timeadaptive 1
#define sdeevolve 0
#define driftterms 1
#define noiseterms 1
#define stratterms 0
#define fieldterms 1
#define jumptol 0.08

// Model selection for evolution
#define polar_x 0
#define polar_z 1
#define XYmodel 0

//#define WUsweeps 5000
//#define Nbins 50
//#define MCsweeps 100

// System size definition
#define system "uniform"
#define Lx 12  /***/
#define Ly 12 /***/
#define Lz 1  
#define d 2
#define N (Lx*Ly*Lz)
#define Nb (d*N)
#define nV 8
#define Pi 3.145159265
#define Delta 0.
#define nobs 24

// These are the key parameters
#define no_real 50000          // disorder realizations
#define nsteps 16              // field values
#define Np 1                   // size of batch
#define hbi (1.35 /(2.*d))     // initial field
#define hbf (1.6 /(2.*d))      // final field
#define dh  ((hbf- hbi) / nsteps )   // field step

// ------------------------
// beta doubling parameters
// ------------------------
#define nmax 10
#define Ne 128
#define Nm (2*Ne)
#define samples 500

// Number of time steps
#define nT 3000

// ------------------------
// Simplex parameters
// ------------------------
#define  MMAX  1000
#define  NMAX  1000
#define  REAL  double
#define  EPSMIN 0.2
#define  zerotol 1e-13

typedef REAL MAT[MMAX][NMAX];

extern int  msamples, psamples;
extern double beta, hb, dt, dto, tlast;

#endif

