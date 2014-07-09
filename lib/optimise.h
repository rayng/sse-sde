#ifndef opt_H
#define opt_H

#include "globalparams.hpp"
#include <stdio.h>
#include <math.h>

namespace optimise
{
  void simplx(MAT a,int m,int n,int m1,int m2,int m3,int *icase,int *izrov, int *iposv);
  void simp1(MAT a,int mm,int *ll,int nll,int iabf,int *kp,REAL *bmax);
  void simp2(MAT a, int m,int n,int *l2,int nl2,int *ip,int kp,REAL *q1);
  void simp3(MAT a,int i1,int k1,int ip,int kp);
}

#endif
