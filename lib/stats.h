class bin
{
  // modify size of of bin

public:
  long double obs[nobs];
  long double physobs[nobs];   // hard coded =/
  
  
  bin() 
    {
      zero();
    }
  
  void operator= ( bin &rhs)
    {
      for(int i=0; i< nobs; i++)
	{
	  obs[i] = rhs.obs[i];
	  //physobs[i] = rhs.physobs[i];
	}
    }
  

  void operator+= ( bin &rhs)
  {
    for(int i=0; i< nobs; i++)
      {
	obs[i] += rhs.obs[i];
	//physobs[i] = rhs.physobs[i];
      }
  }

  void dump(long double O, int i )
  {
    obs[i] += O;
  }
  
  void average(int nn)
  {
    for(int i=0; i<nobs; i++)
      obs[i] = obs[i] / (double) nn;
  }

  void parityaverage(int nn)
  {
    for(int i=0; i<nobs; i++)
      {
	//p=+1
	if(i==13 || i==14 || i==17 || i==18 || i==22)
	  obs[i] = obs[i] / (double) psamples ;
	
	// p=-1
	else if(i==15 || i==16 || i==19 || i==20 || i==23 )
	  obs[i] = obs[i] / (double) msamples ;
	
	else
	  obs[i]=obs[i]/ (double) nn;
      }
    
  }
  
  void zero()
  {
    for(int i=0; i<nobs; i++) {
      obs[i] = 0.;
      physobs[i] =0.;
    }
    
  }
  
  // this only works if the average is normalized.
  void generate_physical_obs()
  {
    double hh=2.*d*hb;
    physobs[0] = obs[6];                                                      // stiffness
    physobs[1] = beta*N*( obs[4] - obs[3]*obs[3]  );                          // compressibility    
    physobs[2] = (obs[8]  - obs[7]*obs[7] - obs[7] ) / (beta *hh*hh *N) ;     // disorder suscpetiblity
    physobs[3] = (obs[10] - obs[9]*obs[9] - obs[9] ) / (beta *hh*hh *N) ;     // hopping susceptibility
  }



  void write_to_file(ofstream& fp)
    {
      int i;
      double deltaE=0.;
      // I want to include that additive constant for the energy as well:
      for(int b=0; b<Nb; b++)
	deltaE += epsilonn[b] + hb;
      deltaE=deltaE/ Nb;
      
      for (i=0; i<nobs; i++)  // ignore the fidelity susceptibility- 13 observables
	if (i==1)
	  fp << setprecision(15) << -obs[i]/(beta*Nb) + deltaE << " " ;
	else 
	  fp << setprecision(15) << obs[i] << " " ;
    }
  

};


/* -----------------------------------------------------
// Time bin
----------------------------------------------------- */
class timebin
{
 private:
  
  //member functions
  void zero();
  void dumpSx(complex<double> val, int pos);
  void dumpSy(complex<double> val, int pos);
  void dumpSz(complex<double> val, int pos);
  void average(int num);
  
  
 public:
  // member data
  complex<double> *Sx, *Sy, *Sz;
  int size;
  
  // Constructor
  timebin();
  timebin(int n);
  
  // Destructor
  virtual ~timebin();

};

// Destuctor
timebin::~timebin()
{
}

timebin::timebin(int n)
{
  size=n;
  
  // allocate arrays
  Sx = new complex<double> [size];
  Sy = new complex<double> [size];
  Sz = new complex<double> [size];
  zero();
}


timebin::timebin()
{
  size=1000;
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

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





/* -----------------------------------------------------
// The trajectory class
----------------------------------------------------- */
class trajectory
{
 public:
  complex<double> z[N], zp[N];
  complex<double> y[N], yp[N];
  complex<double> zetafcn[N];
  
  int SS[N], SSp[N];
  
  void evolve();
  void evolve_y();
  void calcobservables();
  
  void init();
  void init_y();
  void write_variables(ofstream& fp, int t, string var);   // var is either z or y 
  void genzfromy(bool write, ofstream& fp);
  
 private:  

  void calc_zetafcn(complex<double>* y, complex<double>* yp);
  
};


// ++++++++++++++++++++++++++++++
// z-variables 
// ++++++++++++++++++++++++++++++
void trajectory::evolve()  // evolves variables forward by dt  
{
  int i,l,m, t, nIter=3;
  complex<double> R[N], c[N], Tn[N], s[N], cp[N], sp[N], zm[N], zmp[N], Rm[N]; 
  complex<double> sdiffp[N], sdiffm[N], cdiffp[N], cdiffm[N], sdiffppr[N], sdiffmpr[N], cdiffppr[N], cdiffmpr[N] ;
  complex<double> eta[N], zeta[N];
  complex<double> drift=0., driftpr=0., strat=0., stratpr=0., field=0., fieldpr=0., noise, noisepr;
  complex<double> I (0,1.);
  double hh=2.*d*hb;
  
  for(i=0;i<N;i++)
    {
      // expensive functions
      R[i]   = 0.5*(z[i]+zp[i]);
      Rm[i] = R[i];
      zm[i]  = z[i];  
      zmp[i] = zp[i];
      Tn[i]  = tanh(R[i]);
      s[i]   = sinh(z[i]); 
      sp[i]  = sinh(zp[i]);
      c[i]   = cosh(z[i]); 
      cp[i]  = cosh(zp[i]);
      
      m=i+1;
      l=i-1;
      
      if(m>=N)
	m=0;
      if(l<0)
	l=N-1;
      
      // For XY Hamiltonian
      sdiffp[i] = sinh(z[i] - z[m]);
      sdiffm[i] = sinh(z[i] - z[l]);
      cdiffp[i] = cosh(z[i] - z[m]);
      cdiffm[i] = cosh(z[i] - z[l]);
      
      sdiffppr[i] = sinh(zp[i] - zp[m]);
      sdiffmpr[i] = sinh(zp[i] - zp[l]);
      cdiffppr[i] = cosh(zp[i] - zp[m]);
      cdiffmpr[i] = cosh(zp[i] - zp[l]);
      
      // Noise terms:
      eta[i].real() = 1./(sqrt(dt))*normdist(r);
      eta[i].imag() = 1./(sqrt(dt))*normdist(r);
      zeta[i].real() = 1./(sqrt(dt))*normdist(r);
      zeta[i].imag() = 1./(sqrt(dt))*normdist(r);
    }
  
  // Iterate for mid point value
  dt=dto; // reset time step
  
  for(t=0;t<nIter;t++)  
    {
      for(i=0;i<N;i++)
	{
	  l=i-1;  m=i+1;
	  
	  if(l<0)
	    l=N-1;
	  if(m>=N) 
	    m=0;
	  
	  
	  // ---------------------
	  // SxiSxj  hamiltonian
	  // ---------------------
#if(polar_x)
	  
#if (driftterms)
	  drift= -0.5*I*s[i]*( c[m] - s[m]*Tn[m] + c[l] - s[l]*Tn[l] );
	  driftpr= 0.5*I*sp[i]*( cp[m] - sp[m]*Tn[m] + cp[l] - sp[l]*Tn[l] );
#endif

#if (noiseterms)	  
	  noise=sqrt(s[i]*s[m])*conj(eta[i])  + sqrt(s[i]*s[l])*(eta[l]); 
	  noisepr=sqrt(sp[i]*sp[m])*conj(zeta[i])  + sqrt(sp[i]*sp[l])*(zeta[l]); 
#endif
	  
#if(stratterms)
	  strat = 0.25*I*s[i]* (c[l] + c[m]);
	  stratpr = -0.25*I*sp[i]* (cp[l] + cp[m]);
#endif
	  
#if (fieldterms)	  
	  field = I*hh;
	  fieldpr = -I*hh;
#endif
	  
#endif
	  
	  // ---------------------
	  // SziSzj  hamiltonian
	  // ---------------------
#if(polar_z)
	  
#if (driftterms)
	  drift= 0.5*I*(Tn[m] + Tn[l]);
	  driftpr= -0.5*I*(Tn[m] + Tn[l]);
#endif
	  
#if (noiseterms)
	  noise=eta[i]  + I*conj(eta[l]);
	  noisepr=zeta[i]  - I*conj(zeta[l]);
#endif	  
	  
#if(fieldterms)
	  field = I*hh*sinh(z[i]);
	  fieldpr = -I*hh*sinh(zp[i]);
#endif
	  
#endif
	  
	  // ---------------------------
	  // XY model  hamiltonian
	  // ---------------------------
	  
#if(XYmodel)
	  
#if (driftterms)
	  drift= -I*0.5*( sdiffp[i] + cdiffp[i]*Tn[m]  + sdiffm[i] + cdiffm[i]*Tn[l] ) ;  //p=plus, m=minus
	  driftpr=  -I*0.5*( sdiffppr[i] + cdiffppr[i]*Tn[m]  + sdiffmpr[i] + cdiffmpr[i]*Tn[l] ) ;  //p=plus, m=minus
#endif
	  
#if(noiseterms)
	  noise= sqrt(cdiffp[i]/2.)*conj(eta[i]) - I* sqrt(cdiffm[i]/2.)*eta[l];
	  noisepr= sqrt(cdiffppr[i]/2.)*conj(zeta[i]) - I* sqrt(cdiffmpr[i]/2.)*zeta[l];
#endif	  
	  
#if(fieldterms)
	  field = 0.;
	  fieldpr = 0.;
#endif
	  
#endif	  
	  
#if(timeadaptive)
	  complex<double> Rm, R;
	  complex<double> Sm, S;
	  bool small_enough=false;
	  while(!small_enough)
	    {
	      zm[i] = z[i] + (drift + strat + field  + noise)*0.5*dt;
	      zmp[i] = zp[i] + (driftpr + stratpr + fieldpr + noisepr )*0.5*dt;
	      
	      // z variables spiking
	      //if (  ((zm[i]-z[i])/ z[i]).real() > jumptol ||
	      //    ((zm[i]-z[i])/ z[i]).imag() > jumptol )
	      // Check R and S variables spking
	      Rm=(zm[i] +zmp[i])/2.;
	      R=(z[i] +zp[i])/2.;
	      Sm=(zm[i]-zmp[i])/(2.*I);
	      S=(z[i]-zp[i])/(2.*I);
	      
	      if (  ((Rm-R)/ R).real() > jumptol ||
		    ((Sm-S)/S).imag() > jumptol )
		dt=dt*0.5;  // halve it
	      else 
		small_enough =true;
	    }
#else
	  zm[i] = z[i] + (drift + strat + field  + noise)*0.5*dt;
	  zmp[i] = zp[i] + (driftpr + stratpr + fieldpr + noisepr )*0.5*dt;
#endif
	  
	}
      //update using zm, zmp values
      for(i=0; i<N; i++)
	{
	  R[i] = 0.5*(zm[i]+zmp[i]);
	  Tn[i] = tanh(R[i]);
	  s[i] = sinh(zm[i]);
	  c[i] = cosh(zm[i]);
	  cp[i] = cosh(zmp[i]);
	  sp[i] = sinh(zmp[i]);
	  
	  m=i+1;
	  l=i-1;
	  
	  if(m>=N)
	    m=0;
	  if(l<0)
	    l=N-1;
	  
	  // For XY Hamiltonian
	  sdiffp[i] = sinh(zm[i] - zm[m]);
	  sdiffm[i] = sinh(zm[i] - zm[l]);
	  cdiffp[i] = cosh(zm[i] - zm[m]);
	  cdiffm[i] = cosh(zm[i] - zm[l]);
	  
	  sdiffppr[i] = sinh(zmp[i] - zmp[m]);
	  sdiffmpr[i] = sinh(zmp[i] - zmp[l]);
	  cdiffppr[i] = cosh(zmp[i] - zmp[m]);
	  cdiffmpr[i] = cosh(zmp[i] - zmp[l]);
	}
      
      
    }      
  for(i=0;i<N;i++)
    {
      l=i-1; 
      m=i+1;
      
      if(l<0) 
	l=N-1;
      if(m>=N) 
	m=0;
      
#if(polar_x)
      
      
#if (driftterms)
      drift= -0.5*I*s[i]*( c[m]  - s[m]*Tn[m] + c[l] - s[l]*Tn[l] );
      driftpr= 0.5*I*sp[i]*( cp[m] - sp[m]*Tn[m] + cp[l] - sp[l]*Tn[l] );
#endif
      
#if (stratterms)      
      strat = 0.25*I*s[i]* (c[l] + c[m]);
      stratpr = -0.25*I*sp[i]* (cp[l] + cp[m]);
#endif
      
#if (noiseterms)      
      noise=sqrt(s[i]*s[m])*conj(eta[i])  + sqrt(s[i]*s[l])*(eta[l]); 
      noisepr=sqrt(sp[i]*sp[m])*conj(zeta[i])  + sqrt(sp[i]*sp[l])*(zeta[l]); 
#endif
      
#if (fieldterms)      
      field=I*hh;
      fieldpr=-I*hh;
#endif
      
#endif
      
#if(polar_z)
      
#if (driftterms)
      drift= 0.5*I*(Tn[m] + Tn[l]);
      driftpr= -0.5*I*(Tn[m] + Tn[l]);
#endif
      
#if (noiseterms)
      noise=eta[i]  + I*conj(eta[l]);
      noisepr=zeta[i]  - I*conj(zeta[l]);
#endif
      
#if(fieldterms)
      field = I*hh*sinh(z[i]);
      fieldpr = -I*hh*sinh(zp[i]);
#endif
      
#endif
      
      // ---------------------------
      // SxiSxj+SyiSyj  hamiltonian
      // ---------------------------

#if(XYmodel)
      
      
#if (driftterms)
      drift= -I*0.5*( sdiffp[i] + cdiffp[i]*Tn[m]  + sdiffm[i] + cdiffm[i]*Tn[l] ) ;  //p=plus, m=minus
      driftpr=  I*0.5*( sdiffppr[i] + cdiffppr[i]*Tn[m]  + sdiffmpr[i] + cdiffmpr[i]*Tn[l] ) ;  //p=plus, m=minus
#endif
      
      
#if(noiseterms)
      noise= sqrt(cdiffp[i]/2.)*conj(eta[i]) - I* sqrt(cdiffm[i]/2.)*eta[l];
      noisepr= sqrt(cdiffppr[i]/2.)*conj(zeta[i]) - I* sqrt(cdiffmpr[i]/2.)*zeta[l];
#endif	  
      
#if(fieldterms)
      field = 0.;
      fieldpr = 0.;
#endif
      
#endif	  
      z[i] = z[i] + (drift + strat + field  + noise)*dt;
      zp[i] = zp[i] + (driftpr + stratpr + fieldpr + noisepr )*dt;
      
    }
  
}  


void trajectory:: genzfromy(bool write, ofstream& fp)
{
  for(int i=0; i<N; i++)
    {
      if(SS[i] == 1)
	z[i] = -log(y[i]);
      else if(SS[i] == -1)
	z[i] =  log( conj(y[i]) ) ;
      
      if(SSp[i] == 1)
	zp[i] = -log(yp[i]);
      else if(SSp[i]==-1)
	zp[i] =  log( conj(yp[i]) ) ;
    }
  
  if(write)
    {
      fp << tlast  ;
      for(int i=0; i<N; i++)
	fp << setprecision(20) << " " <<  z[i].real()  << " " <<  z[i].imag()  
	   << " " <<  zp[i].real() << " " << zp[i].imag() ;
      fp << endl;  
    }
  
  
}


void trajectory::write_variables(ofstream& fp, int t, string var)
{
  int i;
  
  if(var=="z")
    {
      fp << tlast  ;
      for(i=0; i<N; i++)
	fp << setprecision(20) << " " <<  z[i].real()  << " " <<  z[i].imag()  
	   << " " <<  zp[i].real() << " " << zp[i].imag() ;
      fp<<endl;
    }
  
  if(var=="y")
    genzfromy(true, fp);
  
  
  
}


// ++++++++++++++++++++++++++++++
// y-variables 
// ++++++++++++++++++++++++++++++
void trajectory::evolve_y()  // evolves variables forward by dt  
{
  
  int i,l,m;
  double h=hb*2.*d, tol=1.;
  complex<double> drift=0., strat=0., field=0., noise=0.;
  complex<double> driftpr=0., stratpr=0., fieldpr=0., noisepr=0.;
  complex<double> I (0.,1.);
  
  complex<double> eta[N], zeta[N];
  
  h=2.*h;  //doube the field
  
  // Generate noise
  for(i=0; i<N; i++)
    {
      eta[i].real()=1./sqrt(2.*dt)*normdist(r);
      eta[i].imag()=1./sqrt(2.*dt)*normdist(r);

      zeta[i].real()=1./sqrt(2.*dt)*normdist(r);
      zeta[i].imag()=1./sqrt(2.*dt)*normdist(r);
    }
  
  
  calc_zetafcn(y,yp);
  //h=0.;
  for(i=0; i<N; i++)
    {
      l=i-1;
      m=i+1;
      
      if(l<0)
	l=N-1;
      if(m>=N)
	m=0;
      
      // Non-primed variables
      
      if(SS[i]==1)
	{
	  drift=-I*0.5*y[i]*( zetafcn[l] + zetafcn[m] );
	  field=I*h*0.5*(1. -y[i]*y[i]);
	  noise = y[i]*(eta[i]+I*conj(eta[l]));
	}
      
      else if(SS[i]==-1)
	{
	  drift=-I*0.5*y[i]*( zetafcn[l] + zetafcn[m] );
	  field=-I*h*0.5*(1. -y[i]*y[i]);
	  noise = y[i]*(conj(eta[i])-I*eta[l]);
	}
      
      // primed variables
      if(SSp[i]==1)
	{
	  driftpr=I*0.5*yp[i]*( zetafcn[l] + zetafcn[m] );
	  fieldpr=-I*h*0.5*(1. -yp[i]*yp[i]);
	  noisepr = -yp[i]*(zeta[i]-I*conj(zeta[l]));
	}
      else if(SSp[i]==-1)
	{
	  driftpr=I*0.5*yp[i]*( zetafcn[l] + zetafcn[m] );
	  fieldpr=I*h*0.5*(1. -yp[i]*yp[i]);
	  noisepr = yp[i]*(conj(zeta[i])+I*zeta[l]);
	}
      
      // time step
      y[i] = y[i] + (drift + strat + field  + noise)*dt;
      yp[i] = yp[i] + (driftpr + stratpr + fieldpr + noisepr )*dt;
    }


  // Update sign 
  
  for(i=0; i<N; i++)
    {
      //  Update Sign                                                                          
      if( abs(y[i]) > tol  ) 
	{
	  SS[i]= -SS[i];
	  y[i] = 1./ conj(y[i]);
	}
      
      if( abs(yp[i]) > tol  ) 
	{
	  SSp[i] = -SSp[i];
	  yp[i] = 1./ conj(yp[i]);
	}
    }
  
  
}

// ++++++++++++++++++++++++++++++
// compute the zeta-function
// ++++++++++++++++++++++++++++++
void trajectory::calc_zetafcn(complex<double>* y, complex<double>* yp)
{
  int i;
  for(i=0; i<N; i++)
    if(SS[i]==1 && SSp[i]==1)
      zetafcn[i] = ( 1.-y[i]*yp[i] )/ ( 1.+ y[i]*yp[i])   ;
  
    else if(SS[i]==1 && SSp[i]==-1)
      zetafcn[i] = (   ( conj(yp[i]) - y[i]  ) / ( conj(yp[i]) + y[i]  ) );
  
    else if (SS[i]==-1 && SSp[i] ==1)
      zetafcn[i]= (   ( conj(y[i]) - yp[i]  ) / ( conj(y[i]) + yp[i]  ) );
    
    else if (SS[i]==-1 && SSp[i] ==-1)
      zetafcn[i]= ( conj(y[i])*conj(yp[i]) - 1.  )  / ( conj(y[i])*conj(yp[i]) + 1.  ) ;
    else
      cout << "Problem in zetafcn calculation" << " " << SS[i] << " " << SSp[i] << endl;
  
  
  
}
  
  






