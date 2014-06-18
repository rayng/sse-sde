class bin
{
  
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
// The trajectory class
----------------------------------------------------- */

class trajectory
{
 public:
  complex<double> z[N], zp[N];
  
  void evolve();
  void calcobservables();
  void init();
  void write_variables(ofstream& fp, int t);
  
};


void trajectory::evolve()  // evolves variables forward by dt
{
  int i,l,m, t, nIter=3;
  complex<double> R[N], c[N], Tn[N], s[N], cp[N], sp[N], zm[N], zmp[N], Rm[N];
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
      
      // Noise terms:
      eta[i].real() = 1./(sqrt(dt))*normdist(r);
      eta[i].imag() = 1./(sqrt(dt))*normdist(r);
      zeta[i].real() = 1./(sqrt(dt))*normdist(r);
      zeta[i].imag() = 1./(sqrt(dt))*normdist(r);
    }
  
  
  // Iterate for mid point value
  for(t=0;t<nIter;t++)  
    {
      for(i=0;i<N;i++)
	{
	  l=i-1; 
	  m=i+1;
	  
	  if(l<0) l=N-1;
	  if(m>=N) m=0;
	  

	  // SxiSxj
#if(polar_x)
	  drift= -0.5*I*s[i]*( c[m] - s[m]*Tn[m] + c[l] - s[l]*Tn[l] );
	  driftpr= 0.5*I*sp[i]*( cp[m] - sp[m]*Tn[m] + cp[l] - sp[l]*Tn[l] );
	  noise=sqrt(s[i]*s[m])*conj(eta[i])  + sqrt(s[i]*s[l])*(eta[l]); 
	  noisepr=sqrt(sp[i]*sp[m])*conj(zeta[i])  + sqrt(sp[i]*sp[l])*(zeta[l]); 
	  strat=0.; stratpr=0.;
	  field = I*hh;
	  fieldpr = -I*hh;
#endif


#if(polar_z)
	  drift= 0.5*I*(Tn[m] + Tn[l]);
	  driftpr= -0.5*I*(Tn[m] + Tn[l]);
	  noise=eta[i]  + I*conj(eta[l]);
	  noisepr=zeta[i]  - I*conj(zeta[l]);
	  //field = I*hh*sinh(z[i]);
	  //fieldpr = -I*hh*sinh(zp[i]);
#endif
	  
	  zm[i] = z[i] + (drift + strat + field  + noise)*0.5*dt;
	  zmp[i] = zp[i] + (driftpr + stratpr + fieldpr + noisepr )*0.5*dt;
	  
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
      drift= -0.5*I*s[i]*( c[m]  - s[m]*Tn[m] + c[l] - s[l]*Tn[l] );
      driftpr= 0.5*I*sp[i]*( cp[m] - sp[m]*Tn[m] + cp[l] - sp[l]*Tn[l] );
      //strat = 0.25*I*sinh(z[i])*(cosh(z[i]) + cosh(z[i]));
      //stratpr = -0.25*I*sinh(zp[i])*(cosh(zp[i]) + cosh(zp[i]));
      //strat = 0.25*I*s[i]* (c[l] + c[m]);
      //stratpr = -0.25*I*sp[i]* (cp[l] + cp[m]);
      strat=0.; stratpr=0.;
      noise=sqrt(s[i]*s[m])*conj(eta[i])  + sqrt(s[i]*s[l])*(eta[l]); 
      noisepr=sqrt(sp[i]*sp[m])*conj(zeta[i])  + sqrt(sp[i]*sp[l])*(zeta[l]); 
      field=I*hh;
      fieldpr=-I*hh;
#endif
      
#if(polar_z)
      drift= 0.5*I*(Tn[m] + Tn[l]);
      driftpr= -0.5*I*(Tn[m] + Tn[l]);
      noise=eta[i]  + I*conj(eta[l]);
      noisepr=zeta[i]  - I*conj(zeta[l]);
      //field = I*hh*sinh(z[i]);
      //fieldpr = -I*hh*sinh(zp[i]);
#endif
      
      
      z[i] = z[i] + (drift + strat + field  + noise)*dt;
      zp[i] = zp[i] + (driftpr + stratpr + fieldpr + noisepr )*dt;
      
    }
  
  
  
  
  
  
}


void trajectory::write_variables(ofstream& fp, int t)
{
  int i;
  fp << t*dt  ;
  for(i=0; i<N; i++)
    fp << " " <<  z[i].real()  << " " <<  z[i].imag()  
       << " " <<  zp[i].real() << " " << zp[i].imag() ;
  
}





