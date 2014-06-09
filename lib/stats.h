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
	if(i==13 || i==14 || i==17 || i==18)
	  obs[i] = obs[i] / (double) psamples ;
	
	// p=-1
	else if(i==15 || i==16 || i==19 || i==20)
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
      
      for (i=0; i<nobs-1; i++)  // ignore the fidelity susceptibility- 13 observables
	if (i==1)
	  fp << setprecision(15) << -obs[i]/(beta*Nb) + deltaE << " " ;
	else 
	  fp << setprecision(15) << obs[i] << " " ;
    }
  

};

