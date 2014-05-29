
//
// Over here I will modify the nuber of loops in the problem. 
// The average number of vertices visited should be roughly 2 <nH>
//
#include "preprocc6.hpp"


int GLOBAL_REAL, GLOBAL_RANK;
string dir="../";    /********  change this  ********/
int batch, h_iter;

// -----------------------------------------------------
//
//  MAIN PROGRAM: LET THE FUN BEGIN.
//
// -----------------------------------------------------
int main(int argc, char* argv[])
{
  int i,j, g, rr, x, n=0, nr=0, nn,bb, n_per_interval, local_start, local_end;
  lattice alpha, beta_;
  // Set up field arrays.
  h_iter=convertString(argv[1]);     // This shift is just to add the extra field values I'm calculating. h=5 => h=4.75 now
  batch = convertString(argv[2]);
  hb = hbi+h_iter*dh;                  // field values
  n_per_interval=1;
  
  // Set up file stream
  SETUP_FILE_STREAMS();
  int Mo=20, no=0, NlA, NlB;
  r.Reseed();                // Every batch must get a new seed.   
  
  for(rr=0; rr<n_per_interval; rr++) {
    alpha.init_state="hightemp";    // random initial state
    alpha.set_M(Mo);     // Initialise M
    alpha.set_nH(no);    // Initialise no
    GLOBAL_REAL=rr;     // realization number
    GLOBAL_RANK=batch;  // batch number
    
    // ---------------------------------
    // INITIALISE
    // ---------------------------------
    BUILD_LATTICE();          
    GENERATE_PROBTABLE_SIMPLEX();        
    CHECK_PROBTABLE();
    //PRINT_PROBTABLE2(ferr);
    // ---------------------------------
    // WARMUP CONFIGS: (initial state is too far from equilibrium)
    // ---------------------------------
    alpha.init_config("hightemp");
    
    fp1 << rr <<  " " ;
    
    for( nn=0; nn< nmax; nn++)
      {
	beta = pow(2., nn);
	for(bb=0; bb<Nb; bb++)            // Renormalize the weight
	  for(i=0; i<nV; i++)
	    W[bb][i]=beta*Nb*W[bb][i];
	
	alpha.extend_arrays();
	
	
	// +++++++++++++++++++++++++++++++++++++++++++++
	//  Replica - A
	// ----------------------
	//  Warm up segment
	
	n_sum=0; Nl_sum=0;
	for(i=0; i<4*Ne; i++) {
	  alpha.diagonal_update();
	  alpha.adjust_truncate();
	  alpha.loop_update_eqm(i, ferr);
	}
	alpha.set_Nl(Nl_sum/(4*Ne));
	
	//  Measure segment - No adjusting during measurement
	for( i=0; i < Nm; i++) // Most of the sampling is done here. 
	  {
	    alpha.diagonal_update();
	    alpha.loop_update(i,ferr);
	    alpha.measure();
	  } 
	alpha.XPS.average(Nm);
	alpha.write_to_file(fp1);
	alpha.XPS.zero();
	
	// +++++++++++++++++++++++++++++++++++++++++++++
	// Normalise the weight or it will keep compounding
	// ----------------------
	for(bb=0; bb<Nb; bb++)    
	  for(i=0; i<nV; i++)
	    W[bb][i]=W[bb][i]/(Nb*beta);
      }
    
    fp1 << endl;
    fp1.flush();
    
    
    alpha.free_arrays();
    
    free(bsite);
  }
  
  fp1.close();
  fp2.close();
  ferr.close();
  return 0;
}
// *********************************** END OF MAIN PROGRAM  **********************************************




// -------------------------------------------------------------------------------------------------------

void BUILD_LATTICE()

// This subroutine builds a square lattice which is reflected via the bond list.
// The routine is versatile enough to handle d=1 (chain), 2 (square) and 3 (cube)
// The array bsite[2][Nb] basically tells me the site numbers: element 0 and 1, associated with the bond 
// Nb.
// In here I will also, initialise the staggered phase.  
// 
// -------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------
{
  int i, x,  y, z, v, l, s[4], e, vnum, inleg, exitleg;
  double hav=0.;
  
  // Calloc-ing memory.   <-- just a good exercise I suppose.
  bsite = (int**) calloc( 2,sizeof(int*) );       // This is global, all lattices will share this. 
  for(i=0; i<2; i++)
    bsite[i] = (int*) calloc(Nb, sizeof(int) );
  
  for(z=0; z<Lz; z++)
    for(y=0; y<Ly; y++)
      for(x=0; x<Lx; x++) {
	i = x + y*Lx + z*(Lx*Ly);
	if(d==1) {
	  // x bonds
	  bsite[0][i] = i  ;
	  bsite[1][i] = (x+1)%Lx + y*Lx + z*Lx*Ly ; 
	}
	  
	if(d==2) {
	  // x bonds
	  bsite[0][i] = i  ;
	  bsite[1][i] = (x+1)%Lx + y*Lx + z*Lx*Ly ; 
	  // y bonds
	  bsite[0][i+N] = i  ;
	  bsite[1][i+N] = x + ( (y+1)*Lx )%(Lx*Ly) + z*Lx*Ly; 
	}
	
	if(d==3) {
	  // x bonds
	  bsite[0][i] = i  ;
	  bsite[1][i] = (x+1)%Lx + y*Lx + z*Lx*Ly ; 
	  // y bonds
	  bsite[0][i+N] = i  ;
	  bsite[1][i+N] = x + ( (y+1)*Lx )%(Lx*Ly) + z*Lx*Ly; 
	  // z bonds
	  bsite[0][i+2*N] = i  ;
	  bsite[1][i+2*N] = x + y*Lx + ( (z+1)*(Lx*Ly) ) % N; 
	}
      }
  
  for(i=0; i<N; i++) {
    x = i%Lx;
    y = i%(Lx*Ly)/Lx;
    z = i/(Lx*Ly);
    phi[i] = pow(-1., (x+y+z));
  }
  
  // this need to be shortened ermagawdddd
  // -------------------------------------------------
  // INITIALISE LEGSPIN     
  // - OPTYPE (HARDCODED)
  // - GLOBAL WEIGHTS
  // - EPSILON-NESS
  // -------------------------------------------------
  // vertex 0  ->    O  X
  //                 ----
  //                 X  O   
  optype[0] = 1;
  eps_bool[0] = false;
  legspin[0][2] = -1;  legspin[0][3] = 1;
  legspin[0][0] = 1;  legspin[0][1] = -1;
  
  // vertex 1  ->    O  X
  //                 ----
  //                 O  X
  optype[1] = 0;
  eps_bool[1] = true;
  legspin[1][2] = -1;  legspin[1][3] = 1;
  legspin[1][0] = -1;  legspin[1][1] = 1;
  
  // vertex 2  ->    X  X
  //                 ----
  //                 O  O
  optype[2] = 0;
  eps_bool[2] = false;
  legspin[2][2] = 1;  legspin[2][3] = 1;
  legspin[2][0] = -1;  legspin[2][1] = -1;
  
  // vertex 3  ->    O  O
  //                 ----
  //                 O  O
  optype[3] = 0;
  eps_bool[3] = true;
  legspin[3][2] = -1;  legspin[3][3] = -1;
  legspin[3][0] = -1;  legspin[3][1] = -1;
  
  
  /* ------------------------------------------------------------ */
  
  // vertex 4  ->    X  O
  //                 ----
  //                 O  X
  optype[4] = 1;
  eps_bool[4] = false;
  legspin[4][2] = 1;  legspin[4][3] = -1;
  legspin[4][0] = -1;  legspin[4][1] = 1;
  
  // vertex 5  ->    X  O
  //                 ----
  //                 X  O
  optype[5] = 0;
  eps_bool[5] = true;
  legspin[5][2] = 1;  legspin[5][3] = -1;
  legspin[5][0] = 1;  legspin[5][1] = -1;
  
  // vertex 6  ->    O  O
  //                 ----
  //                 X  X
  optype[6] = 0;
  eps_bool[6] = false;
  legspin[6][2] = -1;  legspin[6][3] = -1;
  legspin[6][0] = 1;  legspin[6][1] = 1;
  
  // vertex 7  ->    X  X
  //                 ----
  //                 X  X
  optype[7] = 0;
  eps_bool[7] = true;
  legspin[7][2] = 1;  legspin[7][3] = 1;
  legspin[7][0] = 1;  legspin[7][1] = 1;
  
  // -------------------------------------------------
  //  INITIALISE VERTEX_ID ARRAY
  // -------------------------------------------------
  for(v=0; v<nV; v++) {
    for(l=0; l<4; l++)
      s[l] = legspin[v][l];
    vtx_id[s[0]+1][s[1]+1][s[2]+1][s[3]+1]= v ;
  }
  
  // -------------------------------------------------
  //  INITIALISE VERTEX_FLIPPING ARRAY
  // -------------------------------------------------
  for(v=0; v<nV; v++) {
    for(i=0;i<4;i++)
      s[i]=legspin[v][i];
    
    for(i=0; i<4; i++)
      for(e=0; e<4; e++) {
	//flip
	s[i]=-s[i];
	s[e]=-s[e];
	vtx_flip[v][i][e]=vtx_id[s[0]+1][s[1]+1][s[2]+1][s[3]+1];
	
	//undo
	s[i]=-s[i];
	s[e]=-s[e];
      }
  }
  
  // -------------------------------------------------
  // INITIALISE RANDOM FIELD
  // -------------------------------------------------
  // Staggered case
  hav=0.;
  for(i=0; i<N; i++) {
    
#if($system==staggered)
    h[i] = phi[i]*hb;
#endif
    
#if($system==uniform)
    h[i] = hb;
#endif

#if($system==disorder)
    h[i] = 2.*hb*r.FloatN()-hb;
    hav+=h[i];
  }
  hav = hav / (double) N;
  for(i=0; i<N; i++)
    h[i] -= hav;
#endif
  
    // -------------------------------------------------
  // INITIALISE VERTEX LIST
  // -------------------------------------------------
  for(i=0; i<nV; i++)
    {
      vnum=0;
      for(l=0; l<4; l++)
	s[l] = legspin[i][l];
      inleg = 0; // fixed
      for(exitleg=0; exitleg<4; exitleg++)
	{
	  s[inleg] = -s[inleg];
	  s[exitleg] = -s[exitleg];
	  vlist[i][vnum] = vtx_id[s[0]+1][s[1]+1][s[2]+1][s[3]+1];  
	  // UNDO! Since we're using the same array .. 
	  s[inleg] = -s[inleg];
	  s[exitleg] = -s[exitleg];
	  vnum++;
	}
    }
  
  // -------------------------------------------------
  // INITIALISE EPSILONN ARRAY
  // -------------------------------------------------
  for(i=0; i<Nb; i++)
    epsilonn[i] = 0.;
  
  // -------------------------------------------------
  //  INITIALISE PROBABILITY TABLES
  // -------------------------------------------------
  for(l=0; l<Nb; l++)
    for(i=0; i<4; i++)              
      for(e=0; e<4; e++)
	for(v=0; v<nV; v++)
	  prob[l][i][e][v]=-1;
  
}

// =======================================================================================================



// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// -------------------------------------------------------------------------------------------------------

void GENERATE_PROBTABLE_SIMPLEX()

// This will initialise the weights which are now bond dependent via the fields as well as 
// generate the probabilities for each vertex. 
//
//
//
// ------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------
{
  // -------------------------------------------------
  // INITIALISE GLOBAL WEIGHTS
  // -------------------------------------------------
  double hi, hj, sol[NN];
  int i,j,k, m, pos,conjpos,b, row_count;
  
  MM=M1+M2+M3;  // Total number of constraints
  
  for(b=0; b<Nb; b++)
    {
      /*printf("\n");
      printf(" Number of variables in E.F.: %d \n", NN); 
      printf(" Number of <= inequalities..: %d \n", M1); 
      printf(" Number of >= inequalities..: %d \n", M2); 
      printf(" Number of = equalities.....: %d \n", M3); 
      */
      hi = h[bsite[0][b]];
      hj = h[bsite[1][b]];
      
      // We need to add back epsilon later on.
      W[b][0]= 0.5;
      W[b][1]= 0.5*(-hi+hj) + hb;         // + epsilonn[b];
      W[b][2]= 0.;
      W[b][3]= 0.5*(-hi-hj) + hb;         // + epsilonn[b];
      
      W[b][4]= 0.5;
      W[b][5]= 0.5*(hi-hj) + hb;          // + epsilonn[b];
      W[b][6]= 0.;
      W[b][7]= 0.5*(hi+hj) + hb;          // + epsilonn[b];
      
      /* -----------------------------------
	 Let the madness begin ... 
	 Create the tableau and initialise
	 if !
	 ----------------------------------- */
      for (i=1; i<=MM+2; i++)
	for (j=1; j<=NN+1; j++)
	  A[i][j]=0.0;
      
      for(i=0; i<MMAX; i++)
	IPOSV[i] = 0.;
      
      for(i=0; i<NMAX; i++)
	IZROV[i] = 0.;
      
      /* ---------------------------------------------------
	 Input BOUNCE FUNCTION <-- first row of tableau A
	 The function is now F = C - (\sum_i bounce_i)
	 where C>0.
	 This will drive \sum_i bounce_i to be minimized by
	 being as close to zero as possible.
	 -------------------------------------------------- */
      row_count=1;  // We don't use  row 0.
      for(k=0; k<nV; k++)
	for(i=0; i<4; i++)
	  for(j=0; j<4; j++)
	    {
	      pos = k*16 + i*4 + j + 2;   // 2 is the constant shift  
	      
	      if(i==j)
		A[1][pos] = -1.;
	    }
      A[1][1] = EPSMIN;   // Any positive number here will do
      row_count ++;
      
      /* ---------------------------------------------------
	 Input  epsilon constraint.   epsilon >= C >0
	 This ensures that epsilon is some positive number as
	 opposed to being zero. 
	 -------------------------------------------------- */
      A[row_count][nV*16+2] = -1;     // epsilon coefficient
      A[row_count][1] = EPSMIN;          // epsilon constant: always column 1
      row_count++;
      
      /* ---------------------------------------------------
	 Input WEIGHT constraints                       
	 e.g  m=2: a^k_11 + a^k_12 + a^k_13 + a^k_14  = W^k_1          
	 e.g  m=3: a^k_21 + a^k_22 + a^k_23 + a^k_24  = W^k_2 
	 :
	 :
	 etc                                         
	 --------------------------------------------------- */
      bool sign_flip;
      for(k=0; k<nV; k++)
	for(m=0; m<4; m++) 
	  {
	    // Constant value: not coeffcients
	    if( W[b][vlist[k][m]] <0 ) 
	      sign_flip=true;   // means negative
	    else 
	      sign_flip=false;
	    A[row_count][1] = fabs(W[b][vlist[k][m]]);   // W[k][m] // consider epsilon. 
	    
	    //Coefficients of a^k_ij:
	    for(j=0; j<4; j++) 
	      {
		pos = k*16 + m*4 + j + 2;  
		if(sign_flip)
		  A[row_count][pos] = +1.; 
		else if(!sign_flip)
		  A[row_count][pos] = -1.; 
	      }
	    
	    // epsilon value:
	    if( eps_bool[vlist[k][m]] )
	      if(sign_flip)
		A[row_count][nV*16+2] = -1.;   // epsilon coefficient
	      else if(!sign_flip)
		A[row_count][nV*16+2] = 1.;    // epsilon coefficient
	    
	    row_count ++;
	  }
      
      /* ---------------------------------------------------
	 Input DETAILED BALANCE constraints             
	 
	 i.e.  a^k_ij = a^k_ji
	 
	 ---------------------------------------------------*/
      for(k=0; k<nV; k++)
	for(i=0; i<4; i++)
	  for(j=0; j<4; j++)
	    if(i>j)
	      {
		pos = k*16 + i*4 + j + 2;   // 2 is the constant shift  
		conjpos = k*16 + j*4 + i + 2;   // 2 is the constant shift  
		A[row_count][pos]= 1.;
		A[row_count][conjpos] = -1.;
		A[row_count][1] = 0.;
		row_count++;
	      }

      /*
      printf("\n Input Table:\n");
      for (i=1; i<=MM+1; i++) 
	{
	  for (j=1; j<=NN+1; j++)
	    ftest << A[i][j] << " " ;
	  ftest<< endl;
	}
      */
      
      /* ---------------------------------------------------
	 Simplex and print results                       
	 --------------------------------------------------- */
      simplx(A,MM,NN,M1,M2,M3,&ICASE,IZROV,IPOSV);
      
      for(i=0; i<NN; i++)
	sol[i]=0.;
      
      if (ICASE==0) {  //result ok.
	//printf("\n Minimum of bounces = %f\n", A[1][1]-TOL);
	for (i=1; i<=NN; i++) {
	  for (j=1; j<=MM; j++)
	    if (IPOSV[j] == i)  {
	      sol[i-1] = A[j+1][1];
	      //printf("  X%d = %f\n", i, A[j+1][1]);
	      goto e3;
	    }
	  //printf("  X%d = %f\n", i, 0.0);
	e3:;}
      }
      else
	printf(" No solution (error code = %d).\n", ICASE);
      //printf("\n");
      
      /* ---------------------------------------------------
	 Now assign the probability tables..
	 --------------------------------------------------- */
      epsilonn[b] = sol[NN-1];
      
      // Insert back epsilon into the weights
      for(k=0; k<nV; k++) {
	if(eps_bool[k])
	  W[b][k] += epsilonn[b];
	//if( abs(W[b][k]) < zerotol && W[b][k]!=0.) {
	//W[b][k]=0.;
	//cout << "Negative weights!: " << b << " " << k << " " << W[b][k]  << endl;  
	//}
      }
      for(k=0; k<nV; k++)
	for(i=0; i<4; i++)
	  for(j=0; j<4; j++) {
	    pos = k*16 + i*4 + j;    
	    
	    if( sol[pos] > zerotol && W[b][vlist[k][i]] > zerotol )
	      {
		prob_keynumbers[b][i][j][vlist[k][i]][0] = sol[pos];
                prob_keynumbers[b][i][j][vlist[k][i]][1] = W[b][vlist[k][i]];
		prob[b][i][j][vlist[k][i]] = sol[pos] /  W[b][vlist[k][i]] ;
	      }
	    else 
	      {
		prob_keynumbers[b][i][j][vlist[k][i]][0] = sol[pos];
                prob_keynumbers[b][i][j][vlist[k][i]][1] = W[b][vlist[k][i]];
		prob[b][i][j][vlist[k][i]] = 0.;
	      }
	  }
      
      //for(i=0; i<nV; i++)
      //W[b][i]=beta*Nb*W[b][i];   // renormalize this bond-dependent weight 
    }
  //ftest << "DONE WITH PROBABILY TABLE GENERATION" << endl;
}

// =======================================================================================================





// -------------------------------------------------------------------------------------------------------
//
void WRITE_HEADER(ofstream& fp)
//
//  This will write the header for the data files.
//
// -------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------
{

  string sp= "      ";
  fp << "# ---------------------------------------------------------------------------------"   << endl;
  fp << "# Here are the results for a: " << Lx << " x " << Ly <<" x " << Lz << " lattice"       << endl;
  fp << "# ---------------------------------------------------------------------------------"   << endl;
  fp << "# Delta/J: "  << Delta                                                                 << endl;
  fp << "# hbi/J: "    << hbi*2.*d                                                              << endl;
  fp << "# hbf/J: "    << hbf*2.*d                                                              << endl;
  fp << "# dh/J: "     << dh*2.*d                                                               << endl;
  fp << "# Nbins: "    << Nbins                                                                 << endl;
  fp << "# Nsamples: " << MCsweeps                                                              << endl;
  fp << "# Nrealizations: " << no_real                                                          << endl;
  fp << "# ================================================================================="   << endl;
  fp << "# Cutoff"  << sp << "T"     << sp  <<  "beta"     << sp << "hb"  << sp << "Delta"    << sp 
     << "E/N"       << sp << "<M>/N" << sp  << "<M^2>/N^2" << sp << "Sus" << sp << "stag <M>" << sp 
     << "stiffness" << sp << "comp"  << endl;
  fp << "# ================================================================================="    << endl;
}


// =======================================================================================================

void CHECK_PROBTABLE()


// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// -------------------------------------------------------------------------------------------------------
{
  double sum=0.;
  bool flag=true;
  for( int b=0; b<Nb; b++)
    for( int k=0; k<nV; k++)
      for( int i=0; i<4; i++)
	{
	  sum=0.;
	  for( int j=0; j<4; j++)
	    sum+=prob[b][i][j][k];
	  
	  // The sum should be 1. at the end of this and we check to ensure that this is so.
	  // This is a flagged case, meaning there are round off errors
	  if( fabs(sum-1.)>zerotol && sum !=0.) 
	    {
	      for(int e=0; e<4; e++)               
		ferr << setprecision(15) << "ERROR --> Bond: " << b << " in: " << i <<  " out: " << e << " prob: " << prob_keynumbers[b][i][e][k][0]<< "/" << prob_keynumbers[b][i][e][k][1] <<    " = " << prob[b][i][e][k] << " " << prob_keynumbers[b][i][e][k][0]/ prob_keynumbers[b][i][e][k][1] << endl;
	      flag=false;
	    }
	}
  
  if(flag)
    ferr << "No problems with REALIZATION: " << GLOBAL_REAL << " computed by:  " << GLOBAL_RANK << endl;
  else 
    ferr << "problems with REALIZATION: " << GLOBAL_REAL << " computed by:  " << GLOBAL_RANK << endl;
  ferr.flush();
  
}


void PRINT_PROBTABLE2(ofstream& fp)
{
  int i,e,v, l;
  
  for(l=0; l<Nb; l++)
    {
      fp << "BOND: " << l << endl;
      for(v=0; v<nV; v++)
	{
	  fp << "Vertex: " << v << endl;
	  for(i=0; i<4; i++)
	    {
	      for(e=0; e<4; e++)
		fp <<  prob[l][i][e][v]  << " " ;
	      
	      fp << endl;
	    }
	  fp << endl;
	}
      
      fp << endl;
    }
}


// =======================================================================================================









// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// -------------------------------------------------------------------------------------------------------
// AUXILLARY FUNCTIONS 
//
//  1. is_even:
//     This determines is an integer is even or not
//
//  2. PRINT_SPINLIST
//
//  3. PRINT_OPSTRING
//
//  4. PRINT_BONDLIST
//
//  5. PRINT_VTXLIST
//  
//  6. PRINT_PROBTABLE
//
//  7. PRINT_ARRAYTABLE
// 
//  8. PRINT_LINKLIST
//
//  9. VERTEX_ID
//
//  10. INIT_VERTEX_SPINS
//
// -------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------


bool is_even(int i)
{
  if(i%2==0)
    return true;
  else 
    return false;
}



// =======================================================================================================



void PRINT_SPINLIST(lattice& a, ofstream& fp)
{
  for(int i=0; i<N; i++)
    fp << a.spin[i] << " " ;
  fp << endl;
  
}


// =======================================================================================================

/*
void PRINT_OPSTRING(ofstream& fp)
{
  fp << "Operator String List" << endl;
  for(int i=0; i<M; i++)
    {
      if (sm[i]==-1)
	fp << i << "-> " << sm[i] << " b: " << "-" << endl;
      else
       fp << i << "-> " << sm[i] << " b: " << sm[i]/2 << endl;
    }
  fp << endl;
}

// =======================================================================================================

void PRINT_LINKLIST(ofstream& fp)
{
  fp << "Link List" << nH << " " << M << endl;
  for(int i=0; i<M; i++)
    {
      if(i==nH)
	fp << "-------------------------------------------" << endl;
      
      fp << i << "-->" << 4*i << ":" <<linklist[4*i] << " " << 4*i+1 << ":" << linklist[4*i+1] << " " << 
	4*i+2 <<":"<< linklist[4*i+2] << " " << 4*i+3 <<":" << linklist[4*i+3] << endl;
    }
}

// =======================================================================================================

void PRINT_BONDLIST(ofstream& fp)
{
  for(int dd=0; dd<d; dd++)
    {
      fp << "direction = " << dd+1 <<endl ;  // direction
      for(int b=0; b<N; b++)
	{
	  fp << "bond: " << b+dd*N << " -> "  ;
	  for(int s=0; s<2; s++)
	    fp << bsite[s][b+dd*N] << " " ;
	  fp << endl;
	}
      fp << endl;
    }
}

// =======================================================================================================


void PRINT_VTXLIST(ofstream& fp)
{
  fp << "Vtx List" << " " << nH << " " << M << endl;
  for(int i=0; i<M; i++)
    {
      if(i==nH) 
	fp << "----------------------------------------------" <<endl;
      fp << i << "-> " << vtx[i] << endl;
    }
  fp << endl;
}


// =======================================================================================================


void PRINT_PROBTABLE(ofstream& fp)
{
  int i,e,v, l;
  
  for(v=0; v<nV; v++)
    {
      fp << "Vertex: " << v << endl;
      for(i=0; i<4; i++)
	{
	  for(e=0; e<4; e++)
	    {
	      fp <<  prob[i][e][v]  << " " ;
	    }
	  fp << endl;
	}
      fp << endl;
    }
}


// =======================================================================================================

void PRINT_ARRAYTABLE2(ofstream& fp)
{
  int i,e,l;
    
  for(l=0; l<Nb; l++)
    {
      // PRINT A TABLE
      fp << "A-table" << endl;
      for(i=0; i<4; i++)
	{
	  for(e=0; e<4; e++)
	    fp << a[l][i][e] << " " ;
	  fp << endl;
	}
      fp << endl;
      
      fp << "B-table" << endl;
      
      // PRINT B TABLE
      for(i=0; i<4; i++)
	{
	  for(e=0; e<4; e++)
	    fp << b[l][i][e] << " " ;
	  fp << endl;
	}
      fp << endl;
    }
  
}


// =======================================================================================================

void PRINT_WEIGHTS(ofstream& fp)
{
  for(int i=0; i<nV; i++)
    fp << "W[" <<i<<"]: "<< W[i] <<endl;
}


// =======================================================================================================


void PRINT_WEIGHTS2(ofstream& fp)
{
  for(int l=0; l<Nb; l++)
    {
      for(int i=0; i<nV; i++)
	fp << "W[" <<l<<"]["<< i << "]: " <<  W[l][i] <<endl;
      
      fp<<endl;
    }
  fp <<endl;
}

// =======================================================================================================


void PRINT_STAG(ofstream& fp)
{
  int i;
  for(int z=0; z<Lz; z++)
    for(int y=0; y<Lx; y++)
      {
	for(int x=0; x<Lx; x++)
	  {
	    i = x + y*Lx + z*Lx*Ly;
	    fp << phi[i] << " " ;
	  }
	fp << endl;
      }
}



// =======================================================================================================

void PRINT_VLIST(ofstream &fp)
{
  for(int i=0; i<nV; i++)
    {
      for(int n=0; n<4; n++)
	fp << vlist[i][n] << " " ;
      
      fp << endl;
    }
}

*/


// =======================================================================================================


double CALC_MAG(int* s)
{
  int Stotal=0;
  for(int i=0; i<N; i++)
    Stotal+= s[i];
  
  return (0.5*Stotal);
}

// =======================================================================================================


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



void SETUP_FILE_STREAMS()
{
  string sys = convertInt(Lx)+"x"+convertInt(Ly);
  // Name of file:
  string infile = dir+"data/"+sys+"/DIST/"+sys+"h"+convertInt(h_iter)+"p"+convertInt(batch)+"beta"+convertDouble(beta);
  string errfile= dir+"data/"+sys+"/ERR/ERRh"+convertInt(h_iter)+"p"+convertInt(batch)+"beta"+convertDouble(beta);
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

  
#if saveconfig
  fp2.open  ( (dir+"CONFIG/init"+convertInt(Lx)+"x"+convertInt(Ly)+"h"+convertInt(h_iter)+"p"+convertInt(batch)+"beta"+convertDouble(beta)).c_str(), ios:: out );
  assert(fp2!=NULL);
#endif
  
  assert(ferr!=NULL);
  assert(fp1!=NULL);
  
}


/*---------------------------------------------------------------------
  WARNING: I have no idea what is going on with the simplex method 
  here on out.  Just rip it off from Numerical recipes of C++ and 
  hope for the best
  DO NOT CHANGE ANYTHING BELOW ...
  -------------------------------------------------------------------- */


// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// -------------------------------------------------------------------------------------------------------

void simplx(MAT a,int m,int n,int m1,int m2,int m3,int *icase,int *izrov, int *iposv)  
{ 
  /*---------------------------------------------------------------------------------------- 
    USES simp1,simp2,simp3 
    Simplex method for linear programming. Input parameters a, m, n, mp, np, m1, m2, and m3, 
    and output parameters a, icase, izrov, and iposv are described above (see reference). 
    Parameters: MMAX is the maximum number of constraints expected; NMAX is the maximum number 
    of variables expected; EPS is the absolute precision, which should be adjusted to the 
    scale of your variables.                                    
    -----------------------------------------------------------------------------------------*/
  
  
  // -------------------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------------------
    
  int i,ip,ir,is,k,kh,kp,m12,nl1,nl2,l1[NMAX],l2[MMAX],l3[MMAX]; 
  REAL bmax,q1,EPS=1e-6; 
  if(m != m1+m2+m3) {
    printf(" Bad input constraint counts in simplx.\n");
    return;
  }
  nl1=n; 
  for (k=1; k<=n; k++) { 
    l1[k]=k;     //Initialize index list of columns admissible for exchange.
    izrov[k]=k;  //Initially make all variables right-hand. 
  } 
  nl2=m;
  for (i=1; i<=m; i++) { 
    if (a[i+1][1] < 0.0) {
      printf(" Bad input tableau in simplx, Constants bi must be nonnegative.\n");
      return;
    }
    l2[i]=i;
    iposv[i]=n+i; 
    /*------------------------------------------------------------------------------------------------
      Initial left-hand variables. m1 type constraints are represented by having their slackv ariable 
      initially left-hand, with no artificial variable. m2 type constraints have their slack 
 variable initially left-hand, with a minus sign, and their artificial variable handled implicitly 
 during their first exchange. m3 type constraints have their artificial variable initially 
 left-hand.     
 ------------------------------------------------------------------------------------------------*/ 
  }
  for (i=1; i<=m2; i++) l3[i]=1;
  ir=0;
  if(m2+m3 == 0) goto e30; //The origin is a feasible starting solution. Go to phase two. 
  ir=1;
  for (k=1; k<=n+1; k++) { //Compute the auxiliary objective function. 
    q1=0.0; 
    for (i=m1+1; i<=m; i++) q1 += a[i+1][k]; 
    a[m+2][k]=-q1; 
  } 
 e10: simp1(a,m+1,l1,nl1,0,&kp,&bmax);    //Find max. coeff. of auxiliary objective fn 
  if(bmax <= EPS && a[m+2][1] < -EPS) { 
    *icase=-1;    //Auxiliary objective function is still negative and can.t be improved, 
    return;       //hence no feasible solution exists.
  }
  else if (bmax <= EPS && a[m+2][1] <= EPS) { 
    //Auxiliary objective function is zero and can.t be improved; we have a feasible starting vector. 
    //Clean out the artificial variables corresponding to any remaining equality constraints by 
    //goto 1.s and then move on to phase two by goto 30.
    m12=m1+m2+1;
    if (m12 <= m)
      for (ip=m12; ip<=m; ip++) 
	if(iposv[ip] == ip+n) {   //Found an artificial variable for an equalityconstraint. 
	  simp1(a,ip,l1,nl1,1,&kp,&bmax); 
	  if(bmax > EPS) goto e1; //Exchange with column corresponding to maximum 
	}                         //pivot element in row. 
    ir=0;
    m12=m12-1;
    if (m1+1 > m12) goto e30; 
    for (i=m1+1; i<=m1+m2; i++)     //Change sign of row for any m2 constraints 
      if(l3[i-m1] == 1)             //still present from the initial basis. 
        for (k=1; k<=n+1; k++) 
          a[i+1][k] *= -1.0; 
    goto e30;                       //Go to phase two. 
  } 

  simp2(a,m,n,l2,nl2,&ip,kp,&q1); //Locate a pivot element (phase one). 
                                         
  if(ip == 0) {                         //Maximum of auxiliary objective function is 
    *icase=-1;                          //unbounded, so no feasible solution exists.
    return; 
  } 
 e1: simp3(a,m+1,n,ip,kp); 
  //Exchange a left- and a right-hand variable (phase one), then update lists. 
  if(iposv[ip] >= n+m1+m2+1) { //Exchanged out an artificial variable for an 
                               //equality constraint. Make sure it stays 
                               //out by removing it from the l1 list. 
    for (k=1; k<=nl1; k++) 
      if(l1[k] == kp) goto e2; 
  e2: nl1=nl1-1; 
    for (is=k; is<=nl1; is++)  l1[is]=l1[is+1]; 
  }
  else {
    if(iposv[ip] < n+m1+1) goto e20;
    kh=iposv[ip]-m1-n; 
    if(l3[kh] == 0) goto e20;  //Exchanged out an m2 type constraint. 
    l3[kh]=0;                  //If it.s the first time, correct the pivot column 
                               //or the minus sign and the implicit 
                               //artificial variable. 
  }   
  a[m+2][kp+1] += 1.0; 
  for (i=1; i<=m+2; i++)  a[i][kp+1] *= -1.0; 
 e20: is=izrov[kp];             //Update lists of left- and right-hand variables. 
  izrov[kp]=iposv[ip]; 
  iposv[ip]=is; 
  if (ir != 0) goto e10;       //if still in phase one, go back to 10. 
  //End of phase one code for finding an initial feasible solution. Now, in phase two, optimize it. 
 e30: simp1(a,0,l1,nl1,0,&kp,&bmax); //Test the z-row for doneness. 
  if(bmax <= EPS) {          //Done. Solution found. Return with the good news. 
    *icase=0; 
    return; 
  }
  simp2(a,m,n,l2,nl2,&ip,kp,&q1);   //Locate a pivot element (phase two). 
  if(ip == 0) {             //Objective function is unbounded. Report and return. 
    *icase=1; 
    return; 
  } 
  simp3(a,m,n,ip,kp);       //Exchange a left- and a right-hand variable (phase two), 
  goto e20;                 //update lists of left- and right-hand variables and 
}                           //return for another iteration.

// The preceding routine makes use of the following utility subroutines: 


// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// -------------------------------------------------------------------------------------------------------
void simp1(MAT a,int mm,int *ll,int nll,int iabf,int *kp,REAL *bmax) 
{
  //Determines the maximum of those elements whose index is contained in the supplied list 
  //ll, either with or without taking the absolute value, as flagged by iabf. 


// -------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------
  
  int k; 
  REAL test; 
  *kp=ll[1]; 
  *bmax=a[mm+1][*kp+1];
  if (nll < 2) return; 
  for (k=2; k<=nll; k++) { 
    if(iabf == 0) 
      test=a[mm+1][ll[k]+1]-(*bmax); 
    else
      test=fabs(a[mm+1][ll[k]+1])-fabs(*bmax); 
    if(test > 0.0) { 
      *bmax=a[mm+1][ll[k]+1]; 
      *kp=ll[k]; 
    } 
  } 
  return; 
} 




// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// -------------------------------------------------------------------------------------------------------

void simp2(MAT a, int m,int n,int *l2,int nl2,int *ip,int kp,REAL *q1) {

  
  
// -------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------
  
  REAL EPS=1e-6; 
  //Locate a pivot element, taking degeneracy into account. 
  int i,ii,k; 
  REAL q,q0,qp; 
  *ip=0; 
  if(nl2 < 1) return;
  for (i=1; i<=nl2; i++) 
    if (a[i+1][kp+1] < -EPS) goto e2; 
  return;  //No possible pivots. Return with message. 
 e2: *q1=-a[l2[i]+1][1]/a[l2[i]+1][kp+1]; 
  *ip=l2[i];
  if (i+1 > nl2) return; 
  for (i=i+1; i<=nl2; i++) { 
    ii=l2[i];
    if(a[ii+1][kp+1] < -EPS) { 
      q=-a[ii+1][1]/a[ii+1][kp+1]; 
      if (q <  *q1) { 
        *ip=ii; 
        *q1=q; 
      }
      else if (q == *q1) {  //We have a degeneracy.
	for (k=1; k<=n; k++) { 
          qp=-a[*ip+1][k+1]/a[*ip+1][kp+1]; 
          q0=-a[ii+1][k+1]/a[ii+1][kp+1]; 
          if (q0 != qp) goto e6; 
        } 
      e6:     if (q0 < qp) *ip=ii; 
      } 
    } 
  } 
  return; 
} 




// ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// -------------------------------------------------------------------------------------------------------
//
void simp3(MAT a,int i1,int k1,int ip,int kp) {
//
// -------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------

  //Matrix operations to exchange a left-hand and right-hand variable (see text). 
  int ii,kk; 
  REAL piv; 
  piv=1.0/a[ip+1][kp+1];
  if (i1 >= 0) 
    for (ii=1; ii<=i1+1; ii++) 
      if (ii-1 != ip) { 
        a[ii][kp+1] *= piv; 
        for (kk=1; kk<=k1+1; kk++) 
          if (kk-1 != kp) 
            a[ii][kk] -= a[ip+1][kk]*a[ii][kp+1]; 
      } 
  for (kk=1; kk<=k1+1; kk++) 
    if(kk-1 !=  kp) a[ip+1][kk] =-a[ip+1][kk]*piv; 
  a[ip+1][kp+1]=piv; 
  return; 
} 



