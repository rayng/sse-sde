#include "lattice.hpp"

// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::buildlattice()
// ++++++++++++++++++++++++++++++++++++++++++++++++++
{
  int i, x,  y, z, v, l, s[4], e, vnum, inleg, exitleg;
  double hav=0.;
  
  // Calloc-ing memory.   <-- just a good exercise I suppose.
  bsite = (int**) calloc( 2,sizeof(int*) );       // This is global, all lattices will share this. 
  for(i=0; i<2; i++)
    bsite[i] = (int*) calloc(Nb, sizeof(int) );
  
  
  //--------------------------
  // CUBIC GEOMETRY
  //--------------------------
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
  
  // Staggered field
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
  optype[2] = 1;
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
  optype[6] = 1;
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
  if(system=="staggered")
    for(i=0; i<N; i++) 
      h[i] = phi[i]*hb;
  
  if(system=="uniform")
    for(i=0; i<N; i++) 
      h[i] = hb;
  
  if(system=="disorder")
    {
      hav=0.;
      for(i=0; i<N; i++) 
	h[i] = 2.*hb*r.FloatN()-hb;
      hav+=h[i];          // No more shifting
      hav = hav / (double) N;
      for(i=0; i<N; i++)
	h[i] -= hav;
    }
  
  
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



// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::generate_probtable_simplex()
// ++++++++++++++++++++++++++++++++++++++++++++++++++
{
  MAT  A;
  int  IPOSV[MMAX], IZROV[NMAX];
  int  ICASE;
  int  NN = nV*16+1;         // Number of variables.
  int  MM;                   // Total number of constraints
  int  M1 =  0 ;             // Number of '<=' constraints
  int  M2 =  1 ;             // Number of '>=' constraints
  int  M3 =  nV*6 + nV*4;    // Number of '=' constraints
  REAL R;
  
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
      W[b][0]= 0.25;
      W[b][1]= 0.;         // + epsilonn[b];
      W[b][2]= 0.25;
      W[b][3]= -hb;         // + epsilonn[b];
      
      W[b][4]= 0.25;
      W[b][5]= 0.;          // + epsilonn[b];
      W[b][6]= 0.25;
      W[b][7]= hb;          // + epsilonn[b];
      
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
      //blank();
      
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
	
	if( abs(W[b][k]) < zerotol && W[b][k]!=0.) {
	  W[b][k]=0.;
	  cout << "Negative weights!: " << b << " " << k << " " << W[b][k]  << endl;  
	}
	
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
      
      for(i=0; i<nV; i++)
	W[b][i]=beta*Nb*W[b][i];   // renormalize this bond-dependent weight 
      
    }


}









// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::init_traj_z()
// ++++++++++++++++++++++++++++++++++++++++++++++++++
{
  int i;
  //If z-variables
  for(i=0; i<N; i++)
    {
      if(spin[i]==1)
	{
	  traj.z[i].real()=5.;
	  traj.z[i].imag()=0.;
	  traj.zp[i]=conj(traj.z[i]);
	}
      else if(spin[i]==-1)
	{
	  traj.z[i]=-5.;
	  traj.z[i].imag()=0.;
	  traj.zp[i]=conj(traj.z[i]);
	}
    }
}




// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::init_traj_y()
// ++++++++++++++++++++++++++++++++++++++++++++++++++
{
  int i;
  //If y-variables
  
  for(i=0; i<N; i++)
    {
      traj.y[i].real()=1e-10;
      traj.y[i].imag()=0.;
      traj.yp[i]=conj(traj.y[i]);
      traj.SS[i]=spin[i];
      traj.SSp[i]=traj.SS[i];
    }
}



// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::init_config(string state)
// ++++++++++++++++++++++++++++++++++++++++++++++++++
{
  long int i;
  int l, e, v;
  for(i=0; i<N; i++) 
    if(state=="hightemp"){
      if(r.FixedU<double>() < 0.5)
	spin[i] = -1;
      else
	spin[i] = +1;
    }
    else if(state=="FMup")
      spin[i] = 1;
    else if(state=="FMdown")
      spin[i] = -1;
  
  sm       = new int [M];          // Operator string
  vtx      = new int [M];          // Vertex list
  vb       = new int [M];          // Vb list (bonds?)  
  linklist = new int [4*M];        // Linklist
  
  for(i=0; i<M; i++) 
    {
      sm[i]  = -1;      
      vtx[i] = -2; 
      vb[i]  = -2; 
      linklist[4*i]   = -1;  
      linklist[4*i+1] = -1;  
      linklist[4*i+2] = -1;  
      linklist[4*i+3] = -1;  
    }
}



// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::update_spin_config()
// ++++++++++++++++++++++++++++++++++++++++++++++++++
{ // For after loop update
  int i, p , l;
  for(i=0; i<N; i++)
    if(first[i]!=-1) {
      p = first[i]/4;  // vertex number;
      l = first[i]%4;
      spin[i]  = legspin[vtx[p]][l];
    }
    else {
      //flip with probabiliy half
      if(r.FloatN<double>() > 0.5)
	spin[i] = -spin[i];
    }
}
// **************************************************




// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::diagonal_update()
// ++++++++++++++++++++++++++++++++++++++++++++++++++
{
  long int i, b;
  int v, s0, s1, s2, s3;
  double pinsert, premove, rand;
  
  for(i=0; i<M; i++) 
    {
      rand = r.FloatN<double>();
      
      if( sm[i] == -1 ) 
	{  // identity operator found
	  b = r.IntegerC<long int>(0,Nb-1);
	  
	  //calculate the matrix element:   <alpha(p) | H | alpha(p) >
	  v = vtx_id[spin[bsite[0][b]]+1][spin[bsite[1][b]]+1][spin[bsite[0][b]]+1][spin[bsite[1][b]]+1];  
	  // insert diagonal element  
	  // pinsert = (beta*Nb*W[b][v]) / (double) (M-nH)  ;
	  pinsert = W[b][v] / (double) (M-nH)  ;
	  
	  if (pinsert>=rand) 
	    {
	      sm[i] = 2*b ;
	      nH++; 
	    }
	  else 
	    continue;
	}
      else if ( sm[i]%2 == 0)  // diagonal
	{    
	  b = sm[i]/2;
	  v = vtx_id[spin[bsite[0][b]]+1][spin[bsite[1][b]]+1][spin[bsite[0][b]]+1][spin[bsite[1][b]]+1];  
	  // remove diagonal element  	  
	  //premove =  (double) (M-nH+1) / (beta*Nb * W[b][v]) ; 
	  premove =  (double) (M-nH+1) / W[b][v] ; 
	  if (premove>=rand) 
	    {
	       sm[i]  = -1 ; 
	       nH--; 
	     }
	}
      else  
	{
	  b = sm[i] /2;
	  spin[bsite[0][b]] = -spin[bsite[0][b]];  //update spin
	   spin[bsite[1][b]] = -spin[bsite[1][b]];
	}
      
      
    }
}
// **************************************************



// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::adjust_truncate()
// ++++++++++++++++++++++++++++++++++++++++++++++++++
{
  // This changes the expansion order
  long int Mnew, i, j;
  int *temp;
  
  Mnew = nH+nH/3;
  if (Mnew>M) {
    temp = new int [Mnew];
    assert(temp!=NULL);
    for (i=0; i<M; i++)              // copy sm over
      temp[i]=sm[i];
    for (i=M; i<Mnew; i++)           // set new elements to -1
      temp[i]=-1 ;
    delete sm;
    sm =new int[Mnew];
    for(i=0; i< Mnew; i++)
      sm[i] = temp[i];
    delete temp;
    delete linklist;
    linklist = new int [4*Mnew];
    M = Mnew;                    // New M value;
    for(i=0; i<M; i++) {         // initialise link list after resizing
      linklist[4*i] = -1;
      linklist[4*i+1] = -1;
      linklist[4*i+2] = -1;
      linklist[4*i+3] = -1;      
    }
    delete vtx;
    delete vb;
    vtx = new int [M];
    vb = new int [M];
    assert(vtx!=NULL);
    assert(vb!=NULL);
  } 
  
}
// **************************************************





// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::construct_linklist()
// ++++++++++++++++++++++++++++++++++++++++++++++++++
{
  long int lnk, p, vo, v1, v2, i1, i2, f, l, i;
  for (i=0; i<N; i++) { // Zero list first
    first[i]=-1;      // What is the first vertex leg of spin i
    last[i]=-1;       // What is the last vertex leg of spin i
  }
  
  // Maybe this doesnt matter but we must be careful not to use values that should not be used. 
  for(i=0; i<M; i++) {
    linklist[4*i] = -1;
    linklist[4*i+1] = -1;
    linklist[4*i+2] = -1;
    linklist[4*i+3] = -1;      
    vtx[i] = -2;
    vb[i] = -2;  // what bond is 
  }
  int v=0; // vertex counter: This will go up to nH
  for(p=0; p<M; p++) {
    vo = 4*v;
    lnk =  sm[p]/2;
    
    // Extract spins associated wth bond
    i1 = bsite[0][lnk];
    i2 = bsite[1][lnk];
    v1 = last[i1];
    v2 = last[i2];
    if(sm[p]==-1) // This means that it is a diagonal operator so we  have no legs to number
      continue;
    else {
      if(v1!=-1) {
	linklist[v1] = vo;
	linklist[vo] = v1; 
      }
      else
	first[i1] = vo;
      if(v2!=-1)  {
	linklist[v2] = vo+1;
	linklist[vo+1] = v2;
      }
      else 
	first[i2] = vo+1;
      last[i1] = vo+2;
      last[i2] = vo+3;
      // -------------------------
      // Update vertex list
      // -------------------------
      //if( is_even(sm[p]) )        // i.e. a diagonal operator:
      if( sm[p]%2==0 )        // i.e. a diagonal operator:
	vtx[v] = vtx_id[spin[i1]+1][spin[i2]+1][spin[i1]+1][spin[i2]+1];
      else {
	vtx[v] = vtx_id[spin[i1]+1][spin[i2]+1][-spin[i1]+1][-spin[i2]+1];
	spin[i1] = -spin[i1];     // Flip it because we have an off diagonal operator
	spin[i2] = -spin[i2];
      }
      vb[v] = lnk;                // This stores the bond associated with the vertex
      v++;                          // Increase vertex number by 1
    }
  }
  if( nH!=v )
    cout << "Problem counting vertices in linklist subroutine:" << nH - v  << endl;
  // ----------------------------------------
  // Connect spins across the boundaries
  // ----------------------------------------
  // **This is pretty important. 
  for(i=0; i<N; i++) {
    f= first[i];
    if (f!=-1) {
      l = last[i];
      linklist[f] = l;
      linklist[l] = f;
    }
  }

  
}

// **************************************************



// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::extend_arrays()
// ++++++++++++++++++++++++++++++++++++++++++++++++++
{
  long int i;
  int *temp;
  // +++++++++++++++++++++++++++++++++++++++++
  // sm -arrays
  // ------------
  temp = new int [M];      // old M
  for(i=0; i<M; i++)       // copy into temp arrays
    temp[i] = sm[i];
  delete sm;
  sm = new int[2*M];     // doubled now
  // Now copy over!
  for(i=0; i<M; i++)
    {
      sm[i] = temp[i];
      sm[M+i] = temp[M-1-i];
    }
  // +++++++++++++++++++++++++++++++++++++++++++
  // vtx -arrays  (These dont really have to be copied over 
  // because LOOP_UPDATE algorithm will reinitialise them)
  //  
  // ------------
  for(i=0; i<nH; i++)
    temp[i] = vtx[i];       
  delete vtx;             // This array is over extended for convenience
  vtx = new int[2*M];   // doubled now
  for(i=0; i<nH; i++)     // Now copy over!
    {
      vtx[i] = temp[i];
      vtx[nH+i] = temp[nH-1-i];
    }
  // +++++++++++++++++++++++++++++++++++++++++++
  // vb -arrays
  // ------------
  for(i=0; i<M; i++)
    temp[i] = vb[i];
  delete vb;
  vb = new int[2*M];     // doubled now
  for(i=0; i<nH; i++)     // Now copy over!
    {
      vb[i] = temp[i];
      vb[nH+i] = temp[nH-1-i];
    }
  // +++++++++++++++++++++++++++++++++++++++++++
  // Double M and nH
  // ------------
  set_M(2*M);      // Initialise M
  set_nH(2*nH);    // Initialise no
  delete temp;
  // +++++++++++++++++++++++++++++++++++++++++++
  // Make sure to readjust linklist array
  // ------------
  delete linklist;
  linklist = new int [4*M];
}
// **************************************************



// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::loop_update_eqm(int nsamples, ofstream& ferr)
// ++++++++++++++++++++++++++++++++++++++++++++++++++
{
  long int i,j, j0, e, nn, nvisits=0, localvisits=0, ll=0,  bb, lnk,p,l,  nav, cutoff_visits, maxvisits_per_loop;
  double r1, cum_prob;
  bool found=false, islong=false;
  int* vtx_temp;
  
  n_sum += nH;   // this is global
  nav = n_sum/(nsamples+1);
  cutoff_visits=2.*nav;   // cutoff
  maxvisits_per_loop=50.*nav;
  
  construct_linklist();
  
  //copy over the vtx list. 
  vtx_temp = new int [nH];
  for(i=0; i<nH; i++)
    vtx_temp[i]=vtx[i];
  
  while(1)
    { 
      // Finds number of loops to form based on how many vertices are visited
      if(islong) 
	{
	  update_spin_config();
	  diagonal_update();
	  adjust_truncate();
	  construct_linklist();
	  
	  delete vtx_temp;
	  vtx_temp = new int [nH];
	  for(i=0; i<nH; i++)
	    vtx_temp[i]=vtx[i];   // reset vtx_tmp
	}
      else if(!islong)
	nvisits +=localvisits;  // update number of vertices visited
      
      // -----------------------------------------------------------
      // START OF A NEW LOOP
      // -----------------------------------------------------------
      j0 = r.IntegerC<long int>(0,4*nH-1);   // Pick a random leg. 
      j=j0;
      islong=false;
      found =false;
      localvisits=0;
      
      // -----------------------------------------------------------
      // Now find the exit leg
      // -----------------------------------------------------------
    A: 
      p = j/4;       //vertex number
      l = j%4;       // leg number 
      bb = vb[p];  // this doesnt change in a loop update. It is just being referenced.
      
      if(localvisits>maxvisits_per_loop) {
	ferr << " Long loop in equilibration stage " << localvisits << " vs " <<  maxvisits_per_loop << " at " << beta << endl;
	islong=true;
	goto SKIP;   // start over
      }
      
      //RETRY:
      //cout << "crapped out" << endl;
      r1 = r.FloatN<double>();
      //r1=0.;
      cum_prob =0.;
      found =false;
      
      for(e=0; e<4; e++)  // Find the exit leg. 
	{
	  cum_prob += prob[bb][l][e][vtx_temp[p]];
	  if(r1 <= cum_prob)   
	    {  // we've found the exit leg
	      // Flip the spins at entrance l and exit e of vertex vtx[p]
	      vtx_temp[p] = vtx_flip[vtx_temp[p]][l][e];   
	      //a.vtx[p] = vtx_flip[a.vtx[p]][l][e];   
	      found=true;
	      break;
	    }
	
	  // We should never make it here because we                                                                                                                                                        // are given a chance to break the loop above                                                                                                                                                   	 
	  if(e==3) {
	    ferr << "No exit found for vtx: " << vtx_temp[p] << " located on bond "<< bb <<   " in equilibrim stage!" << endl;
	    ferr << "Vertex: " << vtx_temp[p] << " Entering in: " << l << " exiting via: " << e << endl;
            ferr << setprecision(15)
		 << "Random number: " << r1 << " Probabilities to beat: "
                 << " + " << prob[bb][l][0][vtx_temp[p]] << " + " << prob[bb][l][1][vtx_temp[p]]
                 << " + " << prob[bb][l][2][vtx_temp[p]] << " + " << prob[bb][l][3][vtx_temp[p]]
                 << " = " << cum_prob << endl << endl;
            //goto RETRY;    // in case we crap out again 
	  }
	}
      j = 4*p+e;
      if(j!=j0) { // this means we have not found the original leg. boo	and onward we go.
	j = linklist[j];
	localvisits++;
	if(j!=j0)
	  goto A;
	else 
	  goto B;
      }
      else 
	goto B;
      
      
      // -----------------------------------------------------------
      //  If we make it to B it means that a loop has closed
      //  This is a good point to see if we can safely exit.
      // ----------------------------------------------------------- 
    B:   
      // | ----- Update operator string ------ |
      for(l=0; l<nH; l++)     //update a.vtx.
	vtx[l]=vtx_temp[l];
      
      nn=0;
      for(l=0; l<M; l++)
	if(sm[l]==-1)
	  continue;
	else {
	  lnk = sm[l]/2;
	  sm[l] = 2*lnk+optype[vtx[nn]];   // a.sm is changed here 
	  nn++;
	}
      
      Nl_sum++;               // every time we make it here, it constitutes a loop.
      if(nvisits> cutoff_visits)  // we can stop making loops now
	goto EXIT;
      
    SKIP: // This means we do not implement any of the changes made in the vertex changes. 
      ;
      
    } //end while loop
  
 EXIT:
  // | ----- Update Spin List ----- |
  update_spin_config();
  delete vtx_temp;
  

}

// **************************************************




// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::loop_update(int nsamples, ofstream& ferr)
// ++++++++++++++++++++++++++++++++++++++++++++++++++
{
  long int i,j, j0, e, nn, nvisits=0, localvisits=0, ll=0,  bb, lnk,p,l, max_looplength, max_visits_per_loop, nav, n;
  double r1, cum_prob;
  bool found=false, islong=false;
  int* vtx_temp;
  
  n_sum += nH;   // this is global
  nav = n_sum/(nsamples+1);
  max_visits_per_loop=50.*nav;
  
  construct_linklist();

  n_sum += nH;   // this is global
  nav = n_sum/(nsamples+1);
  max_looplength=50.*nav;   // cutoff
  
  //copy over the vtx list. 
  vtx_temp = new int [nH];
  for(i=0; i<nH; i++)
    vtx_temp[i]=vtx[i];
  
  for(n=0; n<Nl; n++)
    { 
      // Finds number of loops to form based on how many vertices are visited
      if(islong) 
	{
	  update_spin_config();
	  diagonal_update();
	  adjust_truncate();
	  construct_linklist();
	  delete vtx_temp;
	  vtx_temp = new int [nH];
	  for(i=0; i<nH; i++)
	    vtx_temp[i]=vtx[i];   // reset vtx_tmp
	}
      // -----------------------------------------------------------
      // START OF A NEW LOOP
      // -----------------------------------------------------------
      j0 = r.IntegerC<long int>(0,4*nH-1);   // Pick a random leg. 
      j=j0;
      islong=false;
      found =false;
      localvisits=0;
      // -----------------------------------------------------------
      // Now find the exit leg
      // -----------------------------------------------------------
    A: 
      p = j/4;       //vertex number
      l = j%4;       // leg number 
      bb = vb[p];  // this doesnt change in a loop update. It is just being referenced.
      if(localvisits>max_visits_per_loop) {
	ferr << " Long loop in sampling stage " << localvisits << " vs " <<  100.*nav << " at " << beta << endl;
	islong=true;
	goto SKIP;   // start over
      }
      
      r1 = r.FloatN<double>();
      //r1=1.0;
      cum_prob =0.;
      found =false;
      
      for(e=0; e<4; e++)  // Find the exit leg. 
	{
	  cum_prob += prob[bb][l][e][vtx_temp[p]];
	  if(r1 <= cum_prob)   
	    { 
	      vtx_temp[p] = vtx_flip[vtx_temp[p]][l][e];   
	      found=true;
	      break;
	    }
	}
      if(!found)
	ferr << "No exit found for vtx: " << vtx_temp[p] << endl;
      j = 4*p+e;
      if(j!=j0) { 
	// this means we have not found the original leg. 
	// boo. So onward we go.
	j = linklist[j];
	localvisits++;
	if(j!=j0)
	  goto A;
	else 
	  goto B;
      }
      else 
	goto B;
      // -----------------------------------------------------------
      //  If we make it to B it means that a loop has closed
      //  This is a good point to see if we can safely exit.
      // ----------------------------------------------------------- 
    B:   
      // | ----- Update operator string ------ |
      for(l=0; l<nH; l++)     //update a.vtx.
	vtx[l]=vtx_temp[l];
      
      nn=0;
      for(l=0; l<M; l++)
	if(sm[l]==-1)
	  continue;
	else {
	  lnk = sm[l]/2;
	  sm[l] = 2*lnk+optype[vtx[nn]];   // a.sm is changed here 
	  nn++;
	}
    SKIP: // This means we do not implement any of the changes made in the vertex changes. 

      ;
    } //end while loop
  
 EXIT:
  // | ----- Update Spin List ----- |
  
  delete vtx_temp;
  update_spin_config();  
}

// **************************************************


// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::measure()
// ++++++++++++++++++++++++++++++++++++++++++++++++++
{
  long int p, b, i, j, Nxm=0, Nxp=0, Nym=0, Nyp=0, nn=0;
  long int  ndistemp=0, nhoptemp=0;
  long double mtemp=0., stagtemp=0., windnumx=0., windnumy=0., hh=2.*d*hb;
  int parity=1;
  
  // include parity check
#if(parityinclusion)
  long int nhopparp=0, nhopparm=0;
  long int ndisparp=0, ndisparm=0;
#endif
  
#if(FS)
  Array<int, Dynamic,1 > Nhist; 
  vector<int> posn;
#endif
  
  for(i=0; i<N; i++) {
    parity=parity*spin[i];
    mtemp += spin[i];
    stagtemp += phi[i]*spin[i];
  }
	  
  assert(parity==1 || parity ==-1);
  
  if(parity==1)
    psamples++;
  else if (parity==-1)
    msamples++;
  
  mtemp    = 0.5*mtemp / (double) N;
  stagtemp = 0.5*stagtemp / (double) N;
  
  // staggered magnetization
  long double am;
  am = stagtemp;
  for(i=0; i<M; i++)                                               // sum up over operator strings.
    {
      if( sm[i]==-1 )
	continue;
      else if( sm[i]%2 ==1 ) 
	{                                      // off-diagonal 
	  b = sm[i]/2;
	  spin[bsite[0][b]] = -spin[bsite[0][b]];                   // flip spins
	  spin[bsite[1][b]] = -spin[bsite[1][b]];
	  am += 2.0*phi[ bsite[0][b] ]*spin[ bsite[0][b] ];
	  
	  // -------------------------------------------
	  // Calculate the spin stiffness. 
	  // -------------------------------------------
	  // I need to calculate the winding number. 
	  if(b<N) 
	    {
	      // This means it's an x bond. Recall x bonds number from 0 to N  
	      if( spin[bsite[0][b]]==1)
		Nxp++;
	      else if(spin[bsite[1][b]]==1)
		Nxm++;
	    }
	  else if(b>=N && b<2*N) 
	    {                                // This means it's an y bond. Recall y bonds number from N to 2N  
	      if( spin[bsite[0][b]]==1)
		Nyp++;
	      else if( spin[bsite[1][b]]==1)
		Nym++;	      
	    }
	
#if(parityinclusion)
	  if(parity==1)
	    nhopparp++;
	  else
	    nhopparm++;
#endif
	  nhoptemp++;   // ignoring parity hopping terms
	  
#if(FS)
	  posn.push_back(nn);  // this stores the positions of the hopping parameters
#endif
	  nn++;
	}
      else  // non-trivial diagonal operators
	{
	  
#if(parityinclusion)	  
	  if(parity==1)
	    ndisparp++;
	  else
	    ndisparm++;
#endif
	  ndistemp++;   // ignoring parity disorder terms
	  nn++; 
	}
      stagtemp+=am;
    }

  
  assert(nn==nH);
    
  windnumx = (Nxp-Nxm)/ (double) Lx;
  windnumy = (Nyp-Nym)/ (double) Ly;
  
#if FS
  Nhist.resize(nH);
  Nhist.setZero();
  long int m, l;
  // Build histogram.
  for (i=0; i<posn.size(); i++) {
    for (j=0; j<posn.size(); j++) 
      if(i!=j)	{
	m = posn[i]-posn[j];
	if(m>0) 
	  m = m - 1;
	else 
	  m = m + nH - 1;
	//assert(m>=0);
	if(m>nH-2)
	  continue;
	Nhist[m]++;
      }
  }
  
  //Calculate histogram dependent sum
  long double fs=0., X1=0., Amn=0.;
  long int p_,q,r;
  for(m=0; m<=nH-2; m++) {
    p_= m+1;
    q = nH-m-2;
    r = p_+q+1;
    X1 = ( (q-p_) ) / ( sqrt(r) ) ; 
    Amn = 0.5*(1.+erf( X1/sqrt(2.)))*(double) p_/ (double) r;
    fs += Amn*Nhist[m];
  }
  
  fs/=(hh*hh);
  
  
#endif
  

  // ++++++++++++++++++++++++++++++++++++++++++++++++
  // Accumulate the sum of averages
  // ++++++++++++++++++++++++++++++++++++++++++++++++
  long double windnumsq = (windnumx*windnumx + windnumy*windnumy)*0.5;
  XPS.dump( (long double) nH, 0 );                             // Energy
  XPS.dump( (long double) nH, 1 );                             // nop
  XPS.dump( (long double) (nH) * (long double) (nH), 2 );                        // nop^2
  XPS.dump( mtemp,3 );                                    // magnetization
  XPS.dump( mtemp*mtemp,4 );                              // magnetization sq
  XPS.dump( mtemp*mtemp*mtemp*mtemp, 5);                  // mag^4 for binder cumulants
  XPS.dump( stagtemp/(N*nH), 6 );                         // staggered magnetization  <-- this is fishy but we dont need it so whatever.
  XPS.dump( windnumsq*Pi/beta, 7 );                       // stiffness
  XPS.dump( windnumsq*windnumsq, 8);                      // W^4 for binder cumulants
  
  // No parity distinction
  XPS.dump( (long double) ndistemp, 9 );                       // Number disorder
  XPS.dump( (long double) (ndistemp*ndistemp) , 10 );          // Number disorder^2
  XPS.dump( (long double) nhoptemp, 11 );                 // Number hopping operators
  XPS.dump( (long double) (nhoptemp*nhoptemp), 12);       // Number hopping^2
  
  // With parity distinction
  // + sector
#if(parityinclusion)

  if(parity==1)
    {
      XPS.dump( (long double) ndisparp, 13 );                      // Number disorder +
      XPS.dump( (long double) (ndisparp*ndisparp) , 14 );          // 
      
      // + sector
      XPS.dump( (long double) nhopparp, 17 );                      // Number hop +
      XPS.dump( (long double) (nhopparp*nhopparp) , 18 );          // 
    }
  else if(parity==-1) 
    {
      // - sector
      XPS.dump( (long double) ndisparm, 15 );                      // Number disorder +
      XPS.dump( (long double) (ndisparm*ndisparm) , 16 );          // 
      
      // - sector
      XPS.dump( (long double) nhopparm, 19 );                      // Number hop - 
      XPS.dump( (long double) (nhopparm*nhopparm) , 20 );          //
    }

#endif

#if(FS)
  XPS.dump(fs,21);                                   //    0 fs
  if(parity==1)
    XPS.dump(fs,22);                                 //   +1 fs
  else if(parity==-1)
    XPS.dump(fs,23);                                 //   -1 fs
#endif
  
}



// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::write_to_file(ofstream& fp)
// ++++++++++++++++++++++++++++++++++++++++++++++++++
{
  int i;
  double deltaE=0.;
  // I want to include that additive constant for the energy as well:
  for(int b=0; b<Nb; b++)
    deltaE += epsilonn[b] + hb;
  deltaE=deltaE/ Nb;
  
  for (i=0; i<nobs; i++)  // ignore the fidelity susceptibility- 13 observables
    if (i==1)
      fp << setprecision(15) << -XPS.obs[i]/(beta*Nb) + deltaE << " " ;
    else 
      fp << setprecision(15) << XPS.obs[i] << " " ;
  
  //fp << fixed << std::setw( 11 ) << setprecision( 6 ) 
  // << setfill( '0' ) << XPS.obs[i] << " " ;
  
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::free_arrays()
// ++++++++++++++++++++++++++++++++++++++++++++++++++
{
  delete sm;
  delete vtx;
  delete vb;
  delete linklist;
  free(bsite);
  
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::check_probtable(ofstream& ferr)
// ++++++++++++++++++++++++++++++++++++++++++++++++++
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
	  
          // The sum should be 1. at the end of this and we check to ensure that this is so.              // This is a flagged case, meaning there are round off errors                                                                                                                        
          if( fabs(sum-1.)>zerotol && sum !=0.)
            {
              for(int e=0; e<4; e++)
                ferr << setprecision(15) << "ERROR --> Bond: " << b << " in: " << i <<  " out: " << e << " prob: " << prob_keynumbers[b][i][e][k][0]<< "/" << prob_keynumbers[b][i][e][k][1] <<  \
		  " = " << prob[b][i][e][k] << " " << prob_keynumbers[b][i][e][k][0]/ prob_keynumbers[b][i][e][k][1] << endl;
              flag=false;
            }
        }
  
  if(flag)
    ferr << "No problems with REALIZATION: " << GLOBAL_REAL << " computed by:  " << GLOBAL_RANK << endl;
  else
    ferr << "problems with REALIZATION: " << GLOBAL_REAL << " computed by:  " << GLOBAL_RANK << endl;
  ferr.flush();
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++
void lattice::print_probtable(ofstream& fp)
// ++++++++++++++++++++++++++++++++++++++++++++++++++
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



