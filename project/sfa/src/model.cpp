#include "model.h"
#include "files.h"
#include "sfa.h"

ModelnData::ModelnData(void)
{
  randSeed = 0;
  gsl_r = gsl_rng_alloc(gsl_rng_taus);
  rng = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(rng, GetrandSeed());
  rowMeans = 0;
  columnMeans = 0;
  rowCovariance = 0;
  columnCovariance = 0;	
  G = 1; 
  // default parameters
  iters = 20;
  transposep = 0;
  fnOutput.assign("\0");
}

ModelnData::~ModelnData()
{
	vGin.clear();
	fplog.str("\0"); 	
	//cout << "model destructor being called" << endl; 
}

real ModelnData::sfa_model(void)
{
  // read in matrix
  SFA * gs = new SFA(this);
  cout << "Trying to read in matrix with G=" << G << " and N =" << N << "\n";
  if(Getnum_Gin() > 0) {
    if(transposep == 1) {
      Y = FILEIO::read_matrix(vGin[0], N, G);
      Y = transpose(Y);
    } else {
      Y = FILEIO::read_matrix(vGin[0], G, N);
    }
    cout << "read in matrix " << Y->size1 << " by " << Y->size2 << "\n";
    cout << "initializing sfa...\n";
    gs->initialize_sfa();
  } else {
    cout << "Please enter a matrix Y\n";
    exit(1);
  }  
  // set covariance default here                                                                                                                                                      
  if(HasRowCovariance() == 0 && HasColumnCovariance() == 0) {
    SetHasRowCovariance();
  }
  // run some iterations of updates
  gs->single_iteration(GetfnOutput());
  // output the model
  return(1.0);
}

