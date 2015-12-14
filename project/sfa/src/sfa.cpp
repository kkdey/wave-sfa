#include <stdlib.h> 
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include "model.h"
#include "sfa.h"
#include "files.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_statistics_double.h"
#include "fpmath.h"

using namespace std;

SFA::SFA(ModelnData * m)
{
  mnd = m;
  rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng, mnd->GetrandSeed());
  iterations = 0;
  iters = mnd->GetIters(); // default
}

SFA::~SFA()
{
  gsl_matrix_free(F);
  gsl_matrix_free(lambda);
  gsl_vector_free(alpha);
  gsl_vector_free(eta);
  if(mnd->HasRowMeans() == 1) {
    gsl_vector_free(mug);
  }
  if(mnd->HasColumnMeans() == 1) {
    gsl_vector_free(mun);
  }
  gsl_matrix_free(Ytilde);
  for(uint i = 0; i < LLt.size(); i++) {
    gsl_matrix * m = LLt.at(i);
    gsl_matrix_free(m);
  }
  LLt.clear();

  gsl_matrix_free(kbyk);
  gsl_matrix_free(kbyk2);
  gsl_matrix_free(kbyk3);
  gsl_matrix_free(kbyn);
  gsl_matrix_free(gbyk);
  gsl_matrix_free(Omega_g);
  gsl_vector_free(Nvector);

  // free the random number generator
  gsl_rng_free(rng);
}

void SFA::initialize_sfa(void)
{
  iters = mnd->GetIters();

  cout << "G = " << mnd->get_G() << "\n";
  cout << "N = " << mnd->get_N() << "\n";
  cout << "K = " << mnd->get_K() << "\n";

  cout << "Initializing SFA...\n";
  Ytilde = mnd->get_Y();
  Nvector = gsl_vector_alloc(mnd->get_N());
  kbyk = gsl_matrix_alloc(mnd->get_K(), mnd->get_K());
  kbyk2 = gsl_matrix_alloc(mnd->get_K(), mnd->get_K());
  kbyk3 = gsl_matrix_alloc(mnd->get_K(), mnd->get_K());
  Omega_g = gsl_matrix_alloc(mnd->get_K(), mnd->get_N());
  kbyn = gsl_matrix_alloc(mnd->get_K(), mnd->get_N());
  gbyk = gsl_matrix_alloc(mnd->get_G(), mnd->get_K());

  //F = get_K_rows_at_random(mnd->get_Y(), mnd->get_K(), rng);
  F = gsl_matrix_alloc(mnd->get_K(), mnd->get_N());
  //scale_Y(); // only do this when standardizing Y (removing sd)
  for(int k = 0; k < mnd->get_K(); k++) {
    for(int n = 0; n < mnd->get_N(); n++) {
      gsl_matrix_set(F, k, n, gsl_ran_gaussian(rng, 1.0));
    }
  }
  // initialize lambda to all 1.0
  lambda = gsl_matrix_alloc(mnd->get_G(), mnd->get_K());
  gsl_matrix_set_all(lambda, 0.0);
  for(int g = 0; g < mnd->get_G(); g++) {
    LLt.push_back(gsl_matrix_calloc(mnd->get_K(), mnd->get_K()));
  }
  sigma2 = gsl_matrix_calloc(mnd->get_G(), mnd->get_K());

  for(int g = 0; g < mnd->get_G(); g++) {
    gsl_matrix_set(sigma2, g, (gsl_rng_uniform(rng)*mnd->get_K()), double(mnd->get_K())/double(mnd->get_G()));
    //gsl_matrix_set(sigma2, g, (gsl_rng_uniform(rng)*mnd->get_K()), 1.0);
  }
  alpha = gsl_vector_alloc(mnd->get_G());
  eta = gsl_vector_alloc(mnd->get_N());
  gsl_vector_set_all(alpha, 1.0);
  gsl_vector_set_all(eta, 1.0);
  if(mnd->HasRowMeans() == 1) {
    mug = gsl_vector_alloc(mnd->get_G());
    gsl_vector_set_zero(mug);
    update_mu_rows_sfa(); // initialize mu to the means without the factors/loadings; initialize Ytilde
  }
  if(mnd->HasColumnMeans() == 1) {
    mun = gsl_vector_alloc(mnd->get_N());
    gsl_vector_set_zero(mun);
    update_mu_columns_sfa(); // initialize mu to the means without the factors/loadings; initialize Ytilde
  }
  if(mnd->HasColumnMeans() == 1 && mnd->HasRowMeans() == 1) {
    update_mu_additive_sfa();
  }
  cout << "Finished initializing SFA...\n";
}

// Warning: not tested
void SFA::initialize_from_data_sfa(void)
{
  iters = mnd->GetIters();

  cout << "G = " << mnd->get_G() << "\n";
  cout << "N = " << mnd->get_N() << "\n";
  cout << "K = " << mnd->get_K() << "\n";

  cout << "Initializing SFA...\n";
  Ytilde = mnd->get_Y();

  Nvector = gsl_vector_alloc(mnd->get_N());
  kbyk = gsl_matrix_alloc(mnd->get_K(), mnd->get_K());
  kbyk2 = gsl_matrix_alloc(mnd->get_K(), mnd->get_K());
  kbyk3 = gsl_matrix_alloc(mnd->get_K(), mnd->get_K());
  Omega_g = gsl_matrix_alloc(mnd->get_K(), mnd->get_N());
  kbyn = gsl_matrix_alloc(mnd->get_K(), mnd->get_N());
  gbyk = gsl_matrix_alloc(mnd->get_G(), mnd->get_K());
  
  // read in the parameters
  F = FILEIO::read_matrix("input_F.txt", mnd->get_K(), mnd->get_N());
  lambda = FILEIO::read_matrix("input_lambda.txt", mnd->get_G(), mnd->get_K());
  sigma2 = FILEIO::read_matrix("input_sigma2.txt", mnd->get_G(), mnd->get_K());
  if(mnd->HasRowCovariance() == 1) {
    alpha = FILEIO::read_vector("input_alpha.txt", mnd->get_G());
  }
  if(mnd->HasColumnCovariance() == 1) {
    eta = FILEIO::read_vector("input_eta.txt", mnd->get_N());
  }
  if(mnd->HasRowMeans() == 1) {
    mug = gsl_vector_calloc(mnd->get_G()); // don't bother with the means; just update them here
    update_mu_rows_sfa(); // initialize mu to the means without the factors/loadings; initialize Ytilde
  }
  if(mnd->HasColumnMeans() == 1) {
    mun = gsl_vector_calloc(mnd->get_N()); // don't bother with the means; just update them here
    update_mu_columns_sfa(); // initialize mu to the means without the factors/loadings; initialize Ytilde
  }
  if(mnd->HasColumnMeans() == 1 && mnd->HasRowMeans() == 1) {
    update_mu_additive_sfa();
  }
}

void SFA::update_lambda_sfa()
{
  gsl_vector_view thisY;
  gsl_vector_view thislambda;

  for(int g = 0; g < mnd->get_G();  g++) {
    thislambda = gsl_matrix_row(lambda,g);
    if(vector_sum(&gsl_matrix_row(sigma2, g).vector) == 0.0) {
      gsl_vector_set_zero(&gsl_matrix_row(lambda, g).vector);
      gsl_matrix_set_zero(LLt.at(g)); // added
      continue;
    }
    if(mnd->HasRowCovariance() == 1 && mnd->HasColumnCovariance() == 1) {
      gsl_vector_scale(eta, gsl_vector_get(alpha, g));
    }
    thisY = gsl_matrix_row(Ytilde, g);

    // Compute Omega_g, a KxN vector (using kbyk matrix)
    if(mnd->HasRowCovariance() == 1 && mnd->HasColumnCovariance() == 0) {
      woodbury_matrix_inversion_matrix_scaled(1.0/gsl_vector_get(alpha, g), F, &gsl_matrix_const_row(sigma2, g).vector, F, kbyk, kbyk2, kbyn, Omega_g);
    }
    if(mnd->HasColumnCovariance() == 1) {
      woodbury_matrix_inversion_matrix_scaled(eta, F, &gsl_matrix_const_row(sigma2, g).vector, F, kbyk, kbyk2, kbyn, Omega_g);
    }

    gsl_matrix_set_zero(kbyk2);
    //kbyk2 is sigma_g vector
    for(int k = 0; k < mnd->get_K(); k++) {
      gsl_vector_scale(&gsl_matrix_row(Omega_g, k).vector, gsl_matrix_get(sigma2, g, k));
      gsl_matrix_set(kbyk2, k, k, gsl_matrix_get(sigma2, g, k));
    }

    // update lambda
    gsl_blas_dgemv(CblasNoTrans, 1.0, Omega_g, &thisY.vector, 0.0, &thislambda.vector);
    
    outer_product(&thislambda.vector, &thislambda.vector, kbyk);
    gsl_matrix_memcpy(LLt.at(g), kbyk2); // start wtih sigma2
    gsl_matrix_add(LLt.at(g), kbyk); // add outer product
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, Omega_g, F, 0.0, kbyk3);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, kbyk3, kbyk2, 0.0,  kbyk); 

    gsl_matrix_sub(LLt.at(g), kbyk); // bee: if i remove this, works brilliantly

    if(mnd->HasRowCovariance() == 1 && mnd->HasColumnCovariance() == 1) {
      gsl_vector_scale(eta, 1.0/gsl_vector_get(alpha, g));
    }
  }
}

// eta updates cancel, that's why they're not here.
void SFA::update_F_sfa()
{
  gsl_vector * kby1 = &gsl_matrix_row(kbyk2, 0).vector;

  gsl_matrix_memcpy(kbyk3, LLt.at(0));
  gsl_matrix_scale(kbyk3, gsl_vector_get(alpha,0)); 
  //gsl_vector_scale(&(gsl_matrix_diagonal(kbyk3)).vector, gsl_vector_get(alpha,0)); 
  for(int g=1; g < mnd->get_G();  g++) {
    gsl_matrix_memcpy(kbyk, LLt.at(g));
    gsl_matrix_scale(kbyk, gsl_vector_get(alpha, g));
    gsl_matrix_add(kbyk3, kbyk);
  }
  // don't invert a matrix with zeros along the diagonal
  // I've rarely seen this hit.
  while(gsl_vector_min(&gsl_matrix_diagonal(kbyk3).vector) == 0.0) {
    gsl_matrix_set(kbyk3, gsl_vector_min_index(&gsl_matrix_diagonal(kbyk3).vector), gsl_vector_min_index(&gsl_matrix_diagonal(kbyk3).vector), 0.0001);
  }
  invert_square_matrix(kbyk3, kby1, kbyk);
  gsl_matrix_set_zero(Omega_g);
  for(int g = 0; g < mnd->get_G(); g++) {
    outer_product(&(gsl_matrix_row(lambda, g)).vector, &gsl_matrix_row(Ytilde, g).vector, F);
    gsl_matrix_scale(F, gsl_vector_get(alpha, g));
    gsl_matrix_add(Omega_g, F);
  }
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, kbyk, Omega_g, 0.0, F);

}

// Can we combine psi/mu updates for speed?
void SFA::update_psi_sfa()
{
  if(mnd->HasRowCovariance() == 1) {
    double YtY = 0.0;
    double YgFtLg = 0.0;
    gsl_vector * YgFt = gsl_vector_alloc(mnd->get_K());
    gsl_vector_view Yg;
    double FtLLtF = 0.0;
    
    if(mnd->HasColumnCovariance() == 1) {
      for(int k = 0; k < mnd->get_K(); k++) {
	gsl_vector_memcpy(Nvector, &gsl_matrix_row(F, k).vector);
	gsl_vector_mul(Nvector, eta);
	gsl_blas_dgemv(CblasNoTrans, 1.0, F, Nvector, 0.0, &gsl_matrix_column(kbyk, k).vector);
      }
    }
    else {
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, F, F, 0.0, kbyk); 
    }
    for(int g=0; g < mnd->get_G();  g++) {
      Yg = gsl_matrix_row(Ytilde, g);
      if(mnd->HasColumnCovariance() == 1) {
	gsl_vector_memcpy(Nvector, &Yg.vector);
	gsl_vector_mul(Nvector, eta);
	gsl_blas_ddot(&Yg.vector, Nvector, &YtY);
	gsl_blas_dgemv(CblasNoTrans, 1.0, F, Nvector, 0.0, YgFt);
      } else {
	gsl_blas_ddot(&Yg.vector, &Yg.vector, &YtY);
	gsl_blas_dgemv(CblasNoTrans, 1.0, F, &Yg.vector, 0.0, YgFt);
      }
      gsl_blas_ddot(YgFt, &gsl_matrix_const_row(lambda, g).vector, &YgFtLg);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, kbyk, LLt.at(g), 0.0, kbyk2); 
      
      FtLLtF = vector_sum(&gsl_matrix_const_diagonal(kbyk2).vector);
      //gsl_vector_set(alpha, g, float(mnd->get_N())/(YtY - (2.0*YgFtLg) + FtLLtF));  // bee change to 2.0
      gsl_vector_set(alpha, g, (1.0/(((YtY - (2.0*YgFtLg) + FtLLtF)/float(mnd->get_N()))+0.1))); // bee gamma pseudo prior
    }
    gsl_vector_free(YgFt);
  } 
  if (mnd->HasColumnCovariance() == 1) {
    double YtY = 0.0;
    double YnLFn = 0.0;
    gsl_vector * YnL = gsl_vector_alloc(mnd->get_K());
    gsl_vector_view Yn;
    gsl_vector_view Fn;
    double FtLLtF = 0.0;
    gsl_vector * kby1 = &gsl_matrix_row(kbyk, 0).vector;
    gsl_vector * Gvector = &gsl_matrix_column(gbyk, 0).vector;

    //kbyk2 is sum over alphas and LLts
    if(mnd->HasRowCovariance() == 1) {
      gsl_matrix_memcpy(kbyk2, LLt.at(0));
      gsl_matrix_scale(kbyk2, gsl_vector_get(alpha, 0));
      for(int g = 1; g < mnd->get_G(); g++) {
	gsl_matrix_memcpy(kbyk3, LLt.at(g));
	gsl_matrix_scale(kbyk3, gsl_vector_get(alpha, g));
	gsl_matrix_add(kbyk2, kbyk3);
      }
    }
    else {
      gsl_matrix_memcpy(kbyk2, LLt.at(0));
      for(int g = 1; g < mnd->get_G(); g++) {
	gsl_matrix_add(kbyk2, LLt.at(g));
      }
    }
    for(int n=0; n < mnd->get_N();  n++) {
      Yn = gsl_matrix_column(Ytilde, n);
      Fn = gsl_matrix_column(F, n);
      if(mnd->HasRowCovariance() == 1) {
	gsl_vector_memcpy(Gvector, &Yn.vector);
	gsl_vector_mul(Gvector, alpha);
	gsl_blas_ddot(&Yn.vector, Gvector, &YtY);
	gsl_blas_dgemv(CblasTrans, 1.0, lambda, Gvector, 0.0, YnL);
	
      } else {
	gsl_blas_ddot(&Yn.vector, &Yn.vector, &YtY);
	gsl_blas_dgemv(CblasTrans, 1.0, lambda, &Yn.vector, 0.0, YnL);
      }
      gsl_blas_ddot(YnL, &Fn.vector, &YnLFn);
      gsl_blas_dgemv(CblasNoTrans, 1.0, kbyk2, &Fn.vector, 0.0, kby1); 
      gsl_blas_ddot(kby1, &Fn.vector, &FtLLtF);
      
      //gsl_vector_set(eta, n, float(mnd->get_G())/((YtY - (2.0*YnLFn) + FtLLtF))); // bee changed from 2.0 to 1.5
      gsl_vector_set(eta, n, (1.0/(((YtY - (2.0*YnLFn) + FtLLtF)/float(mnd->get_G()))+0.1))); // bee gamma pseudo prior
    }
    gsl_vector_free(YnL);
  }
  if (mnd->HasColumnCovariance() == 1 && mnd->HasRowCovariance() == 1) {
    float alphamax = gsl_vector_max(alpha);
    float alphamin = gsl_vector_min(alpha);
    float etamax = gsl_vector_max(eta);
    float etamin = gsl_vector_min(eta);
    cout << "Alpha and Eta " << alphamax << "," << alphamin << " " << etamax << "," << etamin << '\n';
    if(mnd->get_G() < mnd->get_N()) {
      if(alphamax - alphamin > 3.0) {
	cout << "Rescaling residual variance (alpha) " << (alphamax-alphamin) << '\n';
	gsl_vector_scale(alpha, 3.0/(alphamax-alphamin));
	gsl_vector_scale(eta, (alphamax-alphamin)/3.0);
      }
    } else {
      if(etamax - etamin > 3.0) {
	cout << "Rescaling residual variance (eta) " << (etamax-etamin) << '\n';
	gsl_vector_scale(eta, 3.0/(etamax-etamin));
	gsl_vector_scale(alpha, (etamax-etamin)/3.0);
      }
    }
  }  
}

void SFA::update_sigma2_sfa(int iterations) 
{
  double q2 = 0.0;
  float sigma2gk = 0.0;
  // check to see if any elements of sigma2 are zero.
  // if so, keep around beta matrix to use for all of those sigma2_g,k
  gsl_vector * sigma2p = gsl_vector_alloc(mnd->get_K());

  gsl_vector * Qm = gsl_vector_alloc(mnd->get_K());
  double sm = 0.0;
  double qm = 0.0;
  int iteration_cutoff = 5;

  for(int g = 0; g < mnd->get_G(); g++) {
    gsl_vector_memcpy(sigma2p, &gsl_matrix_const_row(sigma2, g).vector);
    if(mnd->HasRowCovariance() == 1 && mnd->HasColumnCovariance() == 1) {
      gsl_vector_scale(eta, gsl_vector_get(alpha, g));
    }
    
    if(iterations >= iteration_cutoff) {
      if(mnd->HasRowCovariance() == 1 && mnd->HasColumnCovariance() == 0) {
	woodbury_matrix_inversion_matrix_scaled(1.0/gsl_vector_get(alpha, g), F, sigma2p, F, kbyk, kbyk2, kbyn, Omega_g);
      }
      if(mnd->HasColumnCovariance() == 1) {
	woodbury_matrix_inversion_matrix_scaled(eta, F, sigma2p, F, kbyk, kbyk2, kbyn, Omega_g);
      }
      // mutiply by Y/F
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, Omega_g, F, 0.0, kbyk3);
      gsl_blas_dgemv(CblasNoTrans, 1.0, Omega_g, &gsl_matrix_const_row(Ytilde, g).vector, 0.0, Qm);
    }

    for(int k = 0; k < mnd->get_K(); k++) {
      if(iterations < iteration_cutoff) {
	if(mnd->HasRowCovariance() == 1 && mnd->HasColumnCovariance() == 0) {
	  woodbury_matrix_inversion_matrix_scaled(1.0/gsl_vector_get(alpha, g), F, sigma2p, F, kbyk, kbyk2, kbyn, Omega_g);
	}
	if(mnd->HasColumnCovariance() == 1) {
	  woodbury_matrix_inversion_matrix_scaled(eta, F, sigma2p, F, kbyk, kbyk2, kbyn, Omega_g);
	}
	// mutiply by Y/F
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, Omega_g, F, 0.0, kbyk3);
	gsl_blas_dgemv(CblasNoTrans, 1.0, Omega_g, &gsl_matrix_const_row(Ytilde, g).vector, 0.0, Qm);
      }

      sigma2gk = gsl_vector_get(sigma2p, k);
      sm = gsl_matrix_get(kbyk3, k, k);
      qm = gsl_vector_get(Qm, k);
      if(sigma2gk == 0.0) {
	q2 = qm*qm;
      } else {
	q2 = (1.0/sigma2gk)*qm/((1.0/sigma2gk) - sm);
	q2 = q2*q2;
	sm = (1.0/sigma2gk)*sm/((1.0/sigma2gk) - sm);
      }
      sigma2gk = (q2-sm)/(sm*sm);
      if(sigma2gk < 0.0) {
	sigma2gk = 0.0;
      }
      gsl_vector_set(sigma2p, k, sigma2gk);
    }
    gsl_matrix_set_row(sigma2, g, sigma2p);
    if(mnd->HasRowCovariance() == 1 && mnd->HasColumnCovariance() == 1) {
      gsl_vector_scale(eta, 1.0/gsl_vector_get(alpha, g));
    }
  }
  // If a column of sigma2 is all zero, put a small value there so there won't be NANs.
  for(int k = 0; k < mnd->get_K(); k++) {
    if(vector_sum(&gsl_matrix_const_column(sigma2, k).vector) == 0.0) {
      gsl_vector_add_constant(&gsl_matrix_column(sigma2, k).vector, double(mnd->get_K())/double(10*mnd->get_G()));
    }
  }
  gsl_vector_free(sigma2p);
}

double SFA::scale_F_sigma2_sfa()
{
  gsl_vector * kby1 = &gsl_matrix_row(kbyk2, 0).vector;
  compute_sd_rows(F, &gsl_matrix_row(kbyk, 0).vector, kby1);
  FILEIO::print_vector(kby1);
  // Don't let F be divided by zero
  while(gsl_vector_min(kby1) == 0.0) {
    gsl_vector_set(kby1, gsl_vector_min_index(kby1), 1.0);
  }
  for(int i = 0; i < mnd->get_K(); i++) {
    gsl_vector_scale(&gsl_matrix_row(F, i).vector, 1.0/gsl_vector_get(kby1, i));
    gsl_vector_scale(&gsl_matrix_column(sigma2, i).vector, gsl_vector_get(kby1, i)*gsl_vector_get(kby1, i)); 
  }
  double abs_sum = 0.0;
  for(int k = 0; k < mnd->get_K(); k++) {
    abs_sum = (1.0 - gsl_vector_get(kby1, k));
    abs_sum = (abs_sum < 0.0)?-abs_sum:abs_sum;
  }
  return(abs_sum);
}

// remove mean/standard deviation
void SFA::scale_Y()
{
  gsl_vector * temp = &gsl_matrix_row(F,0).vector;
  compute_sd_columns(Ytilde, temp, Nvector);
  FILEIO::print_vector(&gsl_vector_subvector(temp, 0, 100).vector);
  for(int i = 0; i < mnd->get_N(); i++) {
    gsl_vector_add_constant(&gsl_matrix_column(Ytilde,i).vector, -gsl_vector_get(temp,i));
    //gsl_vector_scale(&gsl_matrix_column(Ytilde, i).vector, 1.0/gsl_vector_get(Nvector, i));
  }
  compute_sd_columns(Ytilde, temp, Nvector);
  FILEIO::print_vector(&gsl_vector_subvector(temp, 0, 100).vector);
}

void SFA::update_mu_rows_sfa()
{
  double oldmu = 0.0;
  double newmu = 0.0;
  gsl_vector * thisY = NULL;
  double denom = vector_sum(eta); // added

  for(int g = 0; g < mnd->get_G(); g++) {
    oldmu = gsl_vector_get(mug,g);
    thisY = &gsl_matrix_row(Ytilde, g).vector;
    gsl_blas_dgemv(CblasTrans, 1.0, F, &gsl_matrix_const_row(lambda, g).vector, 0.0, Nvector);
    gsl_vector_add_constant(thisY, oldmu);
    gsl_vector_sub(Nvector, thisY);
    gsl_vector_mul(Nvector, eta); // added
    newmu = (-vector_sum(Nvector)/denom); // changed
    gsl_vector_set(mug, g, newmu);
    
    // update Ytilde
    gsl_vector_add_constant(thisY, -newmu);
  }
}

void SFA::update_mu_columns_sfa()
{
  double oldmu = 0.0;
  double newmu = 0.0;
  gsl_vector * thisY = NULL;
  gsl_vector * Gvector = &gsl_matrix_column(gbyk, 0).vector;
  double denom = vector_sum(alpha);

  for(int n = 0; n < mnd->get_N(); n++) {
    oldmu = gsl_vector_get(mun,n);
    thisY = &gsl_matrix_column(Ytilde, n).vector;
    gsl_blas_dgemv(CblasNoTrans, 1.0, lambda, &gsl_matrix_const_column(F, n).vector, 0.0, Gvector);
    gsl_vector_add_constant(thisY, oldmu);
    gsl_vector_sub(Gvector, thisY);
    gsl_vector_mul(Gvector, alpha);
    newmu = (-vector_sum(Gvector)/denom);
    gsl_vector_set(mun, n, newmu);
    
    // update Ytilde
    gsl_vector_add_constant(thisY, -newmu);
  }
}

// center one mu term, rescale the others.
void SFA::update_mu_additive_sfa()
{
  double rowmm = (1.0/(double)mnd->get_G())*vector_sum(mug);
  gsl_vector_add_constant(mug,-rowmm);
  gsl_vector_add_constant(mun,rowmm);
}


void SFA::write_matrices(string path)
{
  string filename = path+"_F.out";
  FILEIO::write_matrix(filename, F);
  filename = path+"_lambda.out";
  FILEIO::write_matrix(filename, lambda);
  filename = path+"_sigma2.out";
  FILEIO::write_matrix(filename, sigma2);
  if(mnd->HasRowCovariance() == 1) {
    filename = path+"_alpha.out";
    FILEIO::write_vector(filename, alpha);
  }
  if(mnd->HasColumnCovariance() == 1) {
    filename = path+"_eta.out";
    FILEIO::write_vector(filename, eta);
  }
  if(mnd->HasRowMeans() == 1) {
    filename = path+"_mug.out";
    FILEIO::write_vector(filename, mug);
  }
  if(mnd->HasColumnMeans() == 1) {
    filename = path+"_mun.out";
    FILEIO::write_vector(filename, mun);
  }
}

void SFA::single_iteration(string path)
{
  float SMALL = -1e+38;
  float ll = SMALL;
  double scale = 0.0;
  float prev_ll = SMALL;

  for(int i=0; i < iters; i++) {
    // Sample factor loadings: Lambda (sparse)
    cout << "updating lambda\n";
    update_lambda_sfa();
    if(mnd->HasRowCovariance() == 1) {
      cout << "updating rows psi\n";
    }
    if(mnd->HasColumnCovariance() == 1) {
      cout << "updating columns psi\n";
    }
    update_psi_sfa();
    cout << "updating sigma2\n";
    update_sigma2_sfa(iterations);
    cout << "updating F\n";
    update_F_sfa();
    cout << "scaling sigma2 and F\n";
    scale = scale_F_sigma2_sfa();
    if(mnd->HasRowMeans() == 1) {
      cout << "updating mu rows\n";
      update_mu_rows_sfa();
    }
    if(mnd->HasColumnMeans() == 1) {
      cout << "updating mu columns\n";
      update_mu_columns_sfa();
    }
    if(mnd->HasColumnMeans() == 1 && mnd->HasRowMeans() == 1) {
      update_mu_additive_sfa();
    }
    // you can check the marginal log likelihood whenever
    // but it's time consuming, so every 10 iterations seems okay.
    if (iterations % 10 == 0) {
      cout << "computing ll\n";    
      prev_ll = ll;
      ll = marginal_log_likelihood();
      //ll = log_likelihood();
      cout << "Marginal log likelihood at iteration " << iterations << ": " << ll << "\n";
      //cout << "scale: " << (int)(scale * 10000000) << '\n';

      // This is a conservative convergence criterion. If time is critical, change this.
      if((int)(scale*10000000) == 0 && prev_ll >= ll && ll > SMALL) {
	cout << "Converged\n";
      	break;
      }
    }
    iterations += 1;
  }
  cout << "Expected complete Log likelihood at iteration " << iterations << ": " << log_likelihood() << "\n";
  cout << "Marginal log likelihood at iteration " << iterations << ": " << marginal_log_likelihood() << "\n";
  cout << "Residual variance at iteration " << iterations << ": " << residual_variance() << "\n";
  cout << "Residual sum of squares at iteration " << iterations << ": " << residual_sum_squares() << "\n";
  write_matrices(path);
}

// Returns expected complete log likelihood
float SFA::log_likelihood()
{
  float ll = 0.0;
  double onediag = 0.0;
  float psiterm = 0.0;
  gsl_vector * thisY = NULL;
  double betapsirows = 20.0/(double)mnd->get_N();
  double betapsicolumns = 20.0/(double)mnd->get_G();
  double psiprior = 0.0;

  for(int g = 0; g < mnd->get_G(); g++) {
    thisY = &gsl_matrix_row(Ytilde, g).vector;
    gsl_blas_dgemv(CblasTrans, 1.0, F, &gsl_matrix_const_row(lambda, g).vector, 0.0, Nvector);
    gsl_vector_sub(Nvector, thisY);
    gsl_vector_mul(Nvector, Nvector);
    gsl_blas_ddot(Nvector, eta, &onediag);
    ll -= (0.5 * onediag * gsl_vector_get(alpha, g));
    psiterm += (mnd->get_N()) * log(gsl_vector_get(alpha, g));
    if(mnd->HasRowCovariance() == 1) {
      psiprior += (gsl_vector_get(alpha,g)/betapsirows);
    }
  }
  if(mnd->HasColumnCovariance() == 1) {
    for(int n= 0; n < mnd->get_N(); n++) {
      psiterm += (mnd->get_G()) * log(gsl_vector_get(eta,n));
      psiprior += (gsl_vector_get(eta,n) / betapsicolumns);
    }
  }
  return(((-(float)(mnd->get_N()*mnd->get_G())/2.0)*log(2*pi)) + 0.5 * psiterm + ll + psiprior);
}

// Returns marginal log likelihood
float SFA::marginal_log_likelihood()
{
  float total = 0.0;
  gsl_vector * thisY = NULL;
  gsl_vector * thislambda = NULL;
  gsl_vector * nby1 = &gsl_matrix_row(Omega_g, 0).vector;
  double YtBY = 0.0;

  double logdet2 = 1.0;
  double detpsi = 1.0;
  // precompute det(psi_n^-1)                                                                                                                                                                                              
  if(mnd->HasColumnCovariance() == 1) {
    // this is not quite right.
    for(int n = 0; n < mnd->get_N(); n++) {
      detpsi *= gsl_vector_get(eta, n);
    }
  }

  for(int g = 0; g < mnd->get_G(); g++) {
    thislambda = &gsl_matrix_row(lambda,g).vector;
    thisY = &gsl_matrix_row(Ytilde, g).vector;

    if(mnd->HasRowCovariance() == 1 && mnd->HasColumnCovariance() == 1) {
      gsl_vector_scale(eta, gsl_vector_get(alpha, g));
    }
    // compute (F^t\Sigma_gF + \Psi^{-1}_g)^-1
    if(mnd->HasRowCovariance() == 1 && mnd->HasColumnCovariance() == 0) {
      woodbury_matrix_inversion_vector_scaled(1.0/gsl_vector_get(alpha, g), F, &gsl_matrix_const_row(sigma2, g).vector, thisY, kbyk, kbyk2, kbyn, nby1);
    }
    if(mnd->HasColumnCovariance() == 1) {
      woodbury_matrix_inversion_vector_scaled(eta, F, &gsl_matrix_const_row(sigma2, g).vector, thisY, kbyk, kbyk2, kbyn, nby1);
    }
    gsl_blas_ddot(nby1, thisY, &YtBY);

    if(mnd->HasRowCovariance() == 1) {
      detpsi *= (1.0/gsl_vector_get(alpha, g));
    }

    //Sylvester's determinant theorem: det(X + AB) = det(X)det(In + BX^{−1}A),                                                                                                                                             
    gsl_matrix_memcpy(Omega_g, F);
    if(mnd->HasRowCovariance() == 1) {
      gsl_matrix_scale(Omega_g, gsl_vector_get(alpha, g));
    }
    if(mnd->HasColumnCovariance() == 1) {
      for(int k = 0; k < mnd->get_K(); k++) {
        gsl_vector_mul(&gsl_matrix_row(Omega_g, k).vector, eta);
      }
    }
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, Omega_g, F, 0.0, kbyk);
    for(int k = 0; k < mnd->get_K(); k++) {
      gsl_vector_scale(&gsl_matrix_row(kbyk, k).vector, gsl_matrix_get(sigma2, g, k));
      gsl_matrix_set(kbyk, k, k, gsl_matrix_get(kbyk, k, k)+1.0);
    }

    logdet2 = matrix_determinant(kbyk, kbyk2);
    total += log(detpsi*logdet2) + YtBY;

    // put things right
    if(mnd->HasRowCovariance() == 1 && mnd->HasColumnCovariance() == 1) {
      gsl_vector_scale(eta, 1.0/gsl_vector_get(alpha, g));
    }
    if(mnd->HasRowCovariance() == 1) {
      detpsi = detpsi/gsl_vector_get(alpha, g);
    }
  }
  return(-total);
}

double SFA::residual_variance()
{
  double newmu = 0.0;
  gsl_vector * thisY = NULL;
  for(int g = 0; g < mnd->get_G(); g++) {
    thisY = &gsl_matrix_row(Ytilde, g).vector;
    gsl_blas_dgemv(CblasTrans, 1.0, F, &gsl_matrix_const_row(lambda, g).vector, 0.0, Nvector);
    gsl_vector_sub(Nvector, thisY);
    newmu -= (1.0/(float)mnd->get_N())*vector_sum(Nvector);
  }
  newmu = newmu/(double)mnd->get_G();
  double variance = 0.0;
  gsl_matrix_add_constant(Ytilde, -newmu);
  double oneY = 0.0;
  for(int g = 0; g < mnd->get_G(); g++) {
    for(int n = 0; n < mnd->get_N(); n++) {
      oneY = gsl_matrix_get(Ytilde, g, n);
      variance += (oneY * oneY);
    }
  }
  gsl_matrix_add_constant(Ytilde, newmu);

  variance = variance/(double)(mnd->get_G()*mnd->get_N());
  return(variance);
}

double SFA::residual_sum_squares()
{
  double newmu = 0.0;
  double rss = 0.0;
  gsl_vector * thisY = NULL;
  for(int g = 0; g < mnd->get_G(); g++) {
    thisY = &gsl_matrix_row(Ytilde, g).vector;
    gsl_blas_dgemv(CblasTrans, 1.0, F, &gsl_matrix_const_row(lambda, g).vector, 0.0, Nvector);
    gsl_vector_sub(Nvector, thisY);
    gsl_blas_ddot(Nvector, Nvector, &newmu);
    rss += newmu;
  }
  return(rss);
}
