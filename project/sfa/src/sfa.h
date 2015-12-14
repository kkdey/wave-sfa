#ifndef __SFA__
#define __SFA__

#include <string>
#include <map>
#include <vector>
#include "gsl/gsl_rng.h"

using namespace std; 

class SFA
{
 public:
  // data
  int iters;
  gsl_rng * rng;
  
  // initialize data structures
  gsl_matrix * lambda;
  gsl_matrix * F;
  gsl_matrix * sigma2;
  gsl_vector * mug;
  gsl_vector * mun;
  gsl_vector * psi;
  gsl_vector * alpha;
  gsl_vector * eta;
  vector<gsl_matrix*> LLt;

  // temp variables
  gsl_matrix * Ytilde; // G x N matrix
  gsl_vector * Nvector;
  gsl_matrix * kbyk;
  gsl_matrix * kbyk2;
  gsl_matrix * kbyk3;
  gsl_matrix * Omega_g; // K x N matrix
  gsl_matrix * kbyn; // another K x N matrix
  gsl_matrix * gbyk; // G x K matrix
  // model
  ModelnData * mnd;
  int iterations;
  
  // constructor/destructor
  SFA(ModelnData * m);
  ~SFA();

  void initialize_sfa(void);
  void initialize_from_data_sfa(void);
  void update_lambda_sfa(void);
  void update_F_sfa(void);
  void update_sigma2_sfa(int iterations);
  double scale_F_sigma2_sfa();
  void scale_Y();
  void update_psi_sfa(void);
  void update_mu_rows_sfa(void);
  void update_mu_columns_sfa(void);
  void update_mu_additive_sfa(void);
  double objective_function(void);
  void write_matrices(string path);
  void single_iteration(string path);
  float log_likelihood();
  float marginal_log_likelihood();
  void test_matrix_inverse();
  double residual_variance();
  double residual_sum_squares();

};

#endif
