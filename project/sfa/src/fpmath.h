#ifndef __FPMATH_H__
#define __FPMATH_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_multimin.h>
//#include <gsl/gsl_blas.h>
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_randist.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_randist.h"
#include "fp.h"

#define invsr2pi 0.3989423
#define TINY     1e-100
#define pi       3.141592653589


using namespace std;

typedef double real;
//globle random number generator; 
extern const gsl_rng_type * gslType;

class Pair {
	public:	
		real bf; 
		int pos;
		int label; 
};

int imap(int, int, int); 
string getline(streambuf * ); 
void DistinctIntArray(int, int, int, int *); 		// generate arg3 many distinct intergers on [arg1, arg2); 
//void lu_decomp(real ** a, int n, int *indx, int *d);
//void lu_back_sub(real ** a, int n, int *indx, real b[]);
void safe_exit(); 
int compare_pair(const void * a, const void * b); 
int compare(const void * a, const void * b);

void * Allocate1D(size_t, int);
void ** Allocate2D(size_t, int, int);  
void *** Allocate3D(size_t, int, int, int);  
void Free1D(void *);
void Free2D(void **); 
void Free3D(void ***); 

real ** 	Allocate2DMatrix(int, int);
int ** 		Allocate2DIntMatrix(int, int);
short int** Allocate2DShortMatrix(int, int);
char ** 	Allocate2DCharMatrix(int, int);
real *** 	Allocate3DMatrix(int, int, int);
void 	Free2DCharMatrix(char **);
void 	Free2DIntMatrix(int **);
void 	Free2DShortMatrix(short int **);
void 	Free2DMatrix(real **);
void 	Free3DMatrix(real ***);

// helper functions
void matrix_multiply_cov(gsl_matrix * M, gsl_vector * sigma2, gsl_matrix * out);
gsl_matrix * diag(gsl_vector * v);
gsl_vector * invert_vector(gsl_vector * v);
gsl_vector * pow_vector(gsl_vector * v, float p);
void pow_matrix_inplace(gsl_matrix * v, float p);
void mult_rows_matrix_inplace(gsl_matrix * m, gsl_vector * v);
void pow_vector_inplace(gsl_vector * v, float p);
void invert_vector_inplace(gsl_vector * v);
gsl_vector * matrix_row(gsl_matrix * m, int r);
gsl_vector * matrix_column(gsl_matrix * m, int c);
float vector_sum(gsl_vector * v);
float mean(gsl_vector * v);
float var(gsl_vector * v);
void compute_sd_rows(gsl_matrix * m, gsl_vector * v, gsl_vector * out);
void compute_sd_columns(gsl_matrix * m, gsl_vector * v, gsl_vector * out);
//void invert_psd_matrix(gsl_matrix * m);
void invert_square_matrix(gsl_matrix * m, gsl_vector * kby1, gsl_matrix * out);
void sample_multivariate_normal(gsl_vector * mean, gsl_matrix * cov, gsl_rng * rng);
gsl_vector * random_unif_vector(int count, double min, double max, gsl_rng * rng);
gsl_matrix * matrix_correlation(gsl_matrix * m);
gsl_matrix * cholesky_decomposition(gsl_matrix * m);
gsl_matrix * get_eigenvectors(gsl_matrix * m);
gsl_matrix * get_first_k_eigenvectors(gsl_matrix * m, int k);
gsl_vector * new_vector(float x, int size);
gsl_vector * row_average(gsl_matrix * m);
gsl_vector * column_average(gsl_matrix * m);
gsl_matrix * transpose(gsl_matrix * m);
void add_to_matrix(gsl_matrix * m, float n);
bool matrix_has_nan(gsl_matrix * m);
gsl_matrix * get_K_columns_at_random(gsl_matrix * m, int K, gsl_rng * rng);
gsl_matrix * get_K_rows_at_random(gsl_matrix * m, int K, gsl_rng * rng);
int get_max_index(gsl_vector * v);
//void woodbury_matrix_inversion(float A, gsl_matrix * U, gsl_vector * C, gsl_matrix * out);
void woodbury_matrix_inversion_matrix_scaled(float A, gsl_matrix * U, gsl_vector * C, gsl_matrix * F, gsl_matrix * kbyk, gsl_matrix * kbyk2, gsl_matrix * kbyn, gsl_matrix * out);
void woodbury_matrix_inversion_vector_scaled(float A, gsl_matrix * U, gsl_vector * C, gsl_vector * Fk, gsl_matrix * kbyk, gsl_matrix * kbyk2, gsl_matrix * kbyn,  gsl_vector * out);
void woodbury_matrix_inversion_matrix_scaled(gsl_vector * Ainv, gsl_matrix * U, gsl_vector * C, gsl_matrix * F, gsl_matrix * kbyk, gsl_matrix * kbyk2, gsl_matrix * kbyn, gsl_matrix * out);
void woodbury_matrix_inversion_vector_scaled(gsl_vector * Ainv, gsl_matrix * U, gsl_vector * C, gsl_vector * Fk, gsl_matrix * kbyk, gsl_matrix * kbyk2, gsl_matrix * kbyn, gsl_vector * out);
void outer_product(gsl_vector * a, gsl_vector * b, gsl_matrix * out);
void ADAt(gsl_matrix * A, gsl_matrix * A2, gsl_vector * d, gsl_vector * Nby1, gsl_matrix * kbykp);
void dAv(gsl_vector * d, gsl_matrix * A, gsl_vector * v, gsl_vector * Nby1);
void Ad(gsl_matrix * A, gsl_vector * v);
double matrix_determinant(gsl_matrix * m, gsl_matrix * outm);

#endif
