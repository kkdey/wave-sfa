#include "fpmath.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <sys/types.h>
#include "files.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_statistics_double.h"

#if defined (MPI_ENABLED)
#include "mpi.h"
#endif

using namespace std;

string getline(streambuf * pbuf)
{
    char ch;
	string str;
	size_t pos; 
	while((ch = pbuf->sgetc()) != EOF)
	{
		if(ch != '\n' && ch != '\r')
		{
			str.push_back(ch);
			ch = pbuf->snextc();
		}
		else {
			pbuf->sbumpc();  //chomp;
			pos = str.find_first_not_of(";, \t", 0); 
			if(str.empty() || pos  == string::npos || str.at(pos) == '#')
			{
				str.clear();
				continue;
			}
			else
				break;
		}
	}
	return str;
}   //this getline use ;, \t as delimit, and ignore the lines either full of delimit or starting with #. 

int imap(int nK, int x, int y)
{
	int tmax = max(x,y); 
	int tmin = min(x,y); 
	int res = (tmin) * nK - (tmin * tmin - tmin)/2; 
	res += (tmax - tmin);  
	return (res); 
}

void safe_exit() 
{
#if defined (MPI_ENABLED)
	MPI_Barrier(MPI_COMM_WORLD); 
	MPI_Finalize(); 
#endif
	exit(0); 
}

int compare_pair(const void * a, const void * b)
{
	Pair * ta = (Pair *) a;
	Pair * tb = (Pair *) b; 
	if((tb->bf) > (ta->bf)) return 1;
	else if((tb->bf) == (ta->bf)) return 0;
	else return -1;
}

int compare(const void * a, const void * b)
{
	int *ta = (int *) a;
	int *tb = (int *) b;
	if ((*ta) > (*tb)) return 1;
	else if ((*ta) == (*tb)) return 0;
	else return -1;
}

void * Allocate1D(size_t us, int dim)
{
	size_t size = dim * us; 
	char * m = (char *) malloc(size);
	memset(m, 0, size); 
	return (void *) m; 
}

void Free1D(void * m)
{
	if(m == NULL) return; 
	free(m); 
	m = NULL; 
}

void ** Allocate2D(size_t us,  int dim1, int dim2)
{
	char ** m;
	m = (char **) malloc((size_t)(dim1 * us));
	m[0] = (char *) malloc((size_t)(dim1 *dim2 * us));
   	memset(m[0], 0, (size_t) (dim1 * dim2 * us));  
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D real matrix. \n");
		safe_exit();
	}
	for(int i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + (size_t) (dim2 * us);
	}
	return ((void **) (m));
}

void Free2D(void ** m)
{
	if(m == NULL) return; 
	free(m[0]);
	free(m);
	m = NULL; 
}

real ** Allocate2DMatrix(int dim1, int dim2)
{
	int i;
	real ** m;
	
	m = (real **) malloc((size_t)((dim1)*sizeof(real*)));
	m[0] = (real *) malloc((size_t)((dim1*dim2)*sizeof(real)));
	memset(m[0], 0, (dim1*dim2)*sizeof(real)); 
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D real matrix. \n");
		safe_exit();
	}
	for(i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
	}
	return (m);
}

//
//void *** Allocate3DMatrix(size_t us, int dim1, int dim2, int dim3)
//{
//	void *** m; 
//	m = (void ***) malloc((size_t)(dim1 * us));
//	m[0] = (void **) malloc((size_t)(dim1 * dim2 * us));
//	m[0][0] = (void *) malloc((size_t)(dim1 * dim2 * dim3 * us));
//	if (!(m && m[0] &&  m[0][0]))
//	{
//		printf("Error: Problem allocating a 3D real matrix. \n");
//		safe_exit();
//	}
//
//	for (int j = 1; j < dim2; j++)
//	{
//		m[0][j] = m[0][j-1] + dim3;
//	}
//	for (int i = 1; i < dim1; i++)
//	{
//		m[i] = m[i-1] + dim2;
//		m[i][0] = m[i-1][dim2 - 1] + dim3;
//		for(int j = 1; j < dim2; j++)
//		{
//			m[i][j] = m[i][j-1] + dim3;
//		}
//	}
//	return (m);
//}
//
//void Free3DMatrix(void *** m)
//{
//	if(m == NULL) return;
//	free(m[0][0]);
//	free(m[0]);
//	free(m);
//	m = NULL; 
//}

real *** Allocate3DMatrix(int dim1, int dim2, int dim3)
{
	int i, j;
	real *** m; 

	m = (real ***) malloc((size_t)((dim1)*sizeof(real**)));
	m[0] = (real **) malloc((size_t)((dim1) * (dim2) * sizeof(real *)));
	m[0][0] = (real *) malloc((size_t)((dim1) * (dim2) * (dim3) * sizeof(real)));
	memset(m[0][0], 0, (dim1) * (dim2) * (dim3) * sizeof(real));  
	if (!(m && m[0] &&  m[0][0]))
	{
		printf("Error: Problem allocating a 3D real matrix. \n");
		safe_exit();
	}

	for (j = 1; j < dim2; j++)
	{
		m[0][j] = m[0][j-1] + dim3;
	}
	for (i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
		m[i][0] = m[i-1][dim2 - 1] + dim3;
		for(j = 1; j < dim2; j++)
		{
			m[i][j] = m[i][j-1] + dim3;
		}
	}
	return (m);
}

void Free2DMatrix(real ** m)
{
	if(m == NULL) return; 
	free(m[0]);
	free(m);
	m = NULL; 
}

void Free3DMatrix(real *** m)
{
	if(m == NULL) return;
	free(m[0][0]);
	free(m[0]);
	free(m);
	m = NULL; 
}
 

int ** Allocate2DIntMatrix(int dim1, int dim2)
{
	int i;
	int ** m;
	
	m = (int **) malloc((size_t)((dim1)*sizeof(int*)));
	m[0] = (int *) malloc((size_t)((dim1*dim2)*sizeof(int)));
	memset(m[0], 0, (dim1*dim2)*sizeof(int)); 
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D int matrix. \n");
		safe_exit();
	}
	for(i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
	}
	return (m);
}

short int ** Allocate2DShortMatrix(int dim1, int dim2)
{
	int i;
	short int ** m;
	
	m = (short int **) malloc((size_t)((dim1)*sizeof(short int*)));
	m[0] = (short int *) malloc((size_t)((dim1*dim2)*sizeof(short int)));
	memset(m[0], 0, (dim1*dim2)*sizeof(short)); 
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D int matrix. \n");
		safe_exit();
	}
	for(i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
	}
	return (m);
}    

void Free2DShortMatrix(short int ** m)
{
	if(m == NULL) return;
	free(m[0]);
	free(m);
	m = NULL; 
}

void Free2DIntMatrix(int ** m)
{                        
	if(m == NULL) return; 
	free(m[0]);
	free(m);
	m = NULL; 
}  

char  ** Allocate2DCharMatrix(int dim1, int dim2)
{
	int i;
	char ** m;
	
	m = (char **) malloc((size_t)((dim1)*sizeof(char*)));
	m[0] = (char *) malloc((size_t)((dim1*dim2)*sizeof(char)));
	memset(m[0], 0, (dim1*dim2)*sizeof(char)); 
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D char matrix. \n");
		safe_exit();
	}
	for(i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
	}
	return (m);
}    


void Free2DCharMatrix(char ** m)
{
	if(m == NULL) return; 
	free(m[0]);
	free(m);
	m = NULL; 
}  

void DistinctIntArray(int low, int high, int dim, int * A, gsl_rng * rng)
{
  if (dim >= high - low)
		return;  
	
	int i, j, r; 
	int bingle; 
	for (i=0; i<dim; i++)
		A[i] = -1;
	
	int howmany = 0;
	for (i=high - dim; i<high; i++)
	{
		bingle = 0;
		r = low + gsl_rng_uniform_int(rng, i-low);
		for (j = 0; j < howmany; j++)
		{
			if (r == A[j])
			{
				bingle = 1; 
				break;
			}
		}

		if (bingle) A[howmany] = i;
		else A[howmany] = r; 
		howmany++;
	}
}

// computes M^t diag(sigma2) M
// could probably be optimized
void matrix_multiply_cov(gsl_matrix * M, gsl_vector * sigma2, gsl_matrix * out) 
{
  gsl_vector * psi = gsl_vector_alloc(M->size2);
  //cout << "matrix mult: " << M->size1 << " " << M->size2 << "\n";
  //cout << sigma2->size << "\n";
  for(uint i = 0; i < M->size1; i++) {
    gsl_vector_memcpy(psi, sigma2);
    gsl_vector_const_view Mi = gsl_matrix_const_row(M,i);
    gsl_vector_mul(psi, &Mi.vector);
    gsl_blas_dgemv(CblasNoTrans, 1.0, M, psi, 0.0, &gsl_matrix_row(out, i).vector);
  }
}

// Returns a new matrix with vector v in the diagonal
gsl_matrix * diag(gsl_vector * v) 
{
  gsl_matrix * m = gsl_matrix_calloc(v->size, v->size);
  for(uint i = 0; i < v->size; i++) {
    gsl_matrix_set(m, i, i, gsl_vector_get(v,i));
  }
  return(m);
}

gsl_vector * invert_vector(gsl_vector * v)
{
  gsl_vector * vinv = gsl_vector_alloc(v->size);
  for(uint i = 0; i < vinv->size; i++) {
    gsl_vector_set(vinv, i, 1.0/(float)gsl_vector_get(v, i));
  }
  return(vinv);
}

 gsl_vector * pow_vector(gsl_vector * v, float p)
{
  gsl_vector * vinv = gsl_vector_alloc(v->size);
  for(uint i = 0; i < vinv->size; i++) {
    gsl_vector_set(vinv, i, pow(gsl_vector_get(v, i), p));
  }
  return(vinv);
}

void pow_matrix_inplace(gsl_matrix * v, float p)
{
  for(uint i = 0; i < v->size1; i++) {
    for(uint j = 0; j < v->size2; j++) {
      gsl_matrix_set(v, i, j, pow(gsl_matrix_get(v, i, j), p));
    }
  }
}

void mult_rows_matrix_inplace(gsl_matrix * m, gsl_vector * v)
{
  for(uint i = 0; i < m->size1; i++) {
    for(uint j = 0; j < m->size2; j++) {
      gsl_matrix_set(m, i, j, gsl_matrix_get(m, i, j)*gsl_vector_get(v,j));
    }
  }
}

void pow_vector_inplace(gsl_vector * v, float p)
{
  for(uint i = 0; i < v->size; i++) {
    gsl_vector_set(v, i, pow(gsl_vector_get(v, i), p));
  }
}

void invert_vector_inplace(gsl_vector * v)
{
  for(uint i = 0; i < v->size; i++) {
    gsl_vector_set(v, i, 1.0/gsl_vector_get(v, i));
  }
}

gsl_vector * matrix_row(gsl_matrix * m, int r)
{
  gsl_vector * row = gsl_vector_alloc(m->size2);
  gsl_matrix_get_row(row, m, r);
  return(row);
}

gsl_vector * matrix_column(gsl_matrix * m, int c)
{
  gsl_vector * col = gsl_vector_alloc(m->size1);
  gsl_matrix_get_col(col, m, c);
  return(col);
}

float vector_sum(gsl_vector * v)
{
  float total = 0;
  for(uint i = 0; i < v->size; i++) {
    total += gsl_vector_get(v,i);
  }
  return(total);
}

float mean(gsl_vector * v)
{
  float total = 0;
  for(uint i = 0; i < v->size; i++) {
    total += gsl_vector_get(v,i);
  }
  return(total/(double)v->size);
}

float sd(gsl_vector * v)
{
  return(sqrt(var(v)));
}

float var(gsl_vector * v)
{
  float total = 0;
  for(uint i = 0; i < v->size; i++) {
    total += gsl_vector_get(v,i);
  }
  float mean = total/(double)v->size;
  total = 0;
  for(uint i = 0; i < v->size; i++) {
    total += pow(gsl_vector_get(v,i)-mean, 2);
  }
  return(total/(double)v->size);
}

// This computes the inverse based on the Cholesky decomposition 
// in place; input matrix is overwritten.
/*
void invert_psd_matrix(gsl_matrix * m)
{
  //Cholesky decomp
  gsl_linalg_cholesky_decomp(m);
  gsl_linalg_cholesky_invert(m);
}
*/
// This computes the inverse based on the LU decomposition 
// in place; input matrix is overwritten.
void invert_square_matrix(gsl_matrix * m, gsl_vector * kby1, gsl_matrix * out)
{
  //LU decomp
  gsl_permutation * p = gsl_permutation_calloc(m->size1);
  int signum = 1;
  gsl_linalg_LU_decomp(m, p, &signum);
  gsl_vector_set_zero(kby1);
  // DON'T CHANGE THIS TO uint!!!!!!! Trust me!!!
  for(int i = 0; i < m->size1; i++) {
    gsl_vector_set(kby1, i, 1.0);
    if(i - 1 >= 0) {
      gsl_vector_set(kby1, i-1, 0.0);
    }
    gsl_linalg_LU_solve (m, p, kby1, &gsl_matrix_row(out, i).vector);
  }
  gsl_permutation_free(p);
}

// replaces mean term with sample
void sample_multivariate_normal(gsl_vector * mean, gsl_matrix * cov, gsl_rng * rng)
{
  //Cholesky decomp
  FILEIO::write_matrix("covar", cov);
  gsl_matrix * L = cholesky_decomposition(cov);

  // Generate a vector of random normals N(0,1)
  gsl_vector * mvn = gsl_vector_alloc(mean->size);
  for(uint i=0; i < mean->size; i++) {
    gsl_vector_set(mvn, i, gsl_ran_gaussian(rng, 1));
  }
  //FILEIO::print_vector(mvn);
  
  // scale each: Y = LX + mu
  gsl_blas_dgemv(CblasNoTrans, 1.0, L, mvn, 1.0, mean);

  // free new matrices/vectors
  gsl_vector_free(mvn);
  gsl_matrix_free(L);
}

gsl_vector * random_unif_vector(int count, double min, double max, gsl_rng * rng)
{
  gsl_vector * out = gsl_vector_alloc(count);
  for(int i = 0; i < count; i++) {
    gsl_vector_set(out, i, (gsl_rng_uniform(rng)*(max-min))+min);
  }
  return out;
}

gsl_matrix * cholesky_decomposition(gsl_matrix * m)
{
  gsl_matrix * L = gsl_matrix_alloc(m->size1, m->size2);
  gsl_matrix_memcpy(L,m);
  gsl_linalg_cholesky_decomp(L);
  for(uint i=1; i < m->size1; i++) {
    for(uint j = 0; j < i; j++) {
      gsl_matrix_set(L, i, j, 0);
    }
  }
  return(L);
}

gsl_vector * row_average(gsl_matrix * m)
{
  gsl_vector * s = gsl_vector_calloc(m->size1);
  for(uint i = 0; i < m->size2; i++) {
    for(uint j = 0; j < m->size1; j++) {
      gsl_vector_set(s, j, gsl_vector_get(s,j)+gsl_matrix_get(m,j,i));
    }
  }
  for(uint j = 0; j < m->size1; j++) {
    gsl_vector_set(s, j, gsl_vector_get(s,j)/m->size2);
  }
  return(s);
}

gsl_vector * column_average(gsl_matrix * m)
{
  gsl_vector * s = gsl_vector_calloc(m->size2);
  for(uint i = 0; i < m->size1; i++) {
    for(uint j = 0; j < m->size2; j++) {
      gsl_vector_set(s, j, gsl_vector_get(s,j)+gsl_matrix_get(m,i,j));
    }
  }
  for(uint j = 0; j < m->size2; j++) {
    gsl_vector_set(s, j, gsl_vector_get(s,j)/m->size1);
  }
  return(s);
}

 gsl_matrix * get_eigenvectors(gsl_matrix * m)
 {
   // first compute matrix covariance
   gsl_matrix * corm = matrix_correlation(m);
   cout << "Computed correlation...\n";
   gsl_matrix * Z = gsl_matrix_alloc(m->size1, m->size1);
   gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc(m->size1);

   gsl_vector_complex * eval = gsl_vector_complex_alloc(m->size1);
   gsl_matrix_complex * evec = gsl_matrix_complex_alloc(m->size1, m->size1);

   //http://www.gnu.org/software/gsl/manual/html_node/Real-Nonsymmetric-Matrices.html

   gsl_eigen_nonsymmv_Z(corm, eval, evec, Z, w);
   gsl_eigen_nonsymmv_free(w);
   gsl_matrix_free(corm);
   gsl_vector_complex_free(eval);
   gsl_matrix_complex_free(evec);
   cout << "Computed eigenvectors...\n";
   return(Z);
 }

 gsl_matrix * get_first_k_eigenvectors(gsl_matrix * m, int k) {
   gsl_matrix * out = gsl_matrix_calloc(m->size1,k);
   gsl_matrix * eigenv = get_eigenvectors(m);
   k = GSL_MIN(k, eigenv->size2);
   for(uint i=0; i < m->size1; i++) { // rows
     for(uint j=0; j < k; j++) { // columns
       gsl_matrix_set(out, i, j, gsl_matrix_get(eigenv, i, j));
     }
   }
   gsl_matrix_free(eigenv);
   return(out);
 }


gsl_matrix * matrix_correlation(gsl_matrix * m)
{
  gsl_matrix * out = gsl_matrix_alloc(m->size1, m->size1);
  float covr = 1.0;
  gsl_vector * row1 = gsl_vector_alloc(m->size2);
  gsl_vector * row2 = gsl_vector_alloc(m->size2);
  for(uint i = 0; i < m->size1; i++) {
    gsl_matrix_get_row(row1, m,i);
    for(uint j = i; j < m->size1; j++) {
      if(i == j) {
	gsl_matrix_set(out, i, i, 1.0);
      } else {
	gsl_matrix_get_row(row2, m, j);
	covr = gsl_stats_correlation(row1->data, row1->stride, row2->data, row2->stride, m->size2);
	//cout << "Setting " << i << '\t' << j << " to " << covr << '\n';
	gsl_matrix_set(out, i, j, covr);
	gsl_matrix_set(out, j, i, covr);
      }
    }
  }
  return(out);
}

gsl_vector * new_vector(float x, int size)
{
  gsl_vector * v = gsl_vector_alloc(size);
  for(int i = 0; i < size; i++) {
    gsl_vector_set(v,i,x);
  }
  return(v);
}

gsl_matrix * transpose(gsl_matrix * m)
{
  gsl_matrix * mt = gsl_matrix_alloc(m->size2, m->size1);
  gsl_matrix_transpose_memcpy(mt, m);
  //for(uint i = 0; i < m->size1; i++) {
  //gsl_matrix_get_row(&gsl_matrix_column(mt, i).vector, m, i);
  //}
  gsl_matrix_free(m);
  //FILEIO::print_matrix(mt);

  return(mt);
}

void add_to_matrix(gsl_matrix * m, float n)
{
  for(uint i = 0; i < m->size1; i++) {
    for(uint j = 0; j < m->size2; j++) {
      gsl_matrix_set(m, i, j, gsl_matrix_get(m, i, j) + n);
    }
  }
}

bool matrix_has_nan(gsl_matrix * m)
{
  for(uint i = 0; i < m->size1; i++) {
    for(uint j = 0; j < m->size2; j++) {
      float a = gsl_matrix_get(m,i,j);
      if(a != a) return true;
    }
  }
  return false;
}

gsl_matrix * get_K_columns_at_random(gsl_matrix * m, int K, gsl_rng * rng)
{
  gsl_matrix * out = gsl_matrix_alloc(m->size1, K);
  for(int k = 0; k < K; k++) {
    int rn = gsl_rng_uniform_int(rng, m->size2);
    gsl_matrix_get_col(&gsl_matrix_column(out,k).vector, m, rn);
    //gsl_matrix_get_col(&gsl_matrix_column(out,k).vector, m, (k*40+2));
  }
  return(out);
}

void compute_sd_rows(gsl_matrix * m, gsl_vector * v, gsl_vector * out)
{
  gsl_vector_set_zero(v);
  for(uint i = 0; i < m->size2; i++) {
    gsl_vector_add(v, &gsl_matrix_const_column(m, i).vector);
  }
  gsl_vector_scale(v, 1.0/float(m->size2));
  gsl_vector_set_zero(out);
  float sd = 0.0;
  //FILEIO::print_vector(v);
  for(uint i = 0; i < m->size2; i++) {
    for(uint j = 0; j < m->size1; j++) {
      sd = (gsl_vector_get(v, j) - gsl_matrix_get(m, j, i));
      sd = sd*sd;
      gsl_vector_set(out, j, gsl_vector_get(out, j) + sd);
    }
  }
  for(uint j = 0; j < m->size1; j++) {
    gsl_vector_set(out, j, sqrt(gsl_vector_get(out, j)/double(m->size2)));
  }
}

// v out are size2
void compute_sd_columns(gsl_matrix * m, gsl_vector * v, gsl_vector * out)
{
  gsl_vector_set_zero(v);
  for(uint i = 0; i < m->size1; i++) {
    gsl_vector_add(v, &gsl_matrix_const_row(m, i).vector);
  }
  gsl_vector_scale(v, 1.0/float(m->size1));
  gsl_vector_set_zero(out);
  float sd = 0.0;
  for(uint i = 0; i < m->size1; i++) {
    for(uint j = 0; j < m->size2; j++) {
      sd = (gsl_vector_get(v, j) - gsl_matrix_get(m, i, j));
      sd = sd*sd;
      gsl_vector_set(out, j, gsl_vector_get(out, j) + sd);
    }
  }
  for(uint j = 0; j < m->size2; j++) {
    gsl_vector_set(out, j, sqrt(gsl_vector_get(out, j)/double(m->size1)));
  }
}

gsl_matrix * get_K_rows_at_random(gsl_matrix * m, int K, gsl_rng * rng)
{
  gsl_matrix * out = gsl_matrix_alloc(K, m->size2);
  cout << "rows " << m->size1 << '\n';
  cout << "columns " << m->size2 << "\n";
  for(int k = 0; k < K; k++) {
    int rn = gsl_rng_uniform_int(rng, m->size1);
    gsl_matrix_set_row(out, k, &gsl_matrix_row(m, rn).vector);
    //gsl_matrix_get_row(&gsl_matrix_column(out,k).vector, m, (k*40+2));
  }
  return(out);
}

int get_max_index(gsl_vector * v)
{
  int mi = -1;
  double maxn = -100000;
  for(uint i = 0; i < v->size; i++) {
    double vi = gsl_vector_get(v,i);
    if(vi > maxn) {
      maxn = vi;
      mi = i;
    }
  }
  return(mi);
}


void outer_product_like_blas(float a, const gsl_vector * v1, const gsl_vector * v2, float b, gsl_matrix * out)
{
  float ov = 0.0;
  for(uint i = 0; i < v1->size; i++) {
    for(uint j = 0; j < v2->size; j++) {
      ov = gsl_matrix_get(out, i, j);
      gsl_matrix_set(out, i, j, a*gsl_vector_get(v1, i)*gsl_vector_get(v2,j) + b * ov);
    }
  }
}


//(A + U^tCU)^{-1} = A^-1 - A^-1 U^t (C^-1 + U A^-1 U^t)^-1 U A^-1
/**
void woodbury_matrix_inversion(float A, gsl_matrix * U, gsl_vector * C, gsl_matrix * out) 
{
  // pull out KxK matrix from out and compute inverted matrix
  int k = U->size1;
  int n = U->size2;
  gsl_matrix_view subm = gsl_matrix_submatrix(out, 0, 0, k, k);
  gsl_matrix_set_identity(&subm.matrix);
  for(int i = 0; i < k; i++) {
    gsl_matrix_set(&subm.matrix, i, i, 1.0/gsl_vector_get(C, i));
  }
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0/A, U, U, 1.0, &subm.matrix);
  // is the kxk matrix psd?
  invert_psd_matrix(&subm.matrix);
  gsl_matrix_view subm2 = gsl_matrix_submatrix(out, 0, 0, n, k);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, U, &subm.matrix, 0.0, &subm2.matrix);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, (1.0/(A*A)), &subm2.matrix, U, 0.0, out);
  gsl_matrix_scale(out, -1);
  for(int i = 0; i < n; i++) {
    gsl_matrix_set(out, i, i, gsl_matrix_get(out, i, i) + 1.0/A);
  }
}
**/

void woodbury_matrix_inversion_matrix_scaled_oned(float A, gsl_vector_const_view Uk, float Ci, gsl_matrix * F, gsl_matrix * kbyk, gsl_matrix * out)
{
  double cinvvainvu = 0.0;
  gsl_vector_view kby1 = gsl_matrix_row(kbyk, 0);

  gsl_blas_ddot(&Uk.vector, &Uk.vector, &cinvvainvu);
  cinvvainvu = 1.0/(1.0/Ci + (1.0/A) * cinvvainvu);

  gsl_blas_dgemv(CblasNoTrans, -(1.0/A)*cinvvainvu, F, &Uk.vector, 0.0, &kby1.vector);
  gsl_matrix_memcpy(out, F);
  outer_product_like_blas(1.0/A, &kby1.vector, &Uk.vector, 1.0/A, out); 
}

void woodbury_matrix_inversion_matrix_scaled_oned(gsl_vector * Ainv, gsl_vector_const_view Uk, float Ci, gsl_matrix * F, gsl_matrix * kbyk, gsl_matrix * out)
{
  double cinvvainvu = 0.0;
  gsl_vector_view kby1 = gsl_matrix_row(kbyk, 0);
  gsl_vector * outk = &gsl_matrix_row(out, 0).vector;

  gsl_vector_memcpy(outk, &Uk.vector);
  gsl_vector_mul(outk, Ainv);
  gsl_blas_ddot(outk, &Uk.vector, &cinvvainvu);
  cinvvainvu = 1.0/(1.0/Ci + cinvvainvu);

  gsl_blas_dgemv(CblasNoTrans, -cinvvainvu, F, outk, 0.0, &kby1.vector);
  gsl_matrix_memcpy(out, F);
  outer_product_like_blas(1.0, &kby1.vector, &Uk.vector, 1.0, out); 
  Ad(out, Ainv);
}

void woodbury_matrix_inversion_vector_scaled_oned(float A, gsl_vector_const_view Uk, float Ci, gsl_vector * Fk, gsl_vector * nby1, gsl_vector * out)
{
  double cinvvainvu = 0.0;
  double allbutlastv = 0.0;

  gsl_blas_ddot(&Uk.vector, &Uk.vector, &cinvvainvu);
  cinvvainvu = 1.0/(1.0/Ci + (1.0/A) * cinvvainvu);

  gsl_blas_ddot(Fk, &Uk.vector, &allbutlastv);
  allbutlastv = allbutlastv * (-1.0/A)*cinvvainvu;
  gsl_vector_memcpy(out, Fk);
  gsl_vector_memcpy(nby1, &Uk.vector);
  gsl_vector_scale(nby1, allbutlastv);
  gsl_vector_add(out, nby1);
  gsl_vector_scale(out, 1.0/A);
}

void woodbury_matrix_inversion_vector_scaled_oned(gsl_vector * Ainv, gsl_vector_const_view Uk, float Ci, gsl_vector * Fk, gsl_vector * nby1, gsl_vector * out)
{
  double cinvvainvu = 0.0;
  double allbutlastv = 0.0;

  gsl_vector_memcpy(out, &Uk.vector);
  gsl_vector_mul(out, Ainv);
  gsl_blas_ddot(out, &Uk.vector, &cinvvainvu);
  cinvvainvu = 1.0/(1.0/Ci + cinvvainvu);
 
  gsl_vector_memcpy(out, Fk);
  gsl_vector_mul(out, Ainv);
  gsl_blas_ddot(out, &Uk.vector, &allbutlastv);
  allbutlastv = allbutlastv * (-cinvvainvu);

  gsl_vector_memcpy(out, Fk);
  gsl_vector_memcpy(nby1, &Uk.vector);
  gsl_vector_scale(nby1, allbutlastv);
  gsl_vector_add(out, nby1);
  gsl_vector_mul(out, Ainv);
}


//F(A + U^tCU)^{-1} = FA^-1 - FA^-1 U^t (C^-1 + U A^-1 U^t)^-1 U A^-1
// out is KxN
void woodbury_matrix_inversion_matrix_scaled(float A, gsl_matrix * U, gsl_vector * C, gsl_matrix * F, gsl_matrix * kbyk, gsl_matrix * kbyk2, gsl_matrix * kbyn, gsl_matrix * out) 
{
  // First, find zeros in C and remove them (modify Ucopy appropriately)
  uint K = 0;
  int lastIndex = -1;
  gsl_vector * Ci = NULL;
  gsl_matrix * kbykp = NULL;
  gsl_vector * kby1p = NULL;
  gsl_matrix * Up = NULL;
  gsl_matrix * kbykp2 = NULL;
  gsl_matrix * kbykp12 = NULL;
  gsl_matrix * outp = NULL;
  for(uint i = 0; i < C->size; i++) {
    if(gsl_vector_get(C,i) != 0.0) {
      K += 1;
      lastIndex = i;
    }
  }
  //FILEIO::print_vector(C);
  if(K == 0) {
    gsl_matrix_memcpy(out, F);
    gsl_matrix_scale(out, 1.0/A);
    //FILEIO::print_matrix(out);
    return;
  } else if (K == 1) {
    woodbury_matrix_inversion_matrix_scaled_oned(A, gsl_matrix_const_row(U, lastIndex), gsl_vector_get(C, lastIndex), F, kbyk, out); 
    return;
  } else if(K < U->size1) {
    Ci = gsl_vector_alloc(K);
    float ci = 0.0;
    int length = 0;
    for(uint i = 0; i < C->size; i++) {
      ci = gsl_vector_get(C,i);
      if(ci != 0.0) {
	gsl_vector_set(Ci, length, ci);
	gsl_matrix_set_row(kbyn, length, &gsl_matrix_const_row(U, i).vector);
	length += 1;
      }
    }
    Up = &gsl_matrix_const_submatrix(kbyn, 0, 0, K, kbyn->size2).matrix;
    kbykp = &gsl_matrix_const_submatrix(kbyk, 0, 0, K, K).matrix;
    kbykp2 = &gsl_matrix_const_submatrix(kbyk2, 0, 0, K, K).matrix;
    kbykp12 = &gsl_matrix_const_submatrix(kbyk, 0, 0, kbyn->size1, K).matrix;
    outp = &gsl_matrix_const_submatrix(out, 0, 0, K, kbyn->size2).matrix;
    kby1p = &gsl_matrix_column(outp, 0).vector;
  } else { // complete matrix
    Ci = C;
    Up = U;
    kbykp = kbyk;
    kbykp2 = kbyk2;
    kbykp12 = kbyk;
    kby1p = &gsl_matrix_column(out, 0).vector;
    outp = out;
  }
  //cout << "finished pruning\n";
  //FILEIO::print_vector(Ci);
  //FILEIO::print_matrix(Up);

  // pull out KxK matrix from out and compute inverted matrix
  int k = Up->size1;
  
  gsl_matrix_set_identity(kbykp);
  for(uint i = 0; i < k; i++) {
    gsl_matrix_set(kbykp, i, i, 1.0/gsl_vector_get(Ci, i));
  }
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0/A, Up, Up, 1.0, kbykp);
  // is the kxk matrix psd?
  //invert_psd_matrix(kbykp);
  invert_square_matrix(kbykp, kby1p, kbykp2);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, kbykp2, Up, 0.0, outp);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0/A, F, outp, 0.0, kbykp12);
  gsl_matrix_memcpy(out, F);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0/A, kbykp12, Up, 1.0/A, out); 

  if(K < U->size1) {
    gsl_vector_free(Ci);
  }
}


//F(A + U^tCU)^{-1} = FA^-1 - FA^-1 U^t (C^-1 + U A^-1 U^t)^-1 U A^-1
// out is KxN
void woodbury_matrix_inversion_matrix_scaled(gsl_vector * Ainv, gsl_matrix * U, gsl_vector * C, gsl_matrix * F, gsl_matrix * kbyk, gsl_matrix * kbyk2, gsl_matrix * kbyn, gsl_matrix * out) 
{
  // First, find zeros in C and remove them (modify Ucopy appropriately)
  uint K = 0;
  int lastIndex = -1;
  gsl_vector * Ci = NULL;
  gsl_matrix * kbykp = NULL;
  gsl_vector * kby1p = NULL;
  gsl_matrix * Up = NULL;
  gsl_matrix * kbykp2 = NULL;
  gsl_matrix * kbykp12 = NULL;
  gsl_matrix * outp = NULL;
  gsl_vector * Nby1 = NULL;
  for(uint i = 0; i < C->size; i++) {
    if(gsl_vector_get(C,i) != 0.0) {
      K += 1;
      lastIndex = i;
    }
  }
  //cout << "K = " << K << '\n';
  //FILEIO::print_vector(C);
  if(K == 0) {
    gsl_matrix_memcpy(out, F);
    Ad(out, Ainv);
    //FILEIO::print_matrix(out);
    return;
  } else if (K == 1) {
    woodbury_matrix_inversion_matrix_scaled_oned(Ainv, gsl_matrix_const_row(U, lastIndex), gsl_vector_get(C, lastIndex), F, kbyk, out); 
    return;
  } else if(K < U->size1) {
    Ci = gsl_vector_alloc(K);
    float ci = 0.0;
    int length = 0;
    for(uint i = 0; i < C->size; i++) {
      ci = gsl_vector_get(C,i);
      if(ci != 0.0) {
	gsl_vector_set(Ci, length, ci);
	gsl_matrix_set_row(kbyn, length, &gsl_matrix_const_row(U, i).vector);
	length += 1;
      }
    }
    Up = &gsl_matrix_const_submatrix(kbyn, 0, 0, K, kbyn->size2).matrix;
    kbykp = &gsl_matrix_const_submatrix(kbyk, 0, 0, K, K).matrix;
    kbykp2 = &gsl_matrix_const_submatrix(kbyk2, 0, 0, K, K).matrix;
    kbykp12 = &gsl_matrix_const_submatrix(kbyk, 0, 0, kbyn->size1, K).matrix;
    outp = &gsl_matrix_const_submatrix(out, 0, 0, K, kbyn->size2).matrix;
    kby1p = &gsl_matrix_column(outp, 0).vector;
    Nby1 = &gsl_matrix_row(kbyn, K).vector;
  } else { // complete matrix
    Ci = C;
    Up = U;
    Nby1 = &gsl_matrix_row(kbyn, 0).vector;
    kbykp = kbyk;
    kbykp2 = kbyk2;
    kbykp12 = kbyk;
    kby1p = &gsl_matrix_column(out, 0).vector;
    outp = out;
  }
  //cout << "finished pruning\n";
  //FILEIO::print_vector(Ci);
  //FILEIO::print_matrix(Up);

  // pull out KxK matrix from out and compute inverted matrix
  int k = Up->size1;
  
  gsl_matrix_set_identity(kbykp);
  for(int i = 0; i < k; i++) {
    gsl_matrix_set(kbykp, i, i, 1.0/gsl_vector_get(Ci, i));
  }
  //gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0/A, Up, Up, 1.0, kbykp);
  ADAt(Up, Up, Ainv, Nby1, kbykp2);
  gsl_matrix_add(kbykp, kbykp2);
  // is the kxk matrix psd?
  //invert_psd_matrix(kbykp);
  invert_square_matrix(kbykp, kby1p, kbykp2);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, kbykp2, Up, 0.0, outp);
  //gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0/A, F, outp, 0.0, kbykp12);
  //cout << "outp " << kbykp12->size1 << " " << kbykp12->size2 << '\n';
  ADAt(F, outp, Ainv, Nby1, kbykp12);
  gsl_matrix_memcpy(out, F);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, kbykp12, Up, 1.0, out); 
  Ad(out, Ainv);

  if(K < U->size1) {
    gsl_vector_free(Ci);
  }
}


//Fk(A + U^tCU)^{-1} = FkA^-1 - FkA^-1 U^t (C^-1 + U A^-1 U^t)^-1 U A^-1
void woodbury_matrix_inversion_vector_scaled(float A, gsl_matrix * U, gsl_vector * C, gsl_vector * Fk, gsl_matrix * kbyk, gsl_matrix * kbyk2, gsl_matrix * kbyn, gsl_vector * out) 
{
  // First, find zeros in C and remove them (modify Ucopy appropriately)
  uint K = 0;
  int lastIndex = -1;
  gsl_vector * Ci = NULL;
  gsl_matrix * kbykp = NULL;
  gsl_matrix * kbykp2 = NULL;
  gsl_matrix * Up = NULL;
  gsl_vector * kby1p = NULL;
  gsl_vector * kby1p2 = NULL;
  gsl_matrix_set_zero(kbyn);
  for(uint i = 0; i < C->size; i++) {
    if(gsl_vector_get(C,i) != 0.0) {
      K += 1;
      lastIndex = i;
    }
  }
  if(K == 0) {
    gsl_vector_memcpy(out, Fk);
    gsl_vector_scale(out, 1.0/A);
    return;
  } else if (K == 1) {
    woodbury_matrix_inversion_vector_scaled_oned(A, gsl_matrix_const_row(U, lastIndex), gsl_vector_get(C, lastIndex), Fk, &gsl_matrix_row(kbyn, 0).vector, out); 
    return;
  } else if(K < U->size1) {
    Ci = gsl_vector_alloc(K);
    float ci = 0.0;
    int length = 0;
    for(uint i = 0; i < C->size; i++) {
      ci = gsl_vector_get(C,i);
      if(ci != 0.0) {
	gsl_vector_set(Ci, length, ci);
	gsl_matrix_set_row(kbyn, length, &gsl_matrix_const_row(U, i).vector);
	length += 1;
      }
    }
    Up = &gsl_matrix_const_submatrix(kbyn, 0, 0, K, kbyn->size2).matrix;
    kbykp = &gsl_matrix_submatrix(kbyk, 0, 0, K, K).matrix;
    kbykp2 = &gsl_matrix_submatrix(kbyk2, 0, 0, K, K).matrix;
    kby1p = gsl_vector_alloc(K);
    kby1p2 = gsl_vector_alloc(K);
  } else { // complete matrix
    Ci = C;
    Up = U;
    kbykp = kbyk;
    kbykp2 = kbyk2;
    kby1p = gsl_vector_alloc(kbyk->size1);
    kby1p2 = gsl_vector_alloc(kbyk->size1);
  }
  //cout << "finished pruning\n";
  //FILEIO::print_vector(Ci);
  //FILEIO::print_matrix(Up);

  // pull out KxK matrix from out and compute inverted matrix
  int k = Up->size1;
  
  gsl_matrix_set_identity(kbykp);
  for(int i = 0; i < k; i++) {
    gsl_matrix_set(kbykp, i, i, 1.0/gsl_vector_get(Ci, i));
  }
  //cout << "here1\n";
  //FILEIO::print_matrix(kbykp);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0/A, Up, Up, 1.0, kbykp);
  // is the kxk matrix psd?
  //cout << "here2\n";
  //FILEIO::print_matrix(kbykp);
  invert_square_matrix(kbykp, kby1p, kbykp2);
  //cout << "here3\n";
  //FILEIO::print_vector(Fk);
  //FILEIO::print_matrix(Up);
  //FILEIO::print_matrix(kbykp2);
  gsl_blas_dgemv(CblasNoTrans, -1.0/A, Up, Fk, 0.0, kby1p);
  
  //cout << "here4\n";
  //cout << "kby1p " << kby1p->size << "\n";
  //FILEIO::print_vector(kby1p);
  gsl_blas_dgemv(CblasNoTrans, 1.0, kbykp2, kby1p, 0.0, kby1p2);
  //cout << "kby1p2 " << kby1p2->size << "\n";
  //FILEIO::print_vector(kby1p2);
  //cout << "here5\n";
  gsl_vector_memcpy(out, Fk);
  //cout << "here6\n";
  gsl_blas_dgemv(CblasTrans, 1.0/A, Up, kby1p2, 1.0/A, out); 

  if(K < U->size1) {
    gsl_vector_free(Ci);
  }
  gsl_vector_free(kby1p);
  gsl_vector_free(kby1p2);
}


//F(A + U^tCU)^{-1} = FA^-1 - FA^-1 U^t (C^-1 + U A^-1 U^t)^-1 U A^-1
void woodbury_matrix_inversion_vector_scaled(gsl_vector * Ainv, gsl_matrix * U, gsl_vector * C, gsl_vector * Fk, gsl_matrix * kbyk, gsl_matrix * kbyk2, gsl_matrix * kbyn, gsl_vector * out) 
{
  // First, find zeros in C and remove them (modify Ucopy appropriately)
  uint K = 0;
  int lastIndex = -1;
  gsl_vector * Ci = NULL;
  gsl_matrix * kbykp = NULL;
  gsl_matrix * kbykp2 = NULL;
  gsl_matrix * Up = NULL;
  gsl_vector * kby1p = NULL;
  gsl_vector * kby1p2 = NULL;
  gsl_vector * Nby1 = NULL;

  gsl_matrix_set_zero(kbyn);
  for(uint i = 0; i < C->size; i++) {
    if(gsl_vector_get(C,i) != 0.0) {
      K += 1;
      lastIndex = i;
    }
  }
  if(K == 0) {
    gsl_vector_memcpy(out, Fk);
    gsl_vector_mul(out, Ainv);
    return;
  } else if (K == 1) {
    woodbury_matrix_inversion_vector_scaled_oned(Ainv, gsl_matrix_const_row(U, lastIndex), gsl_vector_get(C, lastIndex), Fk, &gsl_matrix_row(kbyn, 0).vector, out); 
    return;
  } else if(K < U->size1) {
    Ci = gsl_vector_alloc(K);
    float ci = 0.0;
    int length = 0;
    for(uint i = 0; i < C->size; i++) {
      ci = gsl_vector_get(C,i);
      if(ci != 0.0) {
	gsl_vector_set(Ci, length, ci);
	gsl_matrix_set_row(kbyn, length, &gsl_matrix_const_row(U, i).vector);
	length += 1;
      }
    }
    Up = &gsl_matrix_const_submatrix(kbyn, 0, 0, K, kbyn->size2).matrix;
    kbykp = &gsl_matrix_submatrix(kbyk, 0, 0, K, K).matrix;
    kbykp2 = &gsl_matrix_submatrix(kbyk2, 0, 0, K, K).matrix;
    kby1p = gsl_vector_alloc(K);
    kby1p2 = gsl_vector_alloc(K);
    Nby1 = &gsl_matrix_row(kbyn, K).vector;
  } else { // complete matrix
    Ci = C;
    Up = U;
    kbykp = kbyk;
    kbykp2 = kbyk2;
    kby1p = gsl_vector_alloc(kbyk->size1);
    kby1p2 = gsl_vector_alloc(kbyk->size1);
    Nby1 = &gsl_matrix_row(kbyn, 0).vector;
  }
  //cout << "finished pruning\n";
  //FILEIO::print_vector(Ci);
  //FILEIO::print_matrix(Up);

  // pull out KxK matrix from out and compute inverted matrix
  int k = Up->size1;
  
  gsl_matrix_set_identity(kbykp);
  for(int i = 0; i < k; i++) {
    gsl_matrix_set(kbykp, i, i, 1.0/gsl_vector_get(Ci, i));
  }
  //cout << "here1\n";
  //FILEIO::print_matrix(kbykp);
  ADAt(Up, Up, Ainv, out, kbykp);
  //gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0/A, Up, Up, 1.0, kbykp);
  // is the kxk matrix psd?
  //cout << "here2\n";
  //FILEIO::print_matrix(kbykp);
  invert_square_matrix(kbykp, kby1p, kbykp2);
  //cout << "here3\n";
  //FILEIO::print_vector(Fk);
  //FILEIO::print_matrix(Up);
  //FILEIO::print_matrix(kbykp2);
  
  gsl_vector_memcpy(out, Fk);
  gsl_vector_mul(out, Ainv);
  gsl_blas_dgemv(CblasNoTrans, -1.0, Up, out, 0.0, kby1p);
  
  //cout << "here4\n";
  //cout << "kby1p " << kby1p->size << "\n";
  //FILEIO::print_vector(kby1p);
  gsl_blas_dgemv(CblasNoTrans, 1.0, kbykp2, kby1p, 0.0, kby1p2);
  //cout << "kby1p2 " << kby1p2->size << "\n";
  //FILEIO::print_vector(kby1p2);
  //cout << "here5\n";
  // bee: here
  dAv(Ainv, Up, kby1p2, Nby1);
  //cout << "here6\n";
  gsl_blas_dgemv(CblasTrans, 1.0, Up, kby1p2, 1.0, out); 
  //cout << "here7\n";

  if(K < U->size1) {
    gsl_vector_free(Ci);
  }
  gsl_vector_free(kby1p);
  gsl_vector_free(kby1p2);
}

// computes A %*% diag(d) %*% A2^t
void ADAt(gsl_matrix * A, gsl_matrix * A2, gsl_vector * d, gsl_vector * Nby1, gsl_matrix * kbykp) {
  int M = A->size1;
  int M2 = A2->size1;
  // kbykp is supposed to be a k x kp matrix
  double scal = 0;
  for(int m = 0; m < M; m++) {
    gsl_matrix_get_row(Nby1, A, m);
    gsl_vector_mul(Nby1, d);
    for(int m2 = 0; m2 < M2; m2++) {
      gsl_blas_ddot(Nby1, &gsl_matrix_row(A2,m2).vector, &scal);
      gsl_matrix_set(kbykp, m, m2, scal);
    }
  }
}


// computes diag(d) %*% A^t %*% v
void dAv(gsl_vector * d, gsl_matrix * A, gsl_vector * v, gsl_vector * Nby1) {
  gsl_blas_dgemv(CblasTrans, 1.0, A, v, 0.0, Nby1);
  gsl_vector_mul(Nby1, d);
}

// computes A %*% diag(d), returns in A
void Ad(gsl_matrix * A, gsl_vector * d) {
  int K = A->size1;
  for(int k = 0; k < K; k++) {
    gsl_vector_mul(&gsl_matrix_row(A, k).vector, d);
  }
}


//void lu_decomp(real ** a, int n, int *indx, int *d)
//{
//	int i, j, k;
//	int imax = -1;
//	real big, dum, sum, temp;
//	real *vv = new real[n];
//	*d = 1;
//
//	for (i = 0; i < n; i++)
//	{
//		big = 0.0; 
//		for (j = 0; j < n; j++)
//			if((temp = fabs(a[i][j])) > big) big = temp;
//		if (big == 0.0)
//		{
//			cout << "singular matrix in routine ludcmp" << endl;
//			exit(0); 
//		}
//		vv[i] = 1.0 / big; 
//	}
//
//	for (j = 0; j < n; j++)
//	{
//		for (i = 0; i < j; i++)
//		{
//			sum = a[i][j];
//			for (k = 0; k < i; k++) sum -= a[i][k] * a[k][j];
//			a[i][j] = sum;
//		}
//
//		big = 0.0;
//		for (i = j; i < n; i++)
//		{
//			sum = a[i][j];
//			for (k = 0; k < j; k++) sum -= a[i][k] * a[k][j];
//			a[i][j] = sum;
//			if ((dum = vv[i] * fabs(sum)) >= big)
//			{
//				big = dum; 
//				imax = i; 
//			}
//		}
//
//		if ( j != imax)
//		{
//			for (k = 0; k < n; k++)
//			{
//				dum = a[imax][k];
//				a[imax][k] = a[j][k];
//				a[j][k] = dum; 
//			}
//			*d = -(*d); 
//			vv[imax] = vv[j];
//		}
//		indx[j] = imax;
//		if(a[j][j] == 0.0) a[j][j] = TINY;
//		if (j != n)
//		{
//			dum = 1.0 / a[j][j];
//			for (i=j+1; i < n; i++) a[i][j] *= dum; 
//		}
//	}
//	delete[] vv; 
//}
//
//void lu_back_sub(real ** a, int n, int *indx, real b[])
//{
//	int i, ii = -1, ip, j; 
//	real sum;
//
//	for (i = 0; i < n; i++)
//	{
//		ip = indx[i];
//		sum = b[ip];
//		b[ip] = b[i];
//		if(ii>=0)
//			for (j=ii;j <i;j++) sum -= a[i][j] * b[j];
//		else if(sum) ii = i; 
//		b[i] = sum; 
//	}
//
//	for (i = n-1; i >= 0; i--) 
//	{
//		sum = b[i];
//		for (j=i+1; j < n; j++) sum -= a[i][j] * b[j];
//		b[i] = sum/a[i][i];
//	}
//} 

void outer_product(gsl_vector * a, gsl_vector * b, gsl_matrix * out)
{
  // check to see whether len(a) == len(b) == both dimensions of out
  int n = a->size;
  int p = b->size;
  if(n != out->size1 || p != out->size2) {
    cout << "Error: incorrect dimensions for outer product\n";
    exit(0);
  }
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < p; j++) {
      gsl_matrix_set(out, i, j, gsl_vector_get(a, i)*gsl_vector_get(b,j));
    }
  }
}

// Computes the log determinant of the matrix m. Let outm be the same size as the (square) matrix m
double matrix_determinant(gsl_matrix * m, gsl_matrix * outm)
{
  gsl_matrix_memcpy(outm, m);
  gsl_permutation * p = gsl_permutation_calloc(m->size1);
  int signum = 1;
  gsl_linalg_LU_decomp(outm, p, &signum);
  double det = gsl_linalg_LU_det (outm, signum);
  gsl_permutation_free(p);
  return(det);
}
