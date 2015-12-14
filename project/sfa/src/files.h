#ifndef __FILES__
#define __FILES__

#include <string>

using namespace std; 

class FILEIO 
{

 public:

  static gsl_matrix * read_matrix(string filename, int rows, int cols);
  static void write_matrix(string outfile, gsl_matrix * m);
  static void write_matrix_scaled(string outfile, gsl_matrix * m, double scale);
  static void print_matrix(gsl_matrix * m);
  static void print_vector(gsl_vector * v);
  static gsl_vector * read_vector(string filename, int size);
  static void write_vector(string outfile, gsl_vector * v);
  static void write_vector_scaled(string outfile, gsl_vector * v, double scale);
  static void write_control_file(string outfile, int K, int iterations, float sigma2Lambda);
  static void read_control_file(string infile, int * K, int * iterations, float * sigma2Lambda);
};

#endif // __FILES__
