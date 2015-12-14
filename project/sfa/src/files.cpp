#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <stdio.h>    
#include <stdlib.h>
#include <sys/types.h>
#include "control.h"
#include "model.h"
#include "files.h"
//#include "indiv.h"
//#include "diploid.h"
//#include "haploid.h"
//#if defined (READLINE)
//#include "readline/readline.h"
//#include "readline/history.h"
//#include "ncurses.h"
//#endif
//#if defined (MPI_ENABLED)
//#include "mpi.h"
//#endif 
//#include "fpmath.h"
//#include "fp.h" 

using namespace std;

gsl_matrix * FILEIO::read_matrix(string filename, int rows, int cols)
{
  fstream infile;
  infile.open(filename.data(), ios::in);
  gsl_matrix * m = gsl_matrix_alloc(rows, cols);
  cout << "Opening matrix " << filename << " with " << rows << " rows and " << cols << " columns\n";
  cout << "Opening matrix " << filename << " with " << m->size1 << " rows and " << m->size2 << " columns\n";
  string delimit(",; :\t");
  string line;
  streambuf * pbuf;
  pbuf = infile.rdbuf();

  for (int i = 0; i < rows; i++) {
    line.assign(getline(pbuf)); 
    int beg = line.find_first_not_of(delimit, 0); 
    int end = -1; 
    for (int np = 0; np < cols; np++) {
      end = line.find_first_of(delimit, beg);
      if(end == beg) break; 
      string sv(line, beg, end-beg);
      //if(sv.compare("NP") == 0)
      //indphval.push_back(NP);
      //else if(sv.compare("NA") == 0)
      //indphval.push_back(NA);
      //else
      gsl_matrix_set(m,i,np,atof(sv.data()));
      beg = line.find_first_not_of(delimit, end); 
    } 
  }
  //  if(gsl_matrix_fscanf(f, m) != 0) {
  //printf("Warning: encountered problem reading in matrix from file %s with %d rows and %d cols", 
  //   filename.data(), rows, cols);
  //}
  infile.close();
  //  fclose(f);

  return(m);
}

void FILEIO::write_matrix(string filename, gsl_matrix * m)
{
  fstream outfile;
  outfile.open(filename.data(), ios::out);
  for(uint i = 0; i < m->size1; i++) {
    for(uint j = 0; j < m->size2-1; j++) {
      outfile << gsl_matrix_get(m,i,j) << '\t';
    }
    outfile << gsl_matrix_get(m, i, m->size2-1) << "\n";
  }
  outfile.close();
}

void FILEIO::write_matrix_scaled(string filename, gsl_matrix * m, double scale)
{
  fstream outfile;
  outfile.open(filename.data(), ios::out);
  for(uint i = 0; i < m->size1; i++) {
    for(uint j = 0; j < m->size2-1; j++) {
      outfile << gsl_matrix_get(m,i,j)/scale << '\t';
    }
    outfile << gsl_matrix_get(m, i, m->size2-1)/scale << "\n";
  }
  outfile.close();
}


void FILEIO::print_matrix(gsl_matrix * m)
{
  for(uint i = 0; i < m->size1; i++) {
    for(uint j = 0; j < m->size2-1; j++) {
      cout << gsl_matrix_get(m,i,j) << '\t';
    }
    cout << gsl_matrix_get(m, i, m->size2-1) << "\n";
  }
}

void FILEIO::print_vector(gsl_vector * v)
{
  for(uint i = 0; i < (v->size-1); i++) {
    cout << gsl_vector_get(v,i) << '\t';
  }
  cout << gsl_vector_get(v, v->size-1) << "\n";
}

gsl_vector * FILEIO::read_vector(string filename, int size)
{
  gsl_vector * v = gsl_vector_alloc(size);
  fstream infile;
  infile.open(filename.data(), ios::in);
  cout << "Opening vector " << filename << " with length " << v->size << "\n";
  string delimit(",; :\t");
  string line;
  streambuf * pbuf;
  pbuf = infile.rdbuf();

  line.assign(getline(pbuf)); 
  int beg = line.find_first_not_of(delimit, 0); 
  int end = -1; 
  for (int i = 0; i < size; i++) {
    end = line.find_first_of(delimit, beg);
    if(end == beg) {
      //line.assign(getline(pbuf)); 
      //beg = line.find_first_not_of(delimit, 0); 
      //end = -1; 
      break; 
    }
    string sv(line, beg, end-beg);
    gsl_vector_set(v,i,atof(sv.data()));
    //cout << " found " << sv.data() << '\n';
    beg = line.find_first_not_of(delimit, end);  
  }
  infile.close();
  return(v);
}

void FILEIO::write_vector(string filename, gsl_vector * v)
{
  fstream outfile;
  outfile.open(filename.data(), ios::out);
  for(uint i = 0; i < v->size-1; i++) {
    outfile << gsl_vector_get(v,i) << '\t';
  }
  outfile << gsl_vector_get(v,v->size-1) << "\n";
  outfile.close();
}

void FILEIO::write_vector_scaled(string filename, gsl_vector * v, double scale)
{
  fstream outfile;
  outfile.open(filename.data(), ios::out);
  for(uint i = 0; i < v->size-1; i++) {
    outfile << gsl_vector_get(v,i)/scale << '\t';
  }
  outfile << gsl_vector_get(v,v->size-1)/scale << "\n";
  outfile.close();
}

void FILEIO::write_control_file(string filename, int K, int iterations, float sigma2Lambda)
{
  fstream outfile;
  outfile.open(filename.data(), ios::out);
  outfile << K << '\n';
  outfile << iterations << '\n';
  outfile << sigma2Lambda << '\n';
  outfile.close();
}

void FILEIO::read_control_file(string infile, int * K, int * iterations, float * sigma2Lambda)
{
  FILE * f = fopen(infile.data(), "rb");
  fscanf(f, "%d\n", K);
  fscanf(f, "%d\n", iterations);
  fscanf(f, "%f\n", sigma2Lambda);
  fclose(f);
}
