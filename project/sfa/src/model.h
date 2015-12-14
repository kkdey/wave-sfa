#ifndef __MODEL_H__                
#define __MODEL_H__

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>

#include "fpmath.h"

using namespace std;

class ModelnData
{
public: 
	stringstream fplog; // log;
	inline short get_G() {return(G);}
	inline double get_N() {return(N);}
	inline short get_K() {return(K);}
	inline gsl_matrix * get_Y() {return(Y);}
private:
	double N;
	short K;
	short G;
	gsl_matrix * Y;
	gsl_rng * rng;
	
	gsl_vector * muN;
	gsl_vector * muG;
	gsl_vector * tauinvN;
	gsl_vector * tauinvG;
	gsl_vector * pk;
	gsl_matrix * Lambda;
	gsl_matrix * F;
	gsl_vector * lambdak;
	gsl_vector * sigma2K;
	float sigma2lambda;

	gsl_rng * gsl_r;

	int randSeed;               //user may specify random seed;
	int iters;
	int transposep;
	fstream logfile; 
	vector<string> vGin;
	string fnOutput;  			//output prefix; 
	
	int rowMeans;
	int columnMeans;
	int rowCovariance;
	int columnCovariance;

public:
	ModelnData();
	~ModelnData();
	
	inline void SetHasRowMeans(void) { rowMeans = 1; }
	inline void SetHasColumnMeans(void) { columnMeans = 1; }
	inline int HasRowMeans(void) { return(rowMeans); }
	inline int HasColumnMeans(void) { return(columnMeans); }

	inline void SetHasRowCovariance(void) { rowCovariance = 1; columnCovariance = 0;}
	inline void SetHasColumnCovariance(void) { columnCovariance = 1; rowCovariance = 0;}
	inline int HasRowCovariance(void) { return(rowCovariance); }
	inline int HasColumnCovariance(void) { return(columnCovariance); }

	inline void SetvGin(string s) {vGin.push_back(s);}
	inline int Getnum_Gin(void) {return (int) vGin.size();}

	inline void SetfnOutput(string s) {fnOutput.assign(s);}
	inline string GetfnOutput(void) {return fnOutput;}

	inline void SetTrans(int s) { transposep = s; }
	inline int GetTrans(void) { return transposep; }

	inline void SetrandSeed(long int s){randSeed = s;}

	inline void SetnPH(int s){cout << " Setting G to be "<< s << "\n"; G = s;}	
	inline void SetnFa(int s){K = s; }	
	inline void SetnIN(int s){cout << " Setting N to be " << s << "\n"; N = s;}	
	inline void SetIters(int s){iters = s;}
	inline long int GetIters(void){return iters;}

	inline long int GetrandSeed(void){return randSeed;}

	//interface
	real sfa_model(void);    // read in Y matrix and run SFA

};    

#endif

