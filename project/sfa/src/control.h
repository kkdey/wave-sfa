// define ctrl class  
#ifndef __CONTROL_H__
#define __CONTROL_H__

#include <stdlib.h>
#include <ctype.h>

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "fpmath.h"
#include "model.h"

using namespace std;

typedef struct _gene_ {
	string name; 
	int chr; 
	long start; 
	long end;
} GENE; 

class CtrlParam
{
private:
	map<string, int> hcom; 
	class ModelnData * pMD; 
	int m_nImputeRepeat; 
	
	int has_read; 
	int has_em; 
	int need_em; 
	int batch_em; 
	gsl_rng * gsl_r;

public: 
	CtrlParam(void); //default construct;
	~CtrlParam(void); //destruct. 

	void read_gene_map(string, vector<GENE>&); 
	void BatchRun(int, char **);
	void Run(void);
	void CommandLine(int, char **);
	void OnHelp(void);
	void PrintHeader(void);
	void OnRead(string); 
	void OnSet(string); 
	void OnEM(string); 
	void OnBF(string); 
};

#endif
