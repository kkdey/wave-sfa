#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <stdio.h>    
#include <stdlib.h>
#include "control.h"
#include "model.h"
#include "fpmath.h"
#include "fp.h" 

using namespace std;

#define com_empty 	0
#define com_geno	1
#define com_pheno	2
#define com_out		3
#define com_pos		4
#define com_struct  5
#define com_nPH		6
#define com_clean   7
#define com_rand 	8
#define com_trans         9
#define com_em 		11
#define com_step 	12
#define com_warm 	13
#define com_cluster     14
#define com_nIN         15
#define com_nFa         16
#define com_iters         17
#define com_thin         18
#define com_burn         19
#define com_impute  20
#define com_multi	21
#define com_level 	23
#define com_pval 	24
#define com_msnp 	25
#define com_meansg      26
#define com_meansn      27
#define com_varg      28
#define com_varn      29
#define com_rem     30
#define com_sem     31
#define com_sigma_a 41
#define com_sigma_d 42
#define com_sort    43

#define com_help    100
#define com_wmg		101
#define com_wbg     102   //write best guess; 
#define com_wgd     103   //write genotype distribution; 
#define com_unix	200
#define com_EXIT	201
#define com_HELP	202
#define com_READ 	203
#define com_EM		204
#define com_BF		205
#define com_SET		206
#define com_PRINT	207
#define com_casectrl       1001
#define com_cluster_method 1002
#define com_mask_impute    1003
#define com_gene    2000
#define com_flank   2001
#define com_ssd 	3000
#define com_psd     3001
#define com_file    3002
#define com_num     3003
#define com_df      3004
#define com_mcmc    3210


// construct function contains the defaults 
CtrlParam::CtrlParam(void)
{
	m_nImputeRepeat = 0; 
	has_read = 0; 
	has_em = 0; 
	need_em = 0; 
	pMD = NULL;
	gsl_r = gsl_rng_alloc(gsl_rng_taus);

	hcom[""] = com_empty;
	hcom["-gen"] = com_geno;
	hcom["-o"] = com_out;
	hcom["-out"] = com_out;
	hcom["-g"] = com_nPH; 
	hcom["-G"] = com_nPH;
	hcom["-k"] = com_nFa; 
	hcom["-K"] = com_nFa;
	hcom["-n"] = com_nIN; 
	hcom["-N"] = com_nIN;
	hcom["-r"] = com_rand;
	hcom["-R"] = com_rand;
	hcom["-rand"] = com_rand;
	hcom["-t"] = com_trans;

	hcom["-iter"] = com_iters;

	hcom["-mg"] = com_meansg; // include G (row) means
	hcom["-mn"] = com_meansn; // include N (column) means
	hcom["-vg"] = com_varg; // include G (row) covariance terms
	hcom["-vn"] = com_varn; // include N (column) covariance terms
//help;     
	hcom["-h"] = com_help;
	hcom["-help"] = com_help; 
} 

CtrlParam::~CtrlParam(void)
{
	;
} 

void CtrlParam::OnHelp(void)
{
	cout << endl; 
	cout << " SFA Version 0.1 released in Dec 2009 by B. Engelhardt and M. Stephens" << endl; 
	cout << " Usage example: ./sfa -gen ../input/genotypedata.sfa -g 123 -k 3 -n 12345 -iter 50 -rand 482 -o out" << endl;
	cout << " All options are case-sensitive." << endl;
	cout << " " << endl;
	cout << " FILE I/O RELATED OPTIONS" << endl;
	cout << " -gen(otype): " << endl;
	cout << " -g: " << "number of individuals (rows) in input matrix" << endl;
	cout << " -n: " << "number of SNPs (columns) in input matrix" << endl;
	cout << " -k: " << "dimension of output matrices (i.e., number of factors)" << endl;
	cout << " \t specify input genotype file (n (individuals) rows, p (loci) columns) with no missing values" << endl;
	cout << " -o or -out(put):" << endl; 
	cout << " \t specify the prefix of all output files, use random seed as a default value" << endl; 
	cout << " -r or -rand:" << endl;
	cout << " \t specify random seed, use system time as default" << endl; 
	cout << " -iter: " << "number of iterations of ECME, default is 20     " << endl; 
	cout << " -t: " << "transpose the matrix after it is read in" << endl;
	cout << " " << endl; 
	cout << " MODEL OPTIONS" << endl;
	cout << " -mg: " << "include a mean vector for the n rows (default is no mean vector)" << endl;
	cout << " -vg: " << "the psi variables are the n row variances (default)" << endl;
	cout << " -mn: " << "include a mean vector for the p columns (default is no mean vector)" << endl;
	cout << " -vn: " << "the psi variables are the p column variances" << endl;
	cout << " " << endl; 
	cout << " OTHER OPTIONS" << endl;
	cout << " -h or -help: " << "print this help" << endl;
	cout << " " << endl; 
}

void CtrlParam::PrintHeader(void)
{
  cout << endl; 
  cout << "This software applies sparse factor analysis to the input matrix." << endl;
  cout << "References: Engelhardt and Stephens." << endl;
  cout << "SFA is a command line program. Open a terminal window, navigate to the directory where " << endl;
  cout << "the executable lives, and type sfa -h for more help. See also the pdf instructions." << endl; 
}

#if defined (READLINE)
void CtrlParam::Run(void)
{
  pMD = new class ModelnData; 
  pMD->open_log(); 
  PrintHeader(); 
  
  int exit = 0; 
  while(!exit) {
    string str;
    string ss; 
    str.assign(rl_gets()); 
    str.append(" \n"); 
    int beg = str.find_first_not_of(" ", 0);
    int end = str.find_first_of(" ", beg);
    if(beg == end)
      ss.assign("");
    else
      ss.assign(str.substr(beg, end-beg));
    if(ss.compare("ls") == 0) str.insert(end, " -G ");
    
    switch(hcom[ss]) {
    case com_EXIT:
      pMD->close_log(); 
      exit = 1; 
      break;
    case com_SET:
      OnSet(str);
      break;
    case com_READ:
      OnRead(str);
      break;
    case com_PRINT:
      //OnPrint(str);
      break;
    case com_BF:
      OnBF(str); 
      break;
    case com_HELP:
      cout << endl; 
      cout << " read: \t read genotype and phenotype data. use read -help for more. " << endl;
      cout << " em: \t run em and set up related parameters. use em -help for mroe. " << endl;
      cout << " bf: \t calculate bf values. use bf -help for more." << endl; 
      cout << " set: \t set and print parameter values. use set -help for more." << endl; 
      cout << " exit: \t write log file and exit." << endl; 
      cout << "       \t support popular unix command. ls, man, top, etc..." << endl;
      cout << endl; 
      break;
    case com_unix:
      system(str.data());
      if(ss.compare("clear") == 0)
	PrintHeader(); 
      break;
    case com_empty:
      break; 
    default:
      cout << "command not found, use \"help\" and try again" << endl;
      break;
    }
  }
}

// Not edited
//print out nEMRuns; nMultiSnp; level; blahblah....
//and can change them. 
void CtrlParam::OnSet(string str)
{
  int beg = str.find_first_of("-", 0);
  int end = str.find_first_of(" ", beg); 
  while (beg != (int) string::npos && end != (int) string::npos) {
    string opt(str.substr(beg, min(5,end-beg))); 
    beg = str.find_first_not_of(" ", end); 
    end = str.find_first_of(" ", beg); 
    string arg(str.substr(beg, end-beg));
    switch(hcom[opt])
      {
      case com_iters:
	pMD->SetIters(atoi(arg.data()));
	break;			
      case com_help:
	cout << endl;
	cout << " -e or -em: " << endl;
	cout << " \t specify number of EM runs, default is 5" << endl; 
	cout << " -s or -step: " << endl; 
	cout << " \t specify steps of each EM run, default is 10" << endl; 
	cout << " -w or -warm: " << endl;
	cout << " \t spcefiy steps of warm up EM run, default is 10" << endl; 
	cout << " -c or -clus(ter): " << endl;
	cout << " \t specify number of clusters in EM algorithm, default is 10." << endl; 
	cout << " \t note: if you rerun em, only use -s option" << endl; 
	cout << " -i or -imp(ute): " << endl;
	cout << " \t specify number of imputations for each EM run, default is 200." << endl;
	cout << " -m or -mul(tiple): " << endl;
	cout << " \t specify number of top snps for multiple snp interrogation, default is 20." << endl;
	cout << " -mg or -mn: " << endl;
	cout << " \t include row means (mg) or column means (mn) in model" << endl;
	cout << " -C or -Com(bination): " << endl;
	cout << " \t maximum number of SNP combinations to be calculated, default is 2000." << endl; 
	cout << " -l or -lev(el): " << endl; 
	cout << " \t specify maximum number of SNPs in all combinations, default is 4." << endl;
	cout << " -pval:" << endl;
	cout << " \t calculate p-value if argument is 1." << endl;
	cout << " -msnp:" << endl;
	cout << " \t calculate multi-snp bf if argument is 1." << endl; 
	cout << endl;
	return;
      default:
	return; 
	
      }
    beg = str.find_first_of("-", end); 
    end = str.find_first_of(" ", beg);
  }
  
}

void CtrlParam::OnRead(string str)
{
	str.append(" \n"); 
	int beg = str.find_first_of("-", 0);
	int end = str.find_first_of(" ", beg); 
	while (beg != (int) string::npos && end != (int) string::npos) {
	  if(str.substr(beg, end-beg).compare("-clean") == 0) {
	    if(pMD) {
	      //	pMD->CleanUp(); 
	      delete pMD; 
	    }
	    pMD = new class ModelnData; 
	    m_nImputeRepeat = 0; 
	    has_read = 0; 
	    cout << "clean up the file" << endl;  
	    return; 
	  }
	  //		string opt(str.substr(beg, min(5,end-beg))); 
	  string opt(str.substr(beg, end-beg)); 
	  beg = str.find_first_not_of(" ", end); 
	  end = str.find_first_of(" ", beg); 
	  string arg(str.substr(beg, end-beg));
	  switch(hcom[opt])
	    {
	    case com_geno:
	      if(arg.empty()) continue;
	      pMD->SetvGin(arg);
	      break;
	    case com_out:
	      if(arg.empty()) continue;
	      pMD->SetfnOutput(arg);
	      break;
	    case com_nPH:
	      pMD->SetnPH(atoi(arg.data())); 
	      break;
	    case com_nFa:
	      pMD->SetnFa(atoi(arg.data())); 
	      break;
	    case com_nIN:
	      pMD->SetnIN(atoi(arg.data())); 
	      break;
	    case com_rand:
	      pMD->SetrandSeed(atoi(arg.data()));
	      break;
	    case com_help:
	      cout << endl; 
	      cout << " -gen(otype): " << endl;
	      cout << " \t specify input genotype files, can be used repeatedly to inpute multiple files." << endl;
	      cout << " -o or -out(put):" << endl; 
	      cout << " \t specify the prefix of all output files, use random seed as a default value." << endl; 
	      cout << " -g or -G:" << endl; 
	      cout << " \t specify how many columns in the input, default is 1." << endl;
	      cout << " -n or -N:" << endl; 
	      cout << " \t specify how many rows in the input, default is 1." << endl;
	      cout << " -k or -K:" << endl; 
	      cout << " \t specify number of factors, default is 1." << endl;
	      cout << " -clean : clean up files and memory." << endl; 
	      cout << endl;
	      return; 
	    default:
	      return; 
	    }
	  beg = str.find_first_of("-", end); 
	  end = str.find_first_of(" ", beg); 
	}
	if (pMD->Getnum_Gin() == 0) {   
	  cout << ":o please provide genotype data by -gen. " << endl;
	  return; 
	}
	else
	  int rs = pMD->GetrandSeed();  // default returns -1;
	if (rs <= 0)
	  rs = (unsigned) time(NULL);
	gsl_rng_set(gsl_r, rs);
	string tstr; 
	tstr.assign(pMD->GetfnOutput());
	if(tstr.compare("\0") == 0) {
	  char tmp[100];
	  memset(tmp, '\0', 100);
	  sprintf(tmp, "%d", rs); 
	  tstr.assign(tmp); 
	  pMD->SetfnOutput(tstr); 
	  cout << ":o prefix of output = " <<  tstr << endl;
	}
	else
	  cout << ":o prefix of output = " <<  tstr << endl;
	//if no provided random seed, use system clock;    
}
#endif 

void CtrlParam::BatchRun(int argc, char ** argv)
{
	pMD = new class ModelnData;

	for(int i = 1; i < argc; i++) {   
	  string str;  
	  string opt; 
	  if(argv[i][0] != '-') 
	    continue;
	  
	  str.assign(argv[i]);
	  opt.assign(str, 0, str.length()); 
	  
	  switch (hcom[opt]) {			
	  case com_geno:
	    if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
	    str.clear();
	    str.assign(argv[i+1]);
	    pMD->SetvGin(str);
	    break;
	  case com_out:
	    if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
	    str.assign(argv[i+1]);
	    pMD->SetfnOutput(str);
	    break;

	  case com_meansg:  
	    pMD->SetHasRowMeans();
	    break;
	  case com_meansn:  
	    pMD->SetHasColumnMeans();
	    break;
	  case com_varg:  
	    pMD->SetHasRowCovariance();
	    break;

	  case com_varn:  
	    pMD->SetHasColumnCovariance();
	    break;

	  case com_rand: // seed 
	    if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
	    if(!isdigit(argv[i+1][0]))
	      {
		cout << "wrong argument after option." << endl; 
		exit(0); 
	      }
	    pMD->SetrandSeed(atoi(argv[i+1]));
	    break;

	  case com_iters:
	    if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
	    if(!isdigit(argv[i+1][0]))
	      {
		cout << "wrong argument after option." << endl; 
		exit(0); 
	      }
	    pMD->SetIters(atoi(argv[i+1]));
	    break;			

	  case com_nPH:
	    if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
	    if(!isdigit(argv[i+1][0])) {
	      cout << "wrong argument after option." << endl; 
	      exit(0); 
	    }
	    pMD->SetnPH(atoi(argv[i+1]));
	    break;
	  case com_nFa:
	    if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
	    if(!isdigit(argv[i+1][0]))
	      {
		cout << "wrong argument after option." << endl; 
		exit(0); 
	      }
	    pMD->SetnFa(atoi(argv[i+1]));
	    break;
	  case com_nIN:
	    if(argv[i+1] == NULL || argv[i+1][0] == '-') continue;
	    if(!isdigit(argv[i+1][0])) {
		cout << "wrong argument after option." << endl; 
		exit(0); 
	      }
	    pMD->SetnIN(atoi(argv[i+1]));
	    break;
	  case com_trans:
	    pMD->SetTrans(1);
	    break;
	  case com_help:
	    OnHelp();
	    exit(0);
	    break;
	    	    
	  default:
	    fprintf(stderr,"Bad option %s\n", argv[i]);
	    OnHelp();
	    exit(0);
	    break; 
	  }              
	}
	
	int rs = 0;
	string str; 
	rs = pMD->GetrandSeed();  // default returns -1;
	if (rs <= 0)
	  rs = (unsigned) time(NULL);
	gsl_rng_set(gsl_r, rs);
	
	str.assign(pMD->GetfnOutput());
	if(str.compare("\0") == 0) {
	  char tmp[100];
	  memset(tmp, '\0', 100);
	  sprintf(tmp, "%d", rs); 
	  str.assign(tmp); 
	  pMD->SetfnOutput(str); 
	}
	//if no provided random seed, use system clock;    
	
	string prefix(pMD->GetfnOutput()); 
	pMD->sfa_model();
	delete pMD; 
}

