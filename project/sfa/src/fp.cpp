#include <iostream>
#include <time.h>
#include "fp.h"
#include "control.h"
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

int main(int argc, char * argv[])
{ 
  
  CtrlParam objCtrl;
  
  //check if dir exist. create if not. 
  std::ifstream check("output/");
  if(!check) {
    mkdir("output", S_IRWXU); 
  } 
  
  if(argc <= 1)
    objCtrl.PrintHeader(); 
    objCtrl.BatchRun(argc, argv); 
    
    return 1;                                                          
}

 
