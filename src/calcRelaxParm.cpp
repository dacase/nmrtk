/*********************************************************************************************  
 *   calcRelaxParm.cpp,                                                                      *
 *     a program to calculate the R1, R2 and NOE values from the overall and internal        *
 *     motions.                                                                              * 
 *                                                                                           *
 *                                                                                           *
 *                                                                                           *
 *********************************************************************************************/
#include <cmath>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <string>

#include <vector>
#include <limits>

#include "relaxParm.h"
 
// argv[1] = ini_file  argv[2] = vec_file  argv[3] = out_file 
int main(int argc, char **argv)
{



  time_t seconds, t_start, t_end;

  t_start = time(NULL);
  
  if(argc<4)
    {
      cerr << " Check the command syntax!" << endl;
      cerr << " It should be 'relaxParmFit ini_file vec_file out_file'" << endl;
      return 0;
    }

  string ini_file, vec_file, out_file; 
  ini_file = argv[1]; vec_file = argv[2]; out_file = argv[3];


  cout << "************************************************************************" << endl;
  cout << "*                    START RELAXATION PARAMETER CALCULATIONS.          *" << endl;
  cout << "************************************************************************" << endl;
 
  relaxParm R1R2Comp=relaxParm(ini_file,vec_file,out_file); 
  // R1R2Comp.read();
  cout << "Finished reading rarameters." << endl;

  R1R2Comp.calcCoeffsCFOVec();
  cout << "Finished calculating coefficients for overall rotation CFs." << endl;
  R1R2Comp.calcCoeffsCFT();
  cout << "Finished calculating coefficients and time constants for total rotation CFs." << endl;
  R1R2Comp.calcRelaxParm();
  cout << "Finished calculating relaxation parameters." << endl;
  R1R2Comp.calcRFactor();
  cout << "Finished calculating R factor." << endl;
  R1R2Comp.write(); 
  cout << "Finished writing results." << endl;
  t_end = time(NULL);

  seconds = t_end - t_start;
  cout << "Finished the job in " << seconds << " seconds." << endl;
}
 





 
 
 








