/***************************************************************************************  
 *   relaxParmFit.h, header for relaxParmFit.cpp to optimize relaxation parameters     *
 *  (R1,R2,and NOE) using genetic algoriithm from GAlib.                               *
 *                                                                                     *  
 *                                                                                     *
 **************************************************************************************/
#ifndef RELAXPARMFIT_H

#define RELAXPARMFIT_H

#include <vector>
#include <string>
#include <cmath>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <ga/ga.h>
#define INSTANTIATE_REAL_GENOME
#include <ga/GARealGenome.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "relaxParm.h"

using namespace std; 

class relaxParmFit 
{
  

 public:

  relaxParmFit( string = "relaxParmFit.inp");
  ~relaxParmFit(); 

  void GAOptimization();
  void write();
  
  unsigned long  gaSteps;                 // total number of GA generations
  unsigned int nPrints;                   // the number of GA steps for printing on screen 
  unsigned int gaSeed;                    // seed for the random number generator 
  double minError;                        // if the minError is reached, then the program stops   
  string difFile;                         // file for initial diffusion tensor and internal motions 
  string vecFile;                         // file for spin-spin vectors 
  string outFile;                         // file for output the final result 

  static bool isOptTauOvr;                // control the optimization of overall rotations
  double minTauOvr, maxTauOvr;            // minimum and maximum values of fitting range for overall rotations (in percentage)
  static bool isOptAngles;                // control the optimization of vector angles
  double minPhi, maxPhi, minPsi, maxPsi;  // minimum and maximum values of fitting range for vector angles (in degree)
  static bool isOptTauInt;                // control the optimization of time constants of internal motions
  double minTauInt, maxTauInt;            // minimum and maximum values of fitting range for time constants of internal motions
  static bool isOptAInt;                  // control the optimization of coefficients of internal motions
  double minAInt, maxAInt;                // minimum and maximum values of fitting range for coefficients of internal motions
    
 private:
  
  relaxParm r1r2ParmExp; 
  static relaxParm r1r2ParmOpt;
  static float Objective(GAGenome &);                                     // calculate the R factor                 
  static void calcRelaxParmFromGenome(relaxParm &, const GARealGenome &); // calculate relaxation parameters from genome 
  static float Comparator(const GAGenome &, const GAGenome &);            // define the comparator 
  

  void setFitRange(GARealAlleleSetArray &);                    // set the fitting ranges for all parameters
  void initGAMethod(GASimpleGA &);                             // initialize the GA method for evolving

  void initPopulations(GASimpleGA &);                         // initialize the populations of first generations
  void printPopulations(ostream &, GARealGenome &, const GASimpleGA &); // print populations to output stream
  void writePopResults(ostream &, GARealGenome &, const GASimpleGA &);  // write R1, R2, and NOE from the whole ensemble 

};

#endif
