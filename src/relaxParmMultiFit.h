/***************************************************************************************  
 *   relaxParmMultiFit.h, header for relaxParmMultiFit.cpp to optimize relaxation      *
 *   parameters (R1,R2,and NOE) using genetic algoriithm from GAlib.                   *
 *                                                                                     *  
 *                                                                                     *
 **************************************************************************************/
#ifndef RELAXPARM_MULTIFIT_H

#define RELAXPARM_MULTIFIT_H

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

class relaxParmMultiFit 
{
  

 public:

  relaxParmMultiFit( string = "relaxParmMultiFit.inp");
  ~relaxParmMultiFit(); 

  void GAOptimization();
  void write();
  
  unsigned long  gaSteps;                 // total number of GA generations
  unsigned int nPrints;                   // the number of GA steps for printing on screen 
  unsigned int gaSeed;                    // seed for the random number generator 
  double minError;                        // if the minError is reached, then the program stops   
  string difFile1,difFile2;               // file for initial diffusion tensor and internal motions 
  string vecFile1,vecFile2;               // file for spin-spin vectors 
  string outFile1,outFile2;               // file for output the final result 
  string prmFile,finFile;                 // the initial parameter and final result files

  static bool isOptTauOvr,isOneTensor;    // control the optimization of overall rotations 
  double minTauOvr, maxTauOvr;            // minimum and maximum values of fitting range for overall rotations (in percentage)
  static bool isOptAngles;                // control the optimization of vector angles
  double minPhi, maxPhi, minPsi, maxPsi;  // minimum and maximum values of fitting range for vector angles (in degree)
  static bool isOptTauInt;                // control the optimization of time constants of internal motions
  double minTauInt, maxTauInt;            // minimum and maximum values of fitting range for time constants of internal motions
  static bool isOptAInt;                  // control the optimization of coefficients of internal motions
  double minAInt, maxAInt;                // minimum and maximum values of fitting range for coefficients of internal motions
    
 private:
  
  relaxParm r1r2ParmExp1,r1r2ParmExp2; 
  static relaxParm r1r2ParmOpt,r1r2ParmOpt1,r1r2ParmOpt2;
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
