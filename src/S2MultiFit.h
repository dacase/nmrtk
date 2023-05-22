/***************************************************************************************  
 *   relaxParmMultiFit.h, header for relaxParmMultiFit.cpp to optimize relaxation      *
 *   parameters (R1,R2,and NOE) using genetic algoriithm from GAlib.                   *
 *                                                                                     *  
 *                                                                                     *
 **************************************************************************************/
#ifndef S2_MULTIFIT_H

#define S2_MULTIFIT_H

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


using namespace std; 

class S2MultiFit 
{
  

 public:

  S2MultiFit( string = "S2MultiFit.inp");
  ~S2MultiFit(); 

  void GAOptimization();
  void write();
  
  unsigned long  gaSteps;                 // total number of GA generations
  unsigned int nPrints;                   // the number of GA steps for printing on screen 
  unsigned int gaSeed;                    // seed for the random number generator 
  double minError;                        // if the minError is reached, then the program stops   
  string iniFile;                         // file for S2 of different state  
  string outFile;                         // file for output the final result 

  vector < int > resID;  
  static vector <double> s2Exp; 
  static vector <double> s2Tot; 
  static vector <double> s2Semi;
  static vector <double> s2Close;
  static vector <double> p2Int;
   
 private:
  
  vector <double> s2TotIni;
  vector <double> s2SemiIni;
  vector <double> s2CloseIni;

  static bool isOptS2Ind; 
  static bool isOptS2EQ; 

  bool isRandomS2; 
  double minS2, maxS2; 

  static float Objective(GAGenome &);                               // calculate the R factor                 
  static void calcS2FromGenome(const GARealGenome &);               // calculate relaxation parameters from genome 
  static float Comparator(const GAGenome &, const GAGenome &);      // define the comparator 
  
  void setFitRange(GARealAlleleSetArray &);                         // set the fitting ranges for all parameters
  void initGAMethod(GASimpleGA &);                                  // initialize the GA method for evolving

  void initPopulations(GASimpleGA &);                         // initialize the populations of first generations
  void printPopulations(ostream &, GARealGenome &, const GASimpleGA &); // print populations to output stream
  void writePopResults(ostream &, GARealGenome &, const GASimpleGA &);  // write R1, R2, and NOE from the whole ensemble 

};

#endif
