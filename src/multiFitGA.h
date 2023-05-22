/*********************************************************************************************  
 *   multiFitGA.hpp, header file of multiFitGA.cpp to                                        *
 *          fit the NMR quantities to multiple structure candidates using Genetic Algorithm. *
 *                                                                                           *  
 *   Authors:                                                                                *
 *          written  by Junchao Xia                                                          *
 *                                                                                           *
 *   Modification Histry:                                                                    *
 *          Mar. 19, 2010: originally written Junchao Xia.                                   *
 *          Mar. 23, 2012: rearrange the program structure and input by Junchao Xia          *
 *                                                                                           *
 *********************************************************************************************/
#ifndef MULTI_FIT_GA_H
#define MULTI_FIT_GA_H


#include <cmath>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <vector>
#include <limits>

#include <stdio.h>

#include <ga/ga.h>
#define INSTANTIATE_REAL_GENOME
#include <ga/GARealGenome.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;


class MultiFitGA
{
  
public:

  MultiFitGA(string = "multiFitGA.inp");  // constructor with reading the initial condition from a file

  ~MultiFitGA();                  // destructor,
  
  void randomNMRData();           // add random error to NMR data
  void GeneticAlgorithmFit();     // perform GA fit
  void printPopulations(ostream &, GARealGenome &, const GASimpleGA &);   // print populations to output stream
  void writeToFile();             // write the result to files

       
private:

  static float Objective(GAGenome &);                             // calculate the objective function (R factor)
  static void  calcNMRFromEnsemble(GARealGenome &);               // calculate NMR data from ensemble average
  static float Comparator(const GAGenome &, const GAGenome &);    // define the comparator 

  static vector< float > nmrExpInp;            // experimental nmr data after random gauss error
  static vector< vector<float> > nmrConfInp;   // NMR data for each conformer
  static vector< float > currCalcNMR;          // the current NMR data set from multiple fit 
  static vector< float > probFactors;          // current probability weight factors

  unsigned long  gaSteps;                 // total number of Monte Carlo steps
  unsigned int nPrints;                   // the number of GA steps for printing on screen 
  unsigned int gaSeed;                    // seed for the random number generator 
  double minError;                        // if the minError is reached, then the program stops      
  
  string nmrExpFile;                      // file for the experimental NMR
  int  numNMRPoints;                      // the total number of NMR points
  vector< double > nmrExpInp0;            // the original experimental nmr data   
  vector< double > nmrExpErr;             // the error for experimental nmr data
  unsigned int isRandomErr;               // if add random error to NMR data
  unsigned int nmrSeed;                   // seed for random errors for NMR data
  double sigmaErr;                        // the variance of error   
  
  string nmrConfFile;                     // NMR data file for all conformer
  int numConfs;                   	  // total number of conformers	 

  char rankCtrl;                          // rank criteria to produce the predection list 
  static int  randomMethod;               // control how to use the random values for weighting 

  vector< float > bestWeights;            // keep the best solution for weighting factors  
  vector<int> rankConfID;                 // the conformer ID for the rank
  vector<float> nmrDevs;                  // nmr deviations of conformers 
  int  bestID;                            // the conformation id for the best solution 
  float bestRfactor;                      // Rfactor for the best solution 
 
  vector< float > bestCalcNMR;            // the best NMR data set from multiple fit

};

#endif
