/*********************************************************************************************  
 *   multiFitMC.hpp, header file of multiFitMC.cpp to                                        *
 *          fit the NMR quantities to multiple structure candidates using Monte Carlo.       *
 *                                                                                           *  
 *   Authors:                                                                                *
 *          written  by Junchao Xia                                                          *
 *                                                                                           *
 *   Modification Histry:                                                                    *
 *          Jan. 29, 2010:    originally written Junchao Xia                                 *
 *          Jul. 15, 2010:    added randomMethod to control the way to add Random values     * 
 *          Mar. 23, 2012:   rearrange the program structures by Junchao Xia                 *
 *                                                                                           *
 *********************************************************************************************/
#ifndef MULTI_FIT_MC_H
#define MULTI_FIT_MC_H


#include <cmath>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <vector>
#include <limits>

using namespace std;


class MultiFitMC
{

public:
  

  MultiFitMC(string ="multiFitMC.inp");  // constructor with reading the initial condition from a file
  ~MultiFitMC();                         // destructor,


  void MonteCarloFit();           // perform Monte Carlo fit
  float calcRfactor();            // calculate the R factor 
  void writeToFile();             // write the result to files
       
private:

  unsigned long  MCSteps;                 // total number of Monte Carlo steps
  unsigned int nprints;                   // the number of MC steps for printing on screen 
  unsigned int MCSeed;                    // seed for the random number generator 
  double minError;                        // if the minError is reached, then the program stops      
  
  string nmrExpFile;                      // file for the experimental NMR
  int  numNMRPoints;                      // the total number of NMR points
  vector< float > nmrExpInp;              // experimental nmr data 
  
  string nmrConfFile;                     // NMR data file for all conformer
  int numConfs;                    	  // total number of conformers	 
  vector< vector<float> > nmrConfInp;     // NMR data for each conformer
  
  char rankCtrl;                          // rank criteria to produce the predection list 
  vector< float > probFactors;            // current probability weight factors
  vector< float > bestWeights;            // keep the best solution for weighting factors  

  int randomMethod;                       // control how to use the random values for weighting 

  vector<int> rankConfID;                 // the conformer ID for the rank
  vector<float> nmrDevs;                  // nmr deviations of conformers 
  int  bestID;                            // the conformation id for the best solution 
  float bestRfactor;                      // Rfactor for the best solution 
  vector< float > currCalcNMR;            // the current NMR data set from multiple fit  
  vector< float > bestCalcNMR;            // the best NMR data set from multiple fit

};

#endif
