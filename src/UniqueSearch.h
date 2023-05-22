/*********************************************************************************************  
 *   UniqueSearch.h, header file of UniqueSearch.cpp to find the unique structures.          *
 *                                                                                           *  
 *   Authors:                                                                                *
 *          originally written  by Junchao Xia                                               *
 *                                                                                           *
 *   Modification Histry:                                                                    *
 *         Dec. 10, 2014:    originally written Junchao Xia.                                 *
 *                                                                                           *
 *********************************************************************************************/
#ifndef UNIQUE_SEARCH_HPP
#define UNIQUE_SEARCH_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;


class UniqueSearch
{

  public:

  UniqueSearch();                   // constructor,
  ~UniqueSearch();                  // destructor,
  UniqueSearch(string & );          // constructor with reading the initial condition from a file
  void findUniqueConfs();
  void findUniqueConfsNoEng();
  void writeUniqueConfs();
  vector <int> uniq_ids;
  vector <float> uniq_engs;
  vector < vector < float > > uniq_angs;

  private:

  void readParaFile(ifstream &, ostream &);
  void readEngFile();  
  void readAngFile();                 // input the configurations such as dihedral angle of conformations from a file
  void sortEng(vector <int> &ids, vector <float> &engs); // sort the data such as energies and ids from low to high
  
  
  string engFile;                       // data file name for the resulting value for each conformer  
  string angFile;                       // allowed phi-psi values or other configurations for each conformer
  unsigned int numAng;                             // total number of degree of freedom
  unsigned int numEng;                     	 // total number of conformers
  
  float minEngDiff;                        // minimum value for the difference in data file to judge the unique structur
  vector< float > confEngs;     
  vector< float > minAngDiff;              // minimum values for the conformation difference to judge the unique structure    
  vector< vector<float> > confAngs;       


};

#endif
