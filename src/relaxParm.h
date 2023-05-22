
/***************************************************************************************  
 *   relaxParm.h, header to define a derived class from NMR class for relaxation       *
 *   parameters (R1,R2,and NOE).                                                       *
 *                                                                                     *  
 *                                                                                     *
 **************************************************************************************/
#ifndef RELAXPARM_H

#define RELAXPARM_H

#include <vector>
#include <string>
#include <cmath>
#include <cstring>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "myNamespace.h"
#include "nmr.h"

using namespace std; 

class relaxParm : public nmr 
{

 public:

  relaxParm();                             
  relaxParm(string, string, string);       // constructor with specified initial files
  // relaxparm(relaxParm & );              // copy constructor  

  ~relaxParm();                            // distructor 

  void read();                             // input initial diffusion tensor and inter motions from file
  void printTensor(ostream &);             // print out the information of diffusion tensor
  void write();                            // output diffusion tensor and relaxation parameters

  void calcVecAngles();                    // calculate the angles between spin vectors and difusion tensor          
  void calcCoeffsCFOVec();                 // calculate coefficients of overall RCF for all spin vectors for vecs   
  void calcCoeffsCFOAng();                 // calculate coefficients of overall RCF for all spin vectors for angles
  void updateDiffTensorEVal();             // calculate eigen values of diffusion tensor from of time constants of overall RCF 
  void calcCoeffsCFT();                    // calculate coefficients and time constants of total RCF for all spin vectors
  void calcRelaxParm();                    // calculate relaxation parameters
  void calcRFactor();                      // calculate the R factor 

  void printRelaxParm(ostream &);          // print out relaxation parameters
  void printRelaxInt(ostream &);           // print out information about internal motions
  void printSpinVecs(ostream &);           // print out the x, y, and z direction of spin vectors
  void printVecAngles(ostream &);          // print out the spherical angles of spin vectors in difussion tensor frame 
 
  string getVecFile() const {return vecFile;}  
  void setVecFile(const string vecFile0) {vecFile = vecFile0;}

  unsigned getNumModesOvr() const {return numModesOvr;}
  void setNumModesOvr(const unsigned numModesOvr0) {numModesOvr = numModesOvr0;}
  vector <double> getTauOvr() const {return tauOvr;}
  void setTauOvr(const vector <double> & tauOvr0) {tauOvr = tauOvr0;}

  unsigned getNumSpinVecs() const {return numSpinVecs;}
  void setNumSpinVecs(const unsigned numSpinVecs0) {numSpinVecs = numSpinVecs0;}
  vector < vector < double > > getVecAngles() const {return vecAngles;}
  void setVecAngles(const  vector < vector < double > > & vecAngles0) {vecAngles = vecAngles0;}
  vector < vector < double > > getAOvr() const {return AOvr;}
  void setAOvr(const  vector < vector < double > > & AOvr0) {AOvr = AOvr0;}

  vector < vector < double > > getTauInt() const {return tauInt;}
  void setTauInt(const  vector < vector < double > > & tauInt0) {tauInt = tauInt0;}
  vector < vector < double > > getAInt() const {return AInt;}
  void setAInt(const  vector < vector < double > > & AInt0) {AInt = AInt0;}
 
  vector < vector < double > > getTauTot() const {return tauTot;}
  void setTauTot(const  vector < vector < double > > & tauTot0) {tauTot = tauTot0;}
  vector < vector < double > > getATot() const {return ATot;}
  void setATot(const  vector < vector < double > > & ATot0) {ATot = ATot0;}

  vector < double > getR1() const {return R1;}
  vector < double > getR2() const {return R2;}
  vector < double > getNOE() const {return NOE;}
   
  double getRFactor() const {return RFactor;}

 private: 

  string vecFile;                          // vector file store the spin vectors 

  double diffTensor[3][3];                 // 3 X 3 diffusion tensor   (x,y,z) X (x,y,z)
  double diffTensorEVal[3];                // eigenvalues of diffusion tensor in decreasing order of Z, Y, X
  double diffTensorEVec[3][3];             // eigenvectors of diffusion tensor in decreasing order of Z, Y, X
  double D, L2;                             
  
  unsigned  nBFields;                      // number of B fields 
  vector < double > BFields;               // B values of magnetic field
 
  unsigned numSpinVecs;                    // total number of spin-spin vectors
  vector < unsigned > resI;                // residue ID for the i spin
  vector < unsigned > resJ;                // residue ID for the j spin 
  vector < string > atmTypeI;              // atom type for the i spin
  vector < string > atmTypeJ;              // atom type for the j spin
  vector < vector < string > > spinTypes;  // the spin pair types such as C-H and N-H 
  vector < vector < double > > spinVecs;   // exponential decay constants for internal motions
  vector < vector < double > > vecAngles;  // the angles between spin vectors and the diffusion tensor

  vector < double > R1In, R1;              // input and calculated R1       
  vector < double > R2In, R2;              // input and calculated R2
  vector < double > NOEIn, NOE;            // input and calculated NOE

  unsigned numModesOvr;                    // total number of modes for the overall rotation
  vector < double >  tauOvr;               // exponential decay constants for the overall rotation 
  vector < vector < double > > AOvr;       // coefficients for the exponential decays of overall rotation 
  vector < unsigned > numModesInt;         // total number of modes for internal rotation
  vector < vector < double > > tauInt;     // exponential decay constants for internal motions
  vector < vector < double > > AInt;       // coefficients for the exponential decays of internal motion 
  vector < vector < double > > tauTot;     // exponential decay constants for total RCF
  vector < vector < double > > ATot;       // coefficients for the exponential decays of total RCF
  
  double RFactor;                          // R factor comparing to experimental values

  double getJw(vector <double> & , vector <double> &, double ); // calculate J(w)  
  void assignSpinTypes();                                       // get the spin types from the atom types
};

#endif
