/***************************************************************************************  
 *   nmr.h, header to define the NMR abstract class.                                   *
 *                                                                                     *  
 *                                                                                     *
 **************************************************************************************/
#ifndef NMR_H

#define NMR_H

#include <iostream>
#include <string>

using namespace std;

class nmr
{

 public:

  nmr();
  nmr(string , string );

  virtual ~nmr();

  void setIniFile( string );
  void setOutFile( string );
  string getIniFile();
  string getOutFile();

  virtual void read()=0; 
  virtual void printTensor(ostream &)=0;
  virtual void write()=0; 

 private: 

  string iniFile;
  string outFile; 

};

#endif
