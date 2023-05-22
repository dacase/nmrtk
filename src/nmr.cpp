#include "nmr.h"



nmr::nmr()
{
}

nmr::nmr(string iniFile0, string outFile0):
  iniFile(iniFile0),outFile(outFile0)
{
  // iniFile = iniFile0; 
  // outFile = outFile0; 
}
 
nmr::~nmr()  
{
}
 
void nmr::setIniFile(string iniFile0)  
{
  iniFile = iniFile0;
}

void nmr::setOutFile(string outFile0)  
{
  outFile = outFile0;
}

string nmr::getIniFile()
{
  return iniFile;
}

string nmr::getOutFile()  
{
  return outFile;
}

