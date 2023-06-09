/*********************************************************************************************  
 *   multiFitMC.cpp,                                                                         *
 *   a program to optimize the populations of multi-structures for the NMR data using        *
 *   Monte Carlo algorithm.                                                                  * 
 *                                                                                           *
 *                                                                                           *  
 *   Authors:                                                                                *
 *          originally written by Junchao Xia                                                *
 *                                                                                           *
 *                                                                                           *
 *   Modification Histry:                                                                    *
 *          Jan. 29, 2010:   originally written by Junchao Xia                               *
 *          Mar. 23, 2012:   rearrange the program structures by Junchao Xia                 *
 *                                                                                           *
 *********************************************************************************************/

#include "multiFitMC.h"

MultiFitMC::MultiFitMC(string inputfile)
{
  ifstream inputparameter;

  inputparameter.open(inputfile.data(),ios::in);

  if(inputparameter.fail())
    { cout<<" Error opening the input parameter file.\n"<< inputfile.data() << " doesn't not exist."<<endl;
      cout<<" CHECK FILENAME !\n";
      cout<<" CHECK FILENAME !\n";
    }
  inputparameter.ignore(256,'\n');
  inputparameter.ignore(256,'\n');
  inputparameter.ignore(256,'\n');
 
  inputparameter >> MCSteps; inputparameter.ignore(120,'\n');
  cout << " The total number Monte Carlo steps is " << MCSteps<< endl;

  inputparameter >> nprints; inputparameter.ignore(120,'\n');
  cout << " The number steps for printing is " << nprints<< endl;

  inputparameter >> MCSeed; inputparameter.ignore(120,'\n');
  cout << " The seed for random number generator is " << MCSeed<< endl;

  inputparameter >> minError; inputparameter.ignore(120,'\n');
  cout << " The minimum R factor to stop simulation is " << minError << endl;
 
  inputparameter >> nmrExpFile; inputparameter.ignore(120,'\n');
  cout << " The nmr data file from the experiments is " << nmrExpFile << endl;

  inputparameter >> numNMRPoints; inputparameter.ignore(120,'\n');
  cout << " The total number of NMR points is " << numNMRPoints << endl;

  inputparameter >> nmrConfFile; inputparameter.ignore(120,'\n');
  cout << " The nmr data file for all conformers is " << nmrConfFile << endl;

  inputparameter >> numConfs; inputparameter.ignore(120,'\n');
  cout << " The total number of conformers is " << numConfs << endl;

  inputparameter >> rankCtrl; inputparameter.ignore(120,'\n');
  cout << " The ranking criteria for the prediction list is " << rankCtrl << endl;

  inputparameter >> randomMethod; inputparameter.ignore(120,'\n');
  cout << " The random method for weighting is " << randomMethod << endl;

  inputparameter.close();

  ifstream nmrConfInput;

  // input the nmr data for all conformations
  nmrConfInput.open(nmrConfFile.data(), ios::in);
  if( nmrConfInput.fail() )
    {
      cout << "Failed to open the NMR data file for all conformations " << nmrConfFile.data() << endl;
      exit(0);
    }
  
  for(int iconfs = 0; iconfs < numConfs; iconfs++)
    {
      int tempi, ranki; 
      nmrConfInput >> tempi;
      nmrConfInput >> ranki;
      rankConfID.push_back(ranki);

      float deviation; 
      nmrConfInput >> deviation;
      nmrDevs.push_back(deviation);

      float nmrvalue; 
      vector <float> temp_vect;
      for(int inmr = 0; inmr < numNMRPoints; inmr++)
	{
	  nmrConfInput >> nmrvalue; 
	  temp_vect.push_back(nmrvalue);
	}
      nmrConfInp.push_back(temp_vect);      
      nmrConfInput.ignore(numeric_limits<streamsize>::max(),'\n');
      
    }
  nmrConfInput.close();

  ifstream nmrExpInput;

  // input the experimental nmr data
  nmrExpInput.open(nmrExpFile.data(), ios::in);
  if( nmrExpInput.fail() )
    {
      cout << "Failed to open the NMR data file for experimental values " << nmrExpFile.data() << endl;
      exit(0);
    }

  float nmrvalue; 
  for(int inmr = 0; inmr < numNMRPoints; inmr++)
    {
      nmrExpInput >> nmrvalue; 
      nmrExpInp.push_back(nmrvalue);
      nmrExpInput.ignore(numeric_limits<streamsize>::max(),'\n');
    }
  nmrExpInput.close();

  for(int iconfs = 0; iconfs < numConfs; iconfs++)
    {
      probFactors.push_back(0.0); 
      bestWeights.push_back(0.0); 
    }

  for(int inmr = 0; inmr < numNMRPoints; inmr++)
    {
      currCalcNMR.push_back(0.0);
      bestCalcNMR.push_back(0.0);
    }

}

MultiFitMC::~MultiFitMC()
{
  // need to add operations to free the memory
}


void MultiFitMC::MonteCarloFit()
{
  
  float rfactor;
  unsigned long istep;
  float sumProb;
  int iconfs; 
  int inmr;

  // srand(static_cast<unsigned>(time(0)));
  srand(MCSeed);
  bestRfactor =numeric_limits<float>::max(); 

  for (istep = 0; istep < MCSteps; istep++)
    {

      if (randomMethod == 0) 
	{
	  sumProb = 0.0;

	  for (inmr = 0; inmr < numNMRPoints; inmr++)
	    {currCalcNMR[inmr]= 0.0;}

	  for(iconfs = 0; iconfs < numConfs; iconfs++)
	    {
	      if (iconfs == (numConfs - 1)) 
		{probFactors[iconfs] = 1.0-sumProb;}
	      else
		probFactors[iconfs] = (static_cast<float>(rand())/RAND_MAX)*(1.0-sumProb);
	      sumProb = sumProb + probFactors[iconfs];
	  
	      
	      for (inmr = 0; inmr < numNMRPoints; inmr++)
		{
		  currCalcNMR[inmr] += probFactors[iconfs]*nmrConfInp[iconfs][inmr];
		}
	    }
	  
	}
      else if (randomMethod == 1)
	{

	  sumProb = 0.0; 
	  for(iconfs = 0; iconfs < numConfs; iconfs++)
	    {
	      probFactors[iconfs] = (static_cast<float>(rand())/RAND_MAX);
	      sumProb = sumProb + probFactors[iconfs];
	    }

	  for(iconfs = 0; iconfs < numConfs; iconfs++)
	    probFactors[iconfs]=probFactors[iconfs]/sumProb;

	  for (inmr = 0; inmr < numNMRPoints; inmr++)
	    {currCalcNMR[inmr]= 0.0;}
	  
	  for(iconfs = 0; iconfs < numConfs; iconfs++)
	    {
	      for (inmr = 0; inmr < numNMRPoints; inmr++)
		{
		  currCalcNMR[inmr] += probFactors[iconfs]*nmrConfInp[iconfs][inmr];
		}
	    }
	}
      else 
	{
	  cout << " Wrong value of randomMethod.\n"; 
	  exit(0);
	}
      
      rfactor = calcRfactor();

      if (rfactor < bestRfactor)
	{
	  bestRfactor = rfactor; 
	  bestWeights = probFactors;
	  bestCalcNMR = currCalcNMR;
	  
	}

      if (istep%nprints == 0 ) cout << " MC step = " << istep << " best R factor = " << bestRfactor << endl;

      if (bestRfactor < minError) 
	{
	  cout << " Find the solution with the best R factor = " << bestRfactor << " in MC step = " << istep << endl; 
	  break; 
	}
      
    }
  
  if (istep == MCSteps) 
    {
      cout << " Maximum MC step reached, the best R factor = " << bestRfactor << " in MC step = " << istep << endl; 
    }
  

}

float MultiFitMC::calcRfactor()
{
  float numerator = 0.0;
  float denominator = 0.0; 
  
  for (int inmr = 0; inmr < numNMRPoints; inmr++)
    {
      numerator += (currCalcNMR[inmr]-nmrExpInp[inmr])*(currCalcNMR[inmr]-nmrExpInp[inmr]);
      denominator += nmrExpInp[inmr]*nmrExpInp[inmr];
    }

  return sqrt(numerator/denominator);

}

void MultiFitMC::writeToFile()
{

  ofstream nmrBestOut;

  nmrBestOut.open("nmrBest.out", ios::out);
  
  for (int inmr=0; inmr<numNMRPoints; inmr++)
    {
     
      nmrBestOut << setiosflags(ios::fixed) << setw(10) << setprecision(5) << nmrExpInp[inmr]
                 << setiosflags(ios::fixed) << setw(10) << setprecision(5) << bestCalcNMR[inmr] << endl;

    }

  nmrBestOut.close();

  ofstream probBestOut;
  probBestOut.open("probBest.out", ios::out);
  for (int iconfs=0;iconfs<numConfs;iconfs++)
    {
      
      probBestOut << setw(10) << iconfs + 1 << setw(10) << rankConfID[iconfs]
      << setiosflags(ios::fixed) << setw(15) << setprecision(5) << nmrDevs[iconfs] 
		  << setiosflags(ios::fixed) << setw(15) << setprecision(5) << bestWeights[iconfs] << endl; 
    }
  probBestOut.close();
  
}


// argv[1] = param_file
int main(int argc, char **argv)
{

  string param_file; 

  time_t seconds, t_start, t_end;

  t_start = time(NULL);
  
  if(argc<2)
    {
      cerr << " Check the command syntax!" << endl;
      cerr << " It should be 'multifit param_file'" << endl;
      return 1;
    }


  cout << "************************************************************************" << endl;
  cout << "*                    START MULTIPLE-STRUCTURE FITTING                  *" << endl;
  cout << "*    FROM AN ENSEMBLE of CONFORMATIONS USING MONTE CARLO ALGORITHM.    *" << endl;
  cout << "************************************************************************" << endl;

  cout << "*************************Reading Parameters:****************************" << endl;
  param_file = argv[1];

  MultiFitMC myMultiFit = MultiFitMC(param_file);    // read in the parameters to start program
  cout << " Done with reading parameters." << endl;
  myMultiFit.MonteCarloFit();
  myMultiFit.writeToFile(); 

  cout << " Done with the searching of the multi-structure solution for the NMR data." << endl; 

  t_end = time(NULL);
  seconds =  t_end - t_start;
  cout << "***********************Output the CPU time information:*****************" << endl;
  cout << " The CPU time for this search is " << seconds  << " seconds." << endl;

  return(0); 
 
}
