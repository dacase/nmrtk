
/*********************************************************************************************  
 *   multiFitGA.cpp,                                                                         *
 *   a program to optimize the population of structure ensemble for NMR data using           *
 *   genetic algorithm.                                                                      * 
 *                                                                                           *
 *                                                                                           *  
 *   Authors:                                                                                *
 *          originally written by Junchao Xia                                                *
 *                                                                                           *
 *                                                                                           *
 *   Modification Histry:                                                                    *
 *          Jan. 29, 2010:   originally written by Junchao Xia                               *
 *          Mar. 23, 2012:   rearrange the program structure by Junchao Xia                  *
 *                                                                                           *
 *********************************************************************************************/

#include "multiFitGA.h"

vector< float > MultiFitGA::nmrExpInp;
vector< vector<float> > MultiFitGA::nmrConfInp;
vector< float > MultiFitGA::currCalcNMR;
vector< float > MultiFitGA::probFactors;
int  MultiFitGA::randomMethod;

MultiFitGA::MultiFitGA(string inputFile)
{

  // some defaults:  (maybe add more later)
  pMutation = 0.001;
  pCrossover = 0.9;

  // input the parameters to initialize fitting 
  ifstream infoFile;
  infoFile.open(inputFile.data(),ios::in); 
  if(infoFile.fail())
    { 
      cerr<<" Error opening the input file to initialize the ensemble optimization using GA. \n";  
      cerr<<" Please check the file: " << inputFile << endl; 
      exit(0);
    }
  cout << "\n Reading the input file for GA fitting: " << inputFile << endl;
 
  string linein;
  string lineHead;
  istringstream ssin;
  while( getline(infoFile, linein))
    {
      cout << linein << endl;
      ssin.clear();
      ssin.str(linein);
      ssin >> lineHead; 
      if (lineHead == "#")
	{continue;}
      else if (lineHead == "populationSize:")
	{ ssin >> populationSize; }
      else if (lineHead == "pMutation:")
	{ ssin >> pMutation; }
      else if (lineHead == "pCrossover:")
	{ ssin >> pCrossover; }
      else if (lineHead == "gaSteps:")
	{ ssin >> gaSteps; }
      else if (lineHead == "nPrints:" )
	{ ssin >> nPrints; }
      else if (lineHead == "gaSeed:" )
	{ ssin >> gaSeed; }
      else if (lineHead == "minError:" )
	{ ssin >> minError; }  
      else if (lineHead == "isRandomErr:" )
	{ ssin >> isRandomErr; }
      else if (lineHead == "nmrSeed:")
	{ ssin >> nmrSeed; }
      else if (lineHead == "sigmaErr:")
	{ ssin >> sigmaErr; }
      else if (lineHead == "nmrExpFile:")
	{ ssin >> nmrExpFile; }
      else if (lineHead == "numNMRPoints:")
	{ ssin >> numNMRPoints; }
      else if (lineHead == "nmrConfFile:")
	{ ssin >> nmrConfFile; }
      else if (lineHead == "numConfs:")
	{ ssin >> numConfs; }
      else if (lineHead == "rankCtrl:")
	{ ssin >> rankCtrl; }
      else if (lineHead == "randomMethod:")
	{ ssin >> randomMethod; }
      else
	{continue;}
    }
  infoFile.close();

  // input the nmr data for all conformations
  nmrConfInp.clear();
  ifstream nmrConfInput;
  nmrConfInput.open(nmrConfFile.data(), ios::in);
  if( nmrConfInput.fail() )
    {
      cerr << "Failed to open the NMR data file for all conformations " << nmrConfFile.data() << endl;
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

  // input the experimental nmr data
  nmrExpInp0.clear(); 
  ifstream nmrExpInput;
  nmrExpInput.open(nmrExpFile.data(), ios::in);
  if( nmrExpInput.fail() )
    {
      cerr << "Failed to open the NMR data file for experimental values " << nmrExpFile.data() << endl;
      exit(0);
    }

  float nmrvalue; 
  for(int inmr = 0; inmr < numNMRPoints; inmr++)
    {
      nmrExpInput >> nmrvalue; 
      nmrExpInp0.push_back(nmrvalue);
      nmrExpInp.push_back(nmrvalue);
      nmrExpInput.ignore(numeric_limits<streamsize>::max(),'\n');
    }
  nmrExpInput.close();
  
  if (isRandomErr)  randomNMRData(); 
  
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

MultiFitGA::~MultiFitGA()
{
  // need to add operations to free the memory
}

void MultiFitGA::randomNMRData()
{
  for (unsigned i = 0; i < nmrExpInp.size(); i++)
    nmrExpErr.push_back(nmrExpInp[i]);

  const gsl_rng_type *T;
  gsl_rng *r; 
  gsl_rng_env_setup();
  T = gsl_rng_default; 
  r = gsl_rng_alloc(T);
  gsl_rng_set(r,nmrSeed);
  
  for(int inmr = 0; inmr < numNMRPoints; inmr++)
    {
      nmrExpErr[inmr] = gsl_ran_gaussian(r, sigmaErr);
      nmrExpInp[inmr] = nmrExpInp[inmr] + nmrExpErr[inmr];
    }

}

void MultiFitGA::GeneticAlgorithmFit()
{

  // initalize Genetic Algorithm
  GARandomSeed(gaSeed);                    
  GARealAlleleSet alleles(0.0,1.0);
  GARealGenome genome(numConfs,alleles,Objective);
  genome.initialize(); 
  GASharing scale(Comparator);    

  GASimpleGA ga(genome);  
  ga.minimize();                    // by default we want to minimize the objective
  ga.scaling(scale);                // set the scaling method to our sharing
  ga.populationSize(populationSize); // how many individuals in the population
  ga.nGenerations(gaSteps);         // number of generations to evolve
  ga.pMutation(pMutation);          // likelihood of mutating new offspring
  ga.pCrossover(pCrossover);        // likelihood of crossing over parents
  ga.scoreFilename("score.out");    // name of file for scores
  ga.scoreFrequency(1000);        // keep the scores of every 1000 generations
  ga.flushFrequency(50);            // specify how often to write the score to disk
  ga.selectScores(GAStatistics::AllScores);
  ga.initialize();

  ofstream outfile;
  float rfactor;
  unsigned long istep;

  cout << " printing initial populations to file: pop_ini.out" << endl;
  outfile.open("pop_ini.out", ios::out);
  printPopulations(outfile, genome, ga);
  outfile.close(); 
  
  bestRfactor =numeric_limits<float>::max(); 
  istep = 0; 
  while(!ga.done())
    {

      genome=ga.statistics().bestIndividual(); 
      rfactor = Objective(genome);

      if (rfactor < bestRfactor)
	{
	  bestRfactor = rfactor; 
	  bestWeights = probFactors;
	  bestCalcNMR = currCalcNMR;	  
	}
 

      ga.step();
      istep++; 

      if (istep%nPrints == 0 ) cout << " GA step = " << istep << " best Q factor = " << bestRfactor << endl;
      if (bestRfactor < minError) 
	{
	  cout << " Find the solution with the best error " << bestRfactor << " in GA step " << istep << endl; 
	  break; 
	}   
 
    }
  
  if (istep == gaSteps) 
    {
      cout << " Maximum GA step reached, the best error is " << bestRfactor << " in MC step " << istep << endl; 
    }
  
  cout << " Printing the final populations to file: pop_fin.out" << endl;
  outfile.open("pop_fin.out", ios::out); 
  printPopulations(outfile,genome,ga);
  outfile.close(); 
 
}

void MultiFitGA::printPopulations(ostream & outfile, GARealGenome & genome, const GASimpleGA & ga) 
{
  vector <float> sumProbs; 
 
  for(int ii=0; ii<ga.population().size(); ii++)
    {
      sumProbs.push_back(0.0);
    }

  if (randomMethod == 0)
    {
  
      for(int ii=0; ii<ga.population().size(); ii++)
	{
	  genome = ga.population().individual(ii);
	  outfile << setiosflags(ios::fixed) << setw(10) << setprecision(5) << genome.score() << " ";
	  sumProbs[ii] = 0.0;        
	}
      outfile << endl; 
      
      float probvalue; 
      for (int iconf = 0; iconf < numConfs; iconf++ )
	{
	  for(int ii=0; ii<ga.population().size(); ii++)
	    {
	      genome = ga.population().individual(ii);
	     
	      if (iconf == (numConfs - 1)) 
		probvalue = 1.0-sumProbs[ii];
	      else
		probvalue = (genome.gene(iconf))*(1.0-sumProbs[ii]);
	      sumProbs[ii] = sumProbs[ii] + probvalue;
	      outfile << setiosflags(ios::fixed) << setw(10) << setprecision(5) << probvalue << " ";
	    }
	  outfile << endl;
	}
    }
  else if (randomMethod == 1)
    {
      for(int ii=0; ii<ga.population().size(); ii++)
	{
	  genome = ga.population().individual(ii);
	  outfile << setiosflags(ios::fixed) << setw(10) << setprecision(5) << genome.score() << " ";
	  
	  for (int iconf = 0; iconf < numConfs; iconf++ )
	    {
	      sumProbs[ii] = sumProbs[ii] + genome.gene(iconf);  
	    }
	}
      outfile << endl; 

      float probvalue; 
      for (int iconf = 0; iconf < numConfs; iconf++ )
	{
	  for(int ii=0; ii<ga.population().size(); ii++)
	    {
	      genome = ga.population().individual(ii);
	      probvalue = genome.gene(iconf)/sumProbs[ii];
	      outfile << setiosflags(ios::fixed) << setw(10) << setprecision(5) << probvalue << " ";
	    }
	  outfile << endl;
	}
    }
  else 
    {
      cerr << " Wrong value of randomMethod.\n"; 
      exit(0);    
    }
}

void MultiFitGA::calcNMRFromEnsemble(GARealGenome& genome)
{
 
  for (unsigned inmr = 0; inmr < currCalcNMR.size(); inmr++)
    {currCalcNMR[inmr]= 0.0;}
  
  if (randomMethod == 0)
    {
      
      float sumProb = 0.0; 
      for(unsigned iconfs = 0; iconfs < nmrConfInp.size(); iconfs++)
	{
	  
	  if (iconfs == (nmrConfInp.size() - 1)) 
	    {probFactors[iconfs] = 1.0-sumProb;}
	  else
	    probFactors[iconfs] = (genome.gene(iconfs))*(1.0-sumProb);
	  sumProb = sumProb + probFactors[iconfs];	    
	  
	  for (unsigned inmr = 0; inmr < currCalcNMR.size(); inmr++)
	    {
	      currCalcNMR[inmr] += probFactors[iconfs]*nmrConfInp[iconfs][inmr];
	    }
	}


    }
  else if (randomMethod ==1)
    {
      float sumProb = 0.0; 
      for(unsigned iconfs = 0; iconfs < nmrConfInp.size(); iconfs++)
	sumProb = sumProb + genome.gene(iconfs);
      for(unsigned iconfs = 0; iconfs < nmrConfInp.size(); iconfs++)
	{
	  probFactors[iconfs] = genome.gene(iconfs)/sumProb;
	  
	  for (unsigned inmr = 0; inmr < currCalcNMR.size(); inmr++)
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
}

float MultiFitGA::Objective(GAGenome& g)
{
  float numerator = 0.0;
  float denominator = 0.0; 
  GARealGenome& a = (GARealGenome &)g;
  
  calcNMRFromEnsemble(a);

  for (unsigned int inmr = 0; inmr < nmrExpInp.size(); inmr++)
    {
      numerator += (currCalcNMR[inmr]-nmrExpInp[inmr])*(currCalcNMR[inmr]-nmrExpInp[inmr]);
      denominator += nmrExpInp[inmr]*nmrExpInp[inmr];
    }
  
  return sqrt(numerator/denominator);

}

float MultiFitGA::Comparator(const GAGenome& g1, const GAGenome& g2) {
  GARealGenome& a = (GARealGenome &)g1;
  GARealGenome& b = (GARealGenome &)g2;
  float sum = 0.0; 
  for (int i = 0; i < a.size(); i++ )
    sum+= (a.gene(i) - b.gene(i)) * (a.gene(i) - b.gene(i)); 
  return exp(-sum/a.size());
}


void MultiFitGA::writeToFile()
{

  ofstream nmrBestOut;

  nmrBestOut.open("nmrBest.out", ios::out);
  
  for (int inmr=0; inmr<numNMRPoints; inmr++)
    {
     
      nmrBestOut << setiosflags(ios::fixed) << setw(10) << setprecision(5) << nmrExpInp0[inmr]
		 << setiosflags(ios::fixed) << setw(10) << setprecision(5) << nmrExpInp[inmr]
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
  cout << "*        FROM AN ENSEMBLE of CONFORMATIONS USING GENETIC ALGORITHM.    *" << endl;
  cout << "************************************************************************" << endl;

  cout << "*************************Reading Parameters:****************************" << endl;
  param_file = argv[1];


  MultiFitGA myMultiFit = MultiFitGA(param_file);    // read in the parameters to start program
  cout << " Done with reading parameters." << endl;
  myMultiFit.GeneticAlgorithmFit();
  myMultiFit.writeToFile(); 

  cout << " Done with the searching of the multi-structure solution for the NMR data." << endl; 

  t_end = time(NULL);
  seconds =  t_end - t_start;
  cout << "***********************Output the CPU time information:*****************" << endl;
  cout << " The CPU time for this search is " << seconds  << " seconds." << endl;

  return(0); 
 
}
