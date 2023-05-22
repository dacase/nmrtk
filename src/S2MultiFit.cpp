/*********************************************************************************************  
 *   S2MultiFit.cpp,                                                                         *
 *          a program to fit the total S2 using two-stat model.                              * 
 *                                                                                           *
 *                                                                                           *  
 *   Authors:                                                                                *
 *          originally written by Junchao Xia                                                *
 *                                                                                           *
 *                                                                                           *
 *                                                                                           *
 *********************************************************************************************/
#include <cmath>
#include <cstdlib>
#include <ctime>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <string>

#include <vector>
#include <limits>

#include "S2MultiFit.h"

vector <double> S2MultiFit::s2Exp;
vector <double> S2MultiFit::s2Tot;
vector <double> S2MultiFit::s2Semi;
vector <double> S2MultiFit::s2Close;
vector <double> S2MultiFit::p2Int;

bool S2MultiFit::isOptS2Ind;
bool S2MultiFit::isOptS2EQ;


S2MultiFit::S2MultiFit(string inputFile)
{ 
 
  /* 
    Reading the initial configuration and control parameters 
  */   
  
  ifstream infoFile;
  infoFile.open(inputFile.data(),ios::in); 
  if(infoFile.fail())
    { 
      cerr<<" Error opening the input file to initialize the relaxation parameter fitting. \n";  
      cerr<<" Please check the file: " << inputFile << endl; 
      exit(0);
    }
  cout << "\n Reading the input file for S2 fitting: " << inputFile << endl;
 
  string linein;
  string lineHead; 
  istringstream ssin;
  while( getline(infoFile, linein))
    {
      ssin.clear();
      ssin.str(linein);
      ssin >> lineHead; 
      if (lineHead == "#")
	{continue;}
      else if (lineHead == "gaSteps:")
	{
	  ssin >> gaSteps;
	  cout << " gaSteps: " << gaSteps << endl;
	}
      else if (lineHead == "nPrints:" )
	{ 
	  ssin >> nPrints;
	  cout << " nPrints: " << nPrints << endl;
	}
      else if (lineHead == "gaSeed:" )
	{ 
	  ssin >> gaSeed;
	  cout << " gaSeed: " << gaSeed << endl;
	}
      else if (lineHead == "minError:" )
	{
	  ssin >> minError;
	  cout << " minError: " << minError << endl;
	}  
      else if (lineHead == "iniFile:" )
	{ 
	  ssin >> iniFile;
	  cout << " iniFile: " << iniFile  << endl;
	}
      else if (lineHead == "outFile:")
	{ 
	  ssin >> outFile;
	  cout << " outFile: " << outFile <<  endl;
	} 
      else if (lineHead == "isOptS2EQ:")
	{ 
	  ssin >> isOptS2EQ;
	  cout << " isOptS2EQ: " << isOptS2EQ <<  endl;
	}

     else if (lineHead == "isOptS2Ind:")
	{ 
	  ssin >> isOptS2Ind;
	  ssin >> isRandomS2;
	  ssin >> minS2; 
	  ssin >> maxS2;
	  cout << " isOptS2Ind: " << isOptS2Ind << " isRandomS2: " << isRandomS2 << " minS2: " << minS2 
	       << " maxS2: " << maxS2 << endl;
	}	
      else
	{continue;}
    }
  infoFile.close(); 

  infoFile.open(iniFile.data(),ios::in); 
  if(infoFile.fail())
    { 
      cerr<<" Error opening the ini file for different states. \n";  
      cerr<<" Please check the file: " << iniFile << endl; 
      exit(0);
    }
  
  while( getline(infoFile, linein))
    {
      ssin.clear();
      ssin.str(linein);
      ssin >> lineHead; 
      if (lineHead == "#")
	{continue;}
      else if (lineHead == "S2:")
	{
	  int resin; 
	  double s2tot,s2semi,s2close, p2int;
	  ssin >> resin; ssin >> s2tot; ssin >> s2semi; ssin >> s2close;  
	  resID.push_back(resin);
	  s2Tot.push_back(s2tot);
	  s2Semi.push_back(s2semi);
	  s2Close.push_back(s2close);
	  if (isOptS2EQ) 
	    {
	      ssin >> p2int;
	      p2Int.push_back(p2int);
	    }
	}
      else
	{continue;}
    }
  infoFile.close();
 
  s2Exp = s2Tot; 
  s2SemiIni = s2Semi;
  s2CloseIni = s2Close;
}

S2MultiFit::~S2MultiFit()
{}

void S2MultiFit::GAOptimization()
{

  // initalize Genetic Algorithm
  // set the seed for random number
  GARandomSeed(gaSeed);

  // set the search ranges for fitting parameters
  GARealAlleleSetArray alleles;                        
  setFitRange(alleles);

  cout << " Finished setting the fitting ranges for all parameters." << endl;
  cout << " The total fitting parameters are " << alleles.size() << endl; 

  // initialize the genome to represent a set of parameters
  GARealGenome realGenome(alleles,Objective);
  realGenome.initialize(); 
  cout << "  Finished realGenome.initialize()" << endl;

  // initialize the GA method for evolving the solution
  GASimpleGA ga(realGenome);
  cout << "  Finished ga(realGenome)" << endl;
  initGAMethod(ga);
  cout << "  Finished initGAMethod(ga)" << endl;
  //initGAMethod(ga);
  //cout << "  Finished initGAMethod(ga)" << endl;

  // initalize the populations 
  // initPopulations(ga); 

  // save the initial populations
  cout << " printing initial populations to file: pop_ini.out" << endl;
  ofstream outfile;
  outfile.open("pop_ini.out", ios::out);
  printPopulations(outfile, realGenome, ga);
  outfile.close(); 
  
  float bestRfactor =numeric_limits<float>::max(); 
  unsigned long istep = 0; 
  while(!ga.done())
    { 
      realGenome = ga.statistics().bestIndividual();
      bestRfactor = Objective(realGenome);
      ga.step();
      istep++; 
      
      if (istep%nPrints == 0 ) cout << " GA step = " << istep << " best R factor = " << bestRfactor << endl;
      
      if (bestRfactor < minError) 
	{
	  cout << " Find the solution with the best R factor = " << bestRfactor << " in GA step " << istep << endl; 
	  break; 
	}      
      
    }

  if (istep == gaSteps) 
    {
      cout << " Maximum GA step reached, the best R factor is " << bestRfactor << " in GA step " << istep << endl; 
    }

  cout << " Printing the final populations to file: pop_fin.out" << endl;
  outfile.open("pop_fin.out", ios::out); 
  printPopulations(outfile,realGenome,ga);
  outfile.close(); 
  cout << " Printing the relaxation parameters from final ensemble to file: relax_pop.out" << endl;
  outfile.open("relax_pop.out", ios::out); 
  writePopResults(outfile,realGenome,ga);
  outfile.close(); 

  realGenome = ga.statistics().bestIndividual();
  calcS2FromGenome(realGenome);
  
}

void S2MultiFit::setFitRange(GARealAlleleSetArray & alleles)
{
 
  if (isOptS2Ind)
    {
      if (isRandomS2)
	{
	  for (unsigned i = 0; i < s2SemiIni.size(); i++)
	    alleles.add(0,1.0);
	  for (unsigned i = 0; i < s2CloseIni.size(); i++)
	    alleles.add(0,1.0);
	}
      else
	{
	  for (unsigned i = 0; i < s2SemiIni.size(); i++)
	    alleles.add(minS2*s2SemiIni[i],maxS2);
	  for (unsigned i = 0; i < s2CloseIni.size(); i++)
	    alleles.add(minS2*s2CloseIni[i],maxS2);
	}
    }
  
  alleles.add(0.0,1.0);
  alleles.add(0.0,1.0);
}

void S2MultiFit::initGAMethod(GASimpleGA & ga)
{
  GASharing scale(Comparator);
  ga.minimize();                      // by default we want to minimize the objective
  ga.scaling(scale);                  // set the scaling method to our sharing
  ga.populationSize(25);              // how many individuals in the population
  ga.nGenerations(gaSteps);           // number of generations to evolve
  ga.pMutation(0.001);                // likelihood of mutating new offspring
  ga.pCrossover(0.9);                 // likelihood of crossing over parents
  ga.scoreFilename("score.out");      // name of file for scores
  ga.scoreFrequency(10);              // keep the scores of every generation
  ga.flushFrequency(50);              // specify how often to write the score to disk
  ga.selectScores(GAStatistics::AllScores);
  //ga.parameters(argc, argv, gaTrue);  // parse commands, complain if bogus args
  ga.initialize();
}


void S2MultiFit::initPopulations(GASimpleGA & ga) 
{
  // need to figure out how to initialize the populations as close to initial values 
  // GAGenome *genome;


}

void S2MultiFit::printPopulations(ostream & outfile, GARealGenome & genome, const GASimpleGA & ga) 
{
  for(int ipop=0; ipop<ga.population().size(); ipop++)
    {
      genome = ga.population().individual(ipop);
      outfile << setiosflags(ios::scientific) << setw(12) << setprecision(5) << genome.score() << " ";       
    }
  outfile << endl; 
  
  float tempvalue; 
  for (int iconf = 0; iconf < genome.size(); iconf++ )
    {
      for(int ipop=0; ipop<ga.population().size(); ipop++)
	{
	  genome = ga.population().individual(ipop);
	  tempvalue = genome.gene(iconf);
	  outfile << setiosflags(ios::scientific) << setw(12) << setprecision(5) << tempvalue << " ";
	}
      outfile << endl;
    }
}

void S2MultiFit::writePopResults(ostream & outfile, GARealGenome & genome, const GASimpleGA & ga) 
{
  vector < vector < double > > S2Pop;
  vector < double > S2PopAv;
  vector < double > rfactors;

  for(int ipop=0; ipop<ga.population().size(); ipop++)
    {
      genome = ga.population().individual(ipop);
      calcS2FromGenome(genome);
      vector < double > S20 = s2Tot;
      vector < double > S2temp; 
      for (unsigned i=0; i< S20.size(); i++)
	{
	  S2temp.push_back(S20[i]);
	  if (ipop == 0)
	    S2PopAv.push_back(S20[i]);
	  else
	    S2PopAv[i] += S20[i];	    
	}

      S2Pop.push_back(S2temp);
      double rfactor0 = Objective(genome);
      rfactors.push_back(rfactor0);
    }
  
  for (unsigned i = 0; i <  S2PopAv.size(); i++)
    S2PopAv[i] = S2PopAv[i]/rfactors.size();

  outfile << "# S2 values from the final ensemble and their averages. \n";
  outfile << setiosflags(ios::fixed) << setw(12) << " average ";
  for (unsigned i = 0; i < S2Pop.size(); i++)
    {
      char confi[10];
      sprintf(confi,"%d",i);
      string pref = "conf_"; 
      string temps = pref + confi;
      outfile << setiosflags(ios::fixed) << setw(12) << temps;
    }
  outfile << endl;

  for (unsigned iconf = 0; iconf < S2PopAv.size(); iconf++)
    {
      outfile << setiosflags(ios::fixed) << setw(12) << setprecision(5) << S2PopAv[iconf];	
      
      for (unsigned ipop = 0; ipop < S2Pop.size(); ipop++)
	{
	  outfile << setiosflags(ios::fixed) << setw(12) << setprecision(5) << S2Pop[ipop][iconf];  
	}
      outfile << endl;
    }

  outfile << setiosflags(ios::fixed) << setw(12) << setprecision(5) << " ";	  
  for (unsigned ipop = 0; ipop < S2Pop.size(); ipop++)
    {
      outfile << setiosflags(ios::fixed) << setw(12) << setprecision(5) << rfactors[ipop];         
    }
  outfile << endl;
  
}


void S2MultiFit::write()
{

  ofstream outfile;
  outfile.open(outFile.data(), ios::out);

  for (unsigned i = 0; i < s2Tot.size(); i++)
    {
      outfile << setiosflags(ios::fixed) << setw(12) << resID[i]
	      << setiosflags(ios::fixed) << setw(12) << setprecision(5) << s2Exp[i]
              << setiosflags(ios::fixed) << setw(12) << setprecision(5) << s2Tot[i]
	      << setiosflags(ios::fixed) << setw(12) << setprecision(5) << s2Semi[i]
	      << setiosflags(ios::fixed) << setw(12) << setprecision(5) << s2Close[i];
      outfile << endl;
    }
 
  outfile.close(); 

}

void S2MultiFit::calcS2FromGenome(const GARealGenome & a)
{

  unsigned int iconf = 0;
  double p1,p2;
 
  if (isOptS2Ind)
    {
      for (unsigned i = 0; i < s2Semi.size(); i++)
	 {
	   s2Semi[i] = a.gene(iconf);
	   iconf++;
	 }
      for (unsigned i = 0; i < s2Close.size(); i++)
	 {
	   s2Close[i] = a.gene(iconf);
	   iconf++;
	 }
    }
  
  p1 = a.gene(iconf);
  iconf++;
  // p2 = 1-p1;
  p2 = a.gene(iconf);
  iconf++; 

  if (isOptS2EQ)
    {

      for (unsigned i=0; i < s2Tot.size(); i++)
	{
	  s2Tot[i] = p1*p1*s2Semi[i] + p2*p2*s2Close[i] + 2*p1*sqrt(s2Semi[i])*p2*sqrt(s2Close[i])*p2Int[i];
	}
    }
  else 
    {
      
      for (unsigned i=0; i < s2Tot.size(); i++)
	{
	  s2Tot[i] = p1*s2Semi[i] + p2*s2Close[i];
	}
    }
}

float S2MultiFit::Objective(GAGenome& g)
{
 
 
  GARealGenome& a = (GARealGenome &)g;

  calcS2FromGenome(a);

  float rfactor = 0.0;
 
  for (unsigned i=0; i < s2Tot.size(); i++)
    {
      float tmpv = s2Tot[i]/s2Exp[i] - 1.0;
      rfactor =  rfactor + tmpv*tmpv; 
    }
  rfactor = rfactor/s2Tot.size();
  rfactor = sqrt(rfactor);
    
  return rfactor;

}

float S2MultiFit::Comparator(const GAGenome& g1, const GAGenome& g2) {
  GARealGenome& a = (GARealGenome &)g1;
  GARealGenome& b = (GARealGenome &)g2;
  float sum = 0.0; 
  unsigned numConfs = a.size(); 
  for (unsigned int i = 0; i < numConfs; i++ )
    sum+= (a.gene(i) - b.gene(i)) * (a.gene(i) - b.gene(i)); 
  return exp(-sum/numConfs);
}


// argv[1] = ini_file
int main(int argc, char **argv)
{

  time_t seconds, t_start, t_end;

  t_start = time(NULL);
  
  if(argc != 2)
    {
      cerr << " Check the command syntax!" << endl;
      cerr << " It should be 'S2MultiFit ini_file'" << endl;
      return 0;
    }

  string ini_file; 
  ini_file = argv[1];


  cout << "************************************************************************" << endl;
  cout << "*          START S2 MULTIPLE FITTING USING GENETIC ALGORITHM.          *" << endl;
  cout << "************************************************************************" << endl;
  S2MultiFit S2Fit = S2MultiFit(ini_file); 
  cout << "Finished reading initial file." << endl;
 
  S2Fit.GAOptimization();

  S2Fit.write(); 
  
  t_end = time(NULL);
  
  seconds = t_end - t_start;
  cout << "Finished the job in " << seconds << " seconds." << endl;
}
 





 
 
 








