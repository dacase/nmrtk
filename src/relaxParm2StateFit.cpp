/*********************************************************************************************  
 *   relaxParm2StateFit.cpp,                                                                 *
 *          a program to find the R1, R2 and NOE data to multi structure unique structures.  * 
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

#include "relaxParm2StateFit.h"

relaxParm relaxParm2StateFit::r1r2ParmOpt;
relaxParm relaxParm2StateFit::r1r2ParmOpt1;
relaxParm relaxParm2StateFit::r1r2ParmOpt2;

bool relaxParm2StateFit::isOptTauOvr;
bool relaxParm2StateFit::isOneTensor;
bool relaxParm2StateFit::isOptAngles;
bool relaxParm2StateFit::isOptTauInt;
bool relaxParm2StateFit::isOptAInt;
bool relaxParm2StateFit::isOptPopA;

double relaxParm2StateFit::popA; 

relaxParm2StateFit::relaxParm2StateFit(string inputFile)
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
  cout << "\n Reading the input file for relaxation parameter fitting: " << inputFile << endl;
 
  prmFile = inputFile; 

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
      else if (lineHead == "difFile:" )
	{ 
	  ssin >> difFile1; ssin >> difFile2;
	  cout << " difFile1: " << difFile1 << " difFile2: " << difFile2  << endl;
	}
      else if (lineHead == "vecFile:")
	{
	  ssin >> vecFile1; ssin >> vecFile2;
	  cout << " vecFile1: " << vecFile1 << " vecFile2: " << vecFile2 << endl;
	}
      else if (lineHead == "outFile:")
	{ 
	  ssin >> outFile1; ssin >> outFile2;
	  cout << " outFile1: " << outFile1 << " outFile2: " << outFile2 <<  endl;
	}	
      else if (lineHead == "finFile:")
	{ 
	  ssin >> finFile;
	  cout << " finFile: " << finFile <<  endl;
	}	
      else if (lineHead == "isOptTauOvr:")
	{ 
	  // bool tempb;
	  // ssin >> tempb;
	  // isOptTauOvr = tempb;
	  ssin >> isOptTauOvr; 
	  ssin >> isOneTensor; 
	  ssin >> minTauOvr; 
	  ssin >> maxTauOvr;
	  cout << " isOptTauOvr: " << isOptTauOvr << " isOneTensor: " << isOneTensor << " minTauOvr: " << minTauOvr 
	       << " maxTauOvr: " << maxTauOvr << endl;
	}	
      else if (lineHead == "isOptAngles:")
	{ 
	  ssin >> isOptAngles; 
	  ssin >> minPhi; ssin >> maxPhi;
	  ssin >> minPsi; ssin >> maxPsi;
	  cout << " isOptAngles: " << isOptAngles << " minPhi: " << minPhi << " maxPhi: " << maxPhi
	       << " minPsi: " << minPsi << " maxPsi: " << maxPsi<< endl;	
	}	
      else if (lineHead == "isOptTauInt:")
	{ 
	  ssin >> isOptTauInt; 
	  ssin >> minTauInt; 
	  ssin >> maxTauInt;
	  cout << " isOptTauInt: " << isOptTauInt << " minTauInt: " << minTauInt 
	       << " maxTauInt: " << maxTauInt << endl;
	}	
      else if (lineHead == "isOptAInt:")
	{ 
	  ssin >> isOptAInt; 
	  ssin >> minAInt; 
	  ssin >> maxAInt;
	  cout << " isOptAInt: " << isOptAInt << " minAInt: " << minAInt 
	       << " maxAInt: " << maxAInt << endl;
	}
      else if (lineHead == "isOptPopA:")
	{ 
	  ssin >> isOptPopA; ssin >> popAIni;
	  ssin >> minPopA; 
	  ssin >> maxPopA;
	  cout << " isOptPopA: " << isOptPopA << " popAIni: " << popAIni  << " minPopA: " << minPopA 
	       << " maxPopA: " << maxPopA << endl;
	}	
      else
	{continue;}
    }
  infoFile.close(); 
  r1r2ParmExp1 = relaxParm(difFile1,vecFile1,outFile1);
  r1r2ParmExp2 = relaxParm(difFile2,vecFile2,outFile2);
  r1r2ParmOpt1 = r1r2ParmExp1;
  // r1r2ParmOpt1 = relaxParm(difFile1,vecFile1,outFile1);
  r1r2ParmOpt2 = r1r2ParmExp2;
  // r1r2ParmOpt2 = relaxParm(difFile2,vecFile2,outFile2);
  r1r2ParmOpt = r1r2ParmExp1;
  // r1r2ParmOpt = relaxParm(difFile1,vecFile1,outFile1);
  r1r2ParmOpt.setOutFile(finFile);
  popA = popAIni;

  calcRelaxParmIni();
  r1r2ParmOpt.calcRFactor();
  cout << "\n The initial calculated values of relaxation parameters from two states are:  \n";
  r1r2ParmOpt.printRelaxParm(cout);
  cout << "\n With the initial R factor: " << ""  << r1r2ParmOpt.getRFactor() << " \n";
}

relaxParm2StateFit::~relaxParm2StateFit()
{}

void relaxParm2StateFit::GAOptimization()
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
  // cout << "  Finished realGenome.initialize()" << endl;

  // initialize the GA method for evolving the solution
  GASimpleGA ga(realGenome);
  // cout << "  Finished ga(realGenome)" << endl;
  initGAMethod(ga);
  // cout << "  Finished initialize GA method." << endl;

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
  calcRelaxParmFromGenome(r1r2ParmOpt,realGenome);
  r1r2ParmOpt.calcRelaxParm();
  r1r2ParmOpt.calcRFactor();
  r1r2ParmOpt1.calcRelaxParm();
  r1r2ParmOpt1.calcRFactor();
  r1r2ParmOpt2.calcRelaxParm();
  r1r2ParmOpt2.calcRFactor();

  
}

void relaxParm2StateFit::setFitRange(GARealAlleleSetArray & alleles)
{
  // fitting ranges for overall decay constants
  // e.g. minTauOvr = 0.75, maxTauOvr = 1.5
  if (isOptTauOvr)
    {
      vector <double> tauOvr = r1r2ParmExp1.getTauOvr();
      if (r1r2ParmExp1.getNumModesOvr() == 1)
	{
	  for (unsigned i = 0; i < tauOvr.size(); i++)  
	    alleles.add(minTauOvr*tauOvr[i], maxTauOvr*tauOvr[i]); 
	}
      else if (r1r2ParmExp1.getNumModesOvr() == 3)
	{
	  for (unsigned i = 0; i < (tauOvr.size()-1); i++)  
	    alleles.add(minTauOvr*tauOvr[i], maxTauOvr*tauOvr[i]); 
	}
      else if (r1r2ParmExp1.getNumModesOvr() == 5)
	{
	  for (unsigned i = 0; i < (tauOvr.size()-2); i++)  
	    alleles.add(minTauOvr*tauOvr[i], maxTauOvr*tauOvr[i]); 
	}
      else
	{
	  cerr << "the number of modes of Overall rotations is not supported. numModesOvr = 1, 3, or 5. \n";
	  exit(0);
	}
    }

  // fitting ranges for phi and chi angles between spin vectors and diffusion tensor
  // e.g. minPhi = 0.0, maxPhi = 180.0, minPsi = -180.0, maxPsi = 180
  if (isOptAngles)
    {
      if (r1r2ParmExp1.getNumModesOvr() == 1)
	{
	}   
      else if (r1r2ParmExp1.getNumModesOvr() == 3)
	{
	  for (unsigned i = 0; i < r1r2ParmExp1.getNumSpinVecs(); i++)
	    alleles.add(minPhi, maxPhi);
	}   
      else if (r1r2ParmExp1.getNumModesOvr() == 5)
	{
	  for (unsigned i = 0; i < r1r2ParmExp1.getNumSpinVecs(); i++)
	    alleles.add(minPhi, maxPhi);
	  for (unsigned i = 0; i < r1r2ParmExp1.getNumSpinVecs(); i++)
	    alleles.add(minPsi, maxPsi);
	}
      else
	{
	  cerr << "the number of modes of Overall rotations is not supported. numModesOvr =1, 3, or 5. \n";
	  exit(0);
	}
    }
  
  // fitting ranges for the time constants of internal motions
  // e.g. minTauInt = 0.1, maxTauInt = 10
  if (isOptTauInt)
    {
      vector < vector <double> > tauInt = r1r2ParmExp1.getTauInt();
      for (unsigned i = 0; i < tauInt.size(); i++)
	for (unsigned it = 1; it < tauInt[i].size(); it++)
	  alleles.add(minTauInt*tauInt[i][it],maxTauInt*tauInt[i][it]);
    }

  // fitting ranges for the time constants of internal motions
  // e.g. minAInt = 0.5, maxAInt = 1.0
  if (isOptAInt)
    {
      vector < vector <double> > AInt = r1r2ParmExp1.getAInt();
      for (unsigned i = 0; i < AInt.size(); i++)
	for (unsigned it = 0; it < (AInt[i].size()-1); it++)
	  alleles.add(minAInt*AInt[i][it],maxAInt);
    }

  // fitting ranges for overall decay constants
  // e.g. minTauOvr = 0.75, maxTauOvr = 1.5
  if (isOptTauOvr)
    {
      if ( !isOneTensor)
	{
	  vector <double> tauOvr = r1r2ParmExp2.getTauOvr();
	  if (r1r2ParmExp2.getNumModesOvr() == 1)
	    {
	      for (unsigned i = 0; i < tauOvr.size(); i++)  
		alleles.add(minTauOvr*tauOvr[i], maxTauOvr*tauOvr[i]); 
	    }
	  else if (r1r2ParmExp2.getNumModesOvr() == 3)
	    {
	      for (unsigned i = 0; i < (tauOvr.size()-1); i++)  
		alleles.add(minTauOvr*tauOvr[i], maxTauOvr*tauOvr[i]); 
	    }
	  else if (r1r2ParmExp2.getNumModesOvr() == 5)
	    {
	      for (unsigned i = 0; i < (tauOvr.size()-2); i++)  
		alleles.add(minTauOvr*tauOvr[i], maxTauOvr*tauOvr[i]); 
	    }
	  else
	    {
	      cerr << "the number of modes of Overall rotations is not supported. numModesOvr = 1, 3, or 5. \n";
	      exit(0);
	    }
	}
    }

  // fitting ranges for phi and chi angles between spin vectors and diffusion tensor
  // e.g. minPhi = 0.0, maxPhi = 180.0, minPsi = -180.0, maxPsi = 180
  if (isOptAngles)
    {
      if (r1r2ParmExp2.getNumModesOvr() == 1)
	{
	}   
      else if (r1r2ParmExp2.getNumModesOvr() == 3)
	{
	  for (unsigned i = 0; i < r1r2ParmExp2.getNumSpinVecs(); i++)
	    alleles.add(minPhi, maxPhi);
	}   
      else if (r1r2ParmExp2.getNumModesOvr() == 5)
	{
	  for (unsigned i = 0; i < r1r2ParmExp2.getNumSpinVecs(); i++)
	    alleles.add(minPhi, maxPhi);
	  for (unsigned i = 0; i < r1r2ParmExp2.getNumSpinVecs(); i++)
	    alleles.add(minPsi, maxPsi);
	}
      else
	{
	  cerr << "the number of modes of Overall rotations is not supported. numModesOvr =1, 3, or 5. \n";
	  exit(0);
	}
    }
  
  // fitting ranges for the time constants of internal motions
  // e.g. minTauInt = 0.1, maxTauInt = 10
  if (isOptTauInt)
    {
      vector < vector <double> > tauInt = r1r2ParmExp2.getTauInt();
      for (unsigned i = 0; i < tauInt.size(); i++)
	for (unsigned it = 1; it < tauInt[i].size(); it++)
	  alleles.add(minTauInt*tauInt[i][it],maxTauInt*tauInt[i][it]);
    }

  // fitting ranges for the time constants of internal motions
  // e.g. minAInt = 0.5, maxAInt = 1.0
  if (isOptAInt)
    {
      vector < vector <double> > AInt = r1r2ParmExp2.getAInt();
      for (unsigned i = 0; i < AInt.size(); i++)
	for (unsigned it = 0; it < (AInt[i].size()-1); it++)
	  alleles.add(minAInt*AInt[i][it],maxAInt);
    }

  if (isOptPopA)
    {
      alleles.add(minPopA,maxPopA);

    }
}

void relaxParm2StateFit::initGAMethod(GASimpleGA & ga)
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

void relaxParm2StateFit::initPopulations(GASimpleGA & ga) 
{
  // need to figure out how to initialize the populations as close to initial values 
  // GAGenome *genome;

}

void relaxParm2StateFit::printPopulations(ostream & outfile, GARealGenome & genome, const GASimpleGA & ga) 
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

void relaxParm2StateFit::writePopResults(ostream & outfile, GARealGenome & genome, const GASimpleGA & ga) 
{
  vector < vector < double > > R1R2Pop;
  vector < double > R1R2PopAv;
  vector < double > rfactors;

  for(int ipop=0; ipop<ga.population().size(); ipop++)
    {
      genome = ga.population().individual(ipop);
      calcRelaxParmFromGenome(r1r2ParmOpt, genome);
      vector < double > R10 = r1r2ParmOpt.getR1();
      vector < double > R20 = r1r2ParmOpt.getR2();
      vector < double > NOE0 = r1r2ParmOpt.getNOE();
      vector < double > R1R2temp; 
      for (unsigned i=0; i< R10.size(); i++)
	{
	  R1R2temp.push_back(R10[i]);
	  if (ipop == 0)
	    R1R2PopAv.push_back(R10[i]);
	  else
	    R1R2PopAv[i] += R10[i];	    
	}

      for (unsigned i=0; i< R20.size(); i++)
	{
	  R1R2temp.push_back(R20[i]); 
	  if (ipop == 0)
	    R1R2PopAv.push_back(R20[i]);
	  else
	    R1R2PopAv[i+R10.size()] += R20[i];
	}
      for (unsigned i=0; i< NOE0.size(); i++)
	{
	  R1R2temp.push_back(NOE0[i]);
	  if (ipop == 0)
	    R1R2PopAv.push_back(NOE0[i]);
	  else
	    R1R2PopAv[i+R10.size()+R20.size()] += NOE0[i];
	}
      R1R2Pop.push_back(R1R2temp);
      r1r2ParmOpt.calcRFactor();
      double rfactor0 =  r1r2ParmOpt.getRFactor();
      rfactors.push_back(rfactor0);
    }
  
  for (unsigned i = 0; i <  R1R2PopAv.size(); i++)
    R1R2PopAv[i] = R1R2PopAv[i]/rfactors.size();

  outfile << "# R1, R2, and NOE values from the final ensemble and their averages. \n";
  outfile << setiosflags(ios::fixed) << setw(12) << " average ";
  for (unsigned i = 0; i < R1R2Pop.size(); i++)
    {
      char confi[10];
      sprintf(confi,"%d",i);
      string pref = "conf_"; 
      string temps = pref + confi;
      outfile << setiosflags(ios::fixed) << setw(12) << temps;
    }
  outfile << endl;

  for (unsigned iconf = 0; iconf < R1R2PopAv.size(); iconf++)
    {
      outfile << setiosflags(ios::fixed) << setw(12) << setprecision(5) << R1R2PopAv[iconf];	
      
      for (unsigned ipop = 0; ipop < R1R2Pop.size(); ipop++)
	{
	  outfile << setiosflags(ios::fixed) << setw(12) << setprecision(5) << R1R2Pop[ipop][iconf];  
	}
      outfile << endl;
    }

  outfile << setiosflags(ios::fixed) << setw(12) << setprecision(5) << " ";	  
  for (unsigned ipop = 0; ipop < R1R2Pop.size(); ipop++)
    {
      outfile << setiosflags(ios::fixed) << setw(12) << setprecision(5) << rfactors[ipop];         
    }
  outfile << endl;
  
}


void relaxParm2StateFit::write()
{
  
  r1r2ParmOpt1.write();
  r1r2ParmOpt2.write();
  if (isOptTauOvr && isOneTensor) r1r2ParmOpt.setTauOvr(r1r2ParmOpt1.getTauOvr());
  // r1r2ParmOpt.write();
  ofstream outInfoFile;
  outInfoFile.open(r1r2ParmOpt.getOutFile().data(),ios::out); 
  outInfoFile << "\n Final result about diffusion tensor: \n";
  r1r2ParmOpt.printTensor(outInfoFile); 
  outInfoFile << "\n Final result about relaxation parameters:  \n";
  r1r2ParmOpt.printRelaxParm(outInfoFile);
  outInfoFile << "\n The R factor comparing with initial input: " << r1r2ParmOpt.getRFactor() << " \n";
}

void relaxParm2StateFit::calcRelaxParmFromGenome(relaxParm & r1r2Parm, const GARealGenome & a)
{

  unsigned int iconf = 0;

  vector <double> tauOvr1;
  if (isOptTauOvr)
    {
      tauOvr1=r1r2ParmOpt1.getTauOvr();
      if (tauOvr1.size() == 1)
	{
	  for (unsigned i = 0; i < tauOvr1.size(); i++)
	    {
	      tauOvr1[i] = a.gene(iconf);
	      iconf++;
	    }
	}
      if (tauOvr1.size() == 3)
	{
	  for (unsigned i = 0; i < (tauOvr1.size()-1); i++)
	    {
	      tauOvr1[i] = a.gene(iconf);
	      iconf++;
	    }
	}
      else if (tauOvr1.size() == 5)
	{
	  for (unsigned i = 0; i < (tauOvr1.size()-2); i++)
	    {
	      tauOvr1[i] = a.gene(iconf);
	      iconf++;
	    }
	}
      else
	{
	  cerr << "the number of modes of Overall rotations is not supported. numModesOvr = 1, 3, or 5. \n";
	  exit(0);
	}  
      r1r2ParmOpt1.setTauOvr(tauOvr1); 
      r1r2ParmOpt1.updateDiffTensorEVal(); 
    }
  
  vector <vector <double> > angles1;
  if (isOptAngles)
    {
      angles1 = r1r2ParmOpt1.getVecAngles();
      if (r1r2ParmOpt1.getNumModesOvr() == 1)
	{
	}
      else if (r1r2ParmOpt1.getNumModesOvr() == 3)
	{
	  for (unsigned i = 0; i < angles1.size(); i++)
	    {  
	      angles1[i][0] = a.gene(iconf);
	      iconf++;
	    }
	}
      else if (r1r2ParmOpt1.getNumModesOvr() == 5)
	{
	  for (unsigned i = 0; i < angles1.size(); i++)
	    {  
	      angles1[i][0] = a.gene(iconf);
	      iconf++;
	    } 
	  for (unsigned i = 0; i < angles1.size(); i++)
	    {  
	      angles1[i][1] = a.gene(iconf);
	      iconf++;
	    }  
	}  
      else
	{
	  cerr << "the number of modes of Overall rotations is not supported. numModesOvr = 1, 3, or 5. \n";
	  exit(0);
	}
      r1r2ParmOpt1.setVecAngles(angles1);
      r1r2ParmOpt1.calcCoeffsCFOAng();
    }

  vector < vector <double> > tauInt1;
  if (isOptTauInt)
    {
      tauInt1 = r1r2ParmOpt1.getTauInt();
      for (unsigned i = 0; i < tauInt1.size(); i++)
	{
	  for (unsigned it = 1; it < tauInt1[i].size(); it++)
	    {
	      tauInt1[i][it] = a.gene(iconf);
	      iconf++;
	    }
	}
      r1r2ParmOpt1.setTauInt(tauInt1);
    }

  vector < vector <double> > AInt1;  
  if (isOptAInt)
    {
      AInt1 = r1r2ParmOpt1.getAInt();
      for (unsigned i = 0; i < AInt1.size(); i++)
	{
	  float sumProbs = 0.0;
	  for (unsigned it = 0; it < AInt1[i].size(); it++)
	    {
	      if (it == (AInt1[i].size()- 1)) 
		{
		  AInt1[i][it] = 1.0-sumProbs;
		}
	      else
		{ 
		  AInt1[i][it] = a.gene(iconf)*(1.0-sumProbs);
		  sumProbs = sumProbs + AInt1[i][it];
		  iconf++;
		}
	    }  
	}
      r1r2ParmOpt1.setAInt(AInt1);
    }
  r1r2ParmOpt1.calcCoeffsCFT();
  vector < vector <double> > ATot1=r1r2ParmOpt1.getATot();
  vector < vector <double> > tauTot1=r1r2ParmOpt1.getTauTot();

  vector <double> tauOvr2;
  if (isOptTauOvr)
    {
      if (isOneTensor)
	{
	  tauOvr2=r1r2ParmOpt1.getTauOvr();
	}
      else 
	{
	  tauOvr2=r1r2ParmOpt2.getTauOvr();
	  if (tauOvr2.size() == 1)
	    {
	      for (unsigned i = 0; i < tauOvr2.size(); i++)
		{
		  tauOvr2[i] = a.gene(iconf);
		  iconf++;
		}
	    }
	  if (tauOvr2.size() == 3)
	    {
	      for (unsigned i = 0; i < (tauOvr2.size()-1); i++)
		{
		  tauOvr2[i] = a.gene(iconf);
		  iconf++;
		}
	    }
	  else if (tauOvr2.size() == 5)
	    {
	      for (unsigned i = 0; i < (tauOvr2.size()-2); i++)
		{
		  tauOvr2[i] = a.gene(iconf);
		  iconf++;
		}
	    }
	  else
	    {
	      cerr << "the number of modes of Overall rotations is not supported. numModesOvr = 1, 3, or 5. \n";
	      exit(0);
	    }
	}
      r1r2ParmOpt2.setTauOvr(tauOvr2); 
      r1r2ParmOpt2.updateDiffTensorEVal(); 
    }
  
  vector <vector <double> > angles2;
  if (isOptAngles)
    {
      angles2 = r1r2ParmOpt2.getVecAngles();
      if (r1r2ParmOpt2.getNumModesOvr() == 1)
	{
	}
      else if (r1r2ParmOpt2.getNumModesOvr() == 3)
	{
	  for (unsigned i = 0; i < angles2.size(); i++)
	    {  
	      angles2[i][0] = a.gene(iconf);
	      iconf++;
	    }
	}
      else if (r1r2ParmOpt2.getNumModesOvr() == 5)
	{
	  for (unsigned i = 0; i < angles2.size(); i++)
	    {  
	      angles2[i][0] = a.gene(iconf);
	      iconf++;
	    } 
	  for (unsigned i = 0; i < angles2.size(); i++)
	    {  
	      angles2[i][1] = a.gene(iconf);
	      iconf++;
	    }  
	}  
      else
	{
	  cerr << "the number of modes of Overall rotations is not supported. numModesOvr = 1, 3, or 5. \n";
	  exit(0);
	}
      r1r2ParmOpt2.setVecAngles(angles2);
      r1r2ParmOpt2.calcCoeffsCFOAng();
    }

  vector < vector <double> > tauInt2;
  if (isOptTauInt)
    {
      tauInt2 = r1r2ParmOpt2.getTauInt();
      for (unsigned i = 0; i < tauInt2.size(); i++)
	{
	  for (unsigned it = 1; it < tauInt2[i].size(); it++)
	    {
	      tauInt2[i][it] = a.gene(iconf);
	      iconf++;
	    }
	}
      r1r2ParmOpt2.setTauInt(tauInt2);
    }

  vector < vector <double> > AInt2;  
  if (isOptAInt)
    {
      AInt2 = r1r2ParmOpt2.getAInt();
      for (unsigned i = 0; i < AInt2.size(); i++)
	{
	  float sumProbs = 0.0;
	  for (unsigned it = 0; it < AInt2[i].size(); it++)
	    {
 	      if (it == (AInt2[i].size()- 1)) 
		{
		  AInt2[i][it] = 1.0-sumProbs;
		}
	      else
		{ 
		  AInt2[i][it] = a.gene(iconf)*(1.0-sumProbs);
		  sumProbs = sumProbs + AInt2[i][it];
		  iconf++;
		}
	    }  
	}
      r1r2ParmOpt2.setAInt(AInt2);
    }  
  r1r2ParmOpt2.calcCoeffsCFT();
  vector < vector <double> > ATot2=r1r2ParmOpt2.getATot();
  vector < vector <double> > tauTot2=r1r2ParmOpt2.getTauTot();

  double percent = popA;
  if (isOptPopA)
    {
      percent = a.gene(iconf);
      iconf++;
    }

  vector < vector <double> > ATot; 
  vector < vector <double> > tauTot; 
  for (unsigned ivec=0; ivec < ATot1.size(); ivec++)
    {
      vector <double> tempA;
      vector <double> tempT;
      for (unsigned icoef = 0; icoef < ATot1[ivec].size(); icoef++)
	{
	  tempA.push_back(percent*ATot1[ivec][icoef]);
	  tempT.push_back(tauTot1[ivec][icoef]);
	}
      for (unsigned icoef = 0; icoef < ATot2[ivec].size(); icoef++)
	{
	  tempA.push_back((1.0-percent)*ATot2[ivec][icoef]);
	  tempT.push_back(tauTot2[ivec][icoef]);
	}
      ATot.push_back(tempA); 
      tauTot.push_back(tempT); 
    }

  r1r2Parm.setATot(ATot);
  r1r2Parm.setTauTot(tauTot);
  r1r2Parm.calcRelaxParm();

}

void relaxParm2StateFit::calcRelaxParmIni()
{
 
  vector < vector <double> > ATot1=r1r2ParmOpt1.getATot();
  vector < vector <double> > tauTot1=r1r2ParmOpt1.getTauTot();

  vector < vector <double> > ATot2=r1r2ParmOpt2.getATot();
  vector < vector <double> > tauTot2=r1r2ParmOpt2.getTauTot();

  vector < vector <double> > ATot; 
  vector < vector <double> > tauTot; 
  for (unsigned ivec=0; ivec < ATot1.size(); ivec++)
    {
      vector <double> tempA;
      vector <double> tempT;
      for (unsigned icoef = 0; icoef < ATot1[ivec].size(); icoef++)
	{
	  tempA.push_back(popA*ATot1[ivec][icoef]);
	  tempT.push_back(tauTot1[ivec][icoef]);
	}
      for (unsigned icoef = 0; icoef < ATot2[ivec].size(); icoef++)
	{
	  tempA.push_back((1.0-popA)*ATot2[ivec][icoef]);
	  tempT.push_back(tauTot2[ivec][icoef]);
	}
      ATot.push_back(tempA); 
      tauTot.push_back(tempT); 
    }

  r1r2ParmOpt.setATot(ATot);
  r1r2ParmOpt.setTauTot(tauTot);
  r1r2ParmOpt.calcRelaxParm();

}


float relaxParm2StateFit::Objective(GAGenome& g)
{
 
  GARealGenome& a = (GARealGenome &)g;

  calcRelaxParmFromGenome(r1r2ParmOpt, a);
  
  r1r2ParmOpt.calcRFactor();
    
  return r1r2ParmOpt.getRFactor();

}

float relaxParm2StateFit::Comparator(const GAGenome& g1, const GAGenome& g2) {
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
      cerr << " It should be 'relaxParm2StateFit ini_file'" << endl;
      return 0;
    }

  string ini_file; 
  ini_file = argv[1];


  cout << "************************************************************************" << endl;
  cout << "*       START RELAXATION PARAMETER FITTING USING GENETIC ALGORITHM.    *" << endl;
  cout << "************************************************************************" << endl;
  relaxParm2StateFit R1R2Fit = relaxParm2StateFit(ini_file); 
  cout << "Finished reading initial file." << endl;
 
  R1R2Fit.GAOptimization();

  R1R2Fit.write(); 
  
  t_end = time(NULL);
  
  seconds = t_end - t_start;
  cout << "Finished the job in " << seconds << " seconds." << endl;
}
 





 
 
 








