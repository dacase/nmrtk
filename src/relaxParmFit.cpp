/*********************************************************************************************  
 *   relaxParmFit.cpp,                                                                       *
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

#include "relaxParmFit.h"

relaxParm relaxParmFit::r1r2ParmOpt;
bool relaxParmFit::isOptTauOvr;
bool relaxParmFit::isOptAngles;
bool relaxParmFit::isOptTauInt;
bool relaxParmFit::isOptAInt;

relaxParmFit::relaxParmFit(string inputFile)
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
	  ssin >> difFile;
	  cout << " difFile: " << difFile << endl;
	}
      else if (lineHead == "vecFile:")
	{
	  ssin >> vecFile;
	  cout << " vecFile: " << vecFile << endl;
	}
      else if (lineHead == "outFile:")
	{ 
	  ssin >> outFile;
	  cout << " outFile: " << outFile << endl;
	}	
      else if (lineHead == "isOptTauOvr:")
	{ 
	  // bool tempb;
	  // ssin >> tempb;
	  // isOptTauOvr = tempb;
	  ssin >> isOptTauOvr; 
	  ssin >> minTauOvr; 
	  ssin >> maxTauOvr;
	  cout << " isOptTauOvr: " << isOptTauOvr << " minTauOvr: " << minTauOvr 
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
      else
	{continue;}
    }
  infoFile.close(); 
  relaxParmFit::r1r2ParmOpt = relaxParm(difFile,vecFile,outFile);
  r1r2ParmExp = r1r2ParmOpt;
  // r1r2ParmExp = relaxParm(difFile,vecFile,outFile);

}

relaxParmFit::~relaxParmFit()
{}

void relaxParmFit::GAOptimization()
{

  // initalize Genetic Algorithm
  // set the seed for random number
  GARandomSeed(gaSeed);

  // set the search ranges for fitting parameters
  GARealAlleleSetArray alleles;
  setFitRange(alleles);

  // initialize the genome to represent a set of parameters
  GARealGenome realGenome(alleles,Objective);
  realGenome.initialize(); 

  // initialize the GA method for evolving the solution
  GASimpleGA ga(realGenome);
  initGAMethod(ga);
  
  // initalize the populations 
  initPopulations(ga); 

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
  r1r2ParmOpt.calcRFactor();

  
}

void relaxParmFit::setFitRange(GARealAlleleSetArray & alleles)
{
  // fitting ranges for overall decay constants
  // e.g. minTauOvr = 0.75, maxTauOvr = 1.5
  if (isOptTauOvr)
    {
      vector <double> tauOvr = r1r2ParmExp.getTauOvr();
      if (r1r2ParmExp.getNumModesOvr() == 1)
	{
	  for (unsigned i = 0; i < tauOvr.size(); i++)  
	    alleles.add(minTauOvr*tauOvr[i], maxTauOvr*tauOvr[i]); 
	}
      else if (r1r2ParmExp.getNumModesOvr() == 3)
	{
	  for (unsigned i = 0; i < (tauOvr.size()-1); i++)  
	    alleles.add(minTauOvr*tauOvr[i], maxTauOvr*tauOvr[i]); 
	}
      else if (r1r2ParmExp.getNumModesOvr() == 5)
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
      if (r1r2ParmExp.getNumModesOvr() == 1)
	{
	}   
      else if (r1r2ParmExp.getNumModesOvr() == 3)
	{
	  for (unsigned i = 0; i < r1r2ParmExp.getNumSpinVecs(); i++)
	    alleles.add(minPhi, maxPhi);
	}   
      else if (r1r2ParmExp.getNumModesOvr() == 5)
	{
	  for (unsigned i = 0; i < r1r2ParmExp.getNumSpinVecs(); i++)
	    alleles.add(minPhi, maxPhi);
	  for (unsigned i = 0; i < r1r2ParmExp.getNumSpinVecs(); i++)
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
      vector < vector <double> > tauInt = r1r2ParmExp.getTauInt();
      for (unsigned i = 0; i < tauInt.size(); i++)
	for (unsigned it = 1; it < tauInt[i].size(); it++)
	  alleles.add(minTauInt*tauInt[i][it],maxTauInt*tauInt[i][it]);
    }

  // fitting ranges for the time constants of internal motions
  // e.g. minAInt = 0.5, maxAInt = 1.0
  if (isOptAInt)
    {
      vector < vector <double> > AInt = r1r2ParmExp.getAInt();
      for (unsigned i = 0; i < AInt.size(); i++)
	for (unsigned it = 0; it < (AInt[i].size()-1); it++)
	  alleles.add(minAInt*AInt[i][it],maxAInt);
    }
}

void relaxParmFit::initGAMethod(GASimpleGA & ga)
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

void relaxParmFit::initPopulations(GASimpleGA & ga) 
{
  // need to figure out how to initialize the populations as close to initial values 
  // GAGenome *genome;

  /* 
  for(int ipop=0; ipop<ga.population().size(); ipop++)
    {
      // genome =  ga.population().individual(ipop);
      
      unsigned int iconf = 0;
      vector <double> tauOvr = r1r2ParmOpt.getTauOvr();
      for (unsigned i = 0; i < tauOvr.size(); i++)
	{
	  // genome.gene(iconf) = tauOvr[i];
	  ga.population().individual(ipop).gene(iconf,tauOvr[i]);
	  iconf++;
	}

      vector <vector <double> > angles = r1r2ParmOpt.getVecAngles();      
      for (unsigned i = 0; i < angles.size(); i++)
	{  
	  // genome.gene(iconf)=angles[i][0];
	  ga.population().individual(ipop).gene(iconf,angles[i][0]);
	  iconf++;
	}
      
      vector < vector <double> > tauInt = r1r2ParmOpt.getTauInt();
      for (unsigned i = 0; i < tauInt.size(); i++)
	{
	  for (unsigned it = 1; it < tauInt[i].size(); it++)
	    {
	      // genome.gene(iconf) = tauInt[i][it];
	      ga.population().individual(ipop).gene(iconf,tauInt[i][it]);
	      iconf++;
	    }
	}
      
      vector < vector <double> > AInt = r1r2ParmOpt.getAInt();
      for (unsigned i = 0; i < AInt.size(); i++)
	{
	  float sumProbs; 
	  for (unsigned it = 0; it < (AInt[i].size()-1); it++)
	    {
	      // genome.gene(iconf) = AInt[i][it]*(1.0-sumProbs);
	      ga.population().individual(ipop).gene(iconf, AInt[i][it]*(1.0-sumProbs));
	      sumProbs = sumProbs + AInt[i][it];
	      iconf++;
	      
	    }  
	}
      
    }
 
  */
}

void relaxParmFit::printPopulations(ostream & outfile, GARealGenome & genome, const GASimpleGA & ga) 
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

void relaxParmFit::writePopResults(ostream & outfile, GARealGenome & genome, const GASimpleGA & ga) 
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


void relaxParmFit::write()
{
  r1r2ParmOpt.write();
}

void relaxParmFit::calcRelaxParmFromGenome(relaxParm & r1r2Parm, const GARealGenome & a)
{

  unsigned int iconf = 0;
  if (isOptTauOvr)
    {
      vector <double> tauOvr=r1r2Parm.getTauOvr();
      if (r1r2Parm.getNumModesOvr() == 1)
	{
	  for (unsigned i = 0; i < tauOvr.size(); i++)
	    {
	      tauOvr[i] = a.gene(iconf);
	      iconf++;
	    }
	}
      if (r1r2Parm.getNumModesOvr() == 3)
	{
	  for (unsigned i = 0; i < (tauOvr.size()-1); i++)
	    {
	      tauOvr[i] = a.gene(iconf);
	      iconf++;
	    }
	}
      else if (r1r2Parm.getNumModesOvr() == 5)
	{
	  for (unsigned i = 0; i < (tauOvr.size()-2); i++)
	    {
	      tauOvr[i] = a.gene(iconf);
	      iconf++;
	    }
	}
      else
	{
	  cerr << "the number of modes of Overall rotations is not supported. numModesOvr = 1, 3, or 5. \n";
	  exit(0);
	}  
      r1r2Parm.setTauOvr(tauOvr); 
      r1r2Parm.updateDiffTensorEVal(); 
    }

  if (isOptAngles)
    {
      vector <vector <double> > angles = r1r2Parm.getVecAngles();
      if (r1r2Parm.getNumModesOvr() == 1)
	{
	}
      else if (r1r2Parm.getNumModesOvr() == 3)
	{
	  for (unsigned i = 0; i < angles.size(); i++)
	    {  
	      angles[i][0] = a.gene(iconf);
	      iconf++;
	    }
	}
      else if (r1r2Parm.getNumModesOvr() == 5)
	{
	  for (unsigned i = 0; i < angles.size(); i++)
	    {  
	      angles[i][0] = a.gene(iconf);
	      iconf++;
	    } 
	  for (unsigned i = 0; i < angles.size(); i++)
	    {  
	      angles[i][1] = a.gene(iconf);
	      iconf++;
	    }  
	}  
      else
	{
	  cerr << "the number of modes of Overall rotations is not supported. numModesOvr = 1, 3, or 5. \n";
	  exit(0);
	}
      r1r2Parm.setVecAngles(angles);
    }

  if (isOptTauInt)
    {
      vector < vector <double> > tauInt = r1r2Parm.getTauInt();
      for (unsigned i = 0; i < tauInt.size(); i++)
	{
	  for (unsigned it = 1; it < tauInt[i].size(); it++)
	    {
	      tauInt[i][it] = a.gene(iconf);
	      iconf++;
	    }
	}
      r1r2Parm.setTauInt(tauInt);
    }
     
  if (isOptAInt)
    {
      vector < vector <double> > AInt = r1r2Parm.getAInt();
      for (unsigned i = 0; i < AInt.size(); i++)
	{
	  float sumProbs = 0.0;
	  for (unsigned it = 0; it < AInt[i].size(); it++)
	    {
	      if (it == (AInt[i].size()- 1)) 
		{
		  AInt[i][it] = 1.0-sumProbs;
		}
	      else
		{ 
		  AInt[i][it] = a.gene(iconf)*(1.0-sumProbs);
		  sumProbs = sumProbs + AInt[i][it];
		  iconf++;
		}
	    }  
	}
      r1r2Parm.setAInt(AInt);
    }
  
  r1r2Parm.calcCoeffsCFOAng();
  r1r2Parm.calcCoeffsCFT();
  r1r2Parm.calcRelaxParm();

}

float relaxParmFit::Objective(GAGenome& g)
{
 
  GARealGenome& a = (GARealGenome &)g;

  calcRelaxParmFromGenome(r1r2ParmOpt, a);
  
  r1r2ParmOpt.calcRFactor();
    
  return r1r2ParmOpt.getRFactor();

}

float relaxParmFit::Comparator(const GAGenome& g1, const GAGenome& g2) {
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
      cerr << " It should be 'relaxParmFit ini_file'" << endl;
      return 0;
    }

  string ini_file; 
  ini_file = argv[1];


  cout << "************************************************************************" << endl;
  cout << "*       START RELAXATION PARAMETER FITTING USING GENETIC ALGORITHM.    *" << endl;
  cout << "************************************************************************" << endl;
  relaxParmFit R1R2Fit = relaxParmFit(ini_file); 
  cout << "Finished reading initial file." << endl;
 
  R1R2Fit.GAOptimization();

  R1R2Fit.write(); 
  
  t_end = time(NULL);
  
  seconds = t_end - t_start;
  cout << "Finished the job in " << seconds << " seconds." << endl;
}
 





 
 
 








