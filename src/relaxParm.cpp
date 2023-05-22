#include "relaxParm.h"
 
using namespace nmrspace; 
relaxParm::relaxParm()
{}

relaxParm::relaxParm(string iniFile0, string  vecFile0, string outFile0):
  nmr(iniFile0,outFile0),vecFile(vecFile0)
{
  // setIniFile(iniFile0); 
  // setOutFile(outFile0); 
  // vectFile=vecFile0;
  
  read();
}
 
relaxParm::~relaxParm()  
{

}

void relaxParm::read()  
{
  /*
    Reading the initial information of diffusion tensor, relaxation parameters
    and internal motions
  */  
  ifstream diffInfoFile;
  diffInfoFile.open(getIniFile().data(),ios::in); 
  if(diffInfoFile.fail())
    { 
      cerr<<" Error opening the input file for diffusion tensor. \n";  
      cerr<<" Please check the file: " << getIniFile().data() << endl; 
      exit(0);
    }
  cout << "\n Reading the initial information of difusion tensor and etc from file: " << getIniFile() << endl;
 
  string linein;
  string lineHead;
  istringstream ssin;
  vector < vector < double > > expData; 

  while( getline(diffInfoFile, linein))
    {
      ssin.clear();
      ssin.str(linein);
      ssin >> lineHead; 
      if (lineHead == "#")
	{continue;}
      else if (lineHead == "RotDiffTensor:")
	{
	  for (int i = 0; i < 9; i++)
	    {
	      ssin >> diffTensor[i/3][i%3];
	    }
	}
      else if (lineHead == "Eigenvalues:" )
	{
	  for (int i = 0; i < 3; i++)
	    {
	      ssin >> diffTensorEVal[i];
	    }
	}
      else if (lineHead == "Eigenvectors:" )
	{
	  for (int i = 0; i < 9; i++)
	    {
	      ssin >> diffTensorEVec[i/3][i%3];
	    }
	}
      else if (lineHead == "RelaxTimeOvr:" )
	{
	  ssin >> numModesOvr;
	  for (unsigned i = 0; i < numModesOvr; i++)
	    {
	      double tau0;
	      ssin >> tau0; tauOvr.push_back(tau0);
	    }
	}  
      else if (lineHead == "BField:" )
	{
	  ssin >> nBFields;
	  double bfield0; 
	  for (unsigned ibf = 0; ibf < nBFields; ibf++)
	    {
	      ssin >> bfield0; 
	      BFields.push_back(bfield0);
	    }
	}
      else if (lineHead == "RelaxTime:")
	{
	  unsigned resI0, resJ0; string atmTypeI0, atmTypeJ0;
	  ssin >> resI0; ssin >> atmTypeI0; 
          ssin >> resJ0; ssin >> atmTypeJ0; 
	  resI.push_back(resI0); atmTypeI.push_back(atmTypeI0);
	  resJ.push_back(resJ0); atmTypeJ.push_back(atmTypeJ0);

	  // double  R10, R20, NOE0; 
	  // ssin >> R10;  ssin >> R20; ssin >> NOE0;
	  // R1In.push_back(1/R10); R2In.push_back(1/R20); NOEIn.push_back(NOE0);
	  double temp0; 
	  vector < double > tempvec; 
	  for (unsigned i = 0; i < 3*nBFields; i++ )
	    {
	      ssin >> temp0; 
	      tempvec.push_back(1/temp0);
	    }
	  expData.push_back(tempvec);
	}
      else if (lineHead == "RelaxRate:")
	{
	  unsigned resI0, resJ0; string atmTypeI0, atmTypeJ0;
	  ssin >> resI0; ssin >> atmTypeI0; 
          ssin >> resJ0; ssin >> atmTypeJ0; 
	  resI.push_back(resI0); atmTypeI.push_back(atmTypeI0);
	  resJ.push_back(resJ0); atmTypeJ.push_back(atmTypeJ0);
	  // double R10, R20, NOE0; 
	  // ssin >> R10;  ssin >> R20; ssin >> NOE0;
	  // R1In.push_back(R10); R2In.push_back(R20); NOEIn.push_back(NOE0);
	  double temp0; 
	  vector < double > tempvec; 
	  for (unsigned i = 0; i < 3*nBFields; i++ )
	    {
	      ssin >> temp0; 
	      tempvec.push_back(temp0);
	    }
	  expData.push_back(tempvec);
	}
      else if (lineHead == "RelaxInt:")
	{
	  unsigned numModesInt0; 
	  ssin >> numModesInt0;
	  numModesInt.push_back(numModesInt0);
	  vector <double> AIntI, tauIntI;
	  for (unsigned i = 0; i <= numModesInt0; i++)
	    {
	      if (i == 0)
		{
		  double AInt0;
		  ssin >> AInt0;
		  AIntI.push_back(AInt0);
		  tauIntI.push_back(0.0);
		}
	      else 
		{
		  double AInt0, tauInt0;
		  ssin >> AInt0; ssin >> tauInt0;
		  AIntI.push_back(AInt0);
		  tauIntI.push_back(tauInt0);
		}
	    }
	  AInt.push_back(AIntI);
	  tauInt.push_back(tauIntI);
	}	
      else
	{continue;}
    }
  diffInfoFile.close(); 

  // numSpinVecs = R1In.size();
  numSpinVecs = expData.size();
  for (unsigned ifl = 0; ifl < nBFields; ifl++)
    {
      for (unsigned iv =0; iv < numSpinVecs; iv++)
	{
	  R1In.push_back(expData[iv][ifl*3]);
	  R2In.push_back(expData[iv][ifl*3+1]);
	  NOEIn.push_back(expData[iv][ifl*3+2]);
	}
    }
  R1 = R1In; 
  R2 = R2In; 
  NOE = NOEIn; 

 
  printTensor(cout);
  cout << "\n The input values of relaxation parameters:  \n";
  printRelaxParm(cout);
  cout << "\n The initial coefficients and time constants of internal motions: \n";
  printRelaxInt(cout); 

  assignSpinTypes();

  // calculate the remaining quantities of diffusion tensor
  if (numModesOvr == 1)
    {
      D = (diffTensorEVal[0]+diffTensorEVal[1]+diffTensorEVal[2])/3.0;
      L2 = (diffTensorEVal[2]*diffTensorEVal[1] + diffTensorEVal[2]*diffTensorEVal[0] +
		       diffTensorEVal[1]*diffTensorEVal[0])/3.0;
    }
  else if (numModesOvr == 3)
    {
      D = (diffTensorEVal[0]+diffTensorEVal[1]+diffTensorEVal[2])/3.0;
      L2 = (diffTensorEVal[2]*diffTensorEVal[1] + diffTensorEVal[2]*diffTensorEVal[0] +
		       diffTensorEVal[1]*diffTensorEVal[0])/3.0;
    }
  else if (numModesOvr == 5)
    { 	  
      D = (diffTensorEVal[0]+diffTensorEVal[1]+diffTensorEVal[2])/3.0;
      L2 = (diffTensorEVal[2]*diffTensorEVal[1] + diffTensorEVal[2]*diffTensorEVal[0] +
		       diffTensorEVal[1]*diffTensorEVal[0])/3.0;
    }
  else
    {
      cerr << "the number of modes of Overall rotations is not supported. numModesOvr = 1, 3, or 5. \n";
      exit(0);
    }
   
  // reading the information about the directions of spin vectors 
  ifstream vecInfoFile;
  vecInfoFile.open(vecFile.data(),ios::in); 
  if(vecInfoFile.fail())
    { 
      cerr<<" Error opening the input file for spin vectors. \n";  
      cerr<<" Please check the file: " << vecFile.data() << endl; 
      exit(0);
    }
  cout << "\n Reading the initial information about the directions of spin vectors from the file: " << vecFile <<endl; 

  while( getline(vecInfoFile, linein))
    {
      ssin.clear();
      ssin.str(linein);
      ssin >> lineHead; 
      if (lineHead == "#")
	{continue;}
      else if (lineHead == "VecData:" )
	{
	  vector <double >  temp_vect; 
	  for (int i = 0; i < 3; i++)
	    {
	      double temp0; 
	      ssin >> temp0;
	      temp_vect.push_back(temp0);
	    } 
	  spinVecs.push_back(temp_vect);
	} 
      else
	{continue;}
    }
  vecInfoFile.close();
  
  cout << "\n The initial information about the directions of spin vectors: " << vecFile <<endl; 
  printSpinVecs(cout);
  calcVecAngles(); 
  cout << "\n The initial sperical angles of spin vectors: " << vecFile <<endl;
  printVecAngles(cout);
  calcCoeffsCFOAng();
  calcCoeffsCFT();
  calcRelaxParm(); 
  calcRFactor(); 

  cout << "\n The initial calculated values of relaxation parameters before the optimization from files " 
       << getIniFile() << " and " << vecFile <<  ":  \n";
  printRelaxParm(cout);
  cout << "\n With the initial R factor: " << ""  << RFactor << " \n";

}

void relaxParm::printTensor(ostream & outStream)
{

  outStream << setiosflags(ios::fixed) << setw(20) << "RotDiffTensor:" ;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      outStream << setiosflags(ios::scientific) << setw(12) << setprecision(5) << diffTensor[i][j];
  outStream << endl;

  outStream << setiosflags(ios::fixed) << setw(20) << "Eigenvalues:";
  for (int i = 0; i < 3; i++)
    outStream << setiosflags(ios::scientific) << setw(12) << setprecision(5)  << diffTensorEVal[i];
  outStream << endl;

  outStream << setiosflags(ios::fixed) << setw(20) << "Eigenvectors:" ;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      outStream << setiosflags(ios::scientific) << setw(12) << setprecision(5) << diffTensorEVec[i][j];
  outStream << endl;

  outStream << setiosflags(ios::fixed) << setw(20) << "RelaxTimeOvr:" 
	    << setiosflags(ios::fixed) << setw(12) << numModesOvr; 
  for (unsigned i = 0; i < numModesOvr; i++)
    outStream << setiosflags(ios::scientific) << setw(12) << setprecision(5)  << tauOvr[i];
  outStream << endl;

}

void relaxParm::printRelaxParm(ostream & outStream)
{
 
  for (unsigned i = 0; i < numSpinVecs; i++ )
    {
      outStream << setiosflags(ios::fixed) << setw(20) << "RelaxRate:" ;
      for (unsigned ifl = 0; ifl < nBFields; ifl++)
	outStream  << setiosflags(ios::scientific) << setw(12) << setprecision(5) << R1[ifl*numSpinVecs+i]
		   << setiosflags(ios::scientific) << setw(12) << setprecision(5) << R2[ifl*numSpinVecs+i]
		   << setiosflags(ios::scientific) << setw(12) << setprecision(5) << NOE[ifl*numSpinVecs+i];
      outStream << endl;
    } 
}

void relaxParm::printRelaxInt(ostream & outStream)
{
  for (unsigned i = 0; i < numSpinVecs; i++ )
    {
      outStream << setiosflags(ios::fixed) << setw(20) << "RelaxTimeInt:"
		<< setiosflags(ios::fixed) << setw(12) << numModesInt[i]; 
      for (unsigned j = 0; j < AInt[i].size(); j++)
	if (j == 0)
	  outStream << setiosflags(ios::scientific) << setw(12) << setprecision(5) << AInt[i][j];
	else
	  outStream << setiosflags(ios::scientific) << setw(12) << setprecision(5) << AInt[i][j]
		    << setiosflags(ios::scientific) << setw(12) << setprecision(5) << tauInt[i][j];
      outStream << endl;
    }
}

void relaxParm::printSpinVecs(ostream & outStream)
{
  for (unsigned i = 0; i < numSpinVecs; i++ )
    {
      outStream << setiosflags(ios::fixed) << setw(20) << "VecData:";
      for (unsigned j = 0; j < spinVecs[i].size(); j++)
	outStream << setiosflags(ios::scientific) << setw(12) << setprecision(5) << spinVecs[i][j];
      outStream << endl;
    }
}

void relaxParm::printVecAngles(ostream & outStream)
{
  for (unsigned i = 0; i < numSpinVecs; i++ )
    {
      outStream << setiosflags(ios::fixed) << setw(20) << "VecAngles:";
      for (unsigned j = 0; j < vecAngles[i].size(); j++)
	outStream << setiosflags(ios::scientific) << setw(12) << setprecision(5) << vecAngles[i][j];
      outStream << endl;
    }
}

void relaxParm::write()  
{
  ofstream outInfoFile;
  outInfoFile.open(getOutFile().data(),ios::out); 

  outInfoFile << "\n Final result about diffusion tensor: \n";
  printTensor(outInfoFile); 
  outInfoFile << "\n Final result about relaxation parameters:  \n";
  printRelaxParm(outInfoFile);
  outInfoFile << "\n Final result about internal motions: \n";
  printRelaxInt(outInfoFile);
  outInfoFile << "\n Final result about spherical angles of spin vectors: \n";
  printVecAngles(outInfoFile);
  outInfoFile << "\n The R factor comparing with initial input: "
              << setiosflags(ios::scientific) << setw(12) << setprecision(8) <<  RFactor << endl;   
  outInfoFile.close();
 
}

void relaxParm::calcVecAngles()
{
  // loop through all spin vectors
  vecAngles.clear();
  for (unsigned ivec = 0; ivec < numSpinVecs; ivec++)
    {
      double normz = sqrt(diffTensorEVec[0][0]*diffTensorEVec[0][0] + 
			  diffTensorEVec[0][1]*diffTensorEVec[0][1] +
			  diffTensorEVec[0][2]*diffTensorEVec[0][2]); 
      
      double normv = sqrt(spinVecs[ivec][0]*spinVecs[ivec][0] +
			  spinVecs[ivec][1]*spinVecs[ivec][1] + 
			  spinVecs[ivec][2]*spinVecs[ivec][2]); 
  
      double cos_theta = (spinVecs[ivec][0]*diffTensorEVec[0][0] + 
			  spinVecs[ivec][1]*diffTensorEVec[0][1] +
			  spinVecs[ivec][2]*diffTensorEVec[0][2])/(normv*normz);
      double theta0 = acos(cos_theta);
      double y0 = spinVecs[ivec][0]*diffTensorEVec[1][0] + 
                  spinVecs[ivec][1]*diffTensorEVec[1][1] +
	          spinVecs[ivec][2]*diffTensorEVec[1][2];
      double x0 = spinVecs[ivec][0]*diffTensorEVec[2][0] + 
                  spinVecs[ivec][1]*diffTensorEVec[2][1] +
	          spinVecs[ivec][2]*diffTensorEVec[2][2];      
      double phi0 = atan2(y0,x0);

      vector <double> angles_temp;
      angles_temp.push_back(theta0*180/PI);
      angles_temp.push_back(phi0*180/PI);
      vecAngles.push_back(angles_temp);  
    }
}

void relaxParm:: calcCoeffsCFOVec()
{
  // loop through all spin vectors
  AOvr.clear();
  for (unsigned ivec = 0; ivec < numSpinVecs; ivec++)
    {
      vector <double> AO;  

      // calculate the coeffs of overall CF for each spin vector 
      if (numModesOvr == 1)
	{
	  double A = 1.0; AO.push_back(A); 
	}
      else if (numModesOvr == 3)
	{
	  double cos_theta, cos_theta2, sin_theta2;
	  double normz = sqrt(diffTensorEVec[0][0]*diffTensorEVec[0][0] + 
			      diffTensorEVec[0][1]*diffTensorEVec[0][1] +
			      diffTensorEVec[0][2]*diffTensorEVec[0][2]); 
	  
	  double normv = sqrt(spinVecs[ivec][0]*spinVecs[ivec][0] +
			      spinVecs[ivec][1]*spinVecs[ivec][1] + 
			      spinVecs[ivec][2]*spinVecs[ivec][2]); 
	  
	  cos_theta = (spinVecs[ivec][0]*diffTensorEVec[0][0] + 
		       spinVecs[ivec][1]*diffTensorEVec[0][1] +
		       spinVecs[ivec][2]*diffTensorEVec[0][2])/(normv*normz);
	  cos_theta2 = cos_theta*cos_theta;
	  sin_theta2 = 1.0 - cos_theta2;
	  
	  double A = (1.5*cos_theta2- 0.5)*(1.5*cos_theta2- 0.5); AO.push_back(A); 
	  double B = 3.0*cos_theta2*sin_theta2;   AO.push_back(B);
	  double C = 0.75*sin_theta2*sin_theta2;  AO.push_back(C);
	}
      else if (numModesOvr == 5)
	{
	  double cos_theta[3], cos_theta2[3], cos_theta4[3], delta[3];
	  double normv = sqrt(spinVecs[ivec][0]*spinVecs[ivec][0] +
			      spinVecs[ivec][1]*spinVecs[ivec][1] + 
			      spinVecs[ivec][2]*spinVecs[ivec][2]); 
	  for (unsigned i=0; i < 3; i++)
	    {
	      double normd = sqrt(diffTensorEVec[i][0]*diffTensorEVec[i][0] + 
				  diffTensorEVec[i][1]*diffTensorEVec[i][1] +
				  diffTensorEVec[i][2]*diffTensorEVec[i][2]); 
	      
	      cos_theta[i] = (spinVecs[ivec][0]*diffTensorEVec[i][0] + 
			      spinVecs[ivec][1]*diffTensorEVec[i][1] +
			      spinVecs[ivec][2]*diffTensorEVec[i][2])/(normv*normd);
	      cos_theta2[i] = cos_theta[i]*cos_theta[i];
	      cos_theta4[i] = cos_theta2[i]*cos_theta2[i];
	    }

	  double tempv = 1.0/sqrt(D*D-L2);
	  for (unsigned i=0; i < 3; i++)
	    {
	      delta[i] = (diffTensorEVal[i]-D)*tempv;
	    } 
	  double A = 6.0*cos_theta2[1]*cos_theta2[0]; AO.push_back(A); 
	  double B = 6.0*cos_theta2[2]*cos_theta2[0]; AO.push_back(B);
	  double C = 6.0*cos_theta2[2]*cos_theta2[1]; AO.push_back(C);
	  double d = (3.0*(cos_theta4[0]+cos_theta4[1]+cos_theta4[0]) -1.0)*0.5;
	  double e = ( delta[2]*(3*cos_theta4[2] + A - 1.0) +  delta[1]*(3*cos_theta4[1] + B  -1.0) +
		       delta[0]*(3*cos_theta4[0] + C - 1.0) )/3.0; 
	  double D1 = d - e;  AO.push_back(D1); 
	  double E1 = d + e;  AO.push_back(E1);
	}
      else
	{
	  cerr << "the number of modes of Overall rotations is not supported. numModesOvr = 1, 3, or 5. \n";
	  exit(0);
	}
      // calcCoeffsCFO(AO, ivec);
      AOvr.push_back(AO);      
    }
}

void relaxParm:: calcCoeffsCFOAng()
{
  // loop through all spin vectors
  AOvr.clear();
  for (unsigned ivec = 0; ivec < numSpinVecs; ivec++)
    {
      vector <double> AO;  

      // calculate the coeffs of overall CF for each spin vector 
      if (numModesOvr == 1)
	{
	  double A = 1.0; AO.push_back(A); 
	}
      else if (numModesOvr == 3)
	{
	  double cos_theta, cos_theta2, sin_theta2;
  
	  cos_theta = cos(vecAngles[ivec][0]*PI/180);
	  cos_theta2 = cos_theta*cos_theta;
	  sin_theta2 = 1.0 - cos_theta2;

	  double A = (1.5*cos_theta2- 0.5)*(1.5*cos_theta2- 0.5); AO.push_back(A); 
	  double B = 3.0*cos_theta2*sin_theta2;   AO.push_back(B);
	  double C = 0.75*sin_theta2*sin_theta2;  AO.push_back(C);
	}
      else if (numModesOvr == 5)
	{


	  double cos_theta[3], cos_theta2[3], cos_theta4[3], delta[3];

	  cos_theta[0] = cos(vecAngles[ivec][0]*PI/180);
	  cos_theta[1] = sin(vecAngles[ivec][0]*PI/180)*sin(vecAngles[ivec][1]*PI/180);
	  cos_theta[2] = sin(vecAngles[ivec][0]*PI/180)*cos(vecAngles[ivec][1]*PI/180);
	  
	  for (unsigned i=0; i < 3; i++)
	    {
	      cos_theta2[i] = cos_theta[i]*cos_theta[i];
	      cos_theta4[i] = cos_theta2[i]*cos_theta2[i];
	    }

	  double tempv = 1.0/sqrt(D*D-L2);
	  for (unsigned i=0; i < 3; i++)
	    {
	      delta[i] = (diffTensorEVal[i]-D)*tempv;
	    } 

	  double A = 6.0*cos_theta2[1]*cos_theta2[0]; AO.push_back(A); 
	  double B = 6.0*cos_theta2[2]*cos_theta2[0]; AO.push_back(B);
	  double C = 6.0*cos_theta2[2]*cos_theta2[1]; AO.push_back(C);
	  double d = (3.0*(cos_theta4[0]+cos_theta4[1]+cos_theta4[0]) -1.0)*0.5;
	  double e = ( delta[2]*(3*cos_theta4[2] + A - 1.0) +  delta[1]*(3*cos_theta4[1] + B - 1.0) +
		       delta[0]*(3*cos_theta4[0] + C - 1.0) )/3.0; 
	  double D = d - e;  AO.push_back(D); 
	  double E = d + e;  AO.push_back(E);

	}
      else
	{
	  cerr << "the number of modes of Overall rotations is not supported. numModesOvr = 1, 3, or 5. \n";
	  exit(0);
	}  
      // calcCoeffsCFO(AO, ivec);
      AOvr.push_back(AO);      
    }
}

void relaxParm:: updateDiffTensorEVal()
{
  if (numModesOvr == 1)
    {
      double tau1_inv = 1.0E12/tauOvr[0];    //  1.0E12 because of tauOvr stored in the unit of ps 
      diffTensorEVal[0] =  tau1_inv/6.0;
      diffTensorEVal[1] = diffTensorEVal[0];
      diffTensorEVal[2] = diffTensorEVal[0];
    }
  else if (numModesOvr == 3)
    {
      double tau1_inv = 1.0E12/tauOvr[0]; 
      double tau2_inv = 1.0E12/tauOvr[1];
      diffTensorEVal[0] = (tau2_inv - 5.0*tau1_inv/6.0);
      diffTensorEVal[1] = tau1_inv/6.0; 
      diffTensorEVal[2] = diffTensorEVal[1];
      tauOvr[2] = 1.0E12/(4.0*diffTensorEVal[0] + 2.0*diffTensorEVal[1]);
      D = (diffTensorEVal[0]+diffTensorEVal[1]+diffTensorEVal[2])/3.0;
    }
  else if (numModesOvr == 5)
    {
      double tau1_inv = 1.0E12/tauOvr[0]; 
      double tau2_inv = 1.0E12/tauOvr[1];
      double tau3_inv = 1.0E12/tauOvr[2];
      double sumTau_inv = (tau1_inv + tau2_inv + tau3_inv)/6.0;
      diffTensorEVal[0] = (tau3_inv - sumTau_inv)/3.0; 
      diffTensorEVal[1] = (tau2_inv - sumTau_inv)/3.0;
      diffTensorEVal[2] = (tau1_inv - sumTau_inv)/3.0;       	  
      D = (diffTensorEVal[0]+diffTensorEVal[1]+diffTensorEVal[2])/3.0;
      L2 = (diffTensorEVal[2]*diffTensorEVal[1] + diffTensorEVal[2]*diffTensorEVal[0] +
		       diffTensorEVal[1]*diffTensorEVal[0])/3.0;
      double tempv = sqrt(D*D-L2);
      tauOvr[3] = 1.0E12/(6*(D + tempv));
      tauOvr[4] = 1.0E12/(6*(D - tempv));
    }
  else
    {
      cerr << "the number of modes of Overall rotations is not supported. numModesOvr = 1, 3, or 5. \n";
      exit(0);
    }
}


void relaxParm:: calcCoeffsCFT()
{
  tauTot.clear();  ATot.clear();
  // loop through all spin vectors
  for (unsigned ivec = 0; ivec < numSpinVecs; ivec++)
    {

      vector <double> AT, TT;  
      double coeff;
      double tau;
      vector <double> AI, AO; 
      vector <double> TI, TO;

      AI = AInt[ivec]; TI = tauInt[ivec]; 
      AO = AOvr[ivec]; TO = tauOvr; 
  
      // calculate the coeffs and time constants of total CF for spin vector ivec    
      for (unsigned i = 0; i< AI.size(); i++)
	{
	  if (i == 0)
	    {
	      for (unsigned j = 0; j< AO.size(); j++)
		{
		  coeff = AI[i]*AO[j];
		  tau =TO[j];
		  AT.push_back(coeff);
		  TT.push_back(tau);
		}
	    }
	  else
	    {
	      for (unsigned j = 0; j< AO.size(); j++)
		{
		  coeff = AI[i]*AO[j];
		  tau = (TI[i]*TO[j])/(TI[i]+TO[j]);
		  AT.push_back(coeff);
		  TT.push_back(tau);
		}
	    }
	}
  
      // change the time unit to second 
      for (unsigned i =0; i< TT.size(); i++)
	TT[i] = TT[i]/1.0e12;
     
      // calcCoeffsCFT(AT,TT, ivec);
      ATot.push_back(AT);
      tauTot.push_back(TT);      
    }
}

void relaxParm::calcRelaxParm()
{
  // Calculate the relaxation parameters R1/R2 from the multi-exponential approximations
  // of C2 rotational correlation function of C13-H1 or N15-H1 vector. The fitting 
  // parameters are passed to get the sectral density function and the R1/R2 by the 
  // method described in the papers:        
  // 1) N. L. Fawzi, A. H. Phillips, J. Z. Ruscio, M. Doucleff, D. E. Wemmer, and T. Head-Gordon,            
  //    J. Am. Chem. Soc. 130, 6145 (2008);                                                                  
  // 2) A. G. Palmer III, M. Rance, and P. E. Wright, J. Am. Chem. Soc. 113 4371 (1991).    

  // loop through all spin vectors 
  for (unsigned ifl = 0; ifl < nBFields;  ifl++)
    {
      for (unsigned ivec = 0; ivec < numSpinVecs; ivec++)
	{
	  
	  vector <double> AT, TT; 
	  
	  // get coefficients and time constants for total correlation function 
	  // getCoeffsCFT(AT,TT,ivec);
	  AT = ATot[ivec]; TT = tauTot[ivec];
  
	  // set default parameters for spin i and spin j
	  double GAMMA_i, GAMMA_j, OMEGA_i, OMEGA_j, DELTA_i, Reff; 
	  if (spinTypes[ivec][0] == "C")
	    {
	      GAMMA_i = Gamma_C;
	      OMEGA_i = 2*PI*LFreq_C*BFields[ifl]; 
	      DELTA_i = Delta_C;
	      Reff = Reff_C;
	    } 
	  else if (spinTypes[ivec][0] == "N")
	    {
	      GAMMA_i = Gamma_N;
	      OMEGA_i = 2*PI*LFreq_N*BFields[ifl]; 
	      DELTA_i = Delta_N;
	      Reff = Reff_N;
	    }
	  else
	    {
	      cerr<<" The i spin type " << spinTypes[ivec][0] << " is not implemented yet.\n";  
	      exit(0);
	    }
	  if (spinTypes[ivec][1] == "H")
	    {
	      GAMMA_j = Gamma_H;
	      OMEGA_j = 2*PI*LFreq_H*BFields[ifl]; 
	    }
	  else 
	    {
	      cerr<<" The j spin type " << spinTypes[ivec][1] << " is not implemented yet.\n";  
	      exit(0);
	    }
	  
	  // Calculated R1DD, R1CSA, R2DD, and R2CSA first, and then R1, R2 and NOE
	  double K, Jw0, JwJmI, JwI, JwJ, JwJpI, R1DD, R2DD, R1CSA, R2CSA, R2a;  
	  K = MU_0*H_BAR*GAMMA_i*GAMMA_j/(4*PI*Reff*Reff*Reff); // note the power of 3 is missing in the paper
	  Jw0   = getJw(AT,TT,0);
	  JwJmI = getJw(AT,TT,OMEGA_j-OMEGA_i);
	  JwI   = getJw(AT,TT,OMEGA_i);
	  JwJ   = getJw(AT,TT,OMEGA_j);
	  JwJpI = getJw(AT,TT,OMEGA_j+OMEGA_i);
	  
	  // note in the  Fawzi's paper an additional 5 is divided for R1DD, R2DD, R1CSA, and R2CSA
	  R1DD = K*K*(JwJmI + 3.0*JwI + 6.0*JwJpI)/20.0;
	  R2DD = K*K*(4.0*Jw0 + JwJmI + 3.0*JwI + 6.0*JwJ +6.0*JwJpI)/40.0; // note 6.0 is 3.0 in the paper (a typo)
	  R1CSA =DELTA_i*DELTA_i*OMEGA_i*OMEGA_i*JwI/5.0;
	  R2CSA =DELTA_i*DELTA_i*OMEGA_i*OMEGA_i*(4.0*Jw0 + 3.0*JwI)/30.0;
	  R2a = 0.0;
	  
	  R1[ifl*numSpinVecs + ivec] = R1DD + R1CSA; 
	  R2[ifl*numSpinVecs + ivec] = R2DD + R2CSA + R2a;     
	  NOE[ifl*numSpinVecs + ivec] = 1 + (GAMMA_j/GAMMA_i)*(6.0*JwJpI - JwJmI)/
	    (JwJmI + (3.0 + 4.0*DELTA_i*DELTA_i*OMEGA_i*OMEGA_i/(K*K))*JwI + 6.0*JwJpI);

	}
    }
   
}

void relaxParm::calcRFactor()
{
  // Calculate the relaxation parameters R1/R2 from the multi-exponential approximations
  
  // loop through all spin vectors
  RFactor = 0.0; 
  double tempv;   
  for (unsigned ivec = 0; ivec < R1.size(); ivec++)
    {
      tempv = R1[ivec]/R1In[ivec]-1;
      RFactor +=  tempv*tempv;
      tempv = R2[ivec]/R2In[ivec]-1;
      RFactor +=  tempv*tempv;
      tempv = NOE[ivec]/NOEIn[ivec]-1;
      RFactor +=  tempv*tempv;
    }
  RFactor = sqrt(RFactor)/(3*R1.size());
}


double relaxParm::getJw(vector <double> & AT0, vector <double> & TT0, double omega0)
{
  double finalvalue = 0.0;
  for (unsigned i= 0; i< TT0.size(); i++ )
    finalvalue = finalvalue + AT0[i]* 2.0 * TT0[i]/(1.0 + omega0*omega0*TT0[i]*TT0[i]);
  return finalvalue;
}

void relaxParm::assignSpinTypes()
{
  // loop through all spin vectors
  spinTypes.clear();
  for (unsigned ivec = 0; ivec < numSpinVecs; ivec++)
    {
      vector <string> spinPair;
      // spinPair.push_back(string(*atmTypeI[ivec].begin()));
      // spinPair.push_back(string(*atmTypeJ[ivec].begin()));
      spinPair.push_back(atmTypeI[ivec]);
      spinPair.push_back(atmTypeJ[ivec]);
      spinTypes.push_back(spinPair);  
    }
}
