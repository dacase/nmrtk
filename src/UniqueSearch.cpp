/*********************************************************************************************  
 *   unique_search.h, header file of unique_search.cpp to find the unique structures.        *
 *                                                                                           *  
 *   Authors:                                                                                *
 *          originally written  by Junchao Xia                                               *
 *                                                                                           *
 *   Modification Histry:                                                                    *
 *         Dec. 10, 2014:    originally written Junchao Xia.                                 *
 *                                                                                           *
 *********************************************************************************************/

#include "UniqueSearch.h"

UniqueSearch::UniqueSearch()
{
  ifstream inputparameter;
  inputparameter.open("unique_search_prm.inp",ios::in);
  if(inputparameter.fail())
    { cout<<" Error opening the input parameter file.\n" 
          << "unique_search_prm.inp doesn't not exist."<<endl;
      cout<<" CHECK FILENAME !\n";
      cout<<" CHECK FILENAME !\n";
    }
  readParaFile(inputparameter,cout);
  inputparameter.close();

  readEngFile();
  readAngFile();

}

UniqueSearch::UniqueSearch(string & inputfile)
{
  ifstream inputparameter;
  inputparameter.open(inputfile.data(),ios::in);
  if(inputparameter.fail())
    { cout<<" Error opening the input parameter file.\n"
          << inputfile.data() << " doesn't not exist."<<endl;
      cout<<" CHECK FILENAME !\n";
      cout<<" CHECK FILENAME !\n";
    }
  readParaFile(inputparameter,cout);
  inputparameter.close();

  readEngFile();
  readAngFile();

}


UniqueSearch::~UniqueSearch()
{
  // need to add operations to free the memory
}

void UniqueSearch::findUniqueConfs()
{

  cout << "*****Searching the unique structures from the allowed conformations:****" << endl;

  bool samestructure =false;
  bool maxangle=false;
  vector <int> conf_ids;         // save IDs after sorting
  vector <float> conf_engs;      // save energies after sorting

  for(unsigned int iconf = 0; iconf < numEng;iconf++) conf_ids.push_back(iconf);
  conf_engs = confEngs;
  sortEng(conf_ids,conf_engs); // sorting the ids and energies ascendingly

  
  ofstream confs_sorted_out;
  confs_sorted_out.open("confs_sorted.out", ios::out);
  for (unsigned int iconf = 0; iconf < conf_ids.size();iconf++)
    {
        confs_sorted_out << setw(10) << iconf + 1 
		     << setw(10) << conf_ids[iconf] 
		     << setw(10) << conf_engs[iconf] << " ";
      for(unsigned int iangle=0; iangle<confAngs[conf_ids[iconf]].size();iangle++)
	{
	  confs_sorted_out << setw(6) << confAngs[conf_ids[iconf]][iangle] << " ";

	} 
      confs_sorted_out << endl;
    }
  confs_sorted_out.close();
  
  cout << " done with sorting energy and outputing conformations sorted." << endl;

  // go through each sorted conformer and find unique structures
  for(unsigned int jconf = 0; jconf < conf_ids.size();jconf++) 
    {
       unsigned int iconf = conf_ids[jconf];

      if (uniq_ids.size()== 0) // no unique structure
	{
	  uniq_ids.push_back(iconf);
	  uniq_engs.push_back(confEngs[iconf]);
	  uniq_angs.push_back(confAngs[iconf]);
	}
      else if (uniq_ids.size()== 1) // one unique structure
	{
	  samestructure = false;
	  maxangle = false;
	  for (unsigned int iangle = 0; iangle < numAng;iangle++)
	    {
	      if (fabs(confAngs[iconf][iangle]-confAngs[uniq_ids[0]][iangle]) > minAngDiff[iangle]) 
		{ maxangle = true;}
	    }

	  if ( fabs(confEngs[iconf]-confEngs[uniq_ids[0]]) < minEngDiff  && !maxangle )
	    {
	      samestructure = true; 
	      if (confEngs[iconf]<confEngs[uniq_ids[0]])
		{
		  uniq_ids[0] = iconf;
		  uniq_engs[0] = confEngs[iconf];
		  uniq_angs[0] = confAngs[iconf]; 
		  
		}
	    }

	  if (!samestructure)
	    {
	      if (confEngs[iconf] > confEngs[uniq_ids[0]]) 
		{
		  uniq_ids.push_back(iconf);
		  uniq_engs.push_back(confEngs[iconf]);
		  uniq_angs.push_back(confAngs[iconf]);
		}
	      else
		{
		  uniq_ids.insert(uniq_ids.begin(),iconf);
		  //cout << "insert to the head of uniq_engs[0]" << endl;
		  uniq_engs.insert(uniq_engs.begin(),confEngs[iconf]);
		  //cout << "insert to the head of uniq_angs[0]" << endl;
		  uniq_angs.insert(uniq_angs.begin(),confAngs[iconf]);
		  //cout << "end of insert to the head of uniq_angs[0]" << endl;
		}
	      
	    }
	}
      else  // more than one unique structure
	{
	  for (unsigned int iunique = 0; iunique < uniq_ids.size();iunique++)
	    {
	      samestructure = false;
	      maxangle = false;
	      for (unsigned int iangle = 0; iangle < numAng;iangle++)
		{
		  if ( fabs(confAngs[iconf][iangle]-confAngs[uniq_ids[iunique]][iangle]) > minAngDiff[iangle] )
		    {maxangle = true;}
		}

	      if (fabs(confEngs[iconf]-confEngs[uniq_ids[iunique]]) < minEngDiff && !maxangle)
		{
		  samestructure = true; 
		  if (confEngs[iconf]<confEngs[uniq_ids[iunique]])
		    {
		      uniq_ids[iunique] = iconf;
		      uniq_engs[iunique] = confEngs[iconf]; 
		      uniq_angs[iunique] = confAngs[iconf];                      
		    }
		  break;
		}   
	    }

	  if (!samestructure)  
	    {
	      //cout << "conformer " << iconf << " is a new structure " << endl;
	      if (confEngs[iconf]>confEngs[uniq_ids[uniq_ids.size()-1]]) 
		{
		  uniq_ids.push_back(iconf);
		  uniq_engs.push_back(confEngs[iconf]);
		  uniq_angs.push_back(confAngs[iconf]);
		  //cout << " add it to the end of vector." << endl;
		}
	      else if (confEngs[iconf]<confEngs[uniq_ids[0]])
		{
		  //cout << "insert to the head of uniq_id." << endl;
		  uniq_ids.insert(uniq_ids.begin(),iconf);
		  //cout << "insert to the head of uniq_engs." << endl;
		  uniq_engs.insert(uniq_engs.begin(),confEngs[iconf]);
		  //cout << "insert to the head of uniq_angs." << endl;
		  uniq_angs.insert(uniq_angs.begin(),confAngs[iconf]);
		  //cout << "end of insert to the head." << endl;
		}
	      else
		{
		  //cout << " need add it in the middle" << endl;
		  int iu;
		  vector< int >::iterator uni_ids_d1itr;
		  vector< float >::iterator uni_engs_d1itr;
		  vector< vector< float > >::iterator uni_angs_d2itr;
		  for (uni_ids_d1itr=uniq_ids.begin(),uni_engs_d1itr=uniq_engs.begin(),uni_angs_d2itr=uniq_angs.begin(),iu =0;
		       uni_ids_d1itr < uniq_ids.end();uni_ids_d1itr++,uni_engs_d1itr++,uni_angs_d2itr++,iu++)
		    {
		      //cout << " the unique structure is " << iu << endl; 
		      if (confEngs[iconf] < uniq_engs[iu])
			{
			  //cout << "insert to the middle of uniq_id." << endl;
			  uniq_ids.insert(uni_ids_d1itr,iconf);
			  uniq_engs.insert(uni_engs_d1itr,confEngs[iconf]);
			  uniq_angs.insert(uni_angs_d2itr,confAngs[iconf]); 
			  //cout << "end of insert to the middle of uniq_id." << endl;
			  break;
			}
		    }
		}
	    }   
	  
	}
      
    } 

  cout << " Done with searching the unique structures from all conformations." << endl;  
 
}

void UniqueSearch::findUniqueConfsNoEng()
{

  cout << "*****Searching the unique structures from the allowed conformations:****" << endl;
  bool samestructure =false;
  bool maxangle=false;
  vector <int> conf_ids;         // save IDs after sorting
  vector <float> conf_engs;      // save energies after sorting

  for(unsigned int iconf = 0; iconf < numEng;iconf++) conf_ids.push_back(iconf);
  conf_engs = confEngs;
  sortEng(conf_ids,conf_engs); // sorting the ids and energies ascendingly

  
  ofstream confs_sorted_out;
  confs_sorted_out.open("confs_sorted.out", ios::out);
  for (unsigned int iconf = 0; iconf < conf_ids.size();iconf++)
    {
        confs_sorted_out << setw(10) << iconf + 1 
		     << setw(10) << conf_ids[iconf] 
		     << setw(10) << conf_engs[iconf] << " ";
      for(unsigned int iangle=0; iangle<confAngs[conf_ids[iconf]].size();iangle++)
	{
	  confs_sorted_out << setw(6) << confAngs[conf_ids[iconf]][iangle] << " ";

	} 
      confs_sorted_out << endl;
    }
  confs_sorted_out.close();
  
  cout << " done with sorting energy and outputing conformations sorted." << endl;

  // go through each sorted conformer and find unique structures
  for(unsigned int jconf = 0; jconf < conf_ids.size();jconf++) 
    {
       unsigned int iconf = conf_ids[jconf];

      if (uniq_ids.size()== 0) // no unique structure
	{
	  uniq_ids.push_back(iconf);
	  uniq_engs.push_back(confEngs[iconf]);
	  uniq_angs.push_back(confAngs[iconf]);
	}
      else if (uniq_ids.size()== 1) // one unique structure
	{
	  samestructure = false;
	  maxangle = false;
	  for (unsigned int iangle = 0; iangle < numAng;iangle++)
	    {
	      if (fabs(confAngs[iconf][iangle]-confAngs[uniq_ids[0]][iangle]) > minAngDiff[iangle]) 
		{ 
                   maxangle = true;
                }
	    }

	  if (!maxangle)
	    {
              samestructure = true;
	      if (confEngs[iconf]<confEngs[uniq_ids[0]])
		{
		  uniq_ids[0] = iconf;
		  uniq_engs[0] = confEngs[iconf];
		  uniq_angs[0] = confAngs[iconf]; 
		  
		}
	    }

	  if (!samestructure)
	    {
	      if (confEngs[iconf] > confEngs[uniq_ids[0]]) 
		{
		  uniq_ids.push_back(iconf);
		  uniq_engs.push_back(confEngs[iconf]);
		  uniq_angs.push_back(confAngs[iconf]);
		}
	      else
		{
		  uniq_ids.insert(uniq_ids.begin(),iconf);
		  //cout << "insert to the head of uniq_engs[0]" << endl;
		  uniq_engs.insert(uniq_engs.begin(),confEngs[iconf]);
		  //cout << "insert to the head of uniq_angs[0]" << endl;
		  uniq_angs.insert(uniq_angs.begin(),confAngs[iconf]);
		  //cout << "end of insert to the head of uniq_angs[0]" << endl;
		}
	      
	    }
	}
      else  // more than one unique structure
	{
	  for (unsigned int iunique = 0; iunique < uniq_ids.size();iunique++)
	    {
	      samestructure = false;
	      maxangle = false;
	      for (unsigned int iangle = 0; iangle < numAng;iangle++)
		{
		  if ( fabs(confAngs[iconf][iangle]-confAngs[uniq_ids[iunique]][iangle]) > minAngDiff[iangle] )
		    {
                       maxangle = true;
                    }
		}

	      if (!maxangle)
		{
		  samestructure = true; 
		  if (confEngs[iconf]<confEngs[uniq_ids[iunique]])
		    {
		      uniq_ids[iunique] = iconf;
		      uniq_engs[iunique] = confEngs[iconf]; 
		      uniq_angs[iunique] = confAngs[iconf];                      
		    }
		  break;
		}   
	    }

	  if (!samestructure)  
	    {
	      //cout << "conformer " << iconf << " is a new structure " << endl;
	      if (confEngs[iconf]>confEngs[uniq_ids[uniq_ids.size()-1]]) 
		{
		  uniq_ids.push_back(iconf);
		  uniq_engs.push_back(confEngs[iconf]);
		  uniq_angs.push_back(confAngs[iconf]);
		  //cout << " add it to the end of vector." << endl;
		}
	      else if (confEngs[iconf]<confEngs[uniq_ids[0]])
		{
		  //cout << "insert to the head of uniq_id." << endl;
		  uniq_ids.insert(uniq_ids.begin(),iconf);
		  //cout << "insert to the head of uniq_engs." << endl;
		  uniq_engs.insert(uniq_engs.begin(),confEngs[iconf]);
		  //cout << "insert to the head of uniq_angs." << endl;
		  uniq_angs.insert(uniq_angs.begin(),confAngs[iconf]);
		  //cout << "end of insert to the head." << endl;
		}
	      else
		{
		  //cout << " need add it in the middle" << endl;
		  int iu;
		  vector< int >::iterator uni_ids_d1itr;
		  vector< float >::iterator uni_engs_d1itr;
		  vector< vector< float > >::iterator uni_angs_d2itr;
		  for (uni_ids_d1itr=uniq_ids.begin(),uni_engs_d1itr=uniq_engs.begin(),uni_angs_d2itr=uniq_angs.begin(),iu =0;
		       uni_ids_d1itr < uniq_ids.end();uni_ids_d1itr++,uni_engs_d1itr++,uni_angs_d2itr++,iu++)
		    {
		      //cout << " the unique structure is " << iu << endl; 
		      if (confEngs[iconf] < uniq_engs[iu])
			{
			  //cout << "insert to the middle of uniq_id." << endl;
			  uniq_ids.insert(uni_ids_d1itr,iconf);
			  uniq_engs.insert(uni_engs_d1itr,confEngs[iconf]);
			  uniq_angs.insert(uni_angs_d2itr,confAngs[iconf]); 
			  //cout << "end of insert to the middle of uniq_id." << endl;
			  break;
			}
		    }
		}
	    }   
	  
	}
      
    } 

  cout << " Done with searching the unique structures from all conformations." << endl;  
 
}

void UniqueSearch::writeUniqueConfs()
{
  ofstream uniqueconfsout;
  uniqueconfsout.open("unique_confs.out", ios::out);
  for (unsigned int iunique = 0; iunique < uniq_ids.size();iunique++)
    {
      uniqueconfsout << setw(6) << iunique + 1 
		     << setw(6) << uniq_ids[iunique] 
		     << setw(10) << uniq_engs[iunique] << " ";
      for(unsigned int iangle=0; iangle<uniq_angs[iunique].size();iangle++)
	{
	  uniqueconfsout << setw(6) << uniq_angs[iunique][iangle] << " ";

	} 

      uniqueconfsout << endl;

    }
  uniqueconfsout.close(); 

  uniqueconfsout.open("unique_angs.out", ios::out);
  for (unsigned int iunique = 0; iunique < uniq_ids.size();iunique++)
    {
      for(unsigned int iangle=0; iangle<uniq_angs[iunique].size();iangle++)
	{
	  uniqueconfsout << setw(6) << uniq_angs[iunique][iangle] << " ";
	} 
      uniqueconfsout << endl;

    }
  uniqueconfsout.close(); 

  uniqueconfsout.open("unique_engs.out", ios::out);
  for (unsigned int iunique = 0; iunique < uniq_ids.size();iunique++)
    {
      uniqueconfsout  << setw(10) << uniq_engs[iunique] << " " << endl;
    }
  uniqueconfsout.close(); 

}


void UniqueSearch::readParaFile(ifstream &inputparameter, ostream & infoout)
{

  string linein;
  string lineHead; 
  istringstream ssin;
  while( getline(inputparameter, linein))
    {
      ssin.clear();
      ssin.str(linein);
      ssin >> lineHead; 
      if (lineHead == "#")
	{continue;}
      else if (lineHead == "EngFile:")
	{
	  ssin >> engFile;
	  cout << " EngFile: " <<  engFile << endl;
	}
      else if (lineHead == "AngFile:" )
	{ 
	  ssin >> angFile;
	  cout << " AngFile: " << angFile << endl;
	}
      else if (lineHead == "NumEng:" )
	{ 
	  ssin >> numEng;
	  cout << " NumEng: " << numEng << endl;
	}
      else if (lineHead == "NumAng:" )
	{
	  ssin >> numAng;
	  cout << " NumAng: " << numAng << endl;
	}  
     else if (lineHead == "MinEngDiff:" )
	{
	  ssin >> minEngDiff; 
	  cout << " MinEngDiff: " << minEngDiff << endl;
	}  
     else if (lineHead == "MinAngDiff:" )
	{
	  float minval; 
	  ssin >> minval;
          minAngDiff.push_back(minval); 
	  cout << " MinAngDiff: " << minval << endl;
	}  
      else
	{continue;}
    }

}

void UniqueSearch::readAngFile()
{
  float angle;
  ifstream angsfilein;

  // input the phi-psi conformations
  angsfilein.open(angFile.data(), ios::in);
  if( angsfilein.fail() )
    {
      cout << "Failed to open the conformations file of angles " << angFile.data() << endl;
      exit(0);
    }
  
  for(unsigned int iconfs = 0; iconfs < numEng; iconfs++)
    {
      vector <float> temp_vect;
      for(unsigned int iangle = 0; iangle < numAng; iangle++)
	{
	  angsfilein >> angle; 
	  temp_vect.push_back(angle);
	}
      confAngs.push_back(temp_vect);      
      angsfilein.ignore(120,'\n');
      
    }
  angsfilein.close();  // input the energy of conformation

}

void UniqueSearch::readEngFile()
{
  float energy;

  ifstream engsfilein;

  // input the energy of conformation
  engsfilein.open(engFile.data(), ios::in);
  if( engsfilein.fail() )
    {
      cout << "Failed to open the energy file of conformations " << engFile.data() << endl;
      exit(0);
    }
  
  for(unsigned int iconfs = 0; iconfs < numEng; iconfs++)
    {
      engsfilein >> energy; 
      confEngs.push_back(energy);
      engsfilein.ignore(120,'\n');

    }
  engsfilein.close();
}


void UniqueSearch::sortEng(vector <int> &ids,vector <float> &engs)
{
  // sort the conformer ids and energies via Cocktail-Shaker Sort Algorithm 
  int tempInt;
  float tempFloat;
  bool sorted;

  sorted = false;

  while (!sorted)
    {
      sorted = true;
      for (unsigned int j = 0; j < engs.size()- 1; j++)
	{

	  if (engs[j]> engs[j+1])
	    {
	      tempFloat = engs[j];
	      tempInt = ids[j];
	      engs[j] = engs[j+1];
	      ids[j] = ids[j+1];
	      engs[j+1] = tempFloat;
	      ids[j+1] = tempInt;
	      sorted = false;
	    }
	}

      if (sorted) return;

      for (int j = engs.size()-3; j >= 0; j--)
	{
	  if (engs[j]> engs[j+1])
	    {
	      tempFloat = engs[j];
	      tempInt = ids[j];
	      engs[j] = engs[j+1];
	      ids[j] = ids[j+1];
	      engs[j+1] = tempFloat;
	      ids[j+1] = tempInt;
	      sorted = false;
	    }

	}

      if (sorted) return;      
      
    }
  
  return;
}



// argv[1] = param_file
int main(int argc, char **argv)
{
  
  string param_file; 
  time_t t_start, t_end;

  t_start = time(NULL);
  
  if(argc<2)
    {
      cerr << " Check the command syntax!" << endl;
      cerr << " It should be 'unique_search param_file'" << endl;
      return 0;
    }


  cout << "************************************************************************" << endl;
  cout << "*                    START UNIQUE CONFORMATION SEARCH                  *" << endl;
  cout << "************************************************************************" << endl;



  cout << "*************************Reading Parameters:****************************" << endl;
  param_file = argv[1];
  UniqueSearch unique_srch = UniqueSearch(param_file);  // read in the parameters to start program
  cout << " Done with reading parameter file and data file." << endl; 

  cout << "*********Searching the unique structures from all conformations:********" << endl;
  // unique_srch.findUniqueConfs(); 
  unique_srch.findUniqueConfsNoEng();  
  cout << " Done with searching the unique structures." << endl;  

  cout << "*********************Outputing the unique structures:*******************" << endl;
  unique_srch.writeUniqueConfs();  
  cout << " Done with outputing the unique structures." << endl;



  cout << "***********************Output the CPU time information:*****************" << endl;
  t_end = time(NULL);
  cout << " The CPU time for this search is " << t_end - t_start << " seconds." << endl;

 
}
