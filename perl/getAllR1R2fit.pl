#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;

#
my $myhead = <<end_of_myhead;
###########################################################################################################
# This a perl program to calcuate the relaxation parameters R1/R2 from the C2 correlation function of     #
# C13-H1 or N15-H1 vector obtained from the AMBER ptraj analysis. C2 is fitted to multi-exponential       #
# functions via Gnuplot. The fitting parameters are passed to get the sectral density function and the    #
# R1/R2 by the method described in the papers:                                                            #
# 1) N. L. Fawzi, A. H. Phillips, J. Z. Ruscio, M. Doucleff, D. E. Wemmer, and T. Head-Gordon,            #
#    J. Am. Chem. Soc. 130, 6145 (2008);                                                                  #
# 2) A. G. Palmer III, M. Rance, and P. E. Wright, J. Am. Chem. Soc. 113 4371 (1991).                     #
#                                                                                                         #
# The total correlation function is fitted to different forms as below (nexp = ?):                        #           
# 0) Ct(t) = (S2+(1-S2)*exp(-t/tau_e))*exp(-t/tau_m)=S2*exp(-t/tau_m)+(1-S2)*exp(-t/tau) mean-field model # 
# 1) Ct(t) = exp(-t/t1)                                                                                   # 
# 2) Ct(t) = a1*exp(-t/t1) + (1-a1)*exp(-t/t2)                                                            #
# 3) Ct(t) = a1*exp(-t/t1) + a2*exp(-t/t2) + (1-a1-a2)*exp(-t/t3)                                         #
# 4) Ct(t) = a1*exp(-t/t1) + a2*exp(-t/t2) + a3*exp(-t/t3) + (1-a1-a2-a3)*exp(-t/t4)                      #  
# 5) Ct(t) = a1*exp(-t/t1) + a2*exp(-t/t2) + a3*exp(-t/t3) + a4*exp(-t/t4) + (1-a0-a1-a2-a3-a4)*exp(-t/t5)#   
#                                                                                                         #
# Junchao Xia, junchao-xia\@biomaps.rutgers.edu, BioMaPs institute, Rutgers University.                    #
###########################################################################################################
end_of_myhead
#
print $myhead; 
print "\n";
print " NOTE: Unlike getR1R2fit.pl, this program fits all data files and calculate all parameters related to all vectors. \n";
print " NOTE: A list of file names,\@CFfiles, need to be specified from input for all C2 correlation functions.           \n";
print " NOTE: After modification, it also can be used for other protons besides CH and NH.\n";
print " NOTE: The final fitting results are saved in allfitted.dat. Please check parameters if are reasonable or not.\n";
print " NOTE: Checking files such as *_gnu.log, *_ind.log, *.ps *_fit.dat, and *_prm.dat might help to understand the errors.\n ";
print "\n";
print "\n";

my $usage = <<end_of_usage;
usage:"should input like this [command] [Bfield] [spin type C or N] [nexp] [parameters from initial guess] [tstart] [tend] [inputfile] [tolError]"
end_of_usage

die $usage unless getopts("");

if ($#ARGV != 7) {die $usage;}

my $Bfield  = $ARGV[0];         # B field in T 
my $SpinType = $ARGV[1];        # Spin type N or C
my $nexp = $ARGV[2];            # the number of exponentials to fit
my $parmGuess = $ARGV[3];       # list of initial values for guess
my $tstart = $ARGV[4];          # start time for fitting
my $tend = $ARGV[5];            # end time for fitting
my $inputfile =$ARGV[6];        # input file for the list of data file 
my $tolError = $ARGV[7];        # acceptable averaged error

my $outputfile = "allfitted.dat";
my $maxIniGuess = 50; 
my $maxIniRefine = 500; 

my (@parmAll,@parmToSave);  
my @errorAll; 

my @CFfiles;   # list of data files

my @parmIni = split(/ /,$parmGuess);

print "The inital parameters for guess are @parmIni \n";


 
open(INPUT, "$inputfile");

while(<INPUT>) #read in data
{
    chomp;
    my @tmp_array=split;
    if (defined $tmp_array[0] )
    {
	push (@CFfiles,$tmp_array[0]);
    }
}

close(INPUT);

if (-f "$outputfile" ) {system "rm $outputfile";}


for (my $ifc = 0; $ifc <= $#CFfiles; $ifc++)
{
  my @file_type = split(/\./,$CFfiles[$ifc]); #splitting file to get extension
  my $base = $file_type[0];
  my $extension = $file_type[$#file_type]; 
  my $datafile  = $CFfiles[$ifc];

  print "-------------------------------------------------------------------------------------\n";
  print "Begin to fit the $ifc th data file, $CFfiles[$ifc].\n";

  system "mv $datafile temp.dat";
  system "grep -A 1000000 \"Correlation functions\" < temp.dat > $datafile";

  #initial fitting 
  system "getR1R2fit.pl $Bfield $SpinType $nexp $parmGuess $tstart $tend $datafile >& $base\_ind.log";

  my $ipara; 

  my $initialSuccess=1; 
  my $ninitial = 1; 
 
  while ($initialSuccess > 0 && $ninitial <= $maxIniGuess)
  {
      # input the parameters from the initial fitting
      open(INPUT, "$base\_prm.dat");
      my @prmfitted;
      while(<INPUT>) #read in data
      {
	  chomp;
	  my @tmp_array=split;
	  if (defined $tmp_array[0] ) # check if the initial fitting is scucessful or not
	  {
	      push (@prmfitted,$tmp_array[2]);
	       
	  }
      }
      close(INPUT);

      # check if the initial fitting is successful or not 
      $initialSuccess = 0;
      if (defined $prmfitted[0])
      {
	  for ($ipara = 0; $ipara<=$#prmfitted; $ipara++)
	  {
	      if ($prmfitted[$ipara] < 0) 
	      {$initialSuccess = 1;}
	  }
      }
      else
      {$initialSuccess = 1;}
      
      # initial fitting is not scucessful and do it again
      if ($initialSuccess > 0)  
      {

	  my @parmsNew; 
	  for ($ipara = 0; $ipara <=$#parmIni; $ipara++)
	  {
	      my $temp =  $parmIni[$ipara]+$parmIni[$ipara]*(rand()-0.5);  
	      push(@parmsNew, $temp); 
	  }
	 
	  my $parmsTry = " "; 
	  for ($ipara = 0; $ipara <=$#parmsNew; $ipara++)
	  {
	      $parmsTry = $parmsTry ." " . $parmsNew[$ipara];
	  
	  }
	  system "getR1R2fit.pl $Bfield $SpinType $nexp $parmsTry $tstart $tend $datafile >& $base\_ind.log";
	  $ninitial = $ninitial + 1; 
	  if ($ninitial % 50 == 0)
	  {
	      print  "Doing the initial search of $ninitial \n";
	  }
      }

  }
  

  if ($initialSuccess== 0) 
  {
      print "An initial fitting is done successfully in trials of $ninitial.\n";
      print "Begin to refine the fitting.\n";
  }
  else 
  {
      print "An initial fitting is not successfully in trials of $maxIniGuess.\n";
      print "The initial guess will be tried again during the refinement.\n"; 
      print "Begin to refine the fitting.\n";
  }


  my $needRefine = 1;
  my $avgError;
  my $nrefine  = 1; 
  my $needReguess = 0;
  my $hasNegative = 0;

  while ($needRefine > 0 && $nrefine <= $maxIniRefine )
  {
      
      # input the refitting data
      open(INPUT, "$base\_prm.dat");
      my @prmfitted;
      my @errorfitted;
      while(<INPUT>) #read in data
      {
	  chomp;
	  my @tmp_array=split;
	  if (defined $tmp_array[2] )
	  {
	      push (@prmfitted,$tmp_array[2]) ;
	      push (@errorfitted,$tmp_array[4]) ;
	  }
      }
      close(INPUT);   
   
      $needRefine = 1;
      $needReguess = 0;
      $hasNegative = 0;

      if (defined $prmfitted[0]) 
      {      
          # calculate the averaged total error
	  $avgError = 0.0; 
	  for ($ipara = 0; $ipara <=$#prmfitted; $ipara++)
	  {
	      $avgError = $avgError + 100*abs($errorfitted[$ipara]/$prmfitted[$ipara]);
	  }
	  $avgError = $avgError/($#prmfitted + 1.0);

	  if ($avgError < $tolError ) 
	  { $needRefine = 0;}

	  for ($ipara = 0; $ipara<=$#prmfitted; $ipara++)
	  {
	      if ($prmfitted[$ipara] < 0) 
	      {
		  $needRefine = 1;
		  $needReguess = 0;
		  $hasNegative = 1;
	      }
	  }

	  # save the current fitting parameters and error
	  @parmToSave = @prmfitted;
	  push @parmAll, [@parmToSave];
	  push (@errorAll, $avgError); 
      }
      else
      { $needReguess = 1;}

    

      if ($needRefine > 0) 
      {
	  my @parmsNew; 

	  if ($hasNegative > 0 || $needReguess > 0)
	  # if ($needReguess > 0)
	  {
	      for ($ipara = 0; $ipara <=$#parmIni; $ipara++)
	      {
		  push(@parmsNew, $parmIni[$ipara]+$parmIni[$ipara]*(rand()-0.5));
	      }
	  }
	  else
	  {
	      for ($ipara = 0; $ipara <=$#prmfitted; $ipara++)
	      {
		  push(@parmsNew, $prmfitted[$ipara]+$prmfitted[$ipara]*(rand()-0.5));
	      }
	  }


	  my $parmsTry = " "; 
	  for ($ipara = 0; $ipara <=$#parmsNew; $ipara++)
	  {
	      $parmsTry = $parmsTry . "$parmsNew[$ipara] ";
	  }
	  system "getR1R2fit.pl $Bfield $SpinType $nexp $parmsTry $tstart $tend $datafile >& $base\_ind.log"; 
	  $nrefine = $nrefine + 1;
	  if ($nrefine % 50 == 0)
	  {
	      print  "Doing the refinement of $nrefine \n";
	  }
      }
      
  }

  # if still need refinement (desired error not found), the output the best one upto now.  
  if ($needRefine > 0)
  {
      print "No solution for the desired error was found during the refinement of $maxIniRefine \n"; 
      # put the error in a hash table 
      my %hashE = ();
      for (my $ierr = 0; $ierr <= $#errorAll; $ierr++ )
      {
	  $hashE{$ierr} = $errorAll[$ierr];
      }

      # Get a list of hash keys sorted by value.
      my @sortedID = sort { $hashE{$a} <=> $hashE{$b} } keys %hashE;

      my @parmBest;
      for ($ipara = 0; $ipara <=$#parmIni; $ipara++)
      {
	  push (@parmBest, $parmAll[$sortedID[0]][$ipara]);
      }    

      print "The best one with error of $errorAll[$sortedID[0]] will be the final fitting.\n";

      my $parmsTry = " "; 
      for ($ipara = 0; $ipara <=$#parmBest; $ipara++)
      {
	  $parmsTry = $parmsTry . "$parmBest[$ipara] ";
      }
      system "getR1R2fit.pl $Bfield $SpinType $nexp $parmsTry $tstart $tend $datafile >& $base\_ind.log"; 
  } 
  else 
  {
      print "Finished the refinement in the total number of trials of $nrefine \n"; 
      print "and the final error is $avgError \% \n";
  }

  print "Fitting to the $ifc th data file, $CFfiles[$ifc], has been finished.\n\n";

  if ($ifc == 0)
  { system "tail -n 2  $base\_fit.dat >> $outputfile";}
  else
  { system "tail -n 1  $base\_fit.dat >> $outputfile";}
  system "mv temp.dat $datafile";
  
}

print "Done! Please check the file, $outputfile, for the final fitting results for all data. \n";


