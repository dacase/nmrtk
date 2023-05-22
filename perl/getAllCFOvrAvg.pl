#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;
my $myhead = <<end_of_myhead;
############################################################################################################
# This a perl program to calculate the averaged correlation function of overall motion from all vectors.   #
#                                                                                                          #
# Junchao Xia, junchao-xia\@biomaps.rutgers.edu, BioMaPs institute, Rutgers University.                     #
############################################################################################################
end_of_myhead
print $myhead;

#
#
my $usage = <<end_of_usage;
usage:"Please input as [command] [input file for name list] [number of points]"
end_of_usage

die $usage unless getopts("");

if($#ARGV==1)
{
    my $inputfile=$ARGV[0];
    my $nPoints = $ARGV[1];

    my @CFfiles; 

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

    my @time;
    my @CFsum;
    for (my $icf = 0; $icf <= $#CFfiles; $icf++)
    { 
 
      my $datafile =  $CFfiles[$icf]; 
      system "tail -n $nPoints $datafile > CorrFunc.temp";
      open(INPUT, "CorrFunc.temp");
      my $ipoints = 0;
      while(<INPUT>) #read in data in th CFfile
      {
	  chomp;
	  my @tmp_array=split;
	  
	  if (defined $tmp_array[1] )
	  {
	      if ($icf == 0)
	      {
		  push (@time,$tmp_array[0]);
		  push (@CFsum,$tmp_array[1]);
	      }
	      else
	      {
		  $CFsum[$ipoints] = $CFsum[$ipoints] + $tmp_array[1];
		  
              }
	      $ipoints ++;
	  }
      }    
      close(INPUT);
       
    }

    my $nfiles = $#CFfiles + 1; 
    for (my $ipoints =0; $ipoints <= $#time; $ipoints++)	
    {
	$CFsum[$ipoints] = $CFsum[$ipoints]/$nfiles;
    }  

    # output the averaged correlation function 
    my $outputfile = "CFOvrAvg.dat";  
    open(OUTPUT, ">$outputfile"); 
    for (my $ipoints = 0; $ipoints <= $#time; $ipoints++ )
    {
	printf OUTPUT " %12.4f %12.4f \n", $time[$ipoints],$CFsum[$ipoints];
	
    }
    close(OUTPUT);
	
}
else 
{
    die $usage
};
