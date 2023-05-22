#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;
#
#
my $usage = <<end_of_usage;
usage:"Please input as [command] [input file for name list] [folder prefix] [Index start] [index end] [number of points]"
end_of_usage

die $usage unless getopts("");

if($#ARGV==4)
{
    my $inputfile=$ARGV[0];
    my $folderPref = $ARGV[1];
    my $idStart = $ARGV[2];
    my $idEnd = $ARGV[3];
    my $nPoints = $ARGV[4];

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


    for (my $icf = 0; $icf <= $#CFfiles; $icf++)
    { 
      my @time;
      my @CFsum;
      my @P2sum;
      my @R6sum;

      for (my $id = $idStart; $id <= $idEnd; $id++)
      {
         my $datafile =  $CFfiles[$icf]; 
         system "tail -n $nPoints $folderPref$id/$datafile > $folderPref$id/CorrFunc.temp";
         open(INPUT, "$folderPref$id/CorrFunc.temp");
         my $ipoints = 0;
         while(<INPUT>) #read in data in th CFfile
         {
           chomp;
           my @tmp_array=split;

           if (defined $tmp_array[3] )
           {
	     if ($id == $idStart)
	      {
                push (@time,$tmp_array[0]);
		push (@CFsum,$tmp_array[1]);
		push (@P2sum,$tmp_array[2]);
	        push (@R6sum,$tmp_array[3]);
	      }
             else
	      {
                $CFsum[$ipoints] = $CFsum[$ipoints] + $tmp_array[1];
                $P2sum[$ipoints] = $P2sum[$ipoints] + $tmp_array[2];
                $R6sum[$ipoints] = $R6sum[$ipoints] + $tmp_array[3];
              }
	     $ipoints ++;
           }
         }    
         close(INPUT);     
      }  

      my $nfiles = $idEnd - $idStart + 1; 
      for (my $ipoints =0; $ipoints <= $#time; $ipoints++)	
      {
	$CFsum[$ipoints] = $CFsum[$ipoints]/$nfiles;
        $P2sum[$ipoints] = $P2sum[$ipoints]/$nfiles;
        $R6sum[$ipoints] = $R6sum[$ipoints]/$nfiles;
      }  

      # output the averaged correlation function 
      my $outputfile = $CFfiles[$icf];  
      open(OUTPUT, ">$outputfile");
      printf OUTPUT " ***** Correlation functions ***** \n";
      printf OUTPUT " %12s %12s %12s %12s \n", "Time", "<C>", "<P2>", "<1/(r^3*r^3)>"; 
      for (my $ipoints = 0; $ipoints <= $#time; $ipoints++ )
       {
	   printf OUTPUT " %12.4f %12.4f %12.4f %12.4f \n", $time[$ipoints],$CFsum[$ipoints],$P2sum[$ipoints],$R6sum[$ipoints];

       }
      close(OUTPUT);
    }
	
}
else 
{
    die $usage
};
