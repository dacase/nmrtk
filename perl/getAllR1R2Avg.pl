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
      my @Resid;
      my @R1;
      my @R2;
      my @R2dR1;
      my @NOE;

      for (my $id = $idStart; $id <= $idEnd; $id++)
      {
         my $datafile =  $CFfiles[$icf]; 
         open(INPUT, "$folderPref$id/$datafile");
         my $ipoints = 0;
         while(<INPUT>) #read in data in th CFfile
         {
           chomp;
           my @tmp_array=split;

           if (defined $tmp_array[3] )
           {
	     if ($id == $idStart)
	      {
                push (@Resid,$tmp_array[0]);
		push (@R1,$tmp_array[1]);
		push (@R2,$tmp_array[2]);
	        push (@R2dR1,$tmp_array[2]/$tmp_array[1]);
	        push (@NOE,$tmp_array[3]);
	      }
             else
	      {
                $R1[$ipoints] = $R1[$ipoints] + $tmp_array[1];
                $R2[$ipoints] = $R2[$ipoints] + $tmp_array[2];
                $R2dR1[$ipoints] = $R2dR1[$ipoints] + $tmp_array[2]/$tmp_array[1];
		$NOE[$ipoints] = $NOE[$ipoints] + $tmp_array[3];
              }
	     $ipoints ++;
           }
         }    
         close(INPUT);     
      }  

      my $nfiles = $idEnd - $idStart + 1; 
      for (my $ipoints =0; $ipoints <= $#Resid; $ipoints++)	
      {
	$R1[$ipoints] = $R1[$ipoints]/$nfiles;
        $R2[$ipoints] = $R2[$ipoints]/$nfiles;
        $R2dR1[$ipoints] = $R2dR1[$ipoints]/$nfiles;
	$NOE[$ipoints] = $NOE[$ipoints]/$nfiles;
      }  

      # output the averaged quantities
      my $outputfile = $CFfiles[$icf];  
      open(OUTPUT, ">avg$outputfile");
      printf OUTPUT " ***** Averaged relaxazation parameters from multi-trajectories ***** \n";
      printf OUTPUT " %12s %12s %12s %12s %12s  \n", "ResID", "<R1>", "<R2>", "<R2/R1>", "<NOE>"; 
      for (my $ipoints = 0; $ipoints <= $#Resid; $ipoints++ )
       {
	   printf OUTPUT " %12d %12.4f %12.4f %12.4f %12.4f \n", $Resid[$ipoints],$R1[$ipoints],$R2[$ipoints],$R2dR1[$ipoints],$NOE[$ipoints];

       }
      close(OUTPUT);
    }
	
}
else 
{
    die $usage
};
