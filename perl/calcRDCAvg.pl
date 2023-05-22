#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;
#
#
my $usage = <<end_of_usage;
usage:"Please input as [command] [input file for Experimental RDCs] [number of RDCs] [input file for RDCs calculated] [number of conformers]"
end_of_usage

    die $usage unless getopts("");

if($#ARGV==3)
{
    my $inputExpRDCs=$ARGV[0];
    my $nrdc=$ARGV[1];
    my $inputCalRDCs=$ARGV[2];
    my $nconform=$ARGV[3];

    my (@RDCs,@rdc);
    my (@RDCsAvg,@rdcAvg);
    my  @rmsd;
    my  @rfactor;
    my  @nfactor;
    my  @exp_rdc;

    # input the experimental RDCs  
    open(INPUT, "$inputExpRDCs");
    @exp_rdc = ();
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[0] )
	{
	    push (@exp_rdc,$tmp_array[0]) ; 
	}
    }
    close(INPUT);
  
    open(INPUT, "$inputCalRDCs");
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	@rdc = ();
        for (my $irdc = 0; $irdc < $nrdc; $irdc++ )
	{
	  push (@rdc,$tmp_array[$irdc+4]);
        }
        push @RDCs,[@rdc];
       
    }
    close(INPUT);

    my @rdcSum; 
    @rdcAvg = ();
    for (my $irdc = 0; $irdc < $nrdc; $irdc++ )
    {
	push (@rdcSum,0.0);
	push (@rdcAvg,0.0);
    }

    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
    {  @rdc = ();

	for (my $irdc = 0; $irdc < $nrdc; $irdc++ )
	{
	   $rdcSum[$irdc] = $rdcSum[$irdc] + $RDCs[$iconform][$irdc];
	}

	for (my $irdc = 0; $irdc < $nrdc; $irdc++ )
	{
	   $rdcAvg[$irdc] = $rdcSum[$irdc]/($iconform+1);
	}
	push @RDCsAvg,[@rdcAvg];

        # calculate the RDC RMSD and R factor 
	my $rmsd_numerator = 0.0;
	my $rmsd_denominator = 0.0;
	my $nfactor_numerator = 0.0;
	for (my $irdc = 0; $irdc < $nrdc;  $irdc++)
	{
	    $rmsd_numerator =  $rmsd_numerator + ($rdcAvg[$irdc] -  $exp_rdc[$irdc])*($rdcAvg[$irdc] - $exp_rdc[$irdc]);
	    $rmsd_denominator = $rmsd_denominator +  $exp_rdc[$irdc]*$exp_rdc[$irdc];
	    $nfactor_numerator =  $nfactor_numerator + ($rdcAvg[$irdc] -  $exp_rdc[$irdc])*($rdcAvg[$irdc] - $exp_rdc[$irdc])/($exp_rdc[$irdc]*$exp_rdc[$irdc]);
	}
	
	push (@rmsd,sqrt($rmsd_numerator/$nrdc)); 
	push (@rfactor,sqrt($rmsd_numerator/$rmsd_denominator));
	push (@nfactor,sqrt($nfactor_numerator/$nrdc));	
    }

    # output all alignment tensors  
    open(OUTPUT, ">RfactorsAll.out");
    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
       { 
	   printf OUTPUT " %10d ", $iconform;

	   printf OUTPUT "%12.6f %12.6f  %12.6f \n",$rmsd[$iconform],$rfactor[$iconform],$nfactor[$iconform];

       }
    close(OUTPUT);

    # output all RDCs 
    open(OUTPUT, ">RDCsAvgAll.out");
    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
       { 
	   printf OUTPUT " %10d ", $iconform;

	   for (my $irdc = 0; $irdc < $nrdc;  $irdc++)
	   {
	       printf OUTPUT "%10.5f ",$RDCsAvg[$iconform][$irdc];
	   }
	   printf OUTPUT "\n";
       }
    close(OUTPUT);

 
	
}
else 
{
    die $usage
};
