#!/usr/bin/perl -w

require 5.002;
use utf8;
use Getopt::Std;
use warnings;
use strict;

#
#
my $head = <<end_of_head;
-----------------------------------------------------------------------------------------------
|   This is a perl program to integrate the probability within a specified range              |
|   by the Bolzmann average. Only one and two dimensions are considered.                      |
|   Note: a modification needs to make when a minimum is close to the boundary.               | 
|                                                                                             |
|      Junchao Xia, BioMaPS Institute, Rutgers University, junchao-xia\@biomaps.rutgers.edu.   |
-----------------------------------------------------------------------------------------------
end_of_head


print $head;

my $usage = <<end_of_usage;

usage:"Please input as [command] [ndim] [minx] [maxx] [miny] [maxy] [input file for the biased potential] [KT]"
end_of_usage


die $usage unless getopts("");

if($#ARGV==6)
{ 
    &Inter2D($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6]);
}
elsif ($#ARGV==4) 
{ 
    &Inter1D($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4]);
} 
else
{ die $usage;}

sub Inter1D{

    my $ndim=$_[0];
    my $minx=$_[1];
    my $maxx=$_[2];
    my $inputBiasPot=$_[3];

    #my $Kb=*8.31451121E-3;      # Kb in KJ/mol
    #my $Kb=8.31451121E-3/4.184;    # Kb in KCal/mol
    #my $temperature=$_[4];
    #my $temperature = 298.15;    # temperature in K 
    #my $KT = $Kb*$temperature;  
    my $KT=$_[4];

    my @xvalue;
    my @biasPot;
    open(INPUT, $inputBiasPot);
    while(<INPUT>) #read in data
    {
        chomp;
        my @tmp_array=split;

        if (defined $tmp_array[1] )
        {
            push (@xvalue,$tmp_array[0]) ;
            push (@biasPot,$tmp_array[1]) ;
        }
    }
    close INPUT;
    my $probMin = 0;   
    my $partFunc = 0.0;
    my $bolzWeight = 0.0;

    for (my $ipoint = 0; $ipoint <=$#xvalue; $ipoint++)
    {
	#$bolzWeight = exp(-$biasPot[$ipoint]/($Kb*$temperature));
	if ( $biasPot[$ipoint] != 0.0 )
	{
	    $bolzWeight = exp(-$biasPot[$ipoint]/$KT);
	    $partFunc =  $partFunc + $bolzWeight;
            
	    if ($xvalue[$ipoint] <= $maxx && $xvalue[$ipoint]>= $minx) 
	    {
		$probMin = $probMin + $bolzWeight;
	    }
	}

    }
 
    # $probMin = $probMin/$partFunc;
    # print "After the integration, the probability within the range is $probMin \n";

    printf "After the integration, the partial function, the integrator, and the probability within the region is \n
           %20.10e %20.10e %20.10e \n", $partFunc, $probMin, $probMin/$partFunc;
   
  
}

sub Inter2D{
    my $ndim=$_[0];
    my $minx=$_[1];
    my $maxx=$_[2];
    my $miny=$_[3];
    my $maxy=$_[4];
    my $inputBiasPot=$_[5];
    #my $Kb=*8.31451121E-3;        # Kb in KJ/mol
    #my $Kb=8.31451121E-3/4.184;   # Kb in KCal/mol
    #my $temperature=$_[6];
    #my $temperature = 298.15;     # temperature in K 
    #my $KT = $Kb*$temperature;  
    my $KT=$_[6];

    my @xvalue;
    my @yvalue;
    my @biasPot;
    open(INPUT, $inputBiasPot);
    while(<INPUT>) #read in data
    {
        chomp;
        my @tmp_array=split;

        if (defined $tmp_array[2] )
        {
            push (@xvalue,$tmp_array[0]) ;
            push (@yvalue,$tmp_array[1]) ;
            push (@biasPot,$tmp_array[2]) ;
        }
    }
    close INPUT;
    print "Total input data points is ", $#xvalue + 1, "\n";
    my $probMin = 0.0;
    my $partFunc = 0.0;
    my $bolzWeight = 0.0;

    
    for (my $ipoint = 0; $ipoint <=$#xvalue; $ipoint++)
    {
	# $bolzWeight = exp(-$biasPot[$ipoint]/($Kb*$temperature));
        if ( $biasPot[$ipoint] != 0.0 )
	{
	    $bolzWeight = exp(-$biasPot[$ipoint]/$KT);
	    $partFunc =  $partFunc + $bolzWeight;
            
	    if ( ($xvalue[$ipoint] <= $maxx) && ($xvalue[$ipoint]>= $minx) && ($yvalue[$ipoint] <= $maxy) && ($yvalue[$ipoint]>= $miny) ) 
	    {
		$probMin = $probMin + $bolzWeight;
	    }
	}
    }   
    #$probMin = $probMin/$partFunc;
    # print "After the integration, the probability within the range is $probMin \n";
    printf "After the integration, the partial function, the integrator, and the probability within the region is \n
           %20.10e %20.10e %20.10e \n", $partFunc, $probMin, $probMin/$partFunc;
   
}
