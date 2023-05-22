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
|   This is a perl program to integrate the probabilitis around energy minima within a radius |
|   by the Bolzmann average. Only one and two dimensions are considered.                      |
|   Note: a modification needs to make when a minimum is close to the boundary.               | 
|                                                                                             |
|      Junchao Xia, BioMaPS Institute, Rutgers University, junchao-xia\@biomaps.rutgers.edu.   |
-----------------------------------------------------------------------------------------------
end_of_head


print $head;

my $usage = <<end_of_usage;

usage:"Please input as [command] [ndim] [input file for minima] [input file for the biased potential] [temperature]"
end_of_usage


die $usage unless getopts("");

if($#ARGV==3)
{

    if ($ARGV[0] == 1)
    {
	&Inter1D($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3]);
    }
    elsif ($ARGV[0] == 2)
    {
	&Inter2D($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3]);
    }
    else
    {
     print "Sorry. ndim does not have a right value!";  
    }

}
else 
{ die $usage;}

sub Inter1D{

    my $ndim = $_[0];
    my $inputMinima=$_[1];
    my $inputBiasPot=$_[2];
    my $temperature=$_[3];
    #  my $Kb=*8.31451121E-3;      # Kb in KJ/mol
    my $Kb=8.31451121E-3/4.184;    # Kb in KCal/mol


    # input the positions of minima and areas to calculate   
    open(INPUT, "$inputMinima");
    my @xvalueMin; 
    my @rvalueMin; 
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[0] )
	{
	    push (@xvalueMin,$tmp_array[0]) ;
 	    push (@rvalueMin,$tmp_array[1]) ;
 
	}
      
    }
    close(INPUT);

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
    my @probMin;
    for (my $imin = 0; $imin <= $#xvalueMin; $imin++)
    {
	push (@probMin,0.0) ;
    }

    my $partFunc = 0.0;
    my $bolzWeight = 0.0;

    for (my $ipoint = 0; $ipoint <=$#xvalue; $ipoint++)
    {
	$bolzWeight = exp(-$biasPot[$ipoint]/($Kb*$temperature));
	$partFunc =  $partFunc + $bolzWeight;

	for (my $imin = 0; $imin <=$#xvalueMin; $imin++)
	{
           my $rdistance = abs($xvalue[$ipoint]-$xvalueMin[$imin]);
            
	    if ($rdistance <= $rvalueMin[$imin]) 
	   {
	      $probMin[$imin] = $probMin[$imin] + $bolzWeight;
           }
	}
    }
    for (my $imin = 0; $imin <=$#xvalueMin; $imin++)
    {
	$probMin[$imin] = $probMin[$imin]/$partFunc;
    }
    print "After the integration, the probabilities around minima are ";
    for (my $imin = 0; $imin <=$#xvalueMin; $imin++)
    {
	print " $probMin[$imin]";
    }
    print "\n";
}

sub Inter2D{
    my $ndim = $_[0];
    my $inputMinima=$_[1];
    my $inputBiasPot=$_[2];
    my $temperature=$_[3];

    #  my $Kb=*8.31451121E-3;      # Kb in KJ/mol
    my $Kb=8.31451121E-3/4.184;    # Kb in KCal/mol


    # input the positions of minima and areas to calculate   
    open(INPUT, "$inputMinima");
    my @xvalueMin; 
    my @yvalueMin;
    my @rvalueMin; 
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[0] )
	{
	    push (@xvalueMin,$tmp_array[0]) ;
 	    push (@yvalueMin,$tmp_array[1]) ;
	    push (@rvalueMin,$tmp_array[2]) ;
	}
      
    }
    close(INPUT);

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
    my @probMin;
    for (my $imin = 0; $imin <= $#xvalueMin; $imin++)
    {
	push (@probMin,0.0) ;
    }

    my $partFunc = 0.0;
    my $bolzWeight = 0.0;

    for (my $ipoint = 0; $ipoint <=$#xvalue; $ipoint++)
    {
	$bolzWeight = exp(-$biasPot[$ipoint]/($Kb*$temperature));
	$partFunc =  $partFunc + $bolzWeight;

	for (my $imin = 0; $imin <=$#xvalueMin; $imin++)
	{
            my $rdistance = sqrt( ($xvalue[$ipoint]-$xvalueMin[$imin])*($xvalue[$ipoint]-$xvalueMin[$imin]) 
				  + ($yvalue[$ipoint]-$yvalueMin[$imin])*($yvalue[$ipoint]-$yvalueMin[$imin]));
            
	    if ($rdistance <= $rvalueMin[$imin]) 
	   {
	      $probMin[$imin] = $probMin[$imin] + $bolzWeight;
           }
	}
    }
    for (my $imin = 0; $imin <=$#xvalueMin; $imin++)
    {
	$probMin[$imin] = $probMin[$imin]/$partFunc;
    }
    print "After the integration, the probabilities around minima are ";
    for (my $imin = 0; $imin <=$#xvalueMin; $imin++)
    {
	print " $probMin[$imin]";
    }
    print "\n";
    
}
