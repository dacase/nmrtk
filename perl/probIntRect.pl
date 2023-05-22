#!/usr/bin/perl -w

require 5.002;
use utf8;
use Getopt::Std;
use warnings;
use strict;

#
#
my $head = <<end_of_head;
-------------------------------------------------------------------------------------------------
|   This is a perl program to obtain the populations around energy minima within a              |
|   specified rectangle using the Bolzmann average. Only one and two dimensions cases           |
|   are implemented.                                                                            |
|                                                                                               |
|   Note: a modification needs to make when a minimum is close to the periodic boundary.        | 
|                                                                                               |
|      Junchao Xia, BioMaPS Institute, Rutgers University, junchao-xia\@biomaps.rutgers.edu.     |
-------------------------------------------------------------------------------------------------
end_of_head


print $head;

my $usage = <<end_of_usage;

usage:"Please input as [command] [ndim] [input file for minima] [input file for the biased potential] [KT]"
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

    #my $Kb=*8.31451121E-3;      # Kb in KJ/mol
    #my $Kb=8.31451121E-3/4.184;    # Kb in KCal/mol
    #my $temperature=$_[3];
    #my $temperature = 298.15;    # temperature in K 
    #my $KT = $Kb*$temperature;  
    my $KT=$_[3];

    # input the specified rectangle region to calculate    
    open(INPUT, "$inputMinima");
    my @xvalueMin; 
    my @xvalueMax; 
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[0] )
	{
	    push (@xvalueMin,$tmp_array[0]) ;
 	    push (@xvalueMax,$tmp_array[1]) ;
 
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
	# $bolzWeight = exp(-$biasPot[$ipoint]/($Kb*$temperature));
	$bolzWeight = exp(-$biasPot[$ipoint]/$KT);

	$partFunc =  $partFunc + $bolzWeight;

	for (my $imin = 0; $imin <=$#xvalueMin; $imin++)
	{
            
	   if ($xvalue[$ipoint] <= $xvalueMax[$imin] && $xvalue[$ipoint]>= $xvalueMin[$imin] )

	   {
	      $probMin[$imin] = $probMin[$imin] + $bolzWeight;
           }
	}
    }
    for (my $imin = 0; $imin <=$#xvalueMin; $imin++)
    {
	$probMin[$imin] = $probMin[$imin]/$partFunc;
    }
    open(OUTPUT, ">pops.out");
    printf "After the integration, the populations around minima are \n ";
    for (my $imin = 0; $imin <=$#xvalueMin; $imin++)
    {
	printf " %20.10e", $probMin[$imin];
	printf OUTPUT " %20.10e", $probMin[$imin];
    }
    printf "\n";
    printf "The populations are also saved to the file, pops.out.\n";
    printf OUTPUT "\n";
    close(OUTPUT);

}

sub Inter2D{
    my $ndim = $_[0];
    my $inputMinima=$_[1];
    my $inputBiasPot=$_[2];

    #my $Kb=*8.31451121E-3;        # Kb in KJ/mol
    #my $Kb=8.31451121E-3/4.184;   # Kb in KCal/mol
    #my $temperature=$_[3];
    #my $temperature = 298.15;     # temperature in K 
    #my $KT = $Kb*$temperature;  
    my $KT=$_[3];


    # input the specified rectangle region to calculate   
    open(INPUT, "$inputMinima");
    my @xvalueMin; 
    my @xvalueMax; 
    my @yvalueMin; 
    my @yvalueMax; 

    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[0] )
	{
	    push (@xvalueMin,$tmp_array[0]) ;
 	    push (@xvalueMax,$tmp_array[1]) ;
	    push (@yvalueMin,$tmp_array[2]) ;
	    push (@yvalueMax,$tmp_array[3]) ;
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
	# $bolzWeight = exp(-$biasPot[$ipoint]/($Kb*$temperature));
	$bolzWeight = exp(-$biasPot[$ipoint]/$KT);
	$partFunc =  $partFunc + $bolzWeight;

	for (my $imin = 0; $imin <=$#xvalueMin; $imin++)
	{
       
            
	   if ($xvalue[$ipoint] <= $xvalueMax[$imin] && $xvalue[$ipoint]>= $xvalueMin[$imin] && $yvalue[$ipoint] <= $yvalueMax[$imin] && $yvalue[$ipoint]>= $yvalueMin[$imin])
	   {
	      $probMin[$imin] = $probMin[$imin] + $bolzWeight;
           }
	}
    }
    for (my $imin = 0; $imin <=$#xvalueMin; $imin++)
    {
	$probMin[$imin] = $probMin[$imin]/$partFunc;
    }
    open(OUTPUT, ">pops.out");
    printf "After the integration, the populations around minima are \n ";
    for (my $imin = 0; $imin <=$#xvalueMin; $imin++)
    {
	printf " %20.10e", $probMin[$imin];
	printf OUTPUT " %20.10e", $probMin[$imin];
    }
    printf "\n";
    printf "The populations are also saved to the file, pops.out.\n";
    printf OUTPUT "\n";
    close(OUTPUT);
    
}
