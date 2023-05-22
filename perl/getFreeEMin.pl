#!/usr/bin/perl -w
require 5.002;
use Getopt::Std;
use warnings;
use strict; 

#
#
my $head = <<end_of_head;
----------------------------------------------------------------------------------------------
|   This is a perl program to find and shift to zero for the global minimum of free energy   |
|   in one dimension or two dimensions.                                                      | 
|                                                                                            |
|      Junchao Xia, BioMaPS Institute, Rutgers University, junchao-xia\@biomaps.rutgers.edu.  |
----------------------------------------------------------------------------------------------
end_of_head

print $head;

my $usage = <<end_of_usage;

usage:
  1D: getFreeEMin.pl nx inputfile outputfile
  2D: getFreeEMin.pl nx ny inputfile outputfile 
end_of_usage
    
die $usage unless getopts("");

if($#ARGV==3)
{
    &freeE2d($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3]);
}
elsif($#ARGV==2)
{
    &freeE1d($ARGV[0],$ARGV[1],$ARGV[2]);
}
else
{ die $usage;}
##############################################################
sub freeE1d {

    my $nx=$_[0]; #number of divisions  
    my $inputfile = $_[1];
    my $outputfile = $_[2];
    my @array1;
    my @array2;

    open(INPUT, $inputfile);
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;

	push (@array1,$tmp_array[0]) ;
	push (@array2,$tmp_array[1]) ;
    }
    close INPUT;

    my $min=$array2[0];
    my $imin= 0;

    my $counter=0;

    foreach (@array2) #searching for minimum and maximum elements
    {

	if($array2[$counter]<$min)
	{
	    $min=$array2[$counter];
	    $imin = $counter;
	}
	 ++$counter;
    }
    print " The minimum is $array1[$imin], $array2[$imin] \n";
    
    open(OUTPUT, ">$outputfile");
    for $counter (0 .. $nx-1)
    {
	my $newvalue = $array2[$counter]- $array2[$imin];
	print OUTPUT "$array1[$counter]  $newvalue \n"
    }
    close OUTPUT;

}
#########################################################################
sub freeE2d {
    
    my $nx=$_[0]; #number of divisions
    my $ny=$_[1]; #number of divisions
    my $inputfile = $_[2];
    my $outputfile = $_[3];
    my @array1;
    my @array2;
    my @array3;

    open(INPUT, $inputfile);
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;

	if (defined $tmp_array[2] )
	{
	    push (@array1,$tmp_array[0]) ;
	    push (@array2,$tmp_array[1]) ;
	    push (@array3,$tmp_array[2]) ;
	}
    }	
    close INPUT;


    my $min=$array3[0];
    my $imin = 0;
 
    my $counter=0;

    foreach (@array3) #searching for minimum and maximum elements
    {

	if($array3[$counter]<$min)
	{
	    $min=$array3[$counter];
	    $imin = $counter;
	}
	 ++$counter;
    }
    print " The minimum is $array1[$imin], $array2[$imin], $array3[$imin] \n";

    $counter=0;

    open(OUTPUT, ">$outputfile");
	
    for my $I (1 .. $nx) #initializing the array
    {	
	for my $J (1 .. $ny)
	{   
	    my $newvalue = $array3[$counter] - $array3[$imin];
	    printf OUTPUT "$array1[$counter]  $array2[$counter]  $newvalue \n";
	     ++$counter; 
	}    
	print OUTPUT "  \n";
 	
    }
    close OUTPUT;
}
