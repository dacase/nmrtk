#!/usr/bin/perl -w
require 5.002;
use Getopt::Std;
use warnings;
use strict;

#
#
my $usage = <<end_of_usage;
usage:"should input like this [command] [nx|nx ny] [nlines] [input file] [output file] "
end_of_usage
    
    
die $usage unless getopts("");


if($#ARGV==4)
{
    &umbrella2d($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4]);
}
elsif($#ARGV==3)
{
    &umbrella1d($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3]);
}
else
{ die $usage;}
##############################################################
sub Umbrella1d {

    my $nx=$_[0]; #number of divisions  
    my $nlines=$_[1]; #number of lines  
    my $inputfile = $_[2];
    my $outputfile = $_[3];

    my $nhead = $nlines-1;

    system "tail -n $nlines $inputfile | head -n $nhead > umbrella.temp";
    system "myreplace.pl \"coeffs =\"  \" \" < umbrella.temp > umbrella.temp1";
    system "myreplace.pl \"\;\" \" \" < umbrella.temp1 > umbrella.temp2";
    system "cp umbrella.temp2 umbrella.in";
    system "rm -f umbrella.temp2 umbrella.temp1 umbrella.temp"; 

    my @coeffs;

    open(INPUT, "<umbrella.temp");
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	for my $i (0 .. $#tmp_array)
	{
	  #push (@coeffs,$tmp_array[$i]);  
	  if (defined $tmp_array[$i] )
   	  {
	    push (@coeffs,$tmp_array[$i]);
    	  }
        }
    }
    close INPUT;

    my $max1=180;
    my $min1=-180;
    my $dx1=($max1-$min1)/$nx;

    open(OUTPUT, ">$outputfile");

    for my $counter (0 .. $nx-1)
    {
	my $xval=$min1+$dx1*($counter+1)-$dx1/2;
 
	print OUTPUT " $xval  $coeffs[$counter] \n";
    }
    close OUTPUT;
    system "myreplace.pl \",\" \" \" < $outputfile > biasPotEng.temp ";
    system "mv biasPotEng.temp $outputfile";
}
#########################################################################
sub umbrella2d {
    

    my $nx=$_[0]; #number of divisions
    my $ny=$_[1]; #number of divisions
    my $nlines=$_[2]; #number of lines  
    my $inputfile = $_[3];
    my $outputfile = $_[4];

    my $nhead = $nlines-1;

    system "tail -n $nlines $inputfile | head -n $nhead > umbrella.temp";
    system "myreplace.pl \"coeffs =\"  \" \" < umbrella.temp > umbrella.temp1";
    system "myreplace.pl \"\;\" \" \" < umbrella.temp1 > umbrella.temp2";
    system "cp umbrella.temp2 umbrella.in";
    system "rm -f umbrella.temp2 umbrella.temp1 umbrella.temp"; 
    
    my @coeffs;

    open(INPUT, "umbrella.in");
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;

	for my $i (0 .. $#tmp_array)
	{	 
	  #push (@coeffs,$tmp_array[$i]);
    	 
	  if (defined $tmp_array[$i] )
   	  {
	    push (@coeffs,$tmp_array[$i]);
    	  }
        }
    }	
    close INPUT;
 
    # open(OUTPUT, ">test.temp");
    # for my $i (0 .. $#coeffs)
    # {
    # 	printf OUTPUT "$i  $coeffs[$i] \n";	
    # }
    # close OUTPUT;

    my $max1=180;
    my $min1=-180;
    my $max2=180;
    my $min2=-180;
    my $dx1=($max1-$min1)/$nx;
    my $dx2=($max2-$min2)/$ny;

    my $counter=0;

    open(OUTPUT, ">$outputfile");
	
    for my $I (1 .. $nx) #initializing the array
    {	
	my $xval=$min1+$dx1*($I)-$dx1/2;

	for my $J (1 .. $ny)
	{   
	    my $yval=$min2+$dx2*($J)-$dx2/2;
	    printf OUTPUT "$xval  $yval $coeffs[$counter] \n";
	    ++$counter;
	    
	}    

	print OUTPUT "  \n";

    }
    close OUTPUT;
    system "myreplace.pl \",\" \" \" < $outputfile > biasPotEng.temp ";
    system "mv biasPotEng.temp $outputfile";
}
