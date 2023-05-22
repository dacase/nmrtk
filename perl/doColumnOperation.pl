#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;
#
#
my $head = <<end_of_head;
***********************************************************************************************
| This is a perl program to perform column operation from two data file.                      |    
|                                                                                             |
|      Junchao Xia, BioMaPS Institute, Rutgers University, junchao-xia\@biomaps.rutgers.edu.   |
***********************************************************************************************
end_of_head
print $head;
#
#
my $usage = <<end_of_usage;
usage:"Please input as [command] [inputfile 1] [columnID of x] [columnID of y1] [inputfile 2] [columnID of y2] [output file of y] [operator]"
end_of_usage

    die $usage unless getopts("");

if($#ARGV==6)
{
    my $infile1=$ARGV[0];
    my $colIDX=$ARGV[1]-1;
    my $colIDY1=$ARGV[2]-1;
    my $infile2=$ARGV[3];
    my $colIDY2=$ARGV[4]-1;
    my $outfile=$ARGV[5];
    my $myOperator = $ARGV[6]; 
    my @x;
    my @y1;
    my @y2; 

    # input the Y1 from the column, $colIDY1, of file, $infile1.
    open(INPUT, "$infile1");
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[$colIDY1] )
	{
	    push (@x,$tmp_array[$colIDX]);
	    push (@y1,$tmp_array[$colIDY1]); 
	}
    }
    close(INPUT);

    # input the total function from the column, $colIDInt, of file, $infileInt.
    open(INPUT, "$infile2");
    while(<INPUT>) #read in data
    {
        chomp;
        my @tmp_array=split;
        if (defined $tmp_array[$colIDY2] )
        {
            push (@y2,$tmp_array[$colIDY2]);         
        }
    }
    close(INPUT);

    open(OUTPUT, ">$outfile");
    for my $ix (0 .. $#x)
    {
	if ( $myOperator eq "+" ) 
	{printf OUTPUT " %20.5f %20.10f \n", $x[$ix], $y1[$ix]+$y2[$ix];}
	elsif ( $myOperator eq "-" )
	{printf OUTPUT " %20.5f %20.10f \n", $x[$ix], $y1[$ix]-$y2[$ix];}
	elsif ( $myOperator eq "*" )
	{printf OUTPUT " %20.5f %20.10f \n", $x[$ix], $y1[$ix]*$y2[$ix];}
	elsif ( $myOperator eq "/" )
	{printf OUTPUT " %20.5f %20.10f \n", $x[$ix], $y1[$ix]/$y2[$ix];}
	elsif ( $myOperator eq "=" ) 
	{printf OUTPUT " %20.5f %20.10f %20.10f \n", $x[$ix], $y1[$ix], $y2[$ix];}
	else
        {print "Operation, $myOperator, not defined yet. \n"}
	
    }
    close(OUTPUT);
    	
    print " done! \n"
}
else 
{
    die $usage
};
