#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;
#
#
my $head = <<end_of_head;
----------------------------------------------------------------------------------------------
|   This is a perl program to extract experimental data for unique conformation.             |  
|                                                                                            |
|      Junchao Xia, BioMaPS Institute, Rutgers University, junchao-xia\@biomaps.rutgers.edu.  |
----------------------------------------------------------------------------------------------
end_of_head


print $head;
#
#
my $usage = <<end_of_usage;
usage:"Please input as [command] [input file for unique confs] [column # of unique confs] [data file to extract] [output file]"
end_of_usage

    die $usage unless getopts("");

if($#ARGV==3)
{
    my $inputfile=$ARGV[0];
    my $ncolumn=$ARGV[1]-1;
    my $datafile=$ARGV[2];
    my $outputfile=$ARGV[3];
    my @rankID;
    my @rankVal;

    # input the rank ID from the $ncolumn of $inputfile

    open(INPUT, "$inputfile");
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[$ncolumn] )
	{
	    push (@rankID,$tmp_array[$ncolumn]) ;
            push (@rankVal,$tmp_array[$ncolumn+1]); 
	}
    }
    close(INPUT);

    open(INPUT, "$datafile");
    my @datatoRank = <INPUT>;
    close(INPUT);

    open(OUTPUT, ">$outputfile");
    for my $iconform (0 .. $#rankID)
    {
	printf OUTPUT " %10d %10d %9.4f %s ",  $iconform, $rankID[$iconform],  $rankVal[$iconform], $datatoRank[$rankID[$iconform]];
    }
    close(OUTPUT);
    	
}
else 
{
    die $usage
};
