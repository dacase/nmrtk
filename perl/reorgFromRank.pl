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
|   This is a perl program to reorganize a data file according to an input rank of one       |  
|   calculated quantity such as NMRs.                                                        | 
|                                                                                            |
|      Junchao Xia, BioMaPS Institute, Rutgers University, junchao-xia\@biomaps.rutgers.edu.  |
----------------------------------------------------------------------------------------------
end_of_head


print $head;
#
#
my $usage = <<end_of_usage;
usage:"Please input as [command] [input file with ranking] [column id of rank] [data file to rank] [output file]"
end_of_usage

    die $usage unless getopts("");

if($#ARGV==3)
{
    my $inputfile=$ARGV[0];
    my $ncolumn=$ARGV[1]-1;
    my $datafile=$ARGV[2];
    my $outputfile=$ARGV[3];
    my  @rankID;

    # input the rank ID from the $ncolumn of $inputfile

    open(INPUT, "$inputfile");
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[$ncolumn] )
	{
	    push (@rankID,$tmp_array[$ncolumn]) ; 
	}
    }
    close(INPUT);

    open(INPUT, "$datafile");
    my @datatoRank = <INPUT>;
    close(INPUT);

    open(OUTPUT, ">$outputfile");
    for my $iconform (0 .. $#rankID)
    {
	print OUTPUT $datatoRank[$rankID[$iconform]];
    }
    close(OUTPUT);
    	
}
else 
{
    die $usage
};
