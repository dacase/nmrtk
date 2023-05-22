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
|   This is a perl program to extract ensemble structures with populations greater than      |
|   a predefined cutoff value.                                                               |  
|                                                                                            | 
|                                                                                            |
|      Junchao Xia, BioMaPS Institute, Rutgers University, junchao-xia\@biomaps.rutgers.edu.  |
----------------------------------------------------------------------------------------------
end_of_head


print $head;
#
#
my $usage = <<end_of_usage;
usage:"Please input as [command] [datafile to extract] [population file] [population column] [pop cutoff] [output file]"
end_of_usage

    die $usage unless getopts("");

if($#ARGV==4)
{
    my $datafile=$ARGV[0];
    my $popfile=$ARGV[1];
    my $ncolumn=$ARGV[2]-1;
    my $pcutoff = $ARGV[3];
    my $outputfile=$ARGV[4];
    my @popValues;

    # input the rank ID from the $ncolumn of $inputfile

    open(INPUT, "$popfile");
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[$ncolumn] )
	{
	    push (@popValues,$tmp_array[$ncolumn]) ; 
	}
    }
    close(INPUT);

    open(INPUT, "$datafile");
    my @datatoExtract = <INPUT>;
    close(INPUT);

    open(OUTPUT, ">$outputfile");
    for my $iconform (0 .. $#popValues)
    {
	if ($popValues[$iconform] > $pcutoff) 
	{
	    print " The conformer $iconform has population, $popValues[$iconform],  greater the cutoff.\n";
	    print OUTPUT $datatoExtract[$iconform];
	} 
	    
    }
    close(OUTPUT);
    	
}
else 
{
    die $usage
};
