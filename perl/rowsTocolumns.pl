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
| This is a perl program to transform the rows to  columns for a txt data file                |    
|                                                                                             |
|      Junchao Xia, BioMaPS Institute, Rutgers University, junchao-xia\@biomaps.rutgers.edu.   |
***********************************************************************************************
end_of_head
print $head;
#
#
my $usage = <<end_of_usage;
usage:"Please input as [command] [inputfile] [output file]"
end_of_usage

    die $usage unless getopts("");

if($#ARGV==1)
{
    my $infile=$ARGV[0];
    my $outfile=$ARGV[1];

    my @rows=();
    my $oldCols; 
    open(INPUT, "<$infile");
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
        push @rows,[@tmp_array];
        $oldCols = $#tmp_array;
    }
    close(INPUT);

    open(OUTPUT, ">$outfile");
    my $i;
    my $j;
    for $j (0 .. $oldCols)
    {
      for $i (0 .. $#rows)
      { printf OUTPUT " %20.5f", $rows[$i][$j];}
      printf OUTPUT " \n";
    }
    close(OUTPUT);
    	
    print " done! \n"
}
else 
{
    die $usage
};
