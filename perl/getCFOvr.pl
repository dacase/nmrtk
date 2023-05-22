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
| This is a perl program to extract the correlation function of overall rotation from the     | 
| total and the internal ones.                                                                |   
|                                                                                             |
|      Junchao Xia, BioMaPS Institute, Rutgers University, junchao-xia\@biomaps.rutgers.edu.   |
***********************************************************************************************
end_of_head
print $head;
#
#
my $usage = <<end_of_usage;
usage:"Please input as [command] [inputfile of total function] [columnID of time] [columnID of total] [inputfile of internal] [columnID of internal] [output file of overall]"
end_of_usage

    die $usage unless getopts("");

if($#ARGV==5)
{
    my $infileTot=$ARGV[0];
    my $colIDTime=$ARGV[1]-1;
    my $colIDTot=$ARGV[2]-1;
    my $infileInt=$ARGV[3];
    my $colIDInt=$ARGV[4]-1;
    my $outfileOvr=$ARGV[5];
    my @time;
    my @total;
    my @internal;
    my @overall; 

    # input the total function from the column, $colIDtot, of file, $infileTot.
    system "mv $infileTot tempTot.dat";
    system "grep -A 1000000 \"0.000\" < tempTot.dat > $infileTot";
    open(INPUT, "$infileTot");
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[$colIDTot] )
	{
	    push (@time,$tmp_array[$colIDTime]);
	    push (@total,$tmp_array[$colIDTot]); 
	}
    }
    close(INPUT);

    # input the total function from the column, $colIDInt, of file, $infileInt.
    system "mv $infileInt tempInt.dat";
    system "grep -A 1000000 \"0.000\" < tempInt.dat > $infileInt";
    open(INPUT, "$infileInt");
    while(<INPUT>) #read in data
    {
        chomp;
        my @tmp_array=split;
        if (defined $tmp_array[$colIDInt] )
        {
            push (@internal,$tmp_array[$colIDInt]);         
        }
    }
    close(INPUT);

    open(OUTPUT, ">$outfileOvr");
    for my $itime (0 .. $#time)
    {
	printf OUTPUT " %20.5f %20.10f \n", $time[$itime], $total[$itime]/$internal[$itime];
    }
    close(OUTPUT);
    	
    print " done! \n"
}
else 
{
    die $usage
};
