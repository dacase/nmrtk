#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;

#
my $myhead = <<end_of_myhead;
###########################################################################################################
# This a perl program to extract serval columns for multiple data files and put them into a single file.  #
#                                                                                                         #
# Junchao Xia, junchao-xia\@biomaps.rutgers.edu, BioMaPs institute, Rutgers University.                    #
###########################################################################################################
end_of_myhead
#
print $myhead; 
print "\n";
print "\n";

my $usage = <<end_of_usage;
usage:"should input like this [command] [input file for fimename list] [output file]"
end_of_usage

die $usage unless getopts("");

if ($#ARGV != 1) {die $usage;}

my $inputfile =$ARGV[0];        # input file for the list of data file 
my $outputfile = $ARGV[1];      # output file
 
 
my @fileList;   # list of data files

my (@colLists,@cols);   # column lists for all data files 
my (@allData,@colData); # save all data


my $ifc, 
my $ival;
my $icol; 
 
open(INPUT, "$inputfile");

while(<INPUT>) #read in data
{
    chomp;
    my @tmp_array=split;
    if (defined $tmp_array[0] )
    {
	push (@fileList,$tmp_array[0]);
	@cols = ();
	for ($icol = 1; $icol <= $#tmp_array; $icol++)
	{
	    push (@cols, $tmp_array[$icol]-1);
	}
	push @colLists, [@cols];
    }
}

close(INPUT);

if (-f "$outputfile" ) {system "rm $outputfile";}


for ($ifc = 0; $ifc <= $#fileList; $ifc++)
{
    #@cols = ();
    #@cols = $colLists[$ifc];
    for ($icol = 0; $icol <= $#cols; $icol++)
    {
        my $colname = $colLists[$ifc][$icol] + 1; 
	print " Get the column $colname from the file $fileList[$ifc]. \n";
	@colData = ();
	my $filename = $fileList[$ifc];
	open(INPUT, "$filename");
	while(<INPUT>) #read in data
	{
	    chomp;
	    my @tmp_array=split;
	    push (@colData, $tmp_array[$colLists[$ifc][$icol]]);
	} 
	close(INPUT);
	push @allData, [@colData];
    }
    
}

open(OUTPUT, ">$outputfile");


for ($ival = 0; $ival <= $#colData; $ival++)
{
    for ($icol = 0; $icol <=$#allData; $icol++)
    {
	print OUTPUT $allData[$icol][$ival]; 
	print OUTPUT " "; 
    }
    print OUTPUT "\n";
}
my $numcols = $#allData + 1; 
my $numpoints = $#colData + 1;  
print " The number of columns and data points extracted are $numcols and $numpoints. \n";  
print " Done! Please check the file, $outputfile, for the extracted column data. \n";


