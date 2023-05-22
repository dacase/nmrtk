#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;

#
my $myhead = <<end_of_myhead;
############################################################################################################
# This a perl program to get all correlation functions of overall motion by dividing the internal motions  #
# from the total correlation functions.                                                                    #
#                                                                                                          #
# Junchao Xia, junchao-xia\@biomaps.rutgers.edu, BioMaPs institute, Rutgers University.                     #
############################################################################################################;

end_of_myhead
#
print $myhead; 
print "\n";
print " NOTE: Unlike getCFOvr.pl, this program get overall correlation functions for all vectors. \n";
print "\n";
print "\n";

my $usage = <<end_of_usage;
usage:"Please input as [command] [columnID of time] [columnID of total] [columnID of internal] [Filename list]"
end_of_usage

die $usage unless getopts("");

if ($#ARGV != 3) {die $usage;}

my $colIDtime = $ARGV[0];      # column ID for time
my $colIDtot = $ARGV[1];       # column ID for total CF
my $colIDint = $ARGV[2];       # column ID for internal CF
my $fileList = $ARGV[3];       # name list for CF files 

my @fileNameInt; 
my @fileNameTot; 
my @fileNameOvr; 
 
open(INPUT, "$fileList");

while(<INPUT>) #read in data
{
    chomp;
    my @tmp_array=split;
    if (defined $tmp_array[0] )
    {
	push (@fileNameInt,$tmp_array[0]);
	push (@fileNameTot,$tmp_array[1]);
	push (@fileNameOvr,$tmp_array[2]);
    }
}
close(INPUT);


for (my $ifc = 0; $ifc <= $#fileNameInt; $ifc++)
{
 
  print " Get the $ifc th file, $fileNameOvr[$ifc].\n";
  system " getCFOvr.pl $fileNameTot[$ifc] $colIDtime $colIDtot $fileNameInt[$ifc] $colIDint $fileNameOvr[$ifc] >& /dev/null";
  
}

print "Done! \n";


