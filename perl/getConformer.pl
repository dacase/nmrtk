#!/usr/bin/perl -w
require 5.002;
use Getopt::Std;
use warnings;
use strict;
#
#
my $usage = <<end_of_usage;
usage:"Please input the total number of atoms, the id number of conformer, the file type, the input file and the output file."
end_of_usage

    die $usage unless getopts("");

if($#ARGV==4)
{
    my $natom=$ARGV[0];
    my $nconform=$ARGV[1];
    my $filetype = $ARGV[2]; 
    my $inputfile = $ARGV[3];
    my $outputfile = $ARGV[4];
    if ($filetype eq "xyz")
    {	
	my $nlines_per_mol = $natom + 2;
	my $nlines_head = $nlines_per_mol * ($nconform - 1);
	my $nlines_tail = $nlines_per_mol * $nconform;
	system "head -n $nlines_tail $inputfile | tail -n $nlines_per_mol  > $outputfile";
    }
    elsif ($filetype eq "pdb_vmd")
    {
	system "head -n 1 $inputfile  >> $outputfile";
	my $nlines_per_mol = $natom + 1;
	my $nlines_head = $nlines_per_mol * ($nconform - 1) + 1;
	my $nlines_tail = $nlines_per_mol * $nconform + 1;
	system "head -n $nlines_tail $inputfile | tail -n $nlines_per_mol  > $outputfile";
    }
   elsif ($filetype eq "pdb_amber")
    {
        my $nlines_per_mol = $natom + 4;
        my $nlines_head = $nlines_per_mol * ($nconform - 1);
        my $nlines_tail = $nlines_per_mol * $nconform;
        system "head -n $nlines_tail $inputfile | tail -n $nlines_per_mol  > $outputfile";
    }

    else
    {
	print " The input file, $inputfile, is not in a right format.\n";
    }
}
else 
{
    die $usage
};
