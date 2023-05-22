#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;

my $myhead = <<end_of_myhead;
############################################################################################################
# This is a perl program to calculate chemical shifts for each snapshot from the AMBER pdb                 #
# trajectory by calling Sparta+.                                                                           #
#                                                                                                          #
# 1) get the trajectory in pdb format from AMBER ptraj.                                                    #
# 2) run this perl script.                                                                                 # 
#                                                                                                          #
#                                                                                                          #
# Junchao Xia, junchao-xia\@biomaps.rutgers.edu, BioMaPs institute, Rutgers University.                     #
############################################################################################################;
end_of_myhead
print $myhead;
print "\n";

#
#
my $usage = <<end_of_usage;
usage:"Please input as [command] [ncs] [Sparta+ Parameters] [input amberpdb file]"
end_of_usage

die $usage unless getopts("");

if($#ARGV==2)
{
    my $ncs=$ARGV[0];
    my $spartaParm = $ARGV[1];
    my $inputpdb = $ARGV[2];

    my $amberpdb = "tempamber.pdb";
    my $csoutfile = "csoutput.tab";

    my (@ChemShifts,@cshift);
    my  @rmsd;
    my  @rfactor;
    my  @exp_cs;
   
    my $iconform = 0; 
    open(INPUTPDB,"<$inputpdb");	
    while(<INPUTPDB>)
    {
	chomp;
	my $line = $_;
	if ($line =~ /^MODEL/)  # begin of new model
	{
	    open(OUTPUTPDB, ">$amberpdb");
	}
	print OUTPUTPDB "$line \n";
	if ($line =~ /^ENDMDL/)  # end of new model
	{
	    close(OUTPUTPDB);

            # call pales and calculate RDCs using the $amberpdb file.
	    system "sparta+ $spartaParm -in $amberpdb -out $csoutfile >& .log";
	
	    # input the calcuated RDCs from Pales  
	    my $ntail = $ncs + 0;
	    system "cat $csoutfile | tail -n $ntail > cs_out.temp";
	    open(INPUT, "cs_out.temp");
	    @cshift = ();
	    while(<INPUT>) #read in data
	    {
		chomp;
		my @tmp_array=split;
		if (defined $tmp_array[4] )
		{
		    push (@cshift,$tmp_array[4]) ; 
		}
	    }
	    close(INPUT);
	    push @ChemShifts, [@cshift];
            
	    
	    $iconform = $iconform + 1; 
	}
    }
    close(INPUTPDB);


    # output all ChemShifts from sparta+
    open(OUTPUT, ">CShiftsAll.out");
    for (my $iconform = 0; $iconform <= $#ChemShifts; $iconform++ )
       { 
	   printf OUTPUT " %10d ", $iconform;

	   for (my $ics = 0; $ics < $ncs;  $ics++)
	   {
	       printf OUTPUT "%9.4f ",$ChemShifts[$iconform][$ics];
	   }
	   printf OUTPUT "\n";
       }
    close(OUTPUT);

	
}
else 
{
    die $usage
};

