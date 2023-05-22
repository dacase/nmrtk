#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;
#
#
my $usage = <<end_of_usage;
usage:"Please input as [command] [number of atoms] [number of conformer] [number of RDCs] [input amberpdb file] [dihedral file]"
end_of_usage

    die $usage unless getopts("");

if($#ARGV==4)
{
    my $natom=$ARGV[0];
    my $nconform=$ARGV[1];
    my $nrdc=$ARGV[2];
    my $inputpdb = $ARGV[3];
    my $dihedralFile = $ARGV[4];

    my $amberpdb = "tempamber.pdb";
    my $ramahpdb = "tempramah.pdb";

    my (@OrdTensors,@orderparm);
    my (@RDCs,@rdc);
    my  @rmsd;
    my  @qfactor;

    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
    {
	# get the ramah pdb file and calculate RDCs via RAMAH
        my $nlines_per_mol = $natom + 4;
        my $nlines_head = $nlines_per_mol * $iconform;
        my $nlines_tail = $nlines_per_mol * ($iconform+1);
        system "head -n $nlines_tail $inputpdb | tail -n $nlines_per_mol  >  $amberpdb ";
	system "ambertoRamahPDB.pl $amberpdb $ramahpdb ";
	system "RAMAH $ramahpdb NONE input.inp output >& ramah.log";

	# input the alignment/order tensor 
	system "cat output.stat | grep -A 6 \"Order Tensor\" | tail -n 5 > OrderTensor.out";
	open(INPUT, "OrderTensor.out");
	@orderparm = ();
	while(<INPUT>) #read in data
	{
	    chomp;
	    my @tmp_array=split;
	    if (defined $tmp_array[2] )
	    {
		push (@orderparm,$tmp_array[2]) ; 
	    }
	}
        close(INPUT);
	push @OrdTensors,[@orderparm];
	
	# input the RDCs  
	my $ntail = $nrdc + 2;
	system "cat output.stat | tail -n $ntail | head -n $nrdc > RDCs.out";
	open(INPUT, "RDCs.out");
	@rdc = ();
	while(<INPUT>) #read in data
	{
	    chomp;
	    my @tmp_array=split;
	    if (defined $tmp_array[6] )
	    {
		push (@rdc,$tmp_array[6]) ; 
	    }
	}
        close(INPUT);
	push @RDCs, [@rdc];


        # input the RDC RMSD 
	system "cat output.stat | grep RDC > rmsd.out ";
	open(INPUT, "rmsd.out");
	while(<INPUT>) #read in data
	{
	    chomp;
	    my @tmp_array=split(/\|/);
	    
	    if (defined $tmp_array[1] )
	    {
		push (@rmsd,$tmp_array[1]) ; 
	    }
	}
        close(INPUT);

	# input the RDC Q factor 
 	system "cat output.stat | grep Q | tail -n 1 > qfactor.out ";
	open(INPUT, "qfactor.out");
	while(<INPUT>) #read in data
	{
	    chomp;
	    my @tmp_array=split;
	    
	    if (defined $tmp_array[2] )
	    {
		push (@qfactor,$tmp_array[2]) ; 
	    }
	}
        close(INPUT);	

    }

    # output all alignment tensors  
    open(OUTPUT, ">OrderTensorsAll.out");
    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
       { 
	   printf OUTPUT " %10d ", $iconform;

	   for (my $iorder = 0; $iorder < 5;  $iorder++)
	   {
	       printf OUTPUT "%10.6f ",$OrdTensors[$iconform][$iorder];
	 
	   }
	   printf OUTPUT "%9.4f %9.3f \n",$rmsd[$iconform],$qfactor[$iconform];

       }
    close(OUTPUT);

    # output all RDCs 
    open(OUTPUT, ">RDCsAll.out");
    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
       { 
	   printf OUTPUT " %10d ", $iconform;

	   for (my $irdc = 0; $irdc < $nrdc;  $irdc++)
	   {
	       printf OUTPUT "%10.5f ",$RDCs[$iconform][$irdc];
	   }
	   printf OUTPUT "\n";
       }
    close(OUTPUT);

    # put and Q factor in a hash table 
    my %hashQ = ();
    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
    {
	$hashQ{$iconform} = $qfactor[$iconform];	
    }

    # Get a list of hash keys sorted by value.
    my @sortedID = sort { $hashQ{$a} <=> $hashQ{$b} } keys %hashQ; 

    
    # output all sorted alignment tensors  
    open(OUTPUT, ">OrderTensorsSorted.out");
    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
       { 
	   printf OUTPUT " %10d %10d ",$iconform, $sortedID[$iconform];

	   for (my $iorder = 0; $iorder < 5;  $iorder++)
	   {
	       printf OUTPUT "%10.6f ",$OrdTensors[$sortedID[$iconform]][$iorder];
	   }
	   printf OUTPUT "%9.4f %9.3f \n",$rmsd[$sortedID[$iconform]],$qfactor[$sortedID[$iconform]];
       }
    close(OUTPUT);

    # output all sorted RDCs 
    open(OUTPUT, ">RDCsSorted.out");
    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
       { 
	   printf OUTPUT " %10d %10d ",$iconform, $sortedID[$iconform];

	   for (my $irdc = 0; $irdc < $nrdc;  $irdc++)
	   {
	       printf OUTPUT "%10.6f ",$RDCs[$sortedID[$iconform]][$irdc];
	   }
	   printf OUTPUT "\n";
       }
    close(OUTPUT);

    open(INPUT, "$dihedralFile");
    my (@dihedrals,@tmp_array);
    my $nangle;

    while(<INPUT>) #read in data
    {
	chomp;
	@tmp_array=split;
	
	if (defined $tmp_array[0] )
	{
	    push @dihedrals,[@tmp_array]; 
	}
	$nangle = $#tmp_array + 1;
    }
    close(INPUT);

    # output all sorted RDCs 
    open(OUTPUT, ">PhiPsiSorted.out");
    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
       { 
	   printf OUTPUT "  %10d %10d ",$iconform, $sortedID[$iconform];

	   for (my $iangle = 0; $iangle < $nangle;  $iangle++)
	   {
	       printf OUTPUT "%10.6f ",$dihedrals[$sortedID[$iconform]][$iangle];
	   }
	   printf OUTPUT "\n";
       }
    close(OUTPUT);
	
}
else 
{
    die $usage
};
