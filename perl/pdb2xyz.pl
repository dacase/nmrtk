#!/usr/bin/perl -w
# ./pdb2xyz.pl amberpdb_file xyz
#   transform Amber PDB to xyz   
require 5.002;
use Getopt::Std;
use warnings;
use strict;

my $usage = <<end_of_usage;
usage:"should input like this [command] [input pdb file] [output xyz file] [number of atoms]"
end_of_usage
    
    die $usage unless getopts("");

if($#ARGV==2)
{
    my $pdbfile= $ARGV[0];
    my $xyzfile = $ARGV[1];
    my $natoms  = $ARGV[2];

    open(PDBIN,$pdbfile);
    my @lines = <PDBIN>;
    close PDBIN;

    open(XYZOUT,">$xyzfile");

    for my $i (0..$#lines) {
	if ($lines[$i] =~ /^ATOM/) {
	    # read PDB ATOM line
	    # my $id = substr($lines[$i],6,5);
	    # my $name = substr($lines[$i],12,4);
	    # my $element = substr($lines[$i],13,1);
	    # my $alt = substr($lines[$i],16,1);
	    # my $resname = substr($lines[$i],17,3);
	    # my $chain = substr($lines[$i],21,1);
	    # my $resid = substr($lines[$i],22,4);
	    # my $x = substr($lines[$i],30,8);
	    # my $y = substr($lines[$i],38,8);
	    # my $z = substr($lines[$i],46,8);
	    # my $occup = substr($lines[$i],54,6);
	    # my $bfactor = substr($lines[$i],60,6);
	    # modify fields
	    # for (($id,$name,$resname,$resid)) {s/^\s+//;s/\s+$//}
	    # $resname =~ s/D([ATGC])/$1/;
	    # $resname =~ s/([ATGC])[35]/$1/;
	    # $name =~ s/([1-9])H(.*)/H$2$1/;

	    my @tmp_array = split(' ', $lines[$i]); 
	    my $element = substr($tmp_array[2],0,1);
	    my $x = $tmp_array[5];
	    my $y = $tmp_array[6];
	    my $z = $tmp_array[7];
	    # print reformatted PDB ATOM line
	    printf XYZOUT "%2s %8.3f %8.3f %8.3f\n",$element,$x,$y,$z;
	}
	elsif ($lines[$i] =~ /^MODEL/){
	    printf XYZOUT "%8d\n", $natoms;
	    printf XYZOUT " xyz file converted by pdb2xyz.pl\n"; 
	}
	else
	{}
    }
    
}
else 
{
    die $usage
};
