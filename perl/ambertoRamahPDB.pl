#!/usr/bin/perl -w
# ./ambertoRamahPDB.pl amberpdb_file ramahpdb_file
#   reformat Amber PDB to be readable by RAMAH
require 5.002;
use Getopt::Std;
use warnings;
use strict;

my $usage = <<end_of_usage;
usage:"should input like this [command] [input amberpdb file] [output ramahpdb file]"
end_of_usage
    
    die $usage unless getopts("");

if($#ARGV==1)
{
    my $amberpdb = $ARGV[0];
    my $ramahpdb = $ARGV[1];

    open(PDBIN,$amberpdb);
    my @lines = <PDBIN>;
    close PDBIN;

    open(PDBOUT,">$ramahpdb");

    for my $i (0..$#lines) {
	if ($lines[$i] =~ /^ATOM/) {
	    # read PDB ATOM line
	    my $id = substr($lines[$i],6,5);
	    my $name = substr($lines[$i],12,4);
	    my $element = substr($lines[$i],13,1);
	    my $alt = substr($lines[$i],16,1);
	    my $resname = substr($lines[$i],17,3);
	    my $chain = substr($lines[$i],21,1);
	    my $resid = substr($lines[$i],22,4);
	    my $x = substr($lines[$i],30,8);
	    my $y = substr($lines[$i],38,8);
	    my $z = substr($lines[$i],46,8);
	    my $occup = substr($lines[$i],54,6);
	    my $bfactor = substr($lines[$i],60,6);
	    # modify fields
	    for (($id,$name,$resname,$resid)) {s/^\s+//;s/\s+$//}
	    $resname =~ s/D([ATGC])/$1/;
	    $resname =~ s/([ATGC])[35]/$1/;
	    $name =~ s/([1-9])H(.*)/H$2$1/;
	    # print reformatted PDB ATOM line
	    printf PDBOUT "ATOM  %5d  %-4s%3s%1s%5d    %8.3f%8.3f%8.3f%6.2f%6.2f %2s\n",$id,$name,$resname,' ',$resid,$x,$y,$z,1,0,$element;
	}
	else {print PDBOUT $lines[$i]}
    }
    
}
else 
{
    die $usage
};
