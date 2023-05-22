#!/usr/bin/perl -w
require 5.002;
use Getopt::Std;
use warnings;
use strict; 

#
#
my $head = <<end_of_head;
------------------------------------------------------------------------------------------------
|   This is a perl program to shift the free energy suface in one dimension or two dimensions. |
|   The global minimum will be shift to FMIN and all value greater than FMAX after shift will  |
|   be set to FMAX.                                                                            | 
|                                                                                              |
|      Junchao Xia, BioMaPS Institute, Rutgers University, junchao-xia\@biomaps.rutgers.edu.    |
------------------------------------------------------------------------------------------------
end_of_head

print $head;

my $usage = <<end_of_usage;

usage:
  1D: getFreeEMin.pl nx fmin fmax inputfile outputfile
  2D: getFreeEMin.pl nx ny fmin fmax inputfile outputfile 
end_of_usage
    
die $usage unless getopts("");

if($#ARGV==5)
{
    &freeE2d($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5]);
}
elsif($#ARGV==4)
{
    &freeE1d($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4]);
}
else
{ die $usage;}
##############################################################
sub freeE1d {

    my $nx=$_[0]; #number of divisions  
    my $fmin = $_[1];
    my $fmax = $_[2];
    my $inputfile = $_[3];
    my $outputfile = $_[4];
    my @array1;
    my @array2;

    open(INPUT, $inputfile);
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;

	push (@array1,$tmp_array[0]) ;
	push (@array2,$tmp_array[1]) ;
    }
    close INPUT;

    my $min=$array2[0];
    my $imin= 0;

    my $counter=0;

    foreach (@array2) #searching for minimum and maximum elements
    {

	if($array2[$counter]<$min)
	{
	    $min=$array2[$counter];
	    $imin = $counter;
	}
	 ++$counter;
    }
    print " The minimum is $array1[$imin], $array2[$imin] \n";
    
    open(OUTPUT, ">$outputfile");
    for $counter (0 .. $nx-1)
    {
	my $newvalue = $array2[$counter]- $array2[$imin];
        if ($newvalue > ($fmax-$fmin))
        {$newvalue = $fmax;}
        else
        {$newvalue = $newvalue+$fmin;}
	print OUTPUT "$array1[$counter] $newvalue \n";
    }
    close OUTPUT;
    print " The free enegy in one dimension has bee shifted and written to $outputfile.\n";

    my @file_type = split(/\./,$outputfile); #splitting file to get the base and extension
    my $plotname = $file_type[0];

    open(OUTPUT, ">template.gnu");
    printf OUTPUT "set pointsize 1.0 \n";
    printf OUTPUT "set data style linespoints \n";
    printf OUTPUT "set xlab \"{/= 24 x }\" \n";
    printf OUTPUT "set ylab \"{/= 24 F(x)}\" \n";
    printf OUTPUT "set yrange\[$fmin:$fmax\] \n";
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set output \"$plotname.ps\" \n";
    printf OUTPUT "plot \"$outputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);
    system "gnuplot template.gnu";
    print " The free energy in one dimension has been plotted into the file, $plotname.ps.\n";

}
#########################################################################
sub freeE2d {
    
    my $nx=$_[0]; #number of divisions
    my $ny=$_[1]; #number of divisions
    my $fmin = $_[2];
    my $fmax = $_[3];
    my $inputfile = $_[4];
    my $outputfile = $_[5];
    my @array1;
    my @array2;
    my @array3;

    open(INPUT, $inputfile);
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;

	if (defined $tmp_array[2] )
	{
	    push (@array1,$tmp_array[0]) ;
	    push (@array2,$tmp_array[1]) ;
	    push (@array3,$tmp_array[2]) ;
	}
    }	
    close INPUT;


    my $min=$array3[0];
    my $imin = 0;
 
    my $counter=0;

    foreach (@array3) #searching for minimum and maximum elements
    {

	if($array3[$counter]<$min)
	{
	    $min=$array3[$counter];
	    $imin = $counter;
	}
	 ++$counter;
    }
    print " The minimum is $array1[$imin], $array2[$imin], $array3[$imin] \n";

    $counter=0;

    open(OUTPUT, ">$outputfile");
	
    for my $I (1 .. $nx) #initializing the array
    {	
	for my $J (1 .. $ny)
	{   
	    my $newvalue = $array3[$counter] - $array3[$imin];
            if ($newvalue > ($fmax-$fmin)) 
            {$newvalue = $fmax;}
	    else 	
            {$newvalue = $newvalue+$fmin;}
            printf OUTPUT "$array1[$counter] $array2[$counter] $newvalue \n";

	     ++$counter; 
	}    
	print OUTPUT "  \n";
 	
    }
    close OUTPUT;

    print " The free enegy in two dimensions has bee shifted and written to $outputfile.\n";

    my $maxCntr = $fmin + 5.0;
    my $dxCntr = 0.2; 
    my $xmin = -180;
    my $xmax = 180;
    my $ymin = -180;
    my $ymax = 180;
    my $zmin = $fmin;
    my $zmax = $fmax;

    my @file_type = split(/\./,$outputfile); #splitting file to get the base and extension
    my $plotname = $file_type[0];

    open(OUTPUT, ">template.gnu");
    printf OUTPUT "set dgrid3d $nx, $ny, 2 \n";
    printf OUTPUT "set view map \n";
    printf OUTPUT "set contour base \n";
    printf OUTPUT "set cntrparam levels incremental $fmin, $dxCntr, $maxCntr \n";
    printf OUTPUT "set pm3d implicit at b \n";
    printf OUTPUT "set style data lines \n";
    printf OUTPUT "unset key \n";
    printf OUTPUT "unset surface \n";
    printf OUTPUT "set xrange \[ $xmin : $xmax \] \n";
    printf OUTPUT "set yrange \[ $ymin : $ymax \] \n";
    printf OUTPUT "set zrange \[ $zmin : $zmax \] \n";
    printf OUTPUT "set xlab \"{/= 24 {/Symbol F}}\" \n";
    printf OUTPUT "set ylab \"{/= 24 {/Symbol Y}}\" \n";
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set output \"$plotname.ps\" \n";
    printf OUTPUT "splot \"$outputfile\" \n";
    close(OUTPUT);
    system "gnuplot template.gnu";
    print " The free energy in two dimensions has been plotted into the file, $plotname.ps.\n";

}
