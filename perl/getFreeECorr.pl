#!/usr/bin/perl -w
require 5.002;
use Getopt::Std;
#
#
$usage = <<end_of_usage;
usage:"should input like this [command] [nx|nx ny] [ABMD file] [corrected file] [output file] "
end_of_usage
    
    
    die $usage unless getopts("");


if($#ARGV==4)
{
    &freeE2d($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4]);
}
elsif($#ARGV==3)
{
    &freeE1d($ARGV[0],$ARGV[1],$ARGV[2],,$ARGV[3]);
}
else
{ die $usage;}
##############################################################
sub freeE1d {

    $nx=$_[0]; #number of divisions  
    my $ABMDfile = $_[1];
    my $Corrfile = $_[2];
    my $outputfile = $_[3];

    open(INPUT, $ABMDfile);
    while(<INPUT>) #read in data
    {
	chomp;
	@tmp_array=split;

	push (@array1,$tmp_array[0]) ;
	push (@array2,$tmp_array[1]) ;
    }
    close INPUT;

    open(INPUT, $Corrfile);
    while(<INPUT>) #read in data
    {
	chomp;
	@tmp_array=split;

	push (@array3,$tmp_array[1]) ;
    }
    close INPUT;
    
    open(OUTPUT, ">$outputfile");
    for $counter (0 .. $nx-1)
    {
	$newvalue = $array2[$counter] + $array3[$counter];
	print OUTPUT "$array1[$counter]  $newvalue \n"
    }
    close OUTPUT;

}
#########################################################################
sub freeE2d {
    

    $nx=$_[0]; #number of divisions
    $ny=$_[1]; #number of divisions
    my $ABMDfile = $_[2];
    my $Corrfile = $_[3];
    my $outputfile = $_[4];

    open(INPUT, $ABMDfile);
    while(<INPUT>) #read in data
    {
	chomp;
	@tmp_array=split;

	if (defined $tmp_array[2] )
	{
	    push (@array1,$tmp_array[0]) ;
	    push (@array2,$tmp_array[1]) ;
	    push (@array3,$tmp_array[2]) ;
	}
    }	
    close INPUT;

    open(INPUT, $Corrfile);
    while(<INPUT>) #read in data
    {
	chomp;
	@tmp_array=split;

	if (defined $tmp_array[2] )
	{
	    push (@array4,$tmp_array[2]) ;
	}
    }	
    close INPUT;

    $counter=0;

    open(OUTPUT, ">$outputfile");
	
    for $I (1 .. $nx) #initializing the array
    {	
	for $J (1 .. $ny)
	{   
	    $ii = $I;
	    $jj = $J;

	    $newvalue = $array3[$counter] + $array4[$counter];
	    printf OUTPUT "$array1[$counter]  $array2[$counter]  $newvalue \n";
	     ++$counter;
	    
	}    
	print OUTPUT "  \n";
 	$ii = $I;
       	$jj = $J;
 	
    }
    close OUTPUT;
}
