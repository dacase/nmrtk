#!/usr/bin/perl -w
require 5.002;
use Getopt::Std;
#
#
$usage = <<end_of_usage;
usage:"should input number of histogram divisions nx or nx ny"
end_of_usage
    
    
    die $usage unless getopts("");


if($#ARGV==5)
{
    &hist2d($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5]);
}
elsif($#ARGV==2)
{
    &hist1d($ARGV[0],$ARGV[1],$ARGV[2]);
}
else
{ die $usage;}
##############################################################
sub hist1d {

    $min=$_[0];
    $max=$_[1];

    $nx=$_[2]; #number of divisions
    while(<STDIN>) #read in data
    {
	chomp;
	push (@array,$_) ;
    }

    $counter=0;
    foreach (@array) 
    {
         ++$counter;
    }

    for $counter (0 .. $nx-1) #initializing the array
    {
        $histo[$counter]=0;
    }
    $dx=($max-$min)/$nx;
    
    
    for $counter (0 .. $#array)
    {
	$tmp=int(($array[$counter]-$min)/$dx);
	if ($tmp==$nx)
	{$tmp=$tmp-1 ;}
	++$histo[$tmp];
    }
    $integral=($#array+1)*$dx; 

    for $counter (0 .. $#histo)
    {
	$histo[$counter]=$histo[$counter]/$integral;
	$tmp=$min+$dx*($counter+1)-$dx/2;
	print"$tmp    $histo[$counter] \n";
    }
}
#########################################################################
sub hist2d {
    
    $min1=$_[0];
    $max1=$_[1];
    $nx=$_[2]; #number of divisions

    $min2=$_[3];
    $max2=$_[4];
    $ny=$_[5]; #number of divisions

    while(<STDIN>) #read in data
    {
	chomp;
	@tmp_array=split;
	push (@array1,$tmp_array[0]) ;
	push (@array2,$tmp_array[1]) ;
    }

    $counter=0;

    foreach (@array1) 
    {
         ++$counter;
    }

    $counter=0;


##ELEMENTS (I,J) WILL BE STORED AS (I,J)=((I-1)*NY+J)-1
## I,J GOING FROM 1 TO NX AND 1 TO NY

    for $I (1 .. $nx) #initializing the array
    {	
	for $J (1 .. $ny)
	{    
	    $histo[(($I-1)*$ny+$J)-1]=0;    
	}    
    }
 
    $dx1=($max1-$min1)/$nx;
    $dx2=($max2-$min2)/$ny;

    for $I (0 .. $#array1)
    {	  
	$tmp1=int(($array1[$I]-$min1)/$dx1); 
	if ($tmp1==$nx)
	{$tmp1=$tmp1-1 ;}
	$tmp2=int(($array2[$I]-$min2)/$dx2); 
	if ($tmp2==$ny)
	{$tmp2=$tmp2-1 ;}
	$x=$tmp1*$ny+$tmp2;
	++$histo[$x];
    }

    $integral=($#array1+1)*$dx1*$dx2; #it's  not +1 cause there's a last ++
    
    for $x (1 .. $nx)
    {  $xval=$min1+$dx1*($x)-$dx1/2;
       for $y (1 .. $ny)
       {
	   $tmp1=(($x-1)*$ny+$y)-1;
	   $histo[$tmp1]=$histo[$tmp1]/($integral);
	   $yval=$min2+$dx2*($y)-$dx2/2;
	   print"$xval    $yval    $histo[$tmp1]\n";
       }
       print"\n";
   }
}
