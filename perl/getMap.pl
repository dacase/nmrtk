#!/usr/bin/perl -w

require 5.002;
use utf8;
use Getopt::Std;
use warnings;
use strict;

#
#
my $head = <<end_of_head;
----------------------------------------------------------------------------------------------
|   This is a perl program to produce 1D or 2D map according to input data.                   |  
|                                                                                             | 
|                                                                                             |
|      Junchao Xia, BioMaPS Institute, Rutgers University, junchao-xia\@biomaps.rutgers.edu.   |
----------------------------------------------------------------------------------------------
end_of_head


print $head;

my $usage = <<end_of_usage;

usage:
  1D: getMap.pl nx minx maxx inputx xcolumn inputy ycolumn output 
  2D: getMap.pl nx minx maxx inputx xcolumn ny miny maxy inputy ycolumn inputz zcolumn output 
end_of_usage
    

die $usage unless getopts("");


if($#ARGV==12)
{
    &hist2d($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8],$ARGV[9],$ARGV[10],$ARGV[11],$ARGV[12]);
}
elsif($#ARGV==7)
{
    &hist1d($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7]);
}
else
{ print "Input options are wrong. Please check the usage below. \n";
  die $usage;
}

##############################################################
sub hist1d {


    my $nx=$_[0];         # number of divisions
    my $minx=$_[1];       # minimum value 
    my $maxx=$_[2];       # maximum value 
    my $inputx=$_[3];     # file name of input data  
    my $xcolumn=$_[4];    # the columm number of input data 
    my $inputy=$_[5];     # file name of input data  
    my $ycolumn=$_[6];    # the columm number of input data 
    my $output=$_[7];     # the file name of output of the histogram 
    my @datax;            # array to store the input data for x
    my @datay;            # array to store the input data for y

    #read in data in X
    open(INPUT, "$inputx");
    while(<INPUT>) 
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[$xcolumn-1] )
	{
	    push (@datax,$tmp_array[$xcolumn-1]) ;
	}
    }
    close(INPUT);

    #read in data in Y
    open(INPUT, "$inputy");
    while(<INPUT>)
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[$ycolumn-1] )
	{
	    push (@datay,$tmp_array[$ycolumn-1]) ;
	}
    }
    close(INPUT);

    # find the maximum and min value in datay 
    my $maxv = -100000000; 
    my $minv =  100000000;
    for my $counter (0 .. $#datay) 
    {
        if ($datay[$counter] > $maxv)
	{$maxv = $datay[$counter];}
	if ($datay[$counter] < $minv)
	{$minv = $datay[$counter];}
    }
 
    #initializing the array
    my @map1D;
    my @storeInd;
    for my $counter (0 .. $nx-1) 
    {
        push (@map1D,$maxv);
        push (@storeInd,0);
    
    }

    my $dx=($maxx-$minx)/$nx;
    for my $counter (0 .. $#datay)
    {
	my $tmp=int(($datax[$counter]-$minx)/$dx);
	if ($tmp>=$nx)
	{$tmp=$nx-1;}
        if ($tmp<=0)
	{$tmp=0;}
        if ($datay[$counter] < $map1D[$tmp])
	{
	    $map1D[$tmp] = $datay[$counter];
	    $storeInd[$tmp] = $counter; 
	}
    } 

    open(OUTPUT, ">$output");
    for my $counter (0 .. $#map1D)
    {
	my $tmp=$minx+$dx*($counter+1)-$dx/2;
	printf OUTPUT "%15.8f %15.8f \n",  $tmp, $map1D[$counter];
    }
    close(OUTPUT);
    open(OUTPUT, ">IndexMap.dat");
    for my $counter (0 .. $#storeInd)
    {
	printf OUTPUT "%15d %15d \n",  $counter, $storeInd[$counter];
    }
    close(OUTPUT);

    print " The map data have been written to the file, $output.\n";
    print " The map indices have been written to the file, IndexMap.dat.\n";
   
    my @file_type = split(/\./,$output); #splitting file to get the base and extension
    my $plotname = $file_type[0];

    open(OUTPUT, ">template.gnu");
    printf OUTPUT "set pointsize 1.0 \n";
    printf OUTPUT "set style data linespoints \n";
    printf OUTPUT "set xlab \"{/= 24 x }\" \n";
    printf OUTPUT "set ylab \"{/= 24 P(x)}\" \n";
    printf OUTPUT "set xrange\[$minx:$maxx\] \n";
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set output \"$plotname.ps\" \n";
    printf OUTPUT "plot \"$output\" u 1:2 t \"\" \n";
    close(OUTPUT);
    system "gnuplot template.gnu";
    print " The map in one dimension has been plotted into the file, $plotname.ps.\n";

}
#########################################################################
sub hist2d {
    
    my $nx=$_[0];         # number of divisions in x
    my $minx=$_[1];       # minimum value in x 
    my $maxx=$_[2];       # maximum value in x
    my $inputx=$_[3];     # file name of input data in x  
    my $xcolumn=$_[4];    # the columm number of input data in x
    my $ny=$_[5];         # number of divisions in y
    my $miny=$_[6];       # minimum value in y 
    my $maxy=$_[7];       # maximum value in y
    my $inputy=$_[8];     # file name of input data in y
    my $ycolumn=$_[9];    # the columm number of input data in y
    my $inputz=$_[10];    # file name of input data in z
    my $zcolumn=$_[11];   # the columm number of input data in z
    my $output=$_[12];    # the file name of output of the histogram 
    my @datax;            # array to store the input data in x
    my @datay;            # array to store the input data in y
    my @dataz;            # array to store the input data in z

    #read in data in X
    open(INPUT, "$inputx");
    while(<INPUT>)
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[$xcolumn-1] )
	{
	    push (@datax,$tmp_array[$xcolumn-1]) ;

	}
    }
    close(INPUT);

    for my $counter (0 .. $#datax) 
    {
        if ($datax[$counter] < 0)
	{
	    $datax[$counter] = $datax[$counter] + 360;
	}
	else 
	{
	    $datax[$counter] = $datax[$counter];  
	}
    }


    #read in data in Y
    open(INPUT, "$inputy");
    while(<INPUT>)
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[$ycolumn-1] )
	{
	    push (@datay,$tmp_array[$ycolumn-1]) ;
	}
    }
    close(INPUT);

    for my $counter (0 .. $#datay) 
    {
        if ($datay[$counter] < 0)
	{
	    $datay[$counter] = $datay[$counter] + 360;
	}
	else 
	{
	    $datay[$counter] = $datay[$counter];  
	}
    }

    #read in data in z
    open(INPUT, "$inputz");
    while(<INPUT>)
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[$zcolumn-1] )
	{
	    push (@dataz,$tmp_array[$zcolumn-1]) ;
	}
    }
    close(INPUT);

    # find the maximum an minimum values in dataz 
    my $maxv = -100000000; 
    my $minv =  100000000;
    for my $counter (0 .. $#dataz) 
    {
        if ($dataz[$counter] > $maxv)
	{$maxv = $dataz[$counter];}
	if ($dataz[$counter] < $minv)
	{$minv = $dataz[$counter];}
    }
 

    ##ELEMENTS (I,J) WILL BE STORED AS (I,J)=((I-1)*NY+J)-1
    ## I,J GOING FROM 1 TO NX AND 1 TO NY
    #initializing the array
    my @map2D;
    my @storeInd;
    for my $I (1 .. $nx) 
    {	
	for my $J (1 .. $ny)
	{   
	    push (@map2D,$maxv);
	    push (@storeInd,0);
	    ## $map2d[(($I-1)*$ny+$J)-1]=$maxv; 
            ## $storeInd[(($I-1)*$ny+$J)-1]=0;
	}    
    }
 
    my $dx=($maxx-$minx)/$nx;
    my $dy=($maxy-$miny)/$ny;

    for my $I (0 .. $#dataz)
    {	  
	my $tmpx=int(($datax[$I]-$minx)/$dx); 
	if ($tmpx>=$nx)
	{$tmpx=$nx-1;}
	if ($tmpx<0)
	{$tmpx=0;}
	my $tmpy=int(($datay[$I]-$miny)/$dy); 
	if ($tmpy>=$ny)
	{$tmpy=$ny-1;}
	if ($tmpy<0)
	{$tmpy=0;}

	my $x=$tmpx*$ny+$tmpy;
	if ($dataz[$I] < $map2D[$x])
	{
	    $map2D[$x] = $dataz[$I];
	    $storeInd[$x] = $I; 
	}
    }

    open(OUTPUT, ">$output");    
    for my $x (1 .. $nx)
    {  my $xval=$minx+$dx*($x)-$dx/2;
       for my $y (1 .. $ny)
       {
	   my $yval=$miny+$dy*($y)-$dy/2;
	   my $tmp1=(($x-1)*$ny+$y)-1;
	   printf OUTPUT "%15.8f %15.8f %15.8f \n", $xval, $yval, $map2D[$tmp1];
       }
       printf OUTPUT "\n";
    }
    close(OUTPUT);

    open(OUTPUT, ">IndexMap.dat");
    for my $counter (0 .. $#storeInd)
    {
	printf OUTPUT "%15d %15d \n",  $counter, $storeInd[$counter];
    }
    close(OUTPUT);

    print " The map data have been written to the file, $output.\n";   
    print " The map indices have been written to the file, IndexMap.dat.\n";

    my @file_type = split(/\./,$output); #splitting file to get the base and extension
    my $plotname = $file_type[0];

    open(OUTPUT, ">template.gnu");	
    printf OUTPUT "set dgrid3d $nx, $ny, 1 \n";
    printf OUTPUT "set view map \n";
    printf OUTPUT "set contour base \n";
    printf OUTPUT "set cntrparam levels incremental $minv, 0.02, $maxv \n";
    printf OUTPUT "set pm3d implicit at b \n";
    # printf OUTPUT "set palette model RGB \n";
    # printf OUTPUT "set palette model CMY \n";
    printf OUTPUT "set palette defined ( 0 0 0 0, 0.25 0 0 1, 0.5 0 1 0, 0.75 1 0 0, 1 1 1 1 ) \n";
    printf OUTPUT "set style data lines \n";
    printf OUTPUT "unset key \n";
    printf OUTPUT "unset surface \n";
    printf OUTPUT "set xlab \"{/= 24 {/Symbol F}}\" \n";
    printf OUTPUT "set ylab \"{/= 24 {/Symbol Y}}\" \n";
    printf OUTPUT "set xtics $minx, 60, $maxx \n";
    printf OUTPUT "set ytics $miny, 60, $maxy \n";
    printf OUTPUT "set xrange [$minx : $maxx] \n";
    printf OUTPUT "set yrange [$miny : $maxy] \n";
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set output \"$plotname.ps\" \n";
    printf OUTPUT "splot \"$output\" \n";
    close(OUTPUT);
    system "gnuplot template.gnu";
    print " The map in two dimensions has been plotted into the file, $plotname.ps.\n";

}
