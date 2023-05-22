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
|   This is a perl program to map PMF to new chemical reaction coordinates from              |
|   the time series from a trajectory file by simple averaging the old PMF.                  |
|                                                                                            |
|   Note: one dimension or two dimensions only.                                              | 
|                                                                                            |
|      Junchao Xia, BioMaPS Institute, Rutgers University, junchao-xia\@biomaps.rutgers.edu.  |
----------------------------------------------------------------------------------------------
end_of_head


print $head;

my $usage = <<end_of_usage;

usage:
  1D: mapPMF2NewCV.pl nx minx maxx inputFE inputHisF xcolumn fcolumn nf minf maxf output
  2D: MapPMF2NewCV.pl nx minx maxx ny miny maxy inputFE inputHisF xcolumn ycolumn fcolumn nf minf maxf output 
end_of_usage
    

die $usage unless getopts("");


if($#ARGV==14)
{
    &FE2DTo1D($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8],$ARGV[9],$ARGV[10],$ARGV[11],$ARGV[12],$ARGV[13],$ARGV[14]);
}
elsif($#ARGV==10)
{
    &FE1DTo1D($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8],$ARGV[9],$ARGV[10]);
}
else
{ print "Input options are wrong. Please check the usage below. \n";
  die $usage;
}

##############################################################
sub FE1DTo1D {


    my $nx=$_[0];         # number of divisions
    my $minx=$_[1];       # minimum value 
    my $maxx=$_[2];       # maximum value 
     my $inputFE=$_[3];    # file name for inputing free energy data 

    my $inputHIS=$_[4];   # file name for inputing histogram data 
    my $xcolumn=$_[5];    # the columm number of x in histogram data   
    my $fcolumn=$_[6];   # the columm number of f in histogram data 
    my $nf=$_[7];         # number of divisions in f (new CV)
    my $minf=$_[8];       # minimum value in f 
    my $maxf=$_[9];       # maximum value in f
    my $output=$_[10];    # the file name for outputing FE in new CV  

    my @arrayX;
    my @arrayF;

    open(INPUT, $inputFE);
    while(<INPUT>) #read in data
    {
        chomp;
        my @tmp_array=split;

        if (defined $tmp_array[1] )
        {
            push (@arrayX,$tmp_array[0]) ;
            push (@arrayF,$tmp_array[1]) ;
        }
    }
    close INPUT;

    my @datax;            # array to store the input data in x
    my @dataf;            # array to store the input data in f

    open(INPUT, "$inputHIS");
    LINEX: while(<INPUT>) #read in HIS data
    {
	chomp;
        next LINEX if ( /^#/ || /^@/ );  # discard comments
	my @tmp_array=split;
	if (defined $tmp_array[$xcolumn-1] )
	{
	    push (@datax,$tmp_array[$xcolumn-1]);
            push (@dataf,$tmp_array[$fcolumn-1]);
	}
    }
    close(INPUT);


    my @histo;
    my @sumf; 
    for my $counter (0 .. $nf-1) #initializing the array
    {
        push (@histo,0);
	push (@sumf,0);
    }

    my $dx=($maxx-$minx)/$nx;
    my $df=($maxf-$minf)/$nf;
   
    for my $counter (0 .. $#datax)
    {
	my $tmpx=int(($datax[$counter]-$minx)/$dx);
	#if ($tmp>=$nx)
	#{$tmp=$nx-1;}
        #if ($tmp<=0)
	#{$tmp=0;}
	my $tmpf=int(($dataf[$counter]-$minf)/$df);
 
	if ($tmpx<$nx && $tmpx>=0 && $tmpf<$nf && $tmpf>=0)
	{

	   ++$histo[$tmpf];
	   $sumf[$tmpf] += $arrayF[$tmpx];
           # no weight is required 
           # $histo[$tmp] += exp(-$arrayF[$tmpx]);
	   # $sumf[$tmp] += $dataf[$counter]*exp(-$arrayF[$tmpx]);
 
	}
    }

    open(OUTPUT, ">$output");
    for my $counter (0 .. $#histo)
    {        
	$sumf[$counter] /= $histo[$counter];
	my $tmpf=$minf+$df*$counter;
	printf OUTPUT " %20.10f %20.10f %20.10f \n", $tmpf, $histo[$counter], $sumf[$counter];
    }
    close(OUTPUT);

    print " The new 1D PMF has been written to the file, $output.\n";
   
    my @file_type = split(/\./,$output); #splitting file to get the base and extension
    my $plotname = $file_type[0];

    open(OUTPUT, ">template.gnu");
    printf OUTPUT "set pointsize 1.0 \n";
    printf OUTPUT "set style data linespoints \n";
    printf OUTPUT "set xlab \"{/= 24 x }\" \n";
    printf OUTPUT "set ylab \"{/= 24 PMF(x)}\" \n";
    printf OUTPUT "set xrange\[$minf:$maxf\] \n";
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set output \"$plotname.ps\" \n";
    printf OUTPUT "plot \"$output\" u 1:3 t \"\" \n";
    close(OUTPUT);
    system "gnuplot template.gnu";
    print " The new 1D PMF has been plotted into the file, $plotname.ps.\n";

}
#########################################################################
sub FE2DTo1D {
    
    my $nx=$_[0];         # number of divisions in x (CV 1)
    my $minx=$_[1];       # minimum value in x 
    my $maxx=$_[2];       # maximum value in x
    my $ny=$_[3];         # number of divisions in y (CV 2)
    my $miny=$_[4];       # minimum value in y 
    my $maxy=$_[5];       # maximum value in y
    my $inputFE=$_[6];    # file name for inputing free energy data 

    my $inputHIS=$_[7];   # file name for inputing histogram data 
    my $xcolumn=$_[8];    # the columm number of x in histogram data  
    my $ycolumn=$_[9];    # the columm number of y in histogram data 
    my $fcolumn=$_[10];   # the columm number of f in histogram data 

    my $nf=$_[11];         # number of divisions in f (new CV)
    my $minf=$_[12];       # minimum value in f 
    my $maxf=$_[13];       # maximum value in f
    my $output=$_[14];    # the file name for outputing FE in new CV  

    my @arrayX;
    my @arrayY;
    my @arrayF;

    open(INPUT, $inputFE);
    while(<INPUT>) #read in data
    {
        chomp;
        my @tmp_array=split;

        if (defined $tmp_array[2] )
        {
            push (@arrayX,$tmp_array[0]) ;
            push (@arrayY,$tmp_array[1]) ;
            push (@arrayF,$tmp_array[2]) ;
        }
    }
    close INPUT;

    my @datax;            # array to store the input data in x
    my @datay;            # array to store the input data in y
    my @dataf;            # array to store the input data in f

    open(INPUT, "$inputHIS");
    LINEX: while(<INPUT>) #read in HIS data
    {
	chomp;
        next LINEX if ( /^#/ || /^@/ );  # discard comments
	my @tmp_array=split;
	if (defined $tmp_array[$xcolumn-1] )
	{
	    push (@datax,$tmp_array[$xcolumn-1]);
            push (@datay,$tmp_array[$ycolumn-1]);
            push (@dataf,$tmp_array[$fcolumn-1]);
	}
    }
    close(INPUT);


##ELEMENTS (I,J) WILL BE STORED AS (I,J)=((I-1)*NY+J)-1
## I,J GOING FROM 1 TO NX AND 1 TO NY

    my @histo;
    my @sumf;
    for my $I ( 0 .. $nf-1) #initializing the array
    {	
         $histo[$I]=0; 
	 $sumf[$I]=0;     
    }
 
    my $dx=($maxx-$minx)/$nx;
    my $dy=($maxy-$miny)/$ny;
    my $df=($maxf-$minf)/$nf;

    for my $I (0 .. $#datax)
    {	  
	my $tmpx=int(($datax[$I]-$minx)/$dx); 
	#if ($tmpx>=$nx)
	#{$tmpx=$nx-1;}
	#if ($tmpx<0)
	#{$tmpx=0;}
	my $tmpy=int(($datay[$I]-$miny)/$dy); 
	#if ($tmpy>=$ny)
	#{$tmpy=$ny-1;}
	#if ($tmpy<0)
	#{$tmpy=0;}
	my $tmpf=int(($dataf[$I]-$minf)/$df);
 
	if ($tmpx<$nx && $tmpx>=0 && $tmpy<$ny && $tmpy>=0 && $tmpf<$nf && $tmpf>=0 )
	{
	    my $x=$tmpx*$ny+$tmpy;
	    ++$histo[$tmpf];
	    $sumf[$tmpf] += $arrayF[$x];
	   # no weight is required 
           # $histo[$x] += exp(-$dataf[$I]);
	   # $sumf[$x] += $dataf[$I]*exp(-$dataf[$I]);
 
	}
    }

    open(OUTPUT, ">$output");    
    for my $i (0 .. $nf-1)
    { 
	my $fval=$minf+$df*$i;
        $sumf[$i] /= $histo[$i]; 
	printf OUTPUT " %20.10f %20.10f %20.10f \n", $fval, $histo[$i], $sumf[$i];
	
    }
    close(OUTPUT);

    print " The new PMF has been written to the file, $output.\n";   
   
    my @file_type = split(/\./,$output); #splitting file to get the base and extension
    my $plotname = $file_type[0];

    open(OUTPUT, ">template.gnu");
    printf OUTPUT "set pointsize 1.0 \n";
    printf OUTPUT "set style data linespoints \n";
    printf OUTPUT "set xlab \"{/= 24 x }\" \n";
    printf OUTPUT "set ylab \"{/= 24 PMF(x)}\" \n";
    printf OUTPUT "set xrange\[$minf:$maxf\] \n";
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set output \"$plotname.ps\" \n";
    printf OUTPUT "plot \"$output\" u 1:3 t \"\" \n";
    close(OUTPUT);
    system "gnuplot template.gnu";
    print " The new 1D PMF has been plotted into the file, $plotname.ps.\n";


}
