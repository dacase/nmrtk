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
|   This is a perl program to obtain the potential energy as the function of chemical        |
|   reaction coordinates from the energy time series from a trajectory file.                 |
|                                                                                            |
|   Note: one dimension or two dimensions only.                                              | 
|                                                                                            |
|      Junchao Xia, BioMaPS Institute, Rutgers University, junchao-xia\@biomaps.rutgers.edu.  |
----------------------------------------------------------------------------------------------
end_of_head


print $head;

my $usage = <<end_of_usage;

usage:
  1D: getPotE-CV.pl nx minx maxx inputx xcolumn inputf fcolumn output
  2D: getPoTE-CV.pl nx minx maxx inputx xcolumn ny miny maxy inputy ycolumn inputf fcolum output 
end_of_usage
    

die $usage unless getopts("");


if($#ARGV==12)
{
    &PotE2d($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8],$ARGV[9],$ARGV[10],$ARGV[11],$ARGV[12]);
}
elsif($#ARGV==7)
{
    &PotE1d($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7]);
}
else
{ print "Input options are wrong. Please check the usage below. \n";
  die $usage;
}

##############################################################
sub PotE1d {


    my $nx=$_[0];         # number of divisions
    my $minx=$_[1];       # minimum value 
    my $maxx=$_[2];       # maximum value 
    my $inputx=$_[3];     # file name of input data for x (cv1)
    my $xcolumn=$_[4];    # the columm number of input data for x (cv1) 
    my $inputf=$_[5];     # file name of input data for potential energy  
    my $fcolumn=$_[6];    # the columm number of input data for potential energy
    my $output=$_[7];     # the file name of output for potentail as function of CVs

    my @datax;            # array to store the input data for x (cv1)
    open(INPUT, "$inputx");
    LINEX: while(<INPUT>) #read in data for x (cv1)
    {
	chomp;
	next LINEX if ( /^#/ || /^@/ );   # discard comments 
	my @tmp_array=split;
	if (defined $tmp_array[$xcolumn-1] )
	{
	    push (@datax,$tmp_array[$xcolumn-1]) ;
	}
    }
    close(INPUT);

    my @dataf;            # array to store the input data for potential energy
    open(INPUT, "$inputf");
    LINEF: while(<INPUT>) #read in data for potential energy
    {
	chomp;
	next LINEF if ( /^#/ || /^@/ );   # discard comments 
	my @tmp_array=split;
	if (defined $tmp_array[$fcolumn-1] )
	{
	    push (@dataf,$tmp_array[$fcolumn-1]) ;
	}
    }
    close(INPUT);

    my @histo;
    my @sumf; 
    for my $counter (0 .. $nx-1) #initializing the array
    {
        push (@histo,0);
	push (@sumf,0);
    }

    my $dx=($maxx-$minx)/$nx;
      
    for my $counter (0 .. $#datax)
    {
	my $tmp=int(($datax[$counter]-$minx)/$dx);
	#if ($tmp>=$nx)
	#{$tmp=$nx-1;}
        #if ($tmp<=0)
	#{$tmp=0;}
	if ($tmp<$nx && $tmp>=0)
	{

	   ++$histo[$tmp];
	   $sumf[$tmp] += $dataf[$counter];
           # no weight is required 
           # $histo[$tmp] += exp(-$dataf[$counter]);
	   # $sumf[$tmp] += $dataf[$counter]*exp(-$dataf[$counter]);
 
	}
    }

    open(OUTPUT, ">$output");
    for my $counter (0 .. $#histo)
    {        
	$sumf[$counter] /= $histo[$counter];
	my $tmp=$minx+$dx*($counter+1)-$dx/2;
	printf OUTPUT "$tmp    $sumf[$counter] \n";
    }
    close(OUTPUT);

    print " The 1D potential function has been written to the file, $output.\n";
   
    my @file_type = split(/\./,$output); #splitting file to get the base and extension
    my $plotname = $file_type[0];

    open(OUTPUT, ">template.gnu");
    printf OUTPUT "set pointsize 1.0 \n";
    printf OUTPUT "set style data linespoints \n";
    printf OUTPUT "set xlab \"{/= 24 x }\" \n";
    printf OUTPUT "set ylab \"{/= 24 PotE(x)}\" \n";
    printf OUTPUT "set xrange\[$minx:$maxx\] \n";
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set output \"$plotname.ps\" \n";
    printf OUTPUT "plot \"$output\" u 1:2 t \"\" \n";
    close(OUTPUT);
    system "gnuplot template.gnu";
    print " The 1D potential function has been plotted into the file, $plotname.ps.\n";

}
#########################################################################
sub PotE2d {
    
    my $nx=$_[0];         # number of divisions in x (CV 1)
    my $minx=$_[1];       # minimum value in x 
    my $maxx=$_[2];       # maximum value in x
    my $inputx=$_[3];     # file name of input data in x  
    my $xcolumn=$_[4];    # the columm number of input data in x 
    my $ny=$_[5];         # number of divisions in y (CV 2)
    my $miny=$_[6];       # minimum value in y 
    my $maxy=$_[7];       # maximum value in y
    my $inputy=$_[8];     # file name of input data in y
    my $ycolumn=$_[9];    # the columm number of input data in y
    my $inputf=$_[10];     # file name of input data for potential energy
    my $fcolumn=$_[11];    # the columm number of input data for potenial energy
    my $output=$_[12];    # the file name of output of the potential energy function 


    my @datax;            # array to store the input data in x
    open(INPUT, "$inputx");
    LINEX: while(<INPUT>) #read in data in x  (CV 1)
    {
	chomp;
        next LINEX if ( /^#/ || /^@/ );  # discard comments
	my @tmp_array=split;
	if (defined $tmp_array[$xcolumn-1] )
	{
	    push (@datax,$tmp_array[$xcolumn-1]) ;
	}
    }
    close(INPUT);

    my @datay;            # array to store the input data in y
    open(INPUT, "$inputy");
    LINEY: while(<INPUT>) #read in data in y  (CV 2)
    {
	chomp;
        next LINEY if ( /^#/ || /^@/ );  # discard comments
	my @tmp_array=split;
	if (defined $tmp_array[$ycolumn-1] )
	{
	    push (@datay,$tmp_array[$ycolumn-1]) ;
	}
    }
    close(INPUT);

   my @dataf;            # array to store the input data for potential energy
    open(INPUT, "$inputf");
    LINEF: while(<INPUT>) #read in data for potential energy
    {
	chomp;
	next LINEF if ( /^#/ || /^@/ );   # discard comments 
	my @tmp_array=split;
	if (defined $tmp_array[$fcolumn-1] )
	{
	    push (@dataf,$tmp_array[$fcolumn-1]) ;
	}
    }
    close(INPUT);


##ELEMENTS (I,J) WILL BE STORED AS (I,J)=((I-1)*NY+J)-1
## I,J GOING FROM 1 TO NX AND 1 TO NY

    my @histo;
    my @sumf;
    for my $I (1 .. $nx) #initializing the array
    {	
	for my $J (1 .. $ny)
	{    
	    $histo[(($I-1)*$ny+$J)-1]=0; 
	    $sumf[(($I-1)*$ny+$J)-1]=0; 
	}    
    }
 
    my $dx=($maxx-$minx)/$nx;
    my $dy=($maxy-$miny)/$ny;

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
	if ($tmpx<$nx && $tmpx>=0 && $tmpy<$ny && $tmpy>=0)
	{
	    my $x=$tmpx*$ny+$tmpy;
	    ++$histo[$x];
	    $sumf[$x] += $dataf[$I];
	   # no weight is required 
           # $histo[$x] += exp(-$dataf[$I]);
	   # $sumf[$x] += $dataf[$I]*exp(-$dataf[$I]);
 
	}
    }

    open(OUTPUT, ">$output");    
    for my $x (1 .. $nx)
    { 
	my $xval=$minx+$dx*($x)-$dx/2;
	for my $y (1 .. $ny)
	{
	    my $tmp1=(($x-1)*$ny+$y)-1;
	    $sumf[$tmp1] /= $histo[$tmp1];
	    my $yval=$miny+$dy*($y)-$dy/2;
	    print OUTPUT "$xval    $yval    $sumf[$tmp1]\n";
	}
	print OUTPUT "\n";
    }
    close(OUTPUT);

    print " The 2D potential function has been written to the file, $output.\n";   

    my @file_type = split(/\./,$output); #splitting file to get the base and extension
    my $plotname = $file_type[0];

    open(OUTPUT, ">template.gnu");	
    printf OUTPUT "set dgrid3d $nx, $ny, 2 \n";
    printf OUTPUT "set view map \n";
    printf OUTPUT "set contour base \n";
    printf OUTPUT "set cntrparam levels auto \n";
    printf OUTPUT "set pm3d implicit at b \n";
    printf OUTPUT "set palette model CMY \n";
    printf OUTPUT "set style data lines \n";
    printf OUTPUT "unset key \n";
    printf OUTPUT "unset surface \n";
    printf OUTPUT "set xlab \"{/= 24 {/Symbol F}}\" \n";
    printf OUTPUT "set ylab \"{/= 24 {/Symbol Y}}\" \n";
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set output \"$plotname.ps\" \n";
    printf OUTPUT "splot \"$output\" \n";
    close(OUTPUT);
    system "gnuplot template.gnu";
    print " The 2D potential function has been plotted into the file, $plotname.ps.\n";

}
