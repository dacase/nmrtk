#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;

my $myhead = <<end_of_myhead;
############################################################################################################
# This is a perl program to rescale the calculated NMR quantities to best match experimantal values.       #
#                                                                                                          #
#                                                                                                          #
# 1) prepare the experimental data file.                                                                   #
# 2) obtain the calculated NMR data and put in a file.                                                     #
# 3) run this perl script.                                                                                 # 
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
usage:"Please input as [command] [experimental data file] [Column No] [calculated data file] [Column No] [output file]"
end_of_usage

die $usage unless getopts("");

my $ax = 0.1; 
my $bx = 0.4; 
my $cx = 10.0;
my $tol = 1e-6;
my $xmin; 
my $fmin;
my @exp_values = (); 
my @cal_values = (); 
my @ids = ();
my @idsToFit = ();

if($#ARGV==4)
{
    my $expDatafile=$ARGV[0];
    my $ncolExp   =$ARGV[1];
    my $calDatafile= $ARGV[2];
    my $ncolCal   =$ARGV[3];
    my $outputfile = $ARGV[4];


    # input the experimental NMRs  
    open(INPUT,"<$expDatafile");
     while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[1] )
	{
	    push (@idsToFit,$tmp_array[0]) ; 
	    push (@exp_values,$tmp_array[$ncolExp-1]) ;  
	}
    }
    close(INPUT);

    # input the calculated NMRs  
    open(INPUT,"<$calDatafile");
     while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[1] )
	{
	    push (@ids,$tmp_array[0]) ; 
	    push (@cal_values,$tmp_array[$ncolCal-1]) ;  
	}
    }
    close(INPUT);

    $ax = 0.1; 
    $bx = 0.4; 
    $cx = 10.0;
    $tol = 1e-6;
    $xmin = 0.0; 
    $fmin = 0.0;
    &Optimize1D();


   # output all NMR quantities scaled by a factor to best fit the expermental values
    open(OUTPUT, ">$outputfile");
   
    for (my $inmr = 0; $inmr <= $#ids;  $inmr++)
    {
      	printf OUTPUT "%10d %9.4f %9.4f \n", $ids[$inmr], $cal_values[$inmr], $xmin*$cal_values[$inmr];
    }
  
    close(OUTPUT);

}
else 
{
    die $usage
};

sub Optimize1D {

    my $R=0.61803399;
    my $C=1.0-$R;
    
    my $f1;
    my $f2;
    my $x0;
    my $x1;
    my $x2;
    my $x3;
    
    $x0=$ax;
    $x3=$cx;
    if (abs($cx-$bx) > abs($bx-$ax)) 
    {
	$x1=$bx;
	$x2=$bx+$C*($cx-$bx);
    }
    else
    {
	$x2=$bx;
	$x1=$bx-$C*($bx-$ax);
    }
    $f1=&Rfactor($x1);
    $f2=&Rfactor($x2);
    while (abs($x3-$x0) > $tol*(abs($x1)+abs($x2))) 
    {
	if ($f2 < $f1) {
	    # &shift3($x0,$x1,$x2,$R*$x2+$C*$x3);
	    $x0=$x1;
	    $x1=$x2;
	    $x2=$R*$x2+$C*$x3;
	    # &shift2($f1,$f2,Rfactor($x2));
	    $f1 = $f2;
	    $f2 = &Rfactor($x2);

	} else {
	    # &shift3($x3,$x2,$x1,$R*$x1+$C*$x0);
	    $x3=$x2;
	    $x2=$x1;
	    $x1=$R*$x1+$C*$x0;
	    # &shift2($f2,$f1,Rfactor($x1));
	    $f2 = $f1;
	    $f1 = &Rfactor($x1);
	}
	
    }
    if ($f1 < $f2)
    {
	$xmin=$x1;
	$fmin=$f1;
    }
    else
    {
	$xmin=$x2;
	$fmin=$f2;
	
    }
	
}

sub Rfactor {

    my $xin=$_[0];  # the input x value
    my $output = 0.0;
    my $value; 

    for (my $i =0; $i <= $#idsToFit; $i++)
    {
	# $value = ($xin*$cal_values[$i]-$exp_values[$i])/$exp_values[$i]; 
	$value = ($xin*$cal_values[$idsToFit[$i]]-$exp_values[$i]); 
        $output = $output + $value*$value;
    }
    return $output;
}

