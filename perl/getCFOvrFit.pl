#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;


my $myhead = <<end_of_myhead;
############################################################################################################
# This a perl program to fit the correlation function of overall motion upto five exponential components   #
# by calling gnuplot.                                                                                      #
#                                                                                                          #
#  1) Co(t) = exp(-t/t1)                                                                                   # 
#  2) Co(t) = a1*exp(-t/t1) + (1-a1)*exp(-t/t2)                                                            #
#  3) Co(t) = a1*exp(-t/t1) + a2*exp(-t/t2) + (1-a1-a2)*exp(-t/t3)                                         #
#  4) Co(t) = a1*exp(-t/t1) + a2*exp(-t/t2) + a3*exp(-t/t3) + (1-a1-a2-a3)*exp(-t/t4)                      #  
#  5) Co(t) = a1*exp(-t/t1) + a2*exp(-t/t2) + a3*exp(-t/t3) + a4*exp(-t/t4) + (1-a0-a1-a2-a3-a4)*exp(-t/t5)# 
#                                                                                                          #
# Junchao Xia, junchao-xia\@biomaps.rutgers.edu, BioMaPs institute, Rutgers University.                     #
############################################################################################################
end_of_myhead
print $myhead;
print "\n";
print " NOTE: This program only fits one data file for the specified correlation function.\n";
print " NOTE: Please input the initial guess of parameters according to the order of appearance above !!! \n";
print "\n";

#
my $usage = <<end_of_usage;
usage:"should input like this [command] [nexp] [parameters from initial guess] [tstart] [tend] [data file]"
end_of_usage

die $usage unless getopts("");

if ($#ARGV < 4 || $#ARGV > 12) {die $usage;}

my $nexp = $ARGV[0];

my @file_type = split(/\./,$ARGV[$#ARGV]); #splitting file to get extension
my $base = $file_type[0];
my $extension = $file_type[$#file_type]; 

my $datafile  = $ARGV[$#ARGV];

 
if ($nexp == 1)
{

    print " Fit the overall rotation correlation function, $datafile, with one exponential. \n";
    &Fit1Exp($ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4]);

    open(INPUT, "fit1ExpOvr.dat");
    my @prmfitted;
    my @errorfitted;
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[2] )
	{
	    push (@prmfitted,$tmp_array[2]) ;
	    push (@errorfitted,$tmp_array[4]) ;
	}
    }
    close(INPUT);

    printf "\n";
    printf "\n";
    printf "\n";
    printf " %15s %15s %15s %15s \n", "Coeff", "error", "Tau", "error";
    printf " %15.6f %15.6f %15.6f %15.6f  \n" ,1.0, 0.0, $prmfitted[0],$errorfitted[0];

    open(OUTPUT, ">$base\_fit.dat");
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" ,1.0, 0.0, $prmfitted[0],$errorfitted[0];
    close(OUTPUT);

    open(OUTPUT, ">$base\_fit_parm.dat");
    printf OUTPUT "%20s %15.6f \n" ,$datafile, $prmfitted[0];
    close(OUTPUT);
    open(OUTPUT, ">$base\_fit_erro.dat");
    printf OUTPUT "%20s %15.6f \n" ,$datafile, $errorfitted[0];

    system "mv fit1ExpOvr.dat $base\_prm.dat";
    system "mv fit1ExpOvr.ps $base.ps" ;
    system "mv fit1ExpOvr.log $base\_gnu.log";

} 
elsif ($nexp == 2)  
{
    print " Fit the overall rotation correlation function, $datafile, with two exponentials. \n"; 
    &Fit2Exp($ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6]);

    open(INPUT, "fit2ExpOvr.dat");
    my @prmfitted;
    my @errorfitted;
    while(<INPUT>) #read in data
    {
        chomp;
        my @tmp_array=split;
        if (defined $tmp_array[2] )
        {
            push (@prmfitted,$tmp_array[2]) ;
            push (@errorfitted,$tmp_array[4]) ;
        }
    }
    close(INPUT);

    printf "\n";
    printf "\n";
    printf "\n";
    printf " %15s %15s %15s %15s \n", "Coeff", "error", "Tau", "error";
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0],$prmfitted[1],$errorfitted[1];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0],$errorfitted[0],$prmfitted[2],$errorfitted[2];

    open(OUTPUT, ">$base\_fit.dat");
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0],$prmfitted[1],$errorfitted[1];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0],$errorfitted[0],$prmfitted[2],$errorfitted[2];
    close(OUTPUT);

    open(OUTPUT, ">$base\_fit_parm.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f \n",$datafile, $prmfitted[0], $prmfitted[1], 1.0-$prmfitted[0],$prmfitted[2];
    close(OUTPUT);
    open(OUTPUT, ">$base\_fit_erro.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f \n" ,$datafile, $errorfitted[0],$errorfitted[1],$errorfitted[0], $errorfitted[2];
    close(OUTPUT);

    system "mv fit2ExpOvr.dat $base\_prm.dat";
    system "mv fit2ExpOvr.ps $base.ps" ;
    system "mv fit2ExpOvr.log $base\_gnu.log";

}
elsif ($nexp == 3)
{
    print " Fit the overall rotation correlation function, $datafile, with three exponentials. \n"; 
    &Fit3Exp($ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8]);

    open(INPUT, "fit3ExpOvr.dat");
    my @prmfitted;
    my @errorfitted;
    while(<INPUT>) #read in data
    {
        chomp;
        my @tmp_array=split;
        if (defined $tmp_array[2] )
        {
            push (@prmfitted,$tmp_array[2]) ;
            push (@errorfitted,$tmp_array[4]) ;
        }
    }
    close(INPUT);

    printf "\n";
    printf "\n";
    printf "\n";
    printf " %15s %15s %15s %15s \n", "Coeff", "error", "Tau", "error";
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0],$prmfitted[1],$errorfitted[1];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[2],$errorfitted[2],$prmfitted[3],$errorfitted[3];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0]-$prmfitted[2],$errorfitted[0],$prmfitted[4],$errorfitted[4];

    open(OUTPUT, ">$base\_fit.dat");
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0],$prmfitted[1],$errorfitted[1];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[2],$errorfitted[2],$prmfitted[3],$errorfitted[3];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0]-$prmfitted[2],$errorfitted[0],$prmfitted[4],$errorfitted[4];
    close(OUTPUT);

    open(OUTPUT, ">$base\_fit_parm.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", $datafile, $prmfitted[0],$prmfitted[1],$prmfitted[2], $prmfitted[3],
                  1.0-$prmfitted[0]-$prmfitted[2] ,$prmfitted[4];
    close(OUTPUT);
    open(OUTPUT, ">$base\_fit_erro.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f \n", $datafile, $errorfitted[0],$errorfitted[1],$errorfitted[2],$errorfitted[3],
                  $errorfitted[0], $errorfitted[4];
    close(OUTPUT);

    system "mv fit3ExpOvr.dat $base\_prm.dat";
    system "mv fit3ExpOvr.ps $base.ps" ;
    system "mv fit3ExpOvr.log $base\_gnu.log";

}
elsif ($nexp == 4)
{
    print " Fit the overall rotation correlation function, $datafile, with four exponentials. \n"; 
    &Fit4Exp($ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8],$ARGV[9],$ARGV[10]);

    open(INPUT, "fit4ExpOvr.dat");
    my @prmfitted;
    my @errorfitted;
    while(<INPUT>) #read in data
    {
        chomp;
        my @tmp_array=split;
        if (defined $tmp_array[2] )
        {
            push (@prmfitted,$tmp_array[2]) ;
            push (@errorfitted,$tmp_array[4]) ;
        }
    }
    close(INPUT);

    printf "\n";
    printf "\n";
    printf "\n";
    printf " %15s %15s %15s %15s \n", "Coeff", "error", "Tau", "error";
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0],$prmfitted[1],$errorfitted[1];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[2],$errorfitted[2],$prmfitted[3],$errorfitted[3];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[4],$errorfitted[4],$prmfitted[5],$errorfitted[5];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0]-$prmfitted[2]-$prmfitted[4],$errorfitted[0],$prmfitted[6],$errorfitted[6];

    open(OUTPUT, ">$base\_fit.dat");
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0],$prmfitted[1],$errorfitted[1];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[2],$errorfitted[2],$prmfitted[3],$errorfitted[3];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[4],$errorfitted[4],$prmfitted[5],$errorfitted[5];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0]-$prmfitted[2]-$prmfitted[4],$errorfitted[0],$prmfitted[6],$errorfitted[6];
    close(OUTPUT);

    open(OUTPUT, ">$base\_fit_parm.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f \n", $datafile, $prmfitted[0],$prmfitted[1],$prmfitted[2], $prmfitted[3],
                  $prmfitted[4],$prmfitted[5],1.0-$prmfitted[0]-$prmfitted[2]-$prmfitted[4],$prmfitted[6];
    close(OUTPUT);
    open(OUTPUT, ">$base\_fit_erro.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f \n", $datafile, $errorfitted[0],$errorfitted[1],$errorfitted[2],$errorfitted[3],
                  $errorfitted[4],$errorfitted[5],$errorfitted[0], $errorfitted[6];
    close(OUTPUT);


    system "mv fit4ExpOvr.dat $base\_prm.dat";
    system "mv fit4ExpOvr.ps $base.ps" ;
    system "mv fit4ExpOvr.log $base\_gnu.log";

}
elsif ($nexp == 5)
{
    print " Fit the overall rotation correlation function, $datafile, with five exponentials. \n"; 
    &Fit5Exp($ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8],$ARGV[9],$ARGV[10],$ARGV[11],$ARGV[12]);

    open(INPUT, "fit5ExpOvr.dat");
    my @prmfitted;
    my @errorfitted;
    while(<INPUT>) #read in data
    {
        chomp;
        my @tmp_array=split;
        if (defined $tmp_array[2] )
        {
            push (@prmfitted,$tmp_array[2]) ;
            push (@errorfitted,$tmp_array[4]) ;
        }
    }
    close(INPUT);

    printf "\n";
    printf "\n";
    printf "\n";

    printf " %15s %15s  %15s %15s \n", "Coeff", "error", "Tau", "error";
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0],$prmfitted[1],$errorfitted[1];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[2],$errorfitted[2],$prmfitted[3],$errorfitted[3];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[4],$errorfitted[4],$prmfitted[5],$errorfitted[5];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[6],$errorfitted[6],$prmfitted[7],$errorfitted[7];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0]-$prmfitted[2]-$prmfitted[4]-$prmfitted[6],$errorfitted[0],$prmfitted[8],$errorfitted[8];

    open(OUTPUT, ">$base\_fit.dat");
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0],$prmfitted[1],$errorfitted[1];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[2],$errorfitted[2],$prmfitted[3],$errorfitted[3];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[4],$errorfitted[4],$prmfitted[5],$errorfitted[5];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[6],$errorfitted[6],$prmfitted[7],$errorfitted[7];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0]-$prmfitted[2]-$prmfitted[4]-$prmfitted[6],$errorfitted[0],$prmfitted[8],$errorfitted[8];
    close(OUTPUT);

    open(OUTPUT, ">$base\_fit_parm.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", $datafile, $prmfitted[0],$prmfitted[1],$prmfitted[2], $prmfitted[3],
                  $prmfitted[4],$prmfitted[5],$prmfitted[6],$prmfitted[7],1.0-$prmfitted[0]-$prmfitted[2]-$prmfitted[4]-$prmfitted[6],$prmfitted[8];
    close(OUTPUT);
    open(OUTPUT, ">$base\_fit_erro.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", $datafile, $errorfitted[0],$errorfitted[1],$errorfitted[2],$errorfitted[3],
                  $errorfitted[4],$errorfitted[5],$errorfitted[6],$errorfitted[7],$errorfitted[0], $errorfitted[8];
    close(OUTPUT);


    system "mv fit5ExpOvr.dat $base\_prm.dat";
    system "mv fit5ExpOvr.ps $base.ps" ;
    system "mv fit5ExpOvr.log $base\_gnu.log";

}
else 
{ die "nexp is not a valid number."}

#################################################################
sub Fit1Exp {

     # C(x) = exp(-x/t1)
   
    my $t1=$_[0];   # initial t1    
    my $tstart=$_[1]; # start time for fitting
    my $tend=$_[2];   # end time for fitting
    my $inputfile = $_[3]; # file to  fit 

    open(OUTPUT, ">fit1ExpOvr.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set stype data points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_O(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "t1 = $t1 \; \n";
    printf OUTPUT "f(x) = exp(-x/t1) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via t1 \n";
    printf OUTPUT "set output \"fit1ExpOvr.ps\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = exp(-x/t1) \", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);
 
    system "rm -rf fit.log";
    system "gnuplot fit1ExpOvr.gnu >& fit1ExpOvr.log";
    print " Please check fit1ExpOvr.log for gnuplot fitting.\n";
    system "cat fit1ExpOvr.log | grep -A 3 \"Final set of parameters\" | tail -n 1 > fit1ExpOvr.dat"; 

}
#################################################################
sub Fit2Exp {
    
    # C(x) = a1*exp(-x/t1) + (1-a1)*exp(-x/t2)

    my $a1=$_[0];   # initial a1 
    my $t1=$_[1];   # initial t1    
    my $t2=$_[2];   # initial t2
    my $tstart=$_[3]; # start time for fitting
    my $tend=$_[4];   # end time for fitting
    my $inputfile = $_[5]; # file to  fit 

    open(OUTPUT, ">fit2ExpOvr.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set style data points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_O(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "a1 = $a1 \; t1 = $t1 \; t2 = $t2\; \n";
    printf OUTPUT "f(x) = a1*exp(-x/t1) + (1-a1)*exp(-x/t2) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via a1, t1, t2 \n";
    printf OUTPUT "set output \"fit2ExpOvr.ps\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = a1*exp(-x/t1) + (1-a1)*exp(-x/t2) \", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);

    system "rm -rf fit.log";
    system "gnuplot fit2ExpOvr.gnu >& fit2ExpOvr.log";
    print " Please check fit2ExpOvr.log for gnuplot fitting.\n";
    system "cat fit2ExpOvr.log | grep -A 5 \"Final set of parameters\" | tail -n 3 > fit2ExpOvr.dat"; 

}

###################################################################
sub Fit3Exp {
    
    # C(x) = a1*exp(-x/t1) + a2*exp(-x/t2) + (1-a1-a2)*exp(-x/t3)
    my $a1=$_[0];   # initial a1 
    my $t1=$_[1];   # initial t1 
    my $a2=$_[2];   # initial a2 
    my $t2=$_[3];   # initial t2
    my $t3=$_[4];   # initial t3
    my $tstart=$_[5]; # start time for fitting
    my $tend=$_[6];   # end time for fitting
    my $inputfile = $_[7]; # file to  fit 

    open(OUTPUT, ">fit3ExpOvr.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set style data style points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_O(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "a1 = $a1 \; t1 = $t1 \; a2 = $a2 \; t2 = $t2 \; t3 = $t3 \; \n";
    printf OUTPUT "f(x) = a1*exp(-x/t1) + a2*exp(-x/t2) + (1-a1-a2)*exp(-x/t3) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via a1, t1, a2, t2, t3 \n";
    printf OUTPUT "set output \"fit3ExpOvr.ps\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = a1*exp(-x/t1) + a2*exp(-x/t2) + (1-a1-a2)*exp(-x/t3)\", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);

    system "rm -rf fit.log";
    system "gnuplot fit3ExpOvr.gnu >& fit3ExpOvr.log";
    print " Please check fit3ExpOvr.log for gnuplot fitting.\n";
    system "cat fit3ExpOvr.log | grep -A 7 \"Final set of parameters\" | tail -n 5 > fit3ExpOvr.dat"; 

}

###################################################################
sub Fit4Exp {
    
    # C(x) = a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + (1--a1-a2-a3)*exp(-x/t4)
    my $a1=$_[0];   # initial a1 
    my $t1=$_[1];   # initial t1 
    my $a2=$_[2];   # initial a2 
    my $t2=$_[3];   # initial t2
    my $a3=$_[4];   # initial a3 
    my $t3=$_[5];   # initial t3
    my $t4=$_[5];   # initial t4
    my $tstart=$_[7]; # start time for fitting
    my $tend=$_[8];   # end time for fitting
    my $inputfile = $_[9]; # file to  fit 

    open(OUTPUT, ">fit4ExpOvr.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set style data points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_O(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "a1 = $a1 \; t1 = $t1 \; a2 = $a2 \; t2 = $t2 \; a3 = $a3 \; t3 = $t3 \; t4 = $t4 \; \n";
    printf OUTPUT "f(x) = a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + (1-a1-a2-a3)*exp(-x/t4) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via a1, t1, a2, t2, a3, t3, t4  \n";
    printf OUTPUT "set output \"fit4ExpOvr.ps\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + (1-a1-a2-a3)*exp(-x/t4) \", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);

    system "rm -rf fit.log";
    system "gnuplot fit4ExpOvr.gnu >& fit4ExpOvr.log";
    print  " Please check fit4ExpOvr.log for gnuplot fitting.\n";
    system "cat fit4ExpOvr.log | grep -A 9 \"Final set of parameters\" | tail -n 7 > fit4ExpOvr.dat"; 

}

###################################################################
sub Fit5Exp {
    
    # C(x) = a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + a4*exp(-x/t4) + (1-a1-a2-a3-a4)*exp(-x/t5)
    my $a1=$_[0];   # initial a1 
    my $t1=$_[1];   # initial t1 
    my $a2=$_[2];   # initial a2 
    my $t2=$_[3];   # initial t2
    my $a3=$_[4];   # initial a3 
    my $t3=$_[5];   # initial t3
    my $a4=$_[6];   # initial a4 
    my $t4=$_[7];   # initial t4
    my $t5=$_[8];   # initial t5
    my $tstart=$_[9]; # start time for fitting
    my $tend=$_[10];   # end time for fitting
    my $inputfile = $_[11]; # file to  fit 

    open(OUTPUT, ">fit5ExpOvr.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set style data points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_O(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "a1 = $a1 \; t1 = $t1 \; a2 = $a2 \; t2 = $t2 \; a3 = $a3 \; t3 = $t3 \; a4 = $a4 \; t4 = $t4 \; t5 = $t5 \; \n";
    printf OUTPUT "f(x) = a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + a4*exp(-x/t4) + (1-a1-a2-a3-a4)*exp(-x/t5) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via a1, t1, a2, t2, a3, t3, a4, t4, t5  \n";
    printf OUTPUT "set output \"fit5ExpOvr.ps\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + a4*exp(-x/t4) + (1-a1-a2-a3-a4)*exp(-x/t5)\", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);

    system "rm -rf fit.log";
    system "gnuplot fit5ExpOvr.gnu >& fit5ExpOvr.log";
    print  " Please check fit5ExpOvr.log for gnuplot fitting.\n";
    system "cat fit5ExpOvr.log | grep -A 11 \"Final set of parameters\" | tail -n 9 > fit5ExpOvr.dat"; 

}
