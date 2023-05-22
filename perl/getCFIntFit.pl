#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;


my $myhead = <<end_of_myhead;
###########################################################################################################
# This a perl program to fit the correlation function of internal motion upto five exponential components.#
#                                                                                                         #
#  1) Ci(t) = a0+(1-a0)*exp(-t/t1)                                                                        # 
#  2) Ci(t) = a0+a1*exp(-t/t1)+(1-a0-a1)*exp(-t/t2)                                                       #
#  3) Ci(t) = a0+a1*exp(-t/t1)+a2*exp(-t/t2)+(1-a0-a1-a2)*exp(-t/t3)                                      #
#  4) Ci(t) = a0+a1*exp(-t/t1)+a2*exp(-t/t2)+a3*exp(-t/t3)+(1-a0-a1-a2-a3)*exp(-t/t4)                     #  
#  5) Ci(t) = a0+a1*exp(-t/t1)+a2*exp(-t/t2)+a3*exp(-t/t3)+a4*exp(-t/t4)+(1-a0-a1-a2-a3-a4)*exp(-t/t5)    # 
#                                                                                                         #
# Junchao Xia, junchao-xia\@biomaps.rutgers.edu, BioMaPs institute, Rutgers University.                    #
###########################################################################################################;
end_of_myhead
print $myhead;
print "\n";
print " NOTE: This program only fits one data file for the specified correlation function.\n";
print " NOTE: Please input the initial guess of parameters according the order of appearance above !!! \n";
print "\n";

#
my $usage = <<end_of_usage;
usage:"should input like this [command] [nexp] [parameters from initial guess] [tstart] [tend] [data file]"
end_of_usage

die $usage unless getopts("");

if ($#ARGV < 5 || $#ARGV > 13) {die $usage;}

my $nexp = $ARGV[0];

my @file_type = split(/\./,$ARGV[$#ARGV]); #splitting file to get extension
my $base = $file_type[0];
my $extension = $file_type[$#file_type]; 

my $datafile  = $ARGV[$#ARGV];
 
if ($nexp == 1)
{
    print " Fit the internal motion correlation function, $datafile, with one exponential. \n";
    &Fit1Exp($ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5]);
  
    open(INPUT, "fit1ExpInt.dat");
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
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0], 8888, 0;
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0],$errorfitted[0],$prmfitted[1],$errorfitted[1];
    open(OUTPUT, ">$base\_fit.dat");
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0],8888, 0;
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0],$errorfitted[0],$prmfitted[1],$errorfitted[1];
    close(OUTPUT);
    open(OUTPUT, ">$base\_fit_parm.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f \n" ,$datafile, $prmfitted[0], 1-$prmfitted[0], $prmfitted[1];
    close(OUTPUT);
    open(OUTPUT, ">$base\_fit_erro.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f \n" ,$datafile, $errorfitted[0],$errorfitted[0],$errorfitted[1];
    close(OUTPUT);
    system "mv fit1ExpInt.dat $base\_prm.dat";
    system "mv fit1ExpInt.ps $base.ps" ;
    system "mv fit1ExpInt.log $base\_gnu.log";

} 
elsif ($nexp == 2)  
{
    print " Fit the internal motion correlation function, $datafile, with two exponentials. \n";
    &Fit2Exp($ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7]);

    open(INPUT, "fit2ExpInt.dat");
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
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0], 8888, 0.0;
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[1],$errorfitted[1],$prmfitted[2],$errorfitted[2];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0]-$prmfitted[1],$errorfitted[0],$prmfitted[3],$errorfitted[3];

    open(OUTPUT, ">$base\_fit.dat");
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0], 8888, 0.0;
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[1],$errorfitted[1],$prmfitted[2],$errorfitted[2];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0]-$prmfitted[1],$errorfitted[0],$prmfitted[3],$errorfitted[3];
    close(OUTPUT); 

    open(OUTPUT, ">$base\_fit_parm.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f %15.6f\n",$datafile, $prmfitted[0], $prmfitted[1],$prmfitted[2],1.0-$prmfitted[0]-$prmfitted[1],$prmfitted[3];
    close(OUTPUT);
    open(OUTPUT, ">$base\_fit_erro.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f %15.6f\n" ,$datafile, $errorfitted[0],$errorfitted[1],$errorfitted[2],$errorfitted[0], $errorfitted[3];
    close(OUTPUT);

    system "mv fit2ExpInt.dat $base\_prm.dat";
    system "mv fit2ExpInt.ps $base.ps" ;
    system "mv fit2ExpInt.log $base\_gnu.log";

}
elsif ($nexp == 3)
{
    print " Fit the internal motion correlation function, $datafile, with three exponentials. \n";    
    &Fit3Exp($ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8],$ARGV[9]);

    open(INPUT, "fit3ExpInt.dat");
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
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0], 8888, 0.0;
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[1],$errorfitted[1],$prmfitted[2],$errorfitted[2];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[3],$errorfitted[3],$prmfitted[4],$errorfitted[4];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0]-$prmfitted[1]-$prmfitted[3],$errorfitted[0],$prmfitted[5],$errorfitted[5];

    open(OUTPUT, ">$base\_fit.dat");
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0], 8888, 0.0;
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[1],$errorfitted[1],$prmfitted[2],$errorfitted[2];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[3],$errorfitted[3],$prmfitted[4],$errorfitted[4];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0]-$prmfitted[1]-$prmfitted[3],$errorfitted[0],$prmfitted[5],$errorfitted[5];
    close(OUTPUT); 

    open(OUTPUT, ">$base\_fit_parm.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n" ,$datafile, $prmfitted[0],$prmfitted[1],$prmfitted[2], $prmfitted[3],
                  $prmfitted[4],1.0-$prmfitted[0]-$prmfitted[1]-$prmfitted[3] ,$prmfitted[5];
    close(OUTPUT);
    open(OUTPUT, ">$base\_fit_erro.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n" ,$datafile, $errorfitted[0],$errorfitted[1],$errorfitted[2],$errorfitted[3],
                  $errorfitted[4],$errorfitted[0], $errorfitted[5];
    close(OUTPUT);

    system "mv fit3ExpInt.dat $base\_prm.dat";
    system "mv fit3ExpInt.ps $base.ps" ;
    system "mv fit3ExpInt.log $base\_gnu.log";
    
}
elsif ($nexp == 4)
{
    print " Fit the internal motion correlation function, $datafile, with four exponentials. \n";    
    &Fit4Exp($ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8],$ARGV[9],$ARGV[10],$ARGV[11]);

    open(INPUT, "fit4ExpInt.dat");
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
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0], 8888, 0.0;
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[1],$errorfitted[1],$prmfitted[2],$errorfitted[2];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[3],$errorfitted[3],$prmfitted[4],$errorfitted[4];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[5],$errorfitted[5],$prmfitted[6],$errorfitted[6];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0]-$prmfitted[1]-$prmfitted[3]-$prmfitted[5],$errorfitted[0],$prmfitted[7],$errorfitted[7];

    open(OUTPUT, ">$base\_fit.dat");
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0], 8888, 0.0;
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[1],$errorfitted[1],$prmfitted[2],$errorfitted[2];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[3],$errorfitted[3],$prmfitted[4],$errorfitted[4];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[5],$errorfitted[5],$prmfitted[6],$errorfitted[6];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0]-$prmfitted[1]-$prmfitted[3]-$prmfitted[5],$errorfitted[0],$prmfitted[7],$errorfitted[7];
    close(OUTPUT);

    open(OUTPUT, ">$base\_fit_parm.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", $datafile, $prmfitted[0],$prmfitted[1],$prmfitted[2], $prmfitted[3],
                  $prmfitted[4],$prmfitted[5],$prmfitted[6],1.0-$prmfitted[0]-$prmfitted[1]-$prmfitted[3]-$prmfitted[5],$prmfitted[7];
    close(OUTPUT);
    open(OUTPUT, ">$base\_fit_erro.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", $datafile, $errorfitted[0],$errorfitted[1],$errorfitted[2],$errorfitted[3],
                  $errorfitted[4],$errorfitted[5],$errorfitted[6],$errorfitted[0], $errorfitted[7];
    close(OUTPUT);

    system "mv fit4ExpInt.dat $base\_prm.dat";
    system "mv fit4ExpInt.ps $base.ps" ;
    system "mv fit4ExpInt.log $base\_gnu.log";
 
}
elsif ($nexp == 5)
{
    print " Fit the internal motion correlation function, $datafile, with five exponentials. \n";   
    &Fit5Exp($ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8],$ARGV[9],$ARGV[10],$ARGV[11],$ARGV[12],$ARGV[13]);

    open(INPUT, "fit5ExpInt.dat");
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
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0], 8888, 0.0;
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[1],$errorfitted[1],$prmfitted[2],$errorfitted[2];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[3],$errorfitted[3],$prmfitted[4],$errorfitted[4];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[5],$errorfitted[5],$prmfitted[6],$errorfitted[6];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[7],$errorfitted[7],$prmfitted[8],$errorfitted[8];
    printf " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0]-$prmfitted[1]-$prmfitted[3]-$prmfitted[5]-$prmfitted[7],$errorfitted[0],$prmfitted[9],$errorfitted[9];

    open(OUTPUT, ">$base\_fit.dat");
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[0],$errorfitted[0], 8888, 0.0;
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[1],$errorfitted[1],$prmfitted[2],$errorfitted[2];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[3],$errorfitted[3],$prmfitted[4],$errorfitted[4];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[5],$errorfitted[5],$prmfitted[6],$errorfitted[6];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , $prmfitted[7],$errorfitted[7],$prmfitted[8],$errorfitted[8];
    printf OUTPUT " %15.6f %15.6f %15.6f %15.6f  \n" , 1.0-$prmfitted[0]-$prmfitted[1]-$prmfitted[3]-$prmfitted[5]-$prmfitted[7],$errorfitted[0],$prmfitted[9],$errorfitted[9];
    close(OUTPUT); 

    open(OUTPUT, ">$base\_fit_parm.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", $datafile, $prmfitted[0],$prmfitted[1],$prmfitted[2], $prmfitted[3],
                  $prmfitted[4],$prmfitted[5],$prmfitted[6],$prmfitted[7],$prmfitted[8],1.0-$prmfitted[0]-$prmfitted[1]-$prmfitted[3]-$prmfitted[5]-$prmfitted[7],$prmfitted[9];
    close(OUTPUT);
    open(OUTPUT, ">$base\_fit_erro.dat");
    printf OUTPUT "%20s %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", $datafile, $errorfitted[0],$errorfitted[1],$errorfitted[2],$errorfitted[3],
                  $errorfitted[4],$errorfitted[5],$errorfitted[6],$errorfitted[7],$errorfitted[8],$errorfitted[0], $errorfitted[9];
    close(OUTPUT);


    system "mv fit5ExpInt.dat $base\_prm.dat";
    system "mv fit5ExpInt.ps $base.ps" ;
    system "mv fit5ExpInt.log $base\_gnu.log";

}
else 
{ die "nexp is not a valid number."}

#################################################################
sub Fit1Exp {

     # C(x) = a0 + (1-a0)*exp(-x/t1)

    my $a0=$_[0];   # initial a0   
    my $t1=$_[1];   # initial t1   
    my $tstart=$_[2]; # start time for fitting
    my $tend=$_[3];   # end time for fitting
    my $inputfile = $_[4]; # file to  fit 

    open(OUTPUT, ">fit1ExpInt.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set style data points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_I(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "a0 = $a0 \; t1 = $t1 \; \n";
    printf OUTPUT "f(x) = a0+(1-a0)*exp(-x/t1) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via a0, t1 \n";
    printf OUTPUT "set output \"fit1ExpInt.ps\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = a0 + (1-a0)*exp(-x/t1) \", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);
 
    system "rm -rf fit.log";
    system "gnuplot fit1ExpInt.gnu >& fit1ExpInt.log";
    print " Please check fit1ExpInt.log for gnuplot fitting.\n";
    system "cat fit1ExpInt.log | grep -A 4 \"Final set of parameters\" | tail -n 2 > fit1ExpInt.dat"; 

}
#################################################################
sub Fit2Exp {
    
    # C(x) = a0 + a1*exp(-x/t1) + (1-a0-a1)*exp(-x/t2)

    my $a0=$_[0];   # initial a0
    my $a1=$_[1];   # initial a1 
    my $t1=$_[2];   # initial t1    
    my $t2=$_[3];   # initial t2
    my $tstart=$_[4]; # start time for fitting
    my $tend=$_[5];   # end time for fitting
    my $inputfile = $_[6]; # file to  fit 

    open(OUTPUT, ">fit2ExpInt.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set style data points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_I(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "a0 = $a0 \; a1 = $a1 \; t1 = $t1 \;  t2 = $t2\; \n";
    printf OUTPUT "f(x) = a0 + a1*exp(-x/t1) + (1-a0-a1)*exp(-x/t2) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via a0, a1, t1, t2 \n";
    printf OUTPUT "set output \"fit2ExpInt.ps\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = a0 + a1*exp(-x/t1) + (1-a0-a1)*exp(-x/t2) \", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);

    system "rm -rf fit.log";
    system "gnuplot fit2ExpInt.gnu >& fit2ExpInt.log";
    print " Please check fit2ExpInt.log for gnuplot fitting.\n";
    system "cat fit2ExpInt.log | grep -A 6 \"Final set of parameters\" | tail -n 4 > fit2ExpInt.dat"; 

}

###################################################################
sub Fit3Exp {
    
    # C(x) = a0 + a1*exp(-x/t1) + a2*exp(-x/t2) + (1-a0-a1-a2)*exp(-x/t3)
    my $a0=$_[0];   # initial a0
    my $a1=$_[1];   # initial a1 
    my $t1=$_[2];   # initial t1 
    my $a2=$_[3];   # initial a2 
    my $t2=$_[4];   # initial t2
    my $t3=$_[5];   # initial t3
    my $tstart=$_[6]; # start time for fitting
    my $tend=$_[7];   # end time for fitting
    my $inputfile = $_[8]; # file to  fit 

    open(OUTPUT, ">fit3ExpInt.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set style data points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_I(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "a0 = $a0 \; a1 = $a1 \; t1 = $t1 \; a2 = $a2 \; t2 = $t2 \; t3 = $t3 \; \n";
    printf OUTPUT "f(x) = a0 + a1*exp(-x/t1) + a2*exp(-x/t2) + (1-a0-a1-a2)*exp(-x/t3) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via a0, a1, t1, a2, t2, t3 \n";
    printf OUTPUT "set output \"fit3ExpInt.ps\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = a0 + a1*exp(-x/t1) + a2*exp(-x/t2) + (1-a0-a1-a2)*exp(-x/t3)\", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);

    system "rm -rf fit.log";
    system "gnuplot fit3ExpInt.gnu >& fit3ExpInt.log";
    print " Please check fit3ExpInt.log for gnuplot fitting.\n";
    system "cat fit3ExpInt.log | grep -A 8 \"Final set of parameters\" | tail -n 6 > fit3ExpInt.dat"; 

}

###################################################################
sub Fit4Exp {
    
    # C(x) = a0 + a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + (1-a0-a1-a2-a3)*exp(-x/t4)
    my $a0=$_[0];   # initial a0
    my $a1=$_[1];   # initial a1 
    my $t1=$_[2];   # initial t1 
    my $a2=$_[3];   # initial a2 
    my $t2=$_[4];   # initial t2
    my $a3=$_[5];   # initial a3 
    my $t3=$_[6];   # initial t3
    my $t4=$_[7];   # initial t4
    my $tstart=$_[8]; # start time for fitting
    my $tend=$_[9];   # end time for fitting
    my $inputfile = $_[10]; # file to  fit 

    open(OUTPUT, ">fit4ExpInt.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set style data points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_I(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "a0 = $a0 \; a1 = $a1 \; t1 = $t1 \; a2 = $a2 \; t2 = $t2 \; a3 = $a3 \; t3 = $t3 \; t4 = $t4 \; \n";
    printf OUTPUT "f(x) = a0 + a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + (1-a0-a1-a2-a3)*exp(-x/t4) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via a0, a1, t1, a2, t2, a3, t3, t4  \n";
    printf OUTPUT "set output \"fit4ExpInt.ps\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = a0 + a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + (1-a0-a1-a2-a3)*exp(-x/t4) \", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);

    system "rm -rf fit.log";
    system "gnuplot fit4ExpInt.gnu >& fit4ExpInt.log";
    print  " Please check fit4ExpInt.log for gnuplot fitting.\n";
    system "cat fit4ExpInt.log | grep -A 10 \"Final set of parameters\" | tail -n 8 > fit4ExpInt.dat"; 

}

###################################################################
sub Fit5Exp {
    
    # C(x) = a0 + a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + a4*exp(-x/t4) + (1-a0-a1-a2-a3-a4)*exp(-x/t5)
    my $a0=$_[0];   # initial a0
    my $a1=$_[1];   # initial a1 
    my $t1=$_[2];   # initial t1 
    my $a2=$_[3];   # initial a2 
    my $t2=$_[4];   # initial t2
    my $a3=$_[5];   # initial a3 
    my $t3=$_[6];   # initial t3
    my $a4=$_[7];   # initial a4 
    my $t4=$_[8];   # initial t4
    my $t5=$_[9];   # initial t5
    my $tstart=$_[10]; # start time for fitting
    my $tend=$_[11];   # end time for fitting
    my $inputfile = $_[12]; # file to  fit 

    open(OUTPUT, ">fit5ExpInt.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set style data points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_I(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "a0 = $a0 \; a1 = $a1 \; t1 = $t1 \; a2 = $a2 \; t2 = $t2 \; a3 = $a3 \; t3 = $t3 \; a4 = $a4 \; t4 = $t4 \; t5 = $t5 \; \n";
    printf OUTPUT "f(x) = a0 + a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + a4*exp(-x/t4) + (1-a0-a1-a2-a3-a4)*exp(-x/t5) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via a0, a1, t1, a2, t2, a3, t3, a4, t4, t5  \n";
    printf OUTPUT "set output \"fit5ExpInt.ps\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = a0 + a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + a4*exp(-x/t4) + (1-a0-a1-a2-a3-a4)*exp(-x/t5)\", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);

    system "rm -rf fit.log";
    system "gnuplot fit5ExpInt.gnu >& fit5ExpInt.log";
    print  " Please check fit5ExpInt.log for gnuplot fitting.\n";
    system "cat  fit5ExpInt.log | grep -A 12 \"Final set of parameters\" | tail -n 10 > fit5ExpInt.dat"; 

}
