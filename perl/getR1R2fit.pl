#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;

#
my $myhead = <<end_of_myhead;
###########################################################################################################
# This a perl program to calcuate the relaxation parameters R1/R2 from the C2 correlation function of     #
# C13-H1 or N15-H1 vector obtained from the AMBER ptraj analysis. C2 is fitted to multi-exponential       #
# functions via Gnuplot. The fitting parameters are passed to get the sectral density function and the    #
# R1/R2 by the method described in the papers:                                                            #
# 1) N. L. Fawzi, A. H. Phillips, J. Z. Ruscio, M. Doucleff, D. E. Wemmer, and T. Head-Gordon,            #
#    J. Am. Chem. Soc. 130, 6145 (2008);                                                                  #
# 2) A. G. Palmer III, M. Rance, and P. E. Wright, J. Am. Chem. Soc. 113 4371 (1991).                     #
#                                                                                                         #
# The total correlation function is fitted to different forms as below (nexp = ?):                        #           
# 0) Ct(t) = (S2+(1-S2)*exp(-t/tau_e))*exp(-t/tau_m)=S2*exp(-t/tau_m)+(1-S2)*exp(-t/tau) mean-field model # 
# 1) Ct(t) = exp(-t/t1)                                                                                   # 
# 2) Ct(t) = a1*exp(-t/t1) + (1-a1)*exp(-t/t2)                                                            #
# 3) Ct(t) = a1*exp(-t/t1) + a2*exp(-t/t2) + (1-a1-a2)*exp(-t/t3)                                         #
# 4) Ct(t) = a1*exp(-t/t1) + a2*exp(-t/t2) + a3*exp(-t/t3) + (1-a1-a2-a3)*exp(-t/t4)                      #  
# 5) Ct(t) = a1*exp(-t/t1) + a2*exp(-t/t2) + a3*exp(-t/t3) + a4*exp(-t/t4) + (1-a0-a1-a2-a3-a4)*exp(-t/t5)#   
#                                                                                                         #
# Junchao Xia, junchao-xia\@biomaps.rutgers.edu, BioMaPs institute, Rutgers University.                    #
###########################################################################################################
end_of_myhead

print $myhead;
print "\n";
print "\n";
print " NOTE: This program only fits one data file for the specified correlation function. \n";
print " NOTE: After some modifications this program also can be used to calculate relaxation \n";
print "       parameters for other protons besides C13 and N15. \n";
print "\n";
print " NOTE: Please check the parameters defined in the program have the right values!!!!!! \n"; 
print "\n";

my $usage = <<end_of_usage;
usage:"should input like this [command] [Bfield] [spin type C or N] [nexp] [parameters from initial guess] [tstart] [tend] [input file]"
end_of_usage

die $usage unless getopts("");

if ($#ARGV < 6 || $#ARGV > 12) {die $usage;}

my $Bfield  = $ARGV[0];         # B field in T  (9.4)
my $SpinType = $ARGV[1];        # Spin type N or C 

my $GAMMA_i;
if ($SpinType eq "N")
{$GAMMA_i = -27.1e6;}       # rad/sT, default value for N15 (found typo in Mar 15, 2012, no negtive before)
elsif  ($SpinType eq "C")
{$GAMMA_i = 67.2274e6;}    # rad/sT, default value for C13
else
{
    print "Spin type is not right.";
    exit(1);
}

my $GAMMA_j = 2.67522128e8;     # rad/sT, default value for H1
my $H_BAR   = 1.05457266e-34;   # Js
my $MU_0    = 12.56e-7;         # Vs/(Am)

my $Reff;
if ($SpinType eq "N")
{$Reff    = 1.01e-10;}        # m,  default value for N-H vector 
elsif  ($SpinType eq "C")
{$Reff    = 1.11e-10;}        # m,  default value for C-H vector
else
{
    print "Spin type is not right.";
    exit(1);
}

my $PI      = 3.14159;        # 
    
my $K  =  $MU_0*$H_BAR*$GAMMA_i*$GAMMA_j/(4*$PI*$Reff*$Reff*$Reff);  # note the power of 3 is missing in the paper  


my $OMEGA_i;
if ($SpinType eq "N")
{$OMEGA_i = 2.0*$PI*4.315e6*$Bfield;}  # Hz, default for Omega_N  (found typo in Mar 15, 2012, set to negtive before)
elsif  ($SpinType eq "C")
{$OMEGA_i = 2.0*$PI*10.705e6*$Bfield;}  # Hz, default for Omega_C
else
{
    print "Spin type is not right.";
    exit(1);
}

my $OMEGA_j = 2.0*$PI*42.577e6*$Bfield;    # Hz, default for Omega_H


my $DELTA_i;              #  Chemical-shift anisotropy  

if ($SpinType eq "N")
{$DELTA_i  = 170/1.0e6;}  #  Chemical-shift anisotropy of N (found typo in Mar 19, 2012, 1.0e6 was 10e6)
elsif  ($SpinType eq "C")
{$DELTA_i  = 25/1.0e6;}  #  Chemical-shift anisotropy of C  Need check the exact value if it is important
else
{
    print "Spin type is not right.";
    exit(1);
}

my $R1DD;
my $R1CSA;
my $R1;
my $T1;
my $R2DD;
my $R2CSA;
my $R2a;
my $R2;
my $T2;
my $NOE;

my $Jw0;
my $JwjSubi;
my $Jwi;
my $Jwj;
my $JwjPlsi;


my $nexp = $ARGV[2];

my @file_type = split(/\./,$ARGV[$#ARGV]); #splitting file to get extension
my $base = $file_type[0];
my $extension = $file_type[$#file_type]; 

my $datafile  = $ARGV[$#ARGV];
 
if ($nexp == 0)
{&FitMF0Exp($ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8]);} 
elsif ($nexp == 1)  
{&Fit1Exp($ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6]);}
elsif ($nexp == 2)
{&Fit2Exp($ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8]);}
elsif ($nexp == 3)
{&Fit3Exp($ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8],$ARGV[9],$ARGV[10]);}
elsif ($nexp == 4)
{&Fit4Exp($ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8],$ARGV[9],$ARGV[10],$ARGV[11],$ARGV[12]);}
elsif ($nexp == 5)
{&Fit5Exp($ARGV[3],$ARGV[4],$ARGV[5],$ARGV[6],$ARGV[7],$ARGV[8],$ARGV[9],$ARGV[10],$ARGV[11],$ARGV[12],$ARGV[13],$ARGV[14]);}
else 
{ die "nexp is not a valid number."}

if ($nexp == 0)
{open(INPUT, "fitMF0ExpTot.dat");}
elsif ($nexp == 1)
{open(INPUT, "fit1ExpTot.dat");}    
elsif ($nexp == 2)
{open(INPUT, "fit2ExpTot.dat");} 
elsif ($nexp == 3)
{open(INPUT, "fit3ExpTot.dat");}
elsif ($nexp == 4)
{open(INPUT, "fit4ExpTot.dat");}
elsif ($nexp == 5)
{open(INPUT, "fit5ExpTot.dat");}
else 
{ die "nexp is not a valid number."}

my $fitSuccess = 0;  

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

if (defined $prmfitted[0])
{print "Gnuplot fitting is sucessful.\n";}
else 
{
    print " Got error when reading the fitting data using Gnuplot.\n"; 
    $fitSuccess = 1; 
}

if ($nexp == 0)
{
    if ($fitSuccess == 0)
    {
	# change the unit of ps to s
	$prmfitted[1] = $prmfitted[1]/1.0e12;
	$prmfitted[2] = $prmfitted[2]/1.0e12;


	$Jw0 = JwMF0Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],0.0);
	$JwjSubi = JwMF0Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$OMEGA_j-$OMEGA_i);
	$Jwi = JwMF0Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$OMEGA_i);
	$Jwj = JwMF0Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$OMEGA_j);
	$JwjPlsi = JwMF0Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$OMEGA_j+$OMEGA_i);

	&calculateR1R2NOE();

	my $S2 = 1.0-$prmfitted[0];
	my $tau_m = $prmfitted[2]*1.0e12; # change the unit of s to ps
	my $tau = $prmfitted[1]*1.0e12;   # change the unit of s to ps

	my $tau_e = 1.0/(1.0/$tau - 1.0/$tau_m);

	print " The normalized time correlation function C2(t) is obtained from the AMBER ptraj calculation. \n";
	print " C2(t) is fitted to a two-exponential function in the model-free approach of isotropic over   \n";
	print " motion, C2(t)=(S2+(1-S2)*exp(-t/tau_e))*exp(-t/tau_m)=S2*exp(-t/tau_m)+(1-S2)*exp(-t/tau),   \n";
	print " where tau^-1 = tau_m^-1 + tau_e^-1.                                                        \n\n";
	print " Final results for the calculations of relaxation parameters paramters as listed follows:     \n";
	
	printf " %10s %10s %10s %10s %10s \n", "Jw0", "JwjSubi", "Jwi", "Jwj", "JwjPlsi";
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f \n\n",$Jw0,$JwjSubi,$Jwi,$Jwj,$JwjPlsi;
	
	printf " %10s %10s %10s %10s %10s %10s %10s \n",
	"S2 value", "error", "Tau_e", " Tau_m", "error", "Tau", "error";
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",  
	$S2,$errorfitted[0],$tau_e,$tau_m,$errorfitted[1],$tau,$errorfitted[2];
       
	printf " %10s %10s %10s %10s %10s %10s %10s %10s \n", "R1DD", "R1CSA", "R2DD", "R2CSA", "R1", "R2", "R2/R1", "NOE";
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n", $R1DD,$R1CSA,$R2DD,$R2CSA,$R1,$R2,$R2/$R1,$NOE;

	open(OUTPUT, ">$base\_fit.dat");
	printf OUTPUT " %20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
          "file name", "S2 value", "error","Tau_e", "Tau_m", "error", "Tau", "error", "R1DD", "R1CSA", "R2DD", "R2CSA", "R1", "R2", "R2/R1", "NOE"; 
	printf OUTPUT " %20s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",  
	$datafile,$S2,$errorfitted[0],$tau_e,$tau_m,$errorfitted[1],$tau,$errorfitted[2],$R1DD,$R1CSA,$R2DD,$R2CSA,$R1,$R2,$R2/$R1,$NOE;
	close(OUTPUT);

	system "mv fitMF0ExpTot.ps $base.ps";
    }
    system "mv fitMF0ExpTot.dat $base\_prm.dat";
    system "mv fitMF0ExpTot.log $base\_gnu.log";
}
elsif ($nexp == 1)
{
    if ($fitSuccess == 0)
    {
	# change the unit of ps to s
	$prmfitted[0] = $prmfitted[0]/1.0e12;

	$Jw0 = Jw1Exp($prmfitted[0],0.0);
	$JwjSubi = Jw1Exp($prmfitted[0],$OMEGA_j-$OMEGA_i);
	$Jwi = Jw1Exp($prmfitted[0],$OMEGA_i);
	$Jwj = Jw1Exp($prmfitted[0],$OMEGA_j);
	$JwjPlsi = Jw1Exp($prmfitted[0],$OMEGA_j+$OMEGA_i);

	&calculateR1R2NOE();

	print " The normalized time correlation function C2(t) is obtained from the AMBER ptraj calculation.\n";
	print " C2(t) is fitted to a one-exponential function, Ct(t) = exp(-t/t1).                        \n\n";
	print " Final results for the calculations of relaxation parameters paramters as listed follows: \n";
	
	printf " %10s %10s %10s %10s %10s \n", "Jw0", "JwjSubi", "Jwi", "Jwj", "JwjPlsi";
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f \n",$Jw0,$JwjSubi,$Jwi,$Jwj,$JwjPlsi;
	
	printf " %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n",
	"Tau1", "error", "R1DD", "R1CSA", "R2DD", "R2CSA", "R1", "R2", "R2/R1", "NOE"; 
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",  
	$prmfitted[0]*1.0e12,$errorfitted[0],$R1DD,$R1CSA,$R2DD,$R2CSA,$R1,$R2,$R2/$R1,$NOE;

	open(OUTPUT, ">$base\_fit.dat");
	printf OUTPUT " %20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n",
	"file name", "Tau 1", "error", "R1DD", "R1CSA", "R2DD", "R2CSA", "R1", "R2", "R2/R1", "NOE"; 
	printf OUTPUT " %20s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",  
	$datafile,$prmfitted[0]*1.0e12,$errorfitted[0],$R1DD,$R1CSA,$R2DD,$R2CSA,$R1,$R2,$R2/$R1,$NOE;
	close(OUTPUT);

	system "mv fit1ExpTot.ps $base.ps" ;
    }
    system "mv fit1ExpTot.dat $base\_prm.dat";
    system "mv fit1ExpTot.log $base\_gnu.log";
}
elsif ($nexp == 2)
{
    if ($fitSuccess == 0)
    {
	# change the unit of ps to s
	$prmfitted[1] = $prmfitted[1]/1.0e12;
	$prmfitted[2] = $prmfitted[2]/1.0e12;

	$Jw0 = Jw2Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],0.0);
	$JwjSubi = Jw2Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$OMEGA_j-$OMEGA_i);
	$Jwi = Jw2Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$OMEGA_i);
	$Jwj = Jw2Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$OMEGA_j);
	$JwjPlsi = Jw2Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$OMEGA_j+$OMEGA_i);

	&calculateR1R2NOE();

	print " The normalized time correlation function C2(t) is obtained from the AMBER ptraj calculation. \n";
	print " C2(t) is fitted to a two-exponential function, Ct(t) = a1*exp(-t/t1) + (1-a1)*exp(-t/t2.   \n\n";
	print " Final results for the calculations of relaxation parameters paramters as listed follows:     \n";

	printf " %10s %10s %10s %10s %10s \n", "Jw0", "JwjSubi", "Jwi", "Jwj", "JwjPlsi";
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f \n",$Jw0,$JwjSubi,$Jwi,$Jwj,$JwjPlsi;

	printf " %10s %10s %10s %10s %10s %10s \n", "A1  ", "error", " Tau1", "error", "Tau2", "error";
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",  
	$prmfitted[0],$errorfitted[0],$prmfitted[1]*1.0e12,$errorfitted[1],$prmfitted[2]*1.0e12, $errorfitted[2];
    
	printf " %10s %10s %10s %10s %10s %10s %10s %10s \n", "R1DD", "R1CSA", "R2DD", "R2CSA", "R1", "R2", "R2/R1", "NOE";
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n", $R1DD,$R1CSA,$R2DD,$R2CSA,$R1,$R2,$R2/$R1,$NOE;

	open(OUTPUT, ">$base\_fit.dat");
	printf OUTPUT " %20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n",
	"file name", "A1", "error", "Tau1", "error","A2", "Tau2", "error", "R1DD", "R1CSA", "R2DD", "R2CSA", "R1", "R2", "R2/R1", "NOE"; 
	printf OUTPUT " %20s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",  
	$datafile, $prmfitted[0],$errorfitted[0],$prmfitted[1]*1.0e12,$errorfitted[1], 1.0-$prmfitted[0],
	$prmfitted[2]*1.0e12, $errorfitted[2],$R1DD,$R1CSA,$R2DD,$R2CSA,$R1,$R2,$R2/$R1,$NOE;
	close(OUTPUT);

	system "mv fit2ExpTot.ps $base.ps" ;
    }

    system "mv fit2ExpTot.dat $base\_prm.dat";
    system "mv fit2ExpTot.log $base\_gnu.log";
}
elsif ($nexp == 3)
{
    if ($fitSuccess == 0)
    {
	# change the unit of ps to s
	$prmfitted[1] = $prmfitted[1]/1.0e12;
	$prmfitted[3] = $prmfitted[3]/1.0e12;
	$prmfitted[4] = $prmfitted[4]/1.0e12;
	
	$Jw0 = Jw3Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$prmfitted[3],$prmfitted[4],0.0);
	$JwjSubi = Jw3Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$prmfitted[3],$prmfitted[4],$OMEGA_j-$OMEGA_i);
	$Jwi = Jw3Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$prmfitted[3],$prmfitted[4],$OMEGA_i);
	$Jwj = Jw3Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$prmfitted[3],$prmfitted[4],$OMEGA_j);
	$JwjPlsi = Jw3Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$prmfitted[3],$prmfitted[4],$OMEGA_j+$OMEGA_i);

	&calculateR1R2NOE();

	print " The normalized time correlation function C2(t) is obtained from the AMBER ptraj calculation.            \n";
	print " C2(t) is fitted to a three-exponential function, Ct(t) = a1*exp(-t/t1) + a2*exp(-t/t2) + (1-a1-a2)*exp(-t/t3). \n\n";
	print " Final results for the calculations of relaxation parameters paramters as listed follows:                \n";

	printf " %10s %10s %10s %10s %10s \n", "Jw0", "JwjSubi", "Jwi", "Jwj", "JwjPlsi";
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f \n",$Jw0,$JwjSubi,$Jwi,$Jwj,$JwjPlsi;
 
	printf " %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n",
	"A1", "error", "Tau1", "error", "A2", "error", "Tau2", "error", "Tau3", "error"; 
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",  
	$prmfitted[0],$errorfitted[0],$prmfitted[1]*1.0e12,$errorfitted[1],$prmfitted[2],$errorfitted[2],$prmfitted[3]*1.0e12, $errorfitted[3],
	$prmfitted[4]*1.0e12,$errorfitted[4];
  
	printf " %10s %10s %10s %10s %10s %10s %10s %10s \n","R1DD", "R1CSA", "R2DD", "R2CSA", "R1", "R2", "R2/R1", "NOE"; 
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",$R1DD,$R1CSA,$R2DD,$R2CSA,$R1,$R2,$R2/$R1,$NOE;

	open(OUTPUT, ">$base\_fit.dat");
	printf OUTPUT " %20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n",
	"file name", "A1", "error", "Tau1", "error", "A2", "error", "Tau2", "error", "A3", "Tau3", "error", "R1DD", "R1CSA", "R2DD", "R2CSA", "R1", "R2", "R2/R1", "NOE"; 
	printf OUTPUT " %20s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",  
	$datafile,$prmfitted[0],$errorfitted[0],$prmfitted[1]*1.0e12,$errorfitted[1],$prmfitted[2],$errorfitted[2],$prmfitted[3]*1.0e12, $errorfitted[3],
	1.0-$prmfitted[0]-$prmfitted[2], $prmfitted[4]*1.0e12,$errorfitted[4],$R1DD,$R1CSA,$R2DD,$R2CSA,$R1,$R2,$R2/$R1,$NOE;
	close(OUTPUT);
	system "mv fit3ExpTot.ps $base.ps" ; 
    }

    system "mv fit3ExpTot.dat $base\_prm.dat";
    system "mv fit3ExpTot.log $base\_gnu.log";
}
elsif ($nexp == 4)
{
    if ($fitSuccess == 0)
    {
	# change the unit of ps to s
	$prmfitted[1] = $prmfitted[1]/1.0e12;
	$prmfitted[3] = $prmfitted[3]/1.0e12;
	$prmfitted[5] = $prmfitted[5]/1.0e12;
	$prmfitted[6] = $prmfitted[6]/1.0e12;

	$Jw0 = Jw4Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$prmfitted[3],$prmfitted[4],$prmfitted[5],$prmfitted[6],0.0);
	$JwjSubi = Jw4Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$prmfitted[3],$prmfitted[4],$prmfitted[5],$prmfitted[6],$OMEGA_j-$OMEGA_i);
	$Jwi = Jw4Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$prmfitted[3],$prmfitted[4],$prmfitted[5],$prmfitted[6],$OMEGA_i);
	$Jwj = Jw4Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$prmfitted[3],$prmfitted[4],$prmfitted[5],$prmfitted[6],$OMEGA_j);
	$JwjPlsi = Jw4Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$prmfitted[3],$prmfitted[4],$prmfitted[5],$prmfitted[6],$OMEGA_j+$OMEGA_i);

	&calculateR1R2NOE();

	print " The normalized time correlation function C2(t) is obtained from the AMBER ptraj calculation. C2(t) is fitted to a \n";
	print " four-exponential function, Ct(t) = a1*exp(-t/t1) + a2*exp(-t/t2) + a3*exp(-t/t3) + (1-a1-a2-a3)*exp(-t/t4).     \n\n";
	print " Final results for the calculations of relaxation parameters paramters as listed follows:                          \n";
 
	printf " %10s %10s %10s %10s %10s \n", "Jw0", "JwjSubi", "Jwi", "Jwj", "JwjPlsi";
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f \n",$Jw0,$JwjSubi,$Jwi,$Jwj,$JwjPlsi;

	printf " %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n",
	"A1", "error", "Tau1", "error", "A2", "error", "Tau2", "error", "A3", "error", "Tau3", "error", "Tau4", "error"; 
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",  
	$prmfitted[0],$errorfitted[0],$prmfitted[1]*1.0e12,$errorfitted[1],$prmfitted[2],$errorfitted[2],$prmfitted[3]*1.0e12, $errorfitted[3],
	$prmfitted[4],$errorfitted[4],$prmfitted[5]*1.0e12,$errorfitted[5],$prmfitted[6]*1.0e12,$errorfitted[6];
       
	printf " %10s %10s %10s %10s %10s %10s %10s %10s \n", "R1DD", "R1CSA", "R2DD", "R2CSA", "R1", "R2", "R2/R1", "NOE"; 
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",$R1DD,$R1CSA,$R2DD,$R2CSA,$R1,$R2,$R2/$R1,$NOE;

	open(OUTPUT, ">$base\_fit.dat");
	printf OUTPUT " %20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s  \n",
	"file name",  "A1", "error", "Tau1", "error", "A2", "error", "Tau2", "error", "A3", "error", "Tau3", "error", "A4", "Tau4", "error",
	"R1DD", "R1CSA", "R2DD", "R2CSA", "R1", "R2", "R2/R1", "NOE"; 
	printf OUTPUT " %20s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",  
	$datafile,$prmfitted[0],$errorfitted[0],$prmfitted[1]*1.0e12,$errorfitted[1],$prmfitted[2],$errorfitted[2],$prmfitted[3]*1.0e12, $errorfitted[3],
	$prmfitted[4],$errorfitted[4],$prmfitted[5]*1.0e12,$errorfitted[5], 1-$prmfitted[0]-$prmfitted[2]-$prmfitted[4],
	$prmfitted[6]*1.0e12,$errorfitted[6],$R1DD,$R1CSA,$R2DD,$R2CSA,$R1,$R2,$R2/$R1,$NOE;
	close(OUTPUT);

	system "mv fit4ExpTot.ps $base.ps" ;
    }
    system "mv fit4ExpTot.dat $base\_prm.dat";
    system "mv fit4ExpTot.log $base\_gnu.log";

}
elsif ($nexp == 5)
{
    if ($fitSuccess == 0)
    {
	# change the unit of ps to s
	$prmfitted[1] = $prmfitted[1]/1.0e12;
	$prmfitted[3] = $prmfitted[3]/1.0e12;
	$prmfitted[5] = $prmfitted[5]/1.0e12;
	$prmfitted[7] = $prmfitted[7]/1.0e12;
	$prmfitted[8] = $prmfitted[8]/1.0e12;

	$Jw0 = Jw5Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$prmfitted[3],$prmfitted[4],$prmfitted[5],$prmfitted[6],$prmfitted[7],$prmfitted[8],0.0);
	$JwjSubi = Jw5Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$prmfitted[3],$prmfitted[4],$prmfitted[5],$prmfitted[6],$prmfitted[7],$prmfitted[8],$OMEGA_j-$OMEGA_i);
	$Jwi = Jw5Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$prmfitted[3],$prmfitted[4],$prmfitted[5],$prmfitted[6],$prmfitted[7],$prmfitted[8],$OMEGA_i);
	$Jwj = Jw5Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$prmfitted[3],$prmfitted[4],$prmfitted[5],$prmfitted[6],$prmfitted[7],$prmfitted[8],$OMEGA_j);
	$JwjPlsi = Jw5Exp($prmfitted[0],$prmfitted[1],$prmfitted[2],$prmfitted[3],$prmfitted[4],$prmfitted[5],$prmfitted[6],$prmfitted[7],$prmfitted[8],$OMEGA_j+$OMEGA_i);

	&calculateR1R2NOE();

	print " The normalized time correlation function C2(t) is obtained from the AMBER ptraj calculation. C2(t) is fitted to a \n";
	print " four-exponential function, Ct(t) = a1*exp(-t/t1) + a2*exp(-t/t2) + a3*exp(-t/t3) + a4*exp(-t/t4) + (1-a0-a1-a2-a3-a4)*exp(-t/t5).\n\n";
	print " Final results for the calculations of relaxation parameters paramters as listed follows:                          \n";
 
	printf " %10s %10s %10s %10s %10s \n", "Jw0", "JwjSubi", "Jwi", "Jwj", "JwjPlsi";
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f \n",$Jw0,$JwjSubi,$Jwi,$Jwj,$JwjPlsi;

	printf " %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n",
	"A1", "error", "Tau1", "error", "A2", "error", "Tau2", "error", "A3", "error", "Tau3", "error", "A4", "error", "Tau4", "error", "Tau5", "error"; 
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",  
	$prmfitted[0],$errorfitted[0],$prmfitted[1]*1.0e12,$errorfitted[1],$prmfitted[2],$errorfitted[2],$prmfitted[3]*1.0e12, $errorfitted[3],
	$prmfitted[4],$errorfitted[4],$prmfitted[5]*1.0e12,$errorfitted[5],$prmfitted[7]*1.0e12,$errorfitted[7],$prmfitted[8]*1.0e12,$errorfitted[8];
       
	printf " %10s %10s %10s %10s %10s %10s %10s %10s \n", "R1DD", "R1CSA", "R2DD", "R2CSA", "R1", "R2", "R2/R1", "NOE"; 
	printf " %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",$R1DD,$R1CSA,$R2DD,$R2CSA,$R1,$R2,$R2/$R1,$NOE;

	open(OUTPUT, ">$base\_fit.dat");
	printf OUTPUT " %20s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n",
	"file name",  "A1", "error", "Tau1", "error", "A2", "error", "Tau2", "error", "A3", "error", "Tau3", "error", "A4","error", "Tau4", "error", "A5", "Tau5", "error",
	"R1DD", "R1CSA", "R2DD", "R2CSA", "R1", "R2", "R2/R1", "NOE"; 
	printf OUTPUT " %20s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n",  
	$datafile,$prmfitted[0],$errorfitted[0],$prmfitted[1]*1.0e12,$errorfitted[1],$prmfitted[2],$errorfitted[2],$prmfitted[3]*1.0e12, $errorfitted[3],
	$prmfitted[4],$errorfitted[4],$prmfitted[5]*1.0e12,$errorfitted[5],$prmfitted[6],$errorfitted[6],$prmfitted[7]*1.0e12,$errorfitted[7],
	1-$prmfitted[0]-$prmfitted[2]-$prmfitted[4]-$prmfitted[6],$prmfitted[8]*1.0e12,$errorfitted[8],$R1DD,$R1CSA,$R2DD,$R2CSA,$R1,$R2,$R2/$R1,$NOE;
	close(OUTPUT);

	system "mv fit5ExpTot.ps $base.ps" ;
    }
    system "mv fit5ExpTot.dat $base\_prm.dat";
    system "mv fit5ExpTot.log $base\_gnu.log";

}

else 
{die "nexp is not a valid number.";}

#################################################################
sub calculateR1R2NOE{    
    # note in the  Fawzi's paper an additional 5 is divided for R1DD, R2DD, R1CSA, and R2CSA
    $R1DD = $K*$K*($JwjSubi + 3.0*$Jwi+6.0*$JwjPlsi)/20.0;
    $R2DD = $K*$K*(4.0*$Jw0 + $JwjSubi + 3.0*$Jwi + 6.0*$Jwj + 6.0*$JwjPlsi)/40.0;  # note 6.0 is 3.0 in the paper
    $R1CSA =$DELTA_i*$DELTA_i*$OMEGA_i*$OMEGA_i*$Jwi/5.0;
    $R2CSA =$DELTA_i*$DELTA_i*$OMEGA_i*$OMEGA_i*(4.0*$Jw0 + 3.0*$Jwi)/30.0;
    $R2a = 0.0;
    $R1 = $R1DD+$R1CSA;
    $R2 = $R2DD+$R2CSA+$R2a;
    $T1 = 1.0/$R1;
    $T2 = 1.0/$R2;
    $NOE = 1 + ($GAMMA_j/$GAMMA_i)*(6.0*$JwjPlsi - $JwjSubi)/
            ($JwjSubi + (3.0 + 4.0*$DELTA_i*$DELTA_i*$OMEGA_i*$OMEGA_i/($K*$K))*$Jwi + 6.0*$JwjPlsi);

    # change the unit of s to ps 
    $Jw0 = $Jw0*1.0e12;
    $JwjSubi = $JwjSubi*1.0e12;
    $Jwi = $Jwi*1.0e12;
    $Jwj = $Jwj*1.0e12;
    $JwjPlsi = $JwjPlsi*1.0e12;
}

#################################################################
sub FitMF0Exp {
    
    # C(x) = a*exp(-x/b) + (1-a)*exp(-x/c)

    my $coeff1=$_[0]; # initial a
    my $tau1=$_[1];   # initial b    
    my $tau2=$_[2];   # initial c
    my $tstart=$_[3]; # start time for fitting
    my $tend=$_[4];   # end time for fitting
    my $inputfile = $_[5]; # file to  fit 

    open(OUTPUT, ">fitMF0ExpTot.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set style data points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_T(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "set yrange\[0.0:1.0\] \n"; 
    printf OUTPUT "a = $coeff1 \; b = $tau1 \; c = $tau2\; \n";
    printf OUTPUT "f(x) = a*exp(-x/b) + (1-a)*exp(-x/c) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via a, b, c \n";
    printf OUTPUT "set output \"fitMF0ExpTot.ps\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = a*exp(-x/b) + (1-a)*exp(-x/c)\", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);

    system "rm -rf fit.log";
    system "gnuplot fitMF0ExpTot.gnu >& fitMF0ExpTot.log";
    print " Please check fitMF0ExpTot.log for gnuplot fitting.\n";
    system "cat fitMF0ExpTot.log | grep -A 5 \"Final set of parameters\" | tail -n 3 > fitMF0ExpTot.dat"; 

}

#################################################################
sub Fit1Exp {

     # C(x) = exp(-x/t1)
   
    my $t1=$_[0];   # initial t1    
    my $tstart=$_[1]; # start time for fitting
    my $tend=$_[2];   # end time for fitting
    my $inputfile = $_[3]; # file to  fit 

    open(OUTPUT, ">fit1ExpTot.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set style data points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_T(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "t1 = $t1 \; \n";
    printf OUTPUT "f(x) = exp(-x/t1) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via t1 \n";
    printf OUTPUT "set output \"fit1ExpTot.ps\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = exp(-x/t1) \", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);
 
    system "rm -rf fit.log";
    system "gnuplot fit1ExpTot.gnu >& fit1ExpTot.log";
    print " Please check fit1ExpTot.log for gnuplot fitting.\n";
    system "cat fit1ExpTot.log | grep -A 3 \"Final set of parameters\" | tail -n 1 > fit1ExpTot.dat"; 

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

    open(OUTPUT, ">fit2ExpTot.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set style data points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_T(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "a1 = $a1 \; t1 = $t1 \; t2 = $t2\; \n";
    printf OUTPUT "f(x) = a1*exp(-x/t1) + (1-a1)*exp(-x/t2) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via a1, t1, t2 \n";
    printf OUTPUT "set output \"fit2ExpTot.ps\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = a1*exp(-x/t1) + (1-a1)*exp(-x/t2) \", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);

    system "rm -rf fit.log";
    system "gnuplot fit2ExpTot.gnu >& fit2ExpTot.log";
    print " Please check fit2ExpTot.log for gnuplot fitting.\n";
    system "cat fit2ExpTot.log | grep -A 5 \"Final set of parameters\" | tail -n 3 > fit2ExpTot.dat"; 

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

    open(OUTPUT, ">fit3ExpTot.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set style data points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_T(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "a1 = $a1 \; t1 = $t1 \; a2 = $a2 \; t2 = $t2 \; t3 = $t3 \; \n";
    printf OUTPUT "f(x) = a1*exp(-x/t1) + a2*exp(-x/t2) + (1-a1-a2)*exp(-x/t3) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via a1, t1, a2, t2, t3 \n";
    printf OUTPUT "set output \"fit3ExpTot.ps\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = a1*exp(-x/t1) + a2*exp(-x/t2) + (1-a1-a2)*exp(-x/t3)\", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);

    system "rm -rf fit.log";
    system "gnuplot fit3ExpTot.gnu >& fit3ExpTot.log";
    print " Please check fit3ExpTot.log for gnuplot fitting.\n";
    system "cat fit3ExpTot.log | grep -A 7 \"Final set of parameters\" | tail -n 5 > fit3ExpTot.dat"; 

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

    open(OUTPUT, ">fit4ExpTot.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set style data points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_T(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "a1 = $a1 \; t1 = $t1 \; a2 = $a2 \; t2 = $t2 \; a3 = $a3 \; t3 = $t3 \; t4 = $t4 \; \n";
    printf OUTPUT "f(x) = a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + (1-a1-a2-a3)*exp(-x/t4) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via a1, t1, a2, t2, a3, t3, t4  \n";
    printf OUTPUT "set output \"fit4ExpTot.ps\" \n";
    #printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + (1-a1-a2-a3)*exp(-x/t4) \", \"$inputfile\" u 1:2 t \"\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \" fitting to four exponential terms \", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);

    system "rm -rf fit.log";
    system "gnuplot fit4ExpTot.gnu >& fit4ExpTot.log";
    print  " Please check fit4ExpTot.log for gnuplot fitting.\n";
    system "cat fit4ExpTot.log | grep -A 9 \"Final set of parameters\" | tail -n 7 > fit4ExpTot.dat"; 

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

    open(OUTPUT, ">fit5ExpTot.gnu");
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set pointsize 0.7 \n";
    printf OUTPUT "set style data points \n";
    printf OUTPUT "set xlab \"{/= 24 Time (ps)}\" \n";
    printf OUTPUT "set ylab \"{/= 24 C_T(t)}\" \n";
    printf OUTPUT "set title \"\" \n";
    printf OUTPUT "set key top right \n";
    printf OUTPUT "set xrange\[$tstart:$tend\] \n"; 
    printf OUTPUT "a1 = $a1 \; t1 = $t1 \; a2 = $a2 \; t2 = $t2 \; a3 = $a3 \; t3 = $t3 \; a4 = $a4 \; t4 = $t4 \; t5 = $t5 \; \n";
    printf OUTPUT "f(x) = a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + a4*exp(-x/t4) + (1-a1-a2-a3-a4)*exp(-x/t5) \n";
    printf OUTPUT "fit \[$tstart : $tend\] f(x) \"$inputfile\" u 1:2 via a1, t1, a2, t2, a3, t3, a4, t4, t5  \n";
    printf OUTPUT "set output \"fit5ExpTot.ps\" \n";
    #printf OUTPUT "plot \[$tstart : $tend\] f(x) t \"f(x) = a1*exp(-x/t1) + a2*exp(-x/t2) + a3*exp(-x/t3) + a4*exp(-x/t4) + (1-a1-a2-a3-a4)*exp(-x/t5)\", \"$inputfile\" u 1:2 t \"\" \n";
    printf OUTPUT "plot \[$tstart : $tend\] f(x) t \" fitting to five exponential forms \", \"$inputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);

    system "rm -rf fit.log";
    system "gnuplot fit5ExpTot.gnu >& fit5ExpTot.log";
    print  " Please check fit5ExpTot.log for gnuplot fitting.\n";
    system "cat fit5ExpTot.log | grep -A 11 \"Final set of parameters\" | tail -n 9 > fit5ExpTot.dat"; 

}

#################################################################
sub JwMF0Exp {
    
    # J(x) = a1*2*t1/(1 + w^2*t1^2)  + (1-a1)*2*t2/(1 + w^2*t2^2)
    my $coeff1=$_[0]; # final a1 
    my $tau1=$_[1];   # final t1    
    my $tau2=$_[2];   # final t2 
    my $omega=$_[3];  # frequency
    my $part1 = $coeff1* 2.0 * $tau1/(1.0 + $omega*$omega*$tau1*$tau1);
    my $part2 = (1-$coeff1) * 2.0 * $tau2/(1.0 + $omega*$omega*$tau2*$tau2);
    my $finalvalue =  $part1 + $part2;
    return $finalvalue;

}
#################################################################
sub Jw5Exp {
    
    # J(x) = a1*2*t1/(1 + w^2*t1^2)  +  a2*2*t2/(1 + w^2*t2^2) +  a3*2*t3/(1 + w^2*t3^2)  + a4*2*t4/(1 + w^2*t4^2) +  (1-a1-a2-a3-a4)*2*t5/(1 + w^2*t5^2)
    my $coeff1=$_[0]; # final a1
    my $tau1=$_[1];   # final t1    
    my $coeff2=$_[2]; # final a2
    my $tau2=$_[3];   # final t2
    my $coeff3=$_[4]; # final a3
    my $tau3=$_[5];   # final t3
    my $coeff4=$_[6]; # final a4
    my $tau4=$_[7];   # final a4
    my $tau5=$_[8];   # final a5
    my $omega=$_[9];  # frequency
    my $part1 = $coeff1* 2.0 * $tau1/(1.0 + $omega*$omega*$tau1*$tau1);
    my $part2 = $coeff2* 2.0 * $tau2/(1.0 + $omega*$omega*$tau2*$tau2);
    my $part3 = $coeff3* 2.0 * $tau3/(1.0 + $omega*$omega*$tau3*$tau3);
    my $part4 = $coeff4* 2.0 * $tau4/(1.0 + $omega*$omega*$tau4*$tau4);
    my $part5 = (1-$coeff1-$coeff2-$coeff3-$coeff4)* 2.0 * $tau5/(1.0 + $omega*$omega*$tau5*$tau5);
    my $finalvalue =  $part1 + $part2 + $part3 + $part4 + $part5;
    return $finalvalue;

}
#################################################################
sub Jw4Exp {
    
    # J(x) = a1*2*t1/(1 + w^2*t1^2)  +  a2*2*t2/(1 + w^2*t2^2) +  a3*2*t3/(1 + w^2*t3^2)  + (1-a1-a2-a3)*2*t4/(1 + w^2*t4^2) 
    my $coeff1=$_[0]; # final a1 
    my $tau1=$_[1];   # final t1    
    my $coeff2=$_[2]; # final a2
    my $tau2=$_[3];   # final t2
    my $coeff3=$_[4]; # final a3
    my $tau3=$_[5];   # final t3
    my $tau4=$_[6];   # final t4 
    my $omega=$_[7];  # frequency
    my $part1 = $coeff1* 2.0 * $tau1/(1.0 + $omega*$omega*$tau1*$tau1);
    my $part2 = $coeff2* 2.0 * $tau2/(1.0 + $omega*$omega*$tau2*$tau2);
    my $part3 = $coeff3* 2.0 * $tau3/(1.0 + $omega*$omega*$tau3*$tau3);
    my $part4 = (1-$coeff1-$coeff2-$coeff3)* 2.0 * $tau4/(1.0 + $omega*$omega*$tau4*$tau4);
    my $finalvalue =  $part1 + $part2 + $part3 + $part4;
    return $finalvalue;

}
#################################################################
sub Jw3Exp {
    
    # J(x) = a1*2*t1/(1 + w^2*t1^2)  +  a2*2*t2/(1 + w^2*t2^2) + (1-a1-a2)*2*t3/(1 + w^2*t3^2) 
    my $coeff1=$_[0]; # final a1 
    my $tau1=$_[1];   # final t1   
    my $coeff2=$_[2]; # final a2
    my $tau2=$_[3];   # final t2
    my $tau3=$_[4];   # final t3
    my $omega=$_[5];  # frequency
    my $part1 = $coeff1* 2.0 * $tau1/(1.0 + $omega*$omega*$tau1*$tau1);
    my $part2 = $coeff2* 2.0 * $tau2/(1.0 + $omega*$omega*$tau2*$tau2);
    my $part3 = (1-$coeff1-$coeff2) * 2.0 * $tau3/(1.0 + $omega*$omega*$tau3*$tau3);
    my $finalvalue =  $part1 + $part2 + $part3;
    return $finalvalue;

}
#################################################################
sub Jw2Exp {
    
    # J(x) = a1*2*t1/(1 + w^2*t1^2)  + (1-a1)*2*t2/(1 + w^2*t2^2)
    my $coeff1=$_[0]; # final a1 
    my $tau1=$_[1];   # final t1   
    my $tau2=$_[2];   # final t2
    my $omega=$_[3];  # frequency
    my $part1 = $coeff1* 2.0 * $tau1/(1.0 + $omega*$omega*$tau1*$tau1);
    my $part2 = (1-$coeff1) * 2.0 * $tau2/(1.0 + $omega*$omega*$tau2*$tau2);
    my $finalvalue =  $part1 + $part2;
    return $finalvalue;

}
#################################################################
sub Jw1Exp {
    
    # J(x) = 2*t1/(1 + w^2*t1^2)
    my $tau1=$_[0]; # final t1      
    my $omega=$_[1];  # frequency
    #print "$_[0] $_[1] \n";
    my $finalvalue =  2.0 * $tau1/(1.0 + $omega*$omega*$tau1*$tau1);
    return $finalvalue;

}
