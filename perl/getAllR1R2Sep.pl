#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;
use Switch;

my $myhead = <<end_of_myhead;
############################################################################################################
# This a perl program to calculate relaxation parameters R1/R2 and NOE from the fitting parameters of      #
# correlation functions of internal and overall motions.                                                   #
#                                                                                                          #
#                                                                                                          #
# Junchao Xia, junchao-xia\@biomaps.rutgers.edu, BioMaPs institute, Rutgers University.                     #
############################################################################################################;
end_of_myhead
print $myhead;
print "\n";

#

my $usage = <<end_of_usage;
usage:"should input like this [command] [Bfield] [spin type C or N] [nvec] [input file from internal motion] [nexpi] [input file from overall motion] [nexpo]"
end_of_usage

die $usage unless getopts("");

if ($#ARGV != 6) {die $usage;}

my $Bfield  = $ARGV[0];          # B field in T (11.7413)
my $SpinType = $ARGV[1];         # Spin type N or C 
my $nvec = $ARGV[2];             # total number of vectors
my $inputInt = $ARGV[3];         # fitting parameters from interal motion
my $nexpInt = $ARGV[4];          # number of exponentials for fitting interal motion
my $inputOvr = $ARGV[5];         # fitting parameters from overall motion
my $nexpOvr = $ARGV[6];          # number of exponentials for fitting overall motion 

my $GAMMA_i;
if ($SpinType eq "N")
{$GAMMA_i = -27.1e6;}             # rad/sT, default value for N15 (found typo in Mar 15, 2012, no negtive before)
elsif  ($SpinType eq "C")
{$GAMMA_i = 67.2274e6;}          # rad/sT, default value for C13
else
{
    print "Spin type is not right.";
    exit(0);
}
	

my $GAMMA_j = 2.67522128e8;     # rad/sT, default value for H1
my $H_BAR   = 1.05457266e-34;   # Js
my $MU_0    = 12.56e-7;         # Vs/(Am)

my $Reff;
if ($SpinType eq "N")
{$Reff    = 1.01e-10;}        # m,  default value for N-H vector 
elsif  ($SpinType eq "C")
{$Reff    = 1.09e-10;}        # m,  default value for C-H vector
else
{
    print "Spin type is not right.";
    exit(0);
}


my $PI      = 3.14159;         
    
my $K  =  $MU_0*$H_BAR*$GAMMA_i*$GAMMA_j/(4*$PI*$Reff*$Reff*$Reff);  # note the power of 3 is missing in the paper  

my $OMEGA_i;
if ($SpinType eq "N")
{$OMEGA_i = 2.0*$PI*4.315e6*$Bfield;}  # Hz, default for Omega_N  (found typo in Mar 15, 2012, set to negtive before)
elsif  ($SpinType eq "C")
{$OMEGA_i = 2.0*$PI*10.705e6*$Bfield;}  # Hz, default for Omega_C
else
{
    print "Spin type is not right.";
    exit(0);
}


my $OMEGA_j = 2.0*$PI*42.577e6*$Bfield;    # Hz, default for Omega_H

my $DELTA_i;              #  Chemical-shift anisotropy  
if ($SpinType eq "N")
{$DELTA_i  = 170/1.0e6;}  #  Chemical-shift anisotropy of N    (found typo in Mar 19, 2012, 1.0e6 was 10e6) 
elsif  ($SpinType eq "C")
{$DELTA_i  = 25/1.0e6;}  #  Chemical-shift anisotropy of C  Need check the exact value if it is important
else
{
    print "Spin type is not right.";
    exit(0);
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


my (@parmTableI,@parmi);
my (@parmTableO,@parmo);

open(INPUT, "$inputInt");
print "The input fitting parameters from internal motions are:\n";
while(<INPUT>) #read in data
{
    chomp;
    @parmi=split;
    if ((defined $parmi[0]) && ($parmi[0] ne "#"))
    {
 	print "@parmi \n"; 
	push @parmTableI, [@parmi];
    }
}	
close(INPUT);

open(INPUT, "$inputOvr");
print "The input fitting parameters from overall motions are:\n";
while(<INPUT>) #read in data 
{
    chomp;
    @parmo=split;
    if ((defined $parmo[0]) && ($parmo[0] ne "#"))
    {
	print "@parmo \n"; 
	push @parmTableO, [@parmo];
    }
}	
close(INPUT);
 
my @AI;
my @TI;
my @AO;
my @TO;
my @AT;
my @TT;


my $npInts = 2*$nexpInt + 1; 
my $npOvrs = 2*$nexpOvr;
 
open(OUTPUT, ">allfittedR1R2.dat");
for (my $ivec = 0; $ivec < $nvec; $ivec++)
{ 
    print " Working on the spin vector $ivec.\n";
    @AI=();
    @TI=();
    push (@AI,$parmTableI[$ivec][0]);
    push (@TI,0.0);

    for (my $inum = 1; $inum < $npInts-1; $inum = $inum + 2)
    {
	push (@AI,$parmTableI[$ivec][$inum]);
	push (@TI,$parmTableI[$ivec][$inum+1]);
    }

    @AO=();
    @TO=();
    for (my $jnum = 0; $jnum < $npOvrs-1; $jnum = $jnum + 2)
    {
	push (@AO,$parmTableO[$ivec][$jnum]);
	push (@TO,$parmTableO[$ivec][$jnum+1]);
    }

    my $coeff;
    my $tau;

    for (my $i = 0; $i<= $#AI; $i++)
    {
	if ($i == 0)
	{
 	   for (my $j = 0; $j<= $#AO; $j++)
	   {
	      $coeff = $AI[$i]*$AO[$j];
	      $tau =$TO[$j];
	    
	      push(@AT,$coeff);
	      push(@TT,$tau);
	   }
        }
	else
	{
	    for (my $j = 0; $j<= $#AO; $j++)
	    {
		$coeff = $AI[$i]*$AO[$j];
		$tau = ($TI[$i]*$TO[$j])/($TI[$i]+$TO[$j]);
		push(@AT,$coeff);
		push(@TT,$tau);
	    }
	}
    }


    # change the time unit to second 
    for (my $i =0; $i<=$#TT; $i++)
    {
	$TT[$i] = $TT[$i]/1.0e12;
    }

    # note in the  Fawzi's paper an additional 5 is divided for R1DD, R2DD, R1CSA, and R2CSA
    $R1DD = $K*$K*(JwExp($OMEGA_j-$OMEGA_i) + 3.0*JwExp($OMEGA_i) + 6.0*JwExp($OMEGA_j+$OMEGA_i))/20.0;

    $R2DD = $K*$K*(4.0*JwExp(0.0) + JwExp($OMEGA_j-$OMEGA_i) + 3.0*JwExp($OMEGA_i) + 6.0*JwExp($OMEGA_j) 
		   +6.0*JwExp($OMEGA_j+$OMEGA_i))/40.0;  # note 6.0 is 3.0 in the paper

    $R1CSA =$DELTA_i*$DELTA_i*$OMEGA_i*$OMEGA_i*JwExp($OMEGA_i)/5.0;
    

    $R2CSA =$DELTA_i*$DELTA_i*$OMEGA_i*$OMEGA_i*(4.0*JwExp(0.0) + 3.0*JwExp($OMEGA_i))/30.0;

    $R2a = 0.0;
    
    $R1 = $R1DD+$R1CSA;
    $R2 = $R2DD+$R2CSA+$R2a;
    $T1 = 1.0/$R1;
    $T2 = 1.0/$R2;

    $NOE = 1 + ($GAMMA_j/$GAMMA_i)*(6.0*JwExp($OMEGA_j+$OMEGA_i) - JwExp($OMEGA_j-$OMEGA_i))/
	(JwExp($OMEGA_j-$OMEGA_i) + (3.0 + 4.0*$DELTA_i*$DELTA_i*$OMEGA_i*$OMEGA_i/($K*$K))*JwExp($OMEGA_i) + 6.0*JwExp($OMEGA_j+$OMEGA_i));

    if ($ivec == 0 ) 
    {
	printf OUTPUT " %10s %10s %10s %10s %10s %10s %10s %10s \n", "#     R1DD", "R1CSA", "R2DD", "R2CSA", "R1", "R2", "R2/R1", "NOE";
        printf OUTPUT " %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n", $R1DD,$R1CSA,$R2DD,$R2CSA,$R1,$R2,$R2/$R1,$NOE;
    }
    else 
    {
    	printf OUTPUT " %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n", $R1DD,$R1CSA,$R2DD,$R2CSA,$R1,$R2,$R2/$R1,$NOE;
    }

}
close(OUTPUT);
print "done! please check the file, allfittedR1R2.dat, for all fitted relaxation parameters.\n"; 

#################################################################
sub JwExp {
    
    my $omega=$_[0]; # the w value
    my $finalvalue = 0.0;
    for (my $i= 0; $i<=$#TT; $i++ )
	{ 
	    $finalvalue = $finalvalue + $AT[$i]* 2.0 * $TT[$i]/(1.0 + $omega*$omega*$TT[$i]*$TT[$i]);
	}
    return $finalvalue;

}
