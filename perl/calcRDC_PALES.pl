#!/usr/bin/perl -w
require 5.002;

use utf8;
use Getopt::Std;
use warnings;
use strict;

my $myhead = <<end_of_myhead;
############################################################################################################
# This is a perl program to calculate Residual Dipolar Coupling values of each snapshot from the AMBER pdb #
# trajectory by calling the PALES.                                                                         #
#                                                                                                          #
# 1) prepare the input file, input.inp, for the PALES RDC calculation.                                     #
# 2) get the trajectory in pdb format from AMBER ptraj.                                                    #
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
usage:"Please input as [command] [number of RDCs] [PALES Parameters] [Scaling or not] [number of atoms] [number of conformer] [input amberpdb file]"
end_of_usage

die $usage unless getopts("");

my $ax = 0.001; 
my $bx = 0.4; 
my $cx = 10.0;
my $tol = 1e-6;
my $xmin; 
my $fmin;
my @exp_values = (); 
my @cal_values = (); 

if($#ARGV==5)
{
    my $nrdc=$ARGV[0];
    my $palesParm = $ARGV[1];
    my $scaleCtrl = $ARGV[2]; 
    my $natom=$ARGV[3];
    my $nconform=$ARGV[4];
    my $inputpdb = $ARGV[5];

    my $amberpdb = "tempamber.pdb";
    my $RDCinpPALES = "input.inp";
    my $RDCoutPALES = "output.tbl";

    my (@OrdTensors,@orderparm);
    my (@RDCs,@rdc);
    my  @rmsd;
    my  @rfactor;
    my  @exp_rdc;
    my  @scaledF;
 

    # input the experimental RDCs  
    my $nhead = $nrdc + 2;
    system "cat input.inp | head -n $nhead | tail -n $nrdc > RDCs_exp.out";
    open(INPUT, "RDCs_exp.out");
    @exp_rdc = ();

    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;
	if (defined $tmp_array[6] )
	{
	    push (@exp_rdc,$tmp_array[6]) ; 
	    push (@exp_values,$tmp_array[6]) ; 
	    push (@cal_values,$tmp_array[6]) ; 
	}
    }
    close(INPUT);
  
    my $iconform = 0; 
    open(INPUTPDB,"<$inputpdb");	
    while(<INPUTPDB>)
    {
	chomp;
	my $line = $_;
	if ($line =~ /^MODEL/)  # begin of new model
	{
	    open(OUTPUTPDB, ">$amberpdb");
	}
	print OUTPUTPDB "$line \n";
	if ($line =~ /^ENDMDL/)  # end of new model
	{
	    close(OUTPUTPDB);

            # call pales and calculate RDCs using the $amberpdb file.
	    system "pales-linux $palesParm -inD $RDCinpPALES -pdb $amberpdb -outD $RDCoutPALES >& pales.log";
	
	    # input the calcuated RDCs from Pales  
	    my $ntail = $nrdc + 0;
	    system "cat output.tbl | tail -n $ntail > RDCs.out";
	    open(INPUT, "RDCs.out");
	    @rdc = ();
	    while(<INPUT>) #read in data
	    {
		chomp;
		my @tmp_array=split;
		if (defined $tmp_array[8] )
		{
		    push (@rdc,$tmp_array[8]) ; 
		}
	    }
	    close(INPUT);
	    push @RDCs, [@rdc];
            my $irdc = 0;
	    for ($irdc = 0; $irdc <= $#rdc;  $irdc++)
	    {
		$cal_values[ $irdc] = $rdc[ $irdc];
	    }
	
	    if ($scaleCtrl eq "Y") 
	    {
		$ax = 0.1; 
		$bx = 0.4; 
		$cx = 10.0;
		$tol = 1e-6;
		$xmin = 0.0; 
		$fmin = 0.0;
		&Optimize1D();
	    }
	    elsif ($scaleCtrl eq "N")
	    {
		$xmin = 1.0;  
	    }
	    else
	    {
		print "\$scaleCtrl has not right value (Y or N)";
		exit(1);
	    }
	    
	    # calculate the RDC RMSD and R factor 
	    my $rmsd_numerator = 0.0;
	    my $rmsd_denominator = 0.0;
	    
	    for ($irdc = 0; $irdc < $nrdc;  $irdc++)
	    {
		$rmsd_numerator =  $rmsd_numerator + ($xmin*$rdc[$irdc] -  $exp_rdc[$irdc])*($xmin*$rdc[$irdc] - $exp_rdc[$irdc]);
		$rmsd_denominator = $rmsd_denominator +  $exp_rdc[$irdc]*$exp_rdc[$irdc];
		
	    }
	    push (@scaledF,$xmin);
	    push (@rmsd,sqrt($rmsd_numerator/$nrdc)); 
	    push (@rfactor,sqrt($rmsd_numerator/$rmsd_denominator)) ;
	
	
	    # input the alignment/order tensor 
	    system "cat output.tbl | grep \"\(Sxx_d,Syy_d,Szz_d\)\" > OrderTensor.out";
	    open(INPUT, "OrderTensor.out");
	    @orderparm = ();
	    while(<INPUT>) #read in data
	    {
		chomp;
		my @tmp_array=split;
		if (defined $tmp_array[3] )
		{
		    push (@orderparm,$xmin*$tmp_array[3]) ;
		    push (@orderparm,$xmin*$tmp_array[4]) ;
		    push (@orderparm,$xmin*$tmp_array[5]) ;
		    
		    
		    my $GDO = $xmin*sqrt(2.0*($tmp_array[3]*$tmp_array[3] + $tmp_array[4]*$tmp_array[4] + $tmp_array[5]*$tmp_array[5])/3.0);
		    my $eta = ($tmp_array[3]- $tmp_array[4])/$tmp_array[5];
		    push (@orderparm,$GDO);
		    push (@orderparm,$eta);
		}
	    }
	    close(INPUT);
	    push @OrdTensors,[@orderparm];

	    $iconform = $iconform + 1; 
	}
    }
    close(INPUTPDB);

    # output all alignment tensors  
    open(OUTPUT, ">OrderTensorsAll.out");
    printf OUTPUT " %10s %10s %12s %12s %12s %12s %12s %10s %10s \n", "Conform ID", "ScaledF", "Sxx","Syy", "Szz", "GDO", "eta", "rmsd", "Rfactor"; 
    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
       { 
	   printf OUTPUT " %10d %10.6f ", $iconform, $scaledF[$iconform];

	   for (my $iorder = 0; $iorder < 5;  $iorder++)
	   {
	       printf OUTPUT "%12.4e ",$OrdTensors[$iconform][$iorder];
	 
	   }
	   printf OUTPUT "%10.5f %10.5f \n",$rmsd[$iconform],$rfactor[$iconform];

       }
    close(OUTPUT);

    # output all RDCs direct from PALES
    open(OUTPUT, ">RDCsAll.out");
    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
       { 
	   printf OUTPUT " %10d %10.5f ", $iconform, $rfactor[$iconform];

	   for (my $irdc = 0; $irdc < $nrdc;  $irdc++)
	   {
	       printf OUTPUT "%9.4f ",$RDCs[$iconform][$irdc];
	   }
	   printf OUTPUT "\n";
       }
    close(OUTPUT);

   # output all RDCs scaled by a factor after PALES to best fit the expermental values
    open(OUTPUT, ">RDCsScaledAll.out");
    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
       { 
	   printf OUTPUT " %10d %10.5f ", $iconform, $rfactor[$iconform];

	   for (my $irdc = 0; $irdc < $nrdc;  $irdc++)
	   {
	       printf OUTPUT "%9.4f ",$scaledF[$iconform]*$RDCs[$iconform][$irdc];
	   }
	   printf OUTPUT "\n";
       }
    close(OUTPUT);

    # put and Q factor in a hash table 
    my %hashQ = ();
    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
    {
	$hashQ{$iconform} = $rfactor[$iconform];	
    }

    # Get a list of hash keys sorted by value.
    my @sortedID = sort { $hashQ{$a} <=> $hashQ{$b} } keys %hashQ; 

    
    # output all sorted alignment tensors  
    open(OUTPUT, ">OrderTensorsSorted.out");
    printf OUTPUT " %10s %10s %10s %12s %12s %12s %12s %12s %10s %10s \n", "Rank ID", "Conform ID", "ScaledF",
                  "Sxx","Syy", "Szz", "GDO", "eta", "rmsd", "Rfactor";
    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
    {    
	printf OUTPUT " %10d %10d %10.6f ",$iconform, $sortedID[$iconform],$scaledF[$sortedID[$iconform]];
	
	for (my $iorder = 0; $iorder < 5;  $iorder++)
	{
	    printf OUTPUT "%12.4e ",$OrdTensors[$sortedID[$iconform]][$iorder];
	}
	printf OUTPUT "%10.5f %10.5f \n",$rmsd[$sortedID[$iconform]],$rfactor[$sortedID[$iconform]];
    }
    close(OUTPUT);

    # output all sorted RDCs 
    open(OUTPUT, ">RDCsSorted.out");
    for (my $iconform = 0; $iconform < $nconform; $iconform++ )
       { 
	   printf OUTPUT " %10d %10d %10.5f ",$iconform, $sortedID[$iconform],$rfactor[$sortedID[$iconform]];

	   for (my $irdc = 0; $irdc < $nrdc;  $irdc++)
	   {
	       printf OUTPUT "%9.4f ",$scaledF[$sortedID[$iconform]]*$RDCs[$sortedID[$iconform]][$irdc];
	   }
	   printf OUTPUT "\n";
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

    for (my $i= 0; $i<=$#exp_values; $i++ )
        {
	    # $value = ($xin*$cal_values[$i]-$exp_values[$i])/$exp_values[$i]; 
	    $value = ($xin*$cal_values[$i]-$exp_values[$i]); 
            $output = $output + $value*$value;
        }
    return $output;
}

