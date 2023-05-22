#!/usr/bin/perl -w
require 5.002;
use Getopt::Std;
use strict;
#
#

my $head = <<end_of_head;
----------------------------------------------------------------------------------------------
|   This is a perl program to calculate the free energy from the histogram or probability    | 
|   distribution in one dimension or two dimensions.                                         | 
|                                                                                            |
|      Junchao Xia, BioMaPS Institute, Rutgers University, junchao-xia\@biomaps.rutgers.edu.  |
----------------------------------------------------------------------------------------------
end_of_head

print $head;

my $usage = <<end_of_usage;

usage:
  1D: getFreeEng.pl nx KT inputfile outputfile
  2D: getFreeEng.pl nx ny KT inputfile outputfile 
end_of_usage
    
    
die $usage unless getopts("");


if($#ARGV==4)
{
    &FreeEng2d($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3],$ARGV[4]);
}
elsif($#ARGV==3)
{
    &FreeEng1d($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3]);
}
else
{ print "Input options are wrong. Please check the usage below. \n";
  die $usage;
}
##############################################################
sub FreeEng1d {

    my $nx=$_[0]; #number of divisions  
    my $KT = $_[1];
    my $inputfile = $_[2];
    my $outputfile = $_[3];

    #  my $tempperature = $_[1]
    #  my $Kb=*8.31451121E-3;      # Kb in KJ/mol
    #  my $Kb=8.31451121E-3/4.184;    # Kb in KCal/mol  
    #  my $KT = $Kb*$temperature;
  
    my $FreeEng = 0.0;

    my @array1;
    my @array2;
    open(INPUT, $inputfile);
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;

	push (@array1,$tmp_array[0]) ;
	push (@array2,$tmp_array[1]) ;
    }
    close INPUT;

    my $min = 1000;
    for my $counter (0 .. $nx-1)
    {
         if ($array2[$counter] > 0)
            {
               if ($array2[$counter] < $min)  
               {$min = $array2[$counter];}
            }
    }

    open(OUTPUT, ">$outputfile");
    for my $counter (0 .. $nx-1)
    {
	if ($array2[$counter] < $min)
	{
	    # $FreeEng = - $temperature*$Kb*log($min);
	    $FreeEng = - $KT*log($min);
	}
	else
	{
	    # $FreeEng = - $temperture*$Kb*log($array2[$counter]);
	    $FreeEng = - $KT*log($array2[$counter]);
	    
	}
	print OUTPUT " $array1[$counter]  $FreeEng \n" ;
    }
    close OUTPUT;
    print " The free energy in one dimension has been written to $outputfile.\n";

    my @file_type = split(/\./,$outputfile); #splitting file to get the base and extension
    my $plotname = $file_type[0];

    open(OUTPUT, ">template.gnu");
    printf OUTPUT "set pointsize 1.0 \n";
    printf OUTPUT "set data style linespoints \n";
    printf OUTPUT "set xlab \"{/= 24 x }\" \n";
    printf OUTPUT "set ylab \"{/= 24 F(x)}\" \n";
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set output \"$plotname.ps\" \n";
    printf OUTPUT "plot \"$outputfile\" u 1:2 t \"\" \n";
    close(OUTPUT);
    system "gnuplot template.gnu";
    print " The free energy in one dimension has been plotted into the file, $plotname.ps.\n";

}
#########################################################################
sub FreeEng2d {
    

    my $nx=$_[0]; #number of divisions
    my $ny=$_[1]; #number of divisions
    my $KT = $_[2];
    my $inputfile = $_[3];
    my $outputfile = $_[4];
    
    #  my $temperature = $_[2];
    #  my $Kb=*8.31451121E-3;        # Kb in KJ/mol
    #  my $Kb = 8.31451121E-3/4.184;    # Kb in KCal/mol
    #  my $KT = $temperature*$Kb; 

    my $FreeEng = 0.0;

    my @array1;
    my @array2;
    my @array3;
    
    my $counter=0;
   
    open(INPUT, $inputfile);
    while(<INPUT>) #read in data
    {
	chomp;
	my @tmp_array=split;

	if (defined $tmp_array[2] )
	{
	    push (@array1,$tmp_array[0]) ;
	    push (@array2,$tmp_array[1]) ;
	    push (@array3,$tmp_array[2]) ;
	}
    }	
    close INPUT;

    my $min = 1000;
    for $counter (0 .. $#array3)
    {
       if ($array3[$counter] > 0)
       {
          if ($array3[$counter] < $min)
          {$min = $array3[$counter];}

       }   

    }

    $counter=0;

    open(OUTPUT, ">$outputfile");
	
    for my $I (1 .. $nx) #initializing the array
    {	
	for my $J (1 .. $ny)
	{   
            if ($array3[$counter] < $min)
	    {
                # $FreeEng = - $temperature*$Kb*log($min);
		$FreeEng = - $KT*log($min);
	    }
            else 
	    {
		#$FreeEng = - $temperature*$Kb*log($array3[$counter]);
		$FreeEng = - $KT*log($array3[$counter]);
	    }	    
	    printf OUTPUT "$array1[$counter]   $array2[$counter]   $FreeEng \n";
	     ++$counter;
	    
	}    
	print OUTPUT "  \n";
    }
    close OUTPUT;
    print " The free energy in two dimensions has been written to $outputfile.\n";

    my @file_type = split(/\./,$outputfile); #splitting file to get the base and extension
    my $plotname = $file_type[0];

    open(OUTPUT, ">template.gnu");
    printf OUTPUT "set dgrid3d $nx, $ny, 2 \n";
    printf OUTPUT "set view map \n";
    printf OUTPUT "set contour base \n";
    printf OUTPUT "set cntrparam levels auto \n";
    printf OUTPUT "set pm3d implicit at b \n";
    printf OUTPUT "set style data lines \n";
    printf OUTPUT "unset key \n";
    printf OUTPUT "unset surface \n";
    printf OUTPUT "set xlab \"{/= 24 {/Symbol F}}\" \n";
    printf OUTPUT "set ylab \"{/= 24 {/Symbol Y}}\" \n";
    printf OUTPUT "set term post enhanced color \"Times-Roman\" 20 \n";
    printf OUTPUT "set output \"$plotname.ps\" \n";
    printf OUTPUT "splot \"$outputfile\" \n";
    close(OUTPUT);
    system "gnuplot template.gnu";
    print " The free energy in two dimensions has been plotted into the file, $plotname.ps.\n";

}
