#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;

die "need 2 arguments: x0 and x1\n"
  unless defined(@ARGV) and scalar @ARGV == 2;

my ($x0, $x1) = ($ARGV[0], $ARGV[1]);

my $nb = 32;
my $delta = 2*&pi()/$nb;

my $origin = -&pi();

my $delta_str = sprintf "%e", $delta;
my $origin_str = sprintf "%e", $origin;

my $x0_str = sprintf "%e", $x0;
my $x1_str = sprintf "%e", $x1;

print <<EOM;
netcdf umbrella {
dimensions:
	row-major = $nb ;
variables:
	double coeffs(row-major) ;

// global attributes:
		:nextents = 1 ;
		:extents = $nb ;
		:periodicity = 1 ;
		:origin = $origin_str ;
		:spacing = $delta_str ;
data:

 coeffs = 
EOM

for (my $j = 0; $j < $nb; ++$j) {
   my $x = $origin + $j*$delta;
   my $E = 100.0;
   $E = 0.0 if $x > $x0 and $x < $x1;
   printf "%e", $E;
   print ',' unless $j == $nb - 1;
   print "\n";
}
print "; }\n";
