#! /usr/bin/perl

use strict;
use warnings;

my $a = <>;
my $b = <>;
my $c = <>;
my $d = <>;

my $prod1 = $a * $b;
my $prod2 = $c * $d;

my $div = 2 * ($prod1 + $prod2);

die "Bad input\n" if ($div == 0);
my $sum_sq2 = sq($d) + sq($c);
my $sum_sq1 = sq($a) + sq($b);
my $cos_th = ($sum_sq2 - $sum_sq1) / $div;

die "cos_th = $cos_th, are you sure you measured it right?\n" if abs($cos_th) >= 1.0;
my $sin_th = sqrt(1 - sq($cos_th));
my $side = sqrt($sum_sq1 + 2 * $prod1 * $cos_th);
my $r = $side/(2*$sin_th);
printf("cos(th)=%.3f side=%.3f r=%.3f\n", $cos_th, $side, $r);

sub sq
{
	my $x = shift;
	return $x * $x;
}
