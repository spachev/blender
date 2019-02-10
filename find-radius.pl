#! /usr/bin/perl

use strict;
use warnings;

my @points = ();

while (my $line = <>)
{
	my $p = parse_line($line);
	next unless $p;
	push(@points, $p);
	last if (scalar(@points) == 3);
}

die "Bad input\n" unless (scalar(@points) == 3);

my $d = 0;
my $sx = 0;
my $sy = 0;
for (my $i = 0; $i < 3; $i++)
{
	my $other_1 = ($i + 4) % 3;
	my $other_2 = ($i + 2) % 3;
	my $dy_other = $points[$other_1]->{y} - $points[$other_2]->{y};
	my $dx_other = $points[$other_1]->{x} - $points[$other_2]->{x};
	$d += $points[$i]->{x} * $dy_other;
	my $sum_sq = sq($points[$i]->{x}) + sq($points[$i]->{y});
	$sx += $sum_sq * $dy_other;
	$sy += $sum_sq * -$dx_other;
}

$d *= 2;
die "Points are on the same line\n" if ($d == 0.0);
my $x = $sx / $d;
my $y = $sy / $d;
my $center = {x => $x, y => $y};

print "Center ($x, $y)\n";

for (my $i = 0; $i < 3; $i++)
{
	my $r = dist($center, $points[$i]);
	printf("r($i) = %.3f\n", $r);
}

sub dist
{
	my $p1 = shift;
	my $p2 = shift;
	return sqrt(sq($p1->{x} - $p2->{x}) + sq($p1->{y} - $p2->{y}));
}

sub sq
{
	my $x = shift;
	return $x * $x;
}

sub parse_line
{
	my $line = shift;
	my $res = {};
	if ($line =~ /G1\s+X([\d\.]+)\s+Y([\d\.]+)/)
	{
		$res->{x} = $1;
		$res->{y} = $2;
		return $res;
	}
	return 0;
}
