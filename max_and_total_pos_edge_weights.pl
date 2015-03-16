#!/usr/bin/perl

use strict;
use warnings;

my $nV = <>;
my $nE = <>;
my $nC = <>;

for (1 .. $nV) {
	<>;
}

my $max;
my $min;
my $totalPos = 0;
my $totalNeg = 0;
my $nPos = 0;
my $nNeg = 0;

while (<>) {
	my ($u, $v, $w) = split;
	$max = $w if !defined $max || $w > $max;
	$min = $w if !defined $min || $w < $min;
	if ($w > 0) {
		$totalPos += $w;
		++$nPos;
	}
	if ($w < 0) {
		$totalNeg += $w;
		++$nNeg;
	}
}

print "max:\t$max\nmin:\t$min\ntotal +ve:\t$totalPos\ntotal -ve:\t$totalNeg\n# +ve:\t$nPos\n# -ve:\t$nNeg\n";
