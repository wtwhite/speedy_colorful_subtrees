#!/usr/bin/perl

use strict;
use warnings;

<>;
<>;
<>;
<>;		#HACK: Don't even care if we miss the first colour or not

my @weights;
while (<>) {
	if (/\A(\d+)\s+(\d+)\s+([-+\d.eE]+)\s*\z/) {
		push @weights, $3;
	}
}

@weights = sort { $a <=> $b } @weights;
#print STDERR map { "$_\n" } @weights;		#DEBUG
my $prev = shift @weights;
my $min = 1e9;		# "Infinity"
foreach (@weights) {
	$min = $_ - $prev if $_ - $prev < $min;
	$prev = $_;
}

print "Minimum difference between any pair of edge weights: $min\n";
