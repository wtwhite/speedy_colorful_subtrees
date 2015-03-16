#!/usr/bin/perl

use strict;
#use warnings;

chomp(my $nV = <>);
chomp(my $nE = <>);
chomp(my $nC = <>);

my @colFor;
for (1 .. $nV) {
	chomp($_ = <>);
	my ($v, $c) = split;
	$colFor[$v] = $c;
}

my @destColours;
while (<>) {
	chomp;
	my ($u, $v, $w) = split;
	++$destColours[$u]{$colFor[$v]};
}

my @freq;
for (my $i = 0; $i < $nV; ++$i) {
	foreach my $c (keys %{$destColours[$i]}) {
		++$freq[$destColours[$i]{$c}];
	}
}

print "# edges\tFreq\n";
for (my $i = 0; $i < @freq; ++$i) {
	print "$i\t$freq[$i]\n";
}
