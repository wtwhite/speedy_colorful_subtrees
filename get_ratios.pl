#!/usr/bin/perl

use strict;
use warnings;

my @nBetter = (0, 0);

<>; <>;
while (<>) {
	last if /\\hline/;
	
	s/\s*&\s*/\t/g;
	my ($baseline) = (split)[1];
	$baseline =~ s/^\\textbf\{(.*)\}$/$1/;# or die "<$baseline>";
	$_ = <>;
	s/\s*&\s*/\t/g;
	my ($suggested) = (split)[5];
	$suggested =~ s/^\\textbf\{(.*)\}$/$1/;# or die "<$suggested>";
	my ($global) = (split)[1];
	$global =~ s/^\\textbf\{(.*)\}$/$1/;# or die "<$suggested>";
	
	print $suggested / $baseline, "\n";
	++$nBetter[$global < $baseline];
}

print "Global better than baseline $nBetter[1] times, other way round $nBetter[0] times.\n";
