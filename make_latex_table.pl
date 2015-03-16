#!/usr/bin/perl

use strict;
use warnings;

my $cols = shift or die "Must specify columns!";
#print "\\begin{table}{$cols}\n";
print "\\begin{longtable}{$cols}\n";

while (<>) {
	chomp;
	s/\r\z//;
	
	my @f;
	while (s/\A([^",]+|"(.*?)"),?//) {		#HACK: Broken regex, but good enough
		my $x = $1;
		$x =~ s/"//g;		#HACK
		$x =~ s/_/\\_/g;
		push @f, $x;
	}
	
	print "\\hline\n" if $. <= 2;
	print join("\t&\t", @f), "\\\\\n";
}

print "\\hline\n";
#print "\\end{table}\n";
print "\\end{longtable}\n";
