#!/usr/bin/perl

use strict;
use warnings;

#   objectiveValue="-120.434617364992"
#   solutionStatusString="integer optimal solution"
#  <variable name="x9041_9346y" index="2074638" value="0"/>

my $objVal;
my $extraInfo = "";
my @onVarLines;
while (<>) {
	if (/objectiveValue="(.*?)"/) {
		$objVal = -$1;		# Flip sign to agree with the flipped signs on the original problems ilp.jar gives to Gurobi.
	} elsif (/solutionStatusString="(.*?)"/) {
		if ($1 ne "integer optimal solution") {
			$extraInfo = " BUT: CPLEX solutionStatusString=\"$1\"!";
		}
	} elsif (m!<variable name="x_(.*?)" index="\d+" value="(.*?)"/>!) {
		#if ($2 + 0 == 1) {
		if ($2 >= 0.5) {		# CPLEX sometimes returns a value like "0.999999999999999" or "-1.11022302462516e-15"
			push @onVarLines, "x$1y 1\n";
		}
	}
}

printf "# Objective value (to 6 d.p.) = \%.6f$extraInfo\n", $objVal;
print sort @onVarLines;
