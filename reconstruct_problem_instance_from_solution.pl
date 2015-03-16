#!/usr/bin/perl

use strict;
use warnings;
use Fatal;

# Syntax:
# reconstruct_problem_instance_from_solution.pl origproblem.txt < solution.grb.sol[.canon] > smallestpossibleproblem.txt
my $origProblemFName = $ARGV[0];
my $origObjFuncScore;

# Read the first line of the solution to get its score
## Objective value = -3.9941259354968331e+01
#x0_1y 0
$_ = <STDIN>;
my ($objFuncScore) = /^# Objective value .*?=\s*(\S+)/ or die "Could not recognise objective function score in line 1, '$_'.";

# Read the rest to learn what edges are in the solution
my %include;
while (<STDIN>) {
	if (/\Ax(\d+)_(\d+)y\s+1\s*\z/) {
		$include{"$1 $2"} = 1;		#HACK
	}
}

# Just copy the "header info" from the original problem file
open my $origProblemF, "<", $origProblemFName;
my $nV = 0 + <$origProblemF>;
my $nE = 0 + <$origProblemF>;
my $nC = 0 + <$origProblemF>;
print "$nV\n";
print scalar(keys %include), "\n";
print "$nC\n";

# A line containing a known solution objective score (which may consist of 2 parts, separated by "+") may optionally follow
$_ = <$origProblemF>;
#if (/\A\s*\S+\s*\z/) {
if (/\A\s*\S+\s*\z/ || /\d.*\+/) {
	$origObjFuncScore = $_;
	$_ = <$origProblemF>;
}

print $objFuncScore, "\n";

print;
for (my $i = 1; $i < $nV; ++$i) {		# We already read 1 line
	print scalar <$origProblemF>;
}

# Finally, read the rest of the original file, and output whichever edges we saw in the solution.
# We need to do it this way because the solution doesn't contain the edge lengths itself.
while (<$origProblemF>) {
	my ($u, $v, $w) = /\A\s*(\d+)\s+(\d+)\s+(\S+)/;
	if (exists $include{"$u $v"}) {
		print;
	}
}

close $origProblemF;
