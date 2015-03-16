#!/usr/bin/perl

use strict;
use warnings;

sub loadGraph($) {
	my ($fn) = @_;
	
	print STDERR "Loading graph from file '$fn'...\n";
	open my $f, "<", $fn or die "Could not open '$fn' for input: $!";
	
	chomp(my $nV = <$f>);
	chomp(my $nE = <$f>);
	chomp(my $nC = <$f>);
	
	# The following line may or may not contain a known optimal objective function value
	my $line = <$f>;
	if ($line =~ /\A\s*\S+\s*\z/) {
		# This line is the objective function value.
		#HACK: We throw this away at the moment.  Hard to do anything better, because the transformation could result in it being invalid!
		print STDERR "Ignoring known objective function value line in file '$fn'...\n";
		$line = <$f>;
	}
	
	my @vColLines;
	push @vColLines, $line;
	for (2 .. $nV) {
		push @vColLines, scalar <$f>;		# Grab newlines and everything (lazy...)
	}
	
	my @edges;
	while (<$f>) {
		chomp;
		my ($u, $v, $w) = split;
		push @edges, [$u, $v, $w];
	}
	
	return ($nV, $nE, $nC, \@vColLines, \@edges);
}

sub saveGraph(@) {
	my ($nV, $nE, $nC, $vColLines, $edges) = @_;
	
	foreach ($nV, $nE, $nC) {
		print "$_\n";
	}
	
	print @$vColLines;
	
	foreach (@$edges) {
		print join("\t", @$_), "\n";
	}
}

sub sortEdges(@) {
	my ($nV, $nE, $nC, $vColLines, $edges) = @_;
	
	return ($nV, $nE, $nC, $vColLines, [sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @$edges]);
}

sub subtractEdges($$) {
	my ($g1, $g2) = @_;
	my @edges;
	
	my $i = 0;
	my $j = 0;
	while ($i < @{$g1->[4]} && $j < @{$g2->[4]}) {
		my $comp = $g1->[4][$i][0] <=> $g2->[4][$j][0] || $g1->[4][$i][1] <=> $g2->[4][$j][1];
		if ($comp == -1) {
			push @edges, [@{$g1->[4][$i++]}];
		} elsif ($comp == 1) {
			++$j;
		} else {
			++$i;
			++$j;
		}
	}
	
	# Arbitrarily take the number of nodes and the colours from graph 1
	return ($g1->[0], scalar @edges, $g1->[2], $g1->[3], \@edges);
}

sub intersectEdges($$) {
	my ($g1, $g2) = @_;
	my @edges;
	
	my $i = 0;
	my $j = 0;
	while ($i < @{$g1->[4]} && $j < @{$g2->[4]}) {
		my $comp = $g1->[4][$i][0] <=> $g2->[4][$j][0] || $g1->[4][$i][1] <=> $g2->[4][$j][1];
		if ($comp == -1) {
			++$i;
		} elsif ($comp == 1) {
			++$j;
		} else {
			push @edges, [@{$g1->[4][$i]}];
			++$i;
			++$j;
		}
	}
	
	# Arbitrarily take the number of nodes and the colours from graph 1
	return ($g1->[0], scalar @edges, $g1->[2], $g1->[3], \@edges);
}

sub unionEdges($$) {
	my ($g1, $g2) = @_;
	my @edges;
	
	my $i = 0;
	my $j = 0;
	while ($i < @{$g1->[4]} && $j < @{$g2->[4]}) {
		my $comp = $g1->[4][$i][0] <=> $g2->[4][$j][0] || $g1->[4][$i][1] <=> $g2->[4][$j][1];
		#print STDERR "i=$i, j=$j: Comparing ($g1->[4][$i][0], $g1->[4][$i][1]) to ($g2->[4][$j][0], $g2->[4][$j][1]): comp = $comp\n";	#DEBUG
		if ($comp == -1) {
			push @edges, [@{$g1->[4][$i]}];
			++$i;
		} elsif ($comp == 1) {
			push @edges, [@{$g2->[4][$j]}];
			++$j;
		} else {
			push @edges, [@{$g1->[4][$i]}];
			++$i;
			++$j;
		}
	}
	
	# There could be leftover edges from one of the two graphs.
	while ($i < @{$g1->[4]}) {
		push @edges, [@{$g1->[4][$i]}];
		++$i;
	}
	while ($j < @{$g2->[4]}) {
		push @edges, [@{$g2->[4][$j]}];
		++$j;
	}
	
	# Arbitrarily take the number of nodes and the colours from graph 1
	return ($g1->[0], scalar @edges, $g1->[2], $g1->[3], \@edges);
}



# Main program
if (@ARGV != 3) {
	die "Must specify 3 arguments: the name of the original graph file, an operator ('-', '&' or '|'), and the name of the second graph file.\n";
}

my %funcForOp = (
	'-' => \&subtractEdges,
	'&' => \&intersectEdges,
	'|' => \&unionEdges
);

my @g1 = sortEdges loadGraph $ARGV[0];
my $op = $funcForOp{$ARGV[1]};
my @g2 = sortEdges loadGraph $ARGV[2];
#my @resultGraph = subtractEdges \@g1, \@g2;
my @resultGraph = &{$op}(\@g1, \@g2);
saveGraph @resultGraph;
