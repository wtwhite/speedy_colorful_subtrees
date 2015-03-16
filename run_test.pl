#!/usr/bin/perl

use strict;
use warnings;
use Fatal;

my $STOPWATCH = "stopwatch -v -s c";
my $FTREDUCE = "ft_reduce";
my $FTREDUCECMD = "$FTREDUCE read write";		# $STOPWATCH gets prepended later.
my $targetDir = "test_results_" . localtime;		# Default to creating a fresh output dir.
$targetDir =~ tr/: /-_/;

#my $S = (? "\\" : '/');		#'HACK: quote just to fix syntax highlighting...
my ($S, $ES);
if (defined $ENV{OS} && $ENV{OS} =~ /^win/i) {
	$S = "\\";
	$ES = "\\\\";		# Ugh.
} else {
	$S = $ES = "/";
}
#print "S=<$S>\n";		#DEBUG

# Fix any slashes.
sub p($) {
	my $x = shift;
#	$x =~ tr|/\\|$S$S|;
	$x =~ s|[/\\]|$S|g;		#HACK: You CANNOT interpolate variables into a tr///!  Jeezus...
	return $x;
}

sub j(@) {
	return join($S, @_);
}

sub bn($) {
	my $x = p shift;
	#$x =~ s|\A.*$S||;
	$x =~ s|\A.*$ES||;
	return $x;
}

sub d($) {
	my $x = p shift;
	#$x =~ s|${S}[^${S}]*||;
	$x =~ s|${ES}[^${ES}]*||;
	$x .= $S if !length $x || ($S eq "\\" && $x =~ /\A[A-Za-z]:\z/);
	return $x;
}

sub safelyQuote(@) {
	if ($S eq '/') {
		return join(' ', map { /[^-A-Za-z0-9_\.=+]/ ? "'$_'" : $_ } @_);			# Unixy: just add single-quotes.
	} else {
		#return join(' ', map { my $x = $_; $x =~ s|([^-A-Za-z0-9_\.=+*?])|^$1|g; $x } @_);		# Windows: carets are safest.
		return join(' ', map { my $x = $_; if (/[^-A-Za-z0-9_\.=+*?() ]/) { $x =~ s|([^-A-Za-z0-9_\.=+*?()])|^$1|g; } else { $x = (/ / ? "\"$_\"" : $_); } $x } @_);		# Windows: carets are safest.
	}
}

# Print to stdout and also to a log file in the target directory.
my $logF;
sub logMsg(@) {
	if (!defined $logF) {
		open $logF, ">>", p "$targetDir/run_test.log";			# Append
	}
	
	foreach (@_) {
		print;
		print $logF $_;
	}
}

sub run($) {
	my ($cmd) = @_;
	logMsg $cmd, "\n";
	system $cmd and die;
}

# Returns them in the same order (otherwise a simpler implementation is possible)
sub uniq(@) {
	my %seen;
	my @x;
	
	foreach (@_) {
		push @x, $_ if !exists $seen{$_};
		$seen{$_} = 1;
	}
	return @x;
}

sub readLine($$) {
	my ($fn, $lineNo) = @_;
	open my $f, "<", $fn;
	for (my $i = 1; $i < $lineNo; ++$i) {
		<$f>;
	}
	
	return scalar <$f>;
}

sub trim($) {
	my ($x) = @_;
	$x =~ s/\A\s+//;
	$x =~ s/\s+\z//;
	return $x;
}

# Works differently on Linux and Windows.
# Nowadays we try to scrape timing info from 
sub scrapeReductionLog($) {
	my ($logFn) = @_;
	
	my ($totalSecs, $nonIoSecs, $maxMemKb, $exitCode);
	open my $f, "<", $logFn;
	local $_;
	while (<$f>) {
		if (/\AIters\tTotSecs\tAvgSecs\t(EdgeDel\t)?Action\s*\z/) {
			my $hasEdgeDelCol = (defined $1);			# Was missing from earlier versions
			while (<$f>) {
				chomp;
				s/\r\z//;
				last if /===== Total =====/;
				#my (undef, $actionTotSecs, undef, $action) = split;
				my ($actionTotSecs, $action);
				if ($hasEdgeDelCol) {
					(undef, $actionTotSecs, undef, undef, $action) = split;
				} else {
					(undef, $actionTotSecs, undef, $action) = split;
				}
				$totalSecs += $actionTotSecs;
				$nonIoSecs += $actionTotSecs if $action ne 'read' && $action ne 'write';
			}
		}
		
		if (/\Astopwatch: M,Process completed with exit code (-?\d+)\.\s*\z/) {
			$exitCode = $1;
		} elsif (/\Astopwatch: T,(.*?)\s*\z/) {
			my @metrics = split /,/, $1;
			$maxMemKb = int($metrics[7] / 1024);	# We use the "Max Working Set Size", not the "Max VM Size".
		}
	}
	
	foreach ($totalSecs, $nonIoSecs, $maxMemKb, $exitCode) {
		die if !defined;		#HACK
	}
	
	close $f;
	return ($totalSecs, $nonIoSecs, $maxMemKb, $exitCode);
}

my $keepOnlySensibleFNames = 1;			# Allows user to specify wildcards much more easily
my $force = 0;							# If 0, only recreate nonexistent or out-of-date files.
my @fileSpecs;							# Input files.

# Process command-line arguments
while (@ARGV) {
	local $_ = shift;
	
	if ($_ eq '--target-dir') {
		$targetDir = shift;
	} elsif ($_ eq '-i' || $_ eq '--input') {
		push @fileSpecs, shift;
	} elsif ($_ eq '-f' || $_ eq '--force') {
		$force = 1;
	} elsif ($_ eq '--no-filter') {
		$keepOnlySensibleFNames = 0;
	} elsif ($_ eq '--') {
		# Remainder is ft_reduce arguments
		$FTREDUCECMD = safelyQuote @ARGV;
		last;
	} else {
		die "Unrecognised command-line argument '$_'.";
	}
}

mkdir $targetDir if !-d $targetDir;
push @fileSpecs, p "test_data/*.txt" if !@fileSpecs;		# Just a default

$FTREDUCECMD = "$STOPWATCH $FTREDUCECMD";
logMsg "# FTREDUCECMD=$FTREDUCECMD\n";
logMsg "# Target dir=$targetDir\n";
logMsg "# Input files: ", join(', ', @fileSpecs), "\n";

@fileSpecs = map { glob $_ } @fileSpecs if $S eq '\\';	# Only expand wildcards on Win
@fileSpecs = grep { /\.txt\z/ && !/\.(?:minimal)\.txt\z/ } @fileSpecs if $keepOnlySensibleFNames;

my %broken;			# Maps filenames to number of missing edges
my $nReductions = 0;

my $summaryFn = p "$targetDir/summary.tdf";
open my $summaryF, ">>", $summaryFn;		# Append
print $summaryF join("\t", qw(path nedges_orig nedges_reduced nedges_missing total_secs nonio_secs mem_kb exit_code)), "\n";

foreach my $inFn (@fileSpecs) {
	my $bfn = bn $inFn;
	my $outFn = j $targetDir, bn $inFn;
	$outFn =~ s/\.txt\z/.reduced.txt/ or die;
	my $logOutFn = "$outFn.ftr.log";
	
	print STDERR "<$inFn> -> <$outFn>, <$logOutFn>\n";		#DEBUG
	if ($force || !-e $outFn || -M $outFn > -M $inFn) {
		# Nonexistent or out-of-date: remake.
		run "$FTREDUCECMD < \"$inFn\" > \"$outFn\" 2> \"$logOutFn\"";
		my $nOrigEdges = trim readLine($inFn, 2);
		my $nReducedEdges = trim readLine($outFn, 2);
		logMsg "# Edges for $bfn: reduced from $nOrigEdges to $nReducedEdges\n";
		++$nReductions;
		
		# Scrape the execution time
		my ($totalSecs, $nonioSecs, $maxMemKb, $exitCode) = scrapeReductionLog $logOutFn;
		
		# Compare with the minimal problem instance.
		my $minimalFn = $inFn;
		$minimalFn =~ s|\.txt\z|.minimal.txt| or die;
		die "Cannot find minimal problem instance '$minimalFn'!" if !-e $minimalFn;
		
		my $missingEdgesFn = $outFn;
		$missingEdgesFn =~ s|\.txt\z|.missing_edges.txt|;
		run "transform_edges.pl \"$minimalFn\" - \"$outFn\" > \"$missingEdgesFn\"";
		
		# See how many edges are missing
		my $nMissingEdges = trim readLine($missingEdgesFn, 2);
		logMsg "# Missing edges in reduced instance of $bfn: $nMissingEdges\n";
		if ($nMissingEdges != 0) {
			$broken{$inFn} = $nMissingEdges;
		}
		
		# Write out a nice summary line
		print $summaryF join("\t", $inFn, $nOrigEdges, $nReducedEdges, $nMissingEdges, $totalSecs, $nonioSecs, $maxMemKb, $exitCode), "\n";
	} else {
		logMsg "# $outFn: already up-to-date.\n";
	}
}

logMsg "\n# Reductions: ", scalar(@fileSpecs), " files examined, $nReductions actually reduced (rest were up-to-date)\n";
logMsg "# SUMMARY of files that produced missing edges:\n";
foreach my $fn (sort keys %broken) {
	logMsg "# $fn: $broken{$fn}\n";
}

if (!keys %broken) {
	logMsg "# ALL FILES CONTAIN ALL EDGES :-)\n";
}

logMsg "# Summary file: $summaryFn\n";
