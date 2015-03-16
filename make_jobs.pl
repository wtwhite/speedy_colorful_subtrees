#!/usr/bin/perl

# /usr/bin/time make combined/mpos43918.ms/__COMBINED__fixed_colours.grb.sol

#HACK: For some reason, qsub chokes if I try to add these lines to the scripts:
#\$ -o $dir
#\$ -e $dir

# For convenient use with -s.
sub cmpCanon($$) {
	if ($_[0] =~ /\.canon\z/) {
		if ($_[1] =~ /\.canon\z/) {
			return 0;
		} else {
			return 1;
		}
	} elsif ($_[1] =~ /\.canon\z/) {
		return -1;
	} else {
		return 0;
	}
}

# Put bigger single-digit reductions *earlier*, since they're more likely to *solve* to completion.
sub gradeReduction($) {
	my ($d) = @_;
	$d = 0 if !length $d;		#HACK: Don't know what to do with it actually
	$d = 9 - $d if length $d == 1;
	return $d;
}

sub cmpReduce($$) {
	if (my ($d1) = $_[0] =~ /\.reduced(\d*)/) {
		if (my ($d2) = $_[1] =~ /\.reduced(\d*)/) {
			return gradeReduction($d1) <=> gradeReduction($d2);
		} else {
			return -1;
		}
	} elsif ($_[1] =~ /\.reduced(\d*)/) {
		return 1;
	} else {
		return 0;
	}
}

#my $transformExpr = "";		# Do nothing
my @transformExprs;
#if (@ARGV && $ARGV[0] eq '-t') {
#	shift;
#	$transformExpr = shift;
#}
my $filterExpr = 1;		# If true, a file will be considered, otherwise it will be skipped.
my $jobNameExpr = 's!([^/]+)\z!make.$1.sh!';		# Used by default to decide which job filename to put a given target file's make command into.  If more than one target maps to the same job filename, they will appear as separate commands in that file (in unspecified order).
my $sortExpr = 'cmpCanon($a, $b) or cmpReduce($a, $b) or $a cmp $b';		# Used by default to put targets ending with .canon last, then targets containing .reduce last, then alphabetically.
#my @files = @ARGV;
my @files;
my $dryRun = 0;
my @makeOptions;		# Extra arguments that will be passed to each make invocation (e.g. variable settings)

# Default to reading filenames from the command line.  Can be changed with -i.
my $get_next_file = sub { $_ = shift @files; };
my $readingFNamesFromStdin = 0;

# Process command-line arguments
while (@ARGV) {
	$_ = shift;
	if ($_ eq '-t') {
		#$transformExpr = shift;
		push @transformExprs, shift;
	} elsif ($_ eq '-f') {
		$filterExpr = shift;
	} elsif ($_ eq '-j') {
		$jobNameExpr = shift;
	} elsif ($_ eq '-n' || $_ eq '--dry-run') {
		$dryRun = 1;
		print STDERR "Dry run only -- no job files will actually be created.\n";
	} elsif ($_ eq '-O' || $_ eq '--make-option') {
		push @makeOptions, shift;
	} elsif ($_ eq '-i' || $_ eq '--from-stdin') {
		$get_next_file = sub { $_ = <>; if (defined) { chomp; s/\r\z//; } $_ };
		$readingFNamesFromStdin = 1;
		print STDERR "Will read files from STDIN instead of the command line.\n";
	} elsif ($_ eq '--') {
		push @files, @ARGV;
		last;
	} elsif (/\A-/) {
		die "Unrecognised option '$_'.";
	} else {
		push @files, $_;
	}
}

push @transformExprs, "" if !@transformExprs;		# If no transformations were specified, default to a single no-op.

my $makeOptionsStr = join(' ', map { my $x = $_; $x =~ s/(\W)/\\$1/g; $x } @makeOptions);
print STDERR "Transform expressions (-t):\n";
foreach (@transformExprs) {
	print STDERR "<$_>\n";
}
print STDERR "Filter expression (-f): <$filterExpr>\n";
print STDERR "Job filename expression (-j): <$jobNameExpr>\n";
print STDERR "Extra options to pass to make: <$makeOptionsStr>\n";
print STDERR "You didn't specify -i but there are no files on the command line: nothing to do!\n" if !@files && !$readingFNamesFromStdin;

sub poor_mans_realpath($) {
	chomp(my $x = `cd '$_[0]' && pwd`);		#HACK: Ugh.
	return $x;
}

# Should be run from /home/xe53kiw/frag_tree_reduction, otherwise the makefile rules won't work.
my $cwd = `pwd`;
chomp $cwd;
die "You should really run this from /home/xe53kiw/frag_tree_reduction -- if you're sure you know what you're doing, disable this die() call and rerun." if $cwd ne '/home/xe53kiw/frag_tree_reduction';	#HACK

@files = glob "combined/mpos*.ms/__COMBINED__.txt" if !@files;
my %targets;		# @{$targets{"some_job_name.sh"}} is a list of targets that will be made by the job script "some_job_name.sh".

#foreach (<combined/mpos*.ms/__COMBINED__.txt>) {
#foreach (@files) {
#foreach (grep { eval $filterExpr } @files) {
while (defined $get_next_file->()) {		# Reads the next filename into $_, from either the command line or STDIN
	next if !eval $filterExpr;
	
	#my ($job) = m!/(mpos.*?)/! or die;
	#my $job = $_;
	#$job =~ s!.*?/!!;
	#$job =~ tr!/!_!;
	my $kludge = /__COMBINED__/ ? "fixed_colours" : "";		#HACK
#	s/\.txt$/$kludge.grb.sol.canon/ or die;
#	eval($transformExpr);		# This is where you need to transform the current filename
	my $savedCopy = $_;
	foreach my $tExpr (@transformExprs) {
		$_ = $savedCopy;
		eval($tExpr);		# This is where you need to transform the current filename
		#print STDERR "<$_>. job=<$job>\n";
		#open my $jobF, ">", "do_orig_$job.sh" or die;
		#open my $jobF, ">", "$dir/do_orig_$job.sh" or die;
		my $jobFn = $_;
	#	$jobFn =~ s!([^/]+)\z!make.$1.sh!;
	#	print STDERR "Writing job script file '$jobFn'...\n";
		foreach ($jobFn) {			# This foreach allows the job name transformation to be conveniently done on $_ :)
			eval($jobNameExpr);
		}
		print STDERR "Adding make command for target '$_' to job script file '$jobFn'...\n";
		push @{$targets{$jobFn}}, $_;
	}
}

# Now actually generate the job scripts
foreach my $jobFn (keys %targets) {
	print STDERR "Writing " . scalar(@{$targets{$jobFn}}) . " make commands to job script file '$jobFn'...\n";
	
	my ($dir, $bn) = $jobFn =~ m!\A(.*)/(.*)\z!;
	if (!defined $dir) {
		$dir = ".";
		$bn = $_;
	}
	
	my $absDir = poor_mans_realpath $dir;		# Necessary because SGE uses $HOME-relative paths for -o and -e!
	
	if (!$dryRun) {
		open my $jobF, ">", $jobFn or die;

		print $jobF <<THE_END;
#!/bin/sh
#\$ -wd $absDir
#\$ -o .
#\$ -e .
echo Host: `hostname`
cd $cwd
export PATH=\$PATH:.
echo Date:
date
echo Environment:
env
trap 'echo "Caught SIGUSR2, will be terminated soon!"' USR2
THE_END
		
		#foreach (@{$targets{$jobFn}}) {
		foreach (sort { eval($sortExpr) } @{$targets{$jobFn}}) {		# $a and $b appear to be visible to the eval'ed expression, thank God
			print $jobF "/usr/bin/time make $makeOptionsStr $_\n";
		}
		
		close $jobF;
	}
}
