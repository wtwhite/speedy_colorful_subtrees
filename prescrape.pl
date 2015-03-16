#!/usr/bin/perl

use strict;
use warnings;
use autodie;

my $FUDGE = 10 / (60 * 60 * 24);		# Allow a fudge factor of 10s when looking for output files matching some given target file.

# Look in a couple of possible places for the ofile corresponding to a target file
sub tryStuff($) {
	my ($targetFn) = @_;
	
	my @poss;
	my $tryPat;
	if ($targetFn =~ /\.sol\z/) {
		$tryPat = $targetFn;
		$tryPat =~ s|([^/\\]+)\z|make.$1.canon.sh.o*| and print STDERR "TRYING: $tryPat\n" and push @poss, glob $tryPat;
		#$tryPat =~ s|([^/\\]+)\z|make.$1.canon.sh.o*| and push @poss, glob $tryPat;
		$tryPat = $targetFn;
		$tryPat =~ s!([^/\\]+?)(?:\.reduced[46])?\.(?:orig|cut)_cplex\.sol\z!make.$1.reduced046.both_cplex.sol.canon.sh.o*! and print STDERR "TRYING: $tryPat\n" and push @poss, glob $tryPat;
		#$tryPat =~ s!([^/\\]+?)(?:\.reduced[46])?\.(?:orig|cut)_cplex\.sol\z!make.$1.reduced046.both_cplex.sol.canon.sh.o*! and push @poss, glob $tryPat;
	} elsif ($targetFn =~ /\.reduced\d*\.txt\z/) {
		$tryPat = $targetFn;
		$tryPat =~ s|([^/\\]+)\z|make.$1.sh.o*| and print STDERR "TRYING: $tryPat\n" and push @poss, glob $tryPat;
		#$tryPat =~ s|([^/\\]+)\z|make.$1.sh.o*| and push @poss, glob $tryPat;
		$tryPat = $targetFn;
		$tryPat =~ s!([^/\\]+?)(?:\.reduced[46])?\.txt\z!make.$1.reduced046.both_cplex.sol.canon.sh.o*! and print STDERR "TRYING: $tryPat\n" and push @poss, glob $tryPat;
		#$tryPat =~ s!([^/\\]+?)(?:\.reduced[46])?\.txt\z!make.$1.reduced046.both_cplex.sol.canon.sh.o*! and push @poss, glob $tryPat;
	} else {
		die "tryStuff(): Don't know what to do with '$targetFn'!";
	}
	
	return @poss;
}

# Attempts to dig up the oldest *.sh.o* file that could correspond to a given target file.
sub find_matching_ofile($) {
	my ($targetFn) = @_;
	
	my @poss;
	push @poss, tryStuff($targetFn);
	my $prebucketJobsFn = $targetFn;
	if ($prebucketJobsFn =~ s!/b\.[0-9a-f]+/!/prebucket_jobs/!) {
		push @poss, tryStuff($prebucketJobsFn);
	}
	
	#print STDERR "poss=<", join(", ", @poss), ">\n";		#DEBUG
	
	my $targetAge = -M $targetFn;
	@poss = sort { -M $b <=> -M $a } grep { -M $_ <= $targetAge + $FUDGE } @poss;
	
	return $poss[0];		# Will be undef if there were no candidates.
}

# Try and scrape the most important values out of the ofile and into a hashref that can be written out as a .run file.
# Note: All environment variable lines are lumped into a single value under the key "env", with newlines separating them.
# This avoids the possibility that some environment variable name will collide with any of the values like command_line that we want to manage specially.
sub read_ofile($$) {
	my ($tfn, $ofn) = @_;
	
	my %rec;
	$rec{scraped_from_ofile} = 1;		# Might be handy for debugging
	$rec{completed_time} = "" . localtime((stat $tfn)[9]);		# Get a human-readable date from $tfn's modification time
	
	open my $f, "<", $ofn or die "$ofn: $!";
	local $_;
	while (<$f>) {
		chomp;
		s/\r\z//;
		if (/\AHost: (.*)/) {
			$rec{hostname} = $1;
		} elsif ($_ eq 'Date:') {
			# Ignore this actually (it's on the following line) -- it's the start date, we're only interested in the end date.
			<$f>;
		} elsif ($_ eq 'Environment:') {
			# Ignore this too
		} elsif (/\A\S+=/) {
			# It's (probably) an environment variable setting
			if (/\AJOB_ID=(.*)/) {
				$rec{job_id} = $1;
			} elsif (/\APWD=(.*)/) {
				$rec{cwd} = $1;
			}
			
			$rec{env} = "" if !exists $rec{env};
			$rec{env} .= "$_\n";
		} else {
			# Must be a command.
			if (($tfn =~ /\.reduced\d*\.txt\z/ && /\bft_reduce\b/) || ($tfn =~ /\.sol\z/ && /\bsolve_(?:cplex|gurobi)\b/)) {
				$rec{command_line} = $_;
			}
		}
	}
	
	return \%rec;
}

my %explains;		# $explains{'path/to/abc.reduced.txt'} is an arrayref whose elements correspond to ofiles -- each contains a hashref with the juicy details.

# The idea here is that we should call this on EVERY *.sh.o* file, and figure out which target files (reductions and solutions) each ofile explains.
sub study_ofile($) {
	my ($ofn) = @_;
	
	my %rec;
	$rec{scraped_from_ofile} = $ofn;		# Might be handy for debugging
	$rec{age} = -M $ofn;		# The finishing time OF THE SCRIPT!  So it won't be used for completed_time for an individual target.  Use the age in seconds for simple comparisons.
	
	open my $f, "<", $ofn or die "$ofn: $!";
	local $_;
	while (<$f>) {
		chomp;
		s/\r\z//;
		if (/\AHost: (.*)/) {
			$rec{hostname} = $1;
		} elsif ($_ eq 'Date:') {
			# Ignore this actually (it's on the following line) -- it's the start date, we're only interested in the end date.
			<$f>;
		} elsif ($_ eq 'Environment:') {
			# Ignore this too
		} elsif (/\A\S+=/) {
			# It's (probably) an environment variable setting
			if (/\AJOB_ID=(.*)/) {
				$rec{job_id} = $1;
			} elsif (/\APWD=(.*)/) {
				$rec{cwd} = $1;
			}
			
			$rec{env} = "" if !exists $rec{env};
			$rec{env} .= "$_\n";
		} else {
			# Must be a command.
			if (/\bft_reduce\s.*?[^2]>\s*(\S+)/) {
				push @{$explains{$1}}, { %rec, command_line => $_ };		# Should be safe to do this immediately as commands always come after all other stuff is loaded into %rec.
				print STDERR "EXPLAINS: $1 explained by $ofn\n";
			} elsif (/\bsolve_(?:cplex|gurobi)\s.*?\s-o\s+(\S+)/) {
				push @{$explains{$1}}, { %rec, command_line => $_ };		# Should be safe to do this immediately as commands always come after all other stuff is loaded into %rec.
				print STDERR "EXPLAINS: $1 explained by $ofn\n";
			}
		}
	}
	
	#return \%rec;
}

# You must have called study_ofile() on a bunch of *.sh.o* files before calling this on any target file.
sub lookup_explanation($) {
	my ($tfn) = @_;
	
	my @explanations;
	if (exists $explains{$tfn}) {
		print STDERR "TRYING: $tfn (", scalar(@{$explains{$tfn}}), " explanations)\n";
		push @explanations, @{$explains{$tfn}};
	} else {
		print STDERR "TRYING: $tfn (no explanations)\n";
	}
	
	# Try back-converting bucket names too.  HACK: Modifies $tfn in the process...
	if ($tfn =~ s!/b\.[0-9a-f]+/!/!) {
		if (exists $explains{$tfn}) {
			print STDERR "TRYING: $tfn (", scalar(@{$explains{$tfn}}), " explanations)\n";
			push @explanations, @{$explains{$tfn}};
		} else {
			print STDERR "TRYING: $tfn (no explanations)\n";
		}
	}
	
	@explanations = sort { $a->{age} <=> $b->{age} } @explanations;		# Get the most recent explanation first
	return $explanations[0];		# Will be undef if no explanations exist.
}

# Given a hashref containing (not necessarily all) fields that
sub write_run_file($$) {
	my ($fn, $rec) = @_;
	
	open my $f, ">", $fn;
	foreach (qw/command_line job_id hostname completed_time cwd/) {
		print $f "$_=", (defined($rec->{$_}) ? $rec->{$_} : "UNKNOWN"), "\n";
	}
	foreach ('scraped_from_ofile') {		# Just to save typing
		print $f "$_=$rec->{$_}\n" if exists $rec->{$_};
	}
	if (exists $rec->{env}) {
		print $f $rec->{env};
	}
	
	close $f;
}

my $dryRun = 0;
if (@ARGV && $ARGV[0] eq '-n') {
	$dryRun = 1;
	print STDERR "Dry run mode.\n";
	shift @ARGV;
}

#my @files = @ARGV;
#foreach my $tfn (@files) {
#	if ($tfn =~ /\.reduced\d*\.txt\z/ || $tfn =~ /\.sol\z/) {
#		my $rfn = "$tfn.run";
#		if (-e $rfn) {
#			print STDERR "ALREADY: $tfn already has a .run file, $rfn.\n";
#		} else {
#			my $ofn = find_matching_ofile $tfn;
#			if (defined $ofn) {
#				print STDERR "MATCHED: $tfn is best-matched with $ofn.\n";		#HACK
#				my $rec = read_ofile $tfn, $ofn;
#				if (!$dryRun) {
#					write_run_file($rfn, $rec);
#				}
#			} else {
#				print STDERR "UNMATCHED: $tfn could not be matched with any ofile!  (Is the fudge factor of $FUDGE too tight?)\n";
#			}
#		}
#	} else {
#		print STDERR "IGNORING: $tfn is not a recognised target filename.\n";
#	}
#}

my @targets;
my %alreadyExists;		# .run files that already exist

# Read files from STDIN
while (<>) {
	chomp;
	s/\r\z//;
	if (/\.reduced\d*\.txt\z/ || /\.sol\z/) {
		print STDERR "TARGET: $_\n";
		push @targets, $_;
	} elsif (/\.sh\.o\d+\z/) {
		# An ofile: study it.
		print STDERR "STUDYING: $_\n";
		study_ofile($_);
	} elsif (/\.run\z/) {
		# An already-existing runfile.
		print STDERR "RUNFILE: $_\n";
		$alreadyExists{$_} = 1;
	} else {
		print STDERR "IGNORING: $_\n";
	}
}

# Now go through all the targets.  Any that don't have a runfile already, we need to look up info about them that we got from studying ofiles.
foreach (@targets) {
	if (exists $alreadyExists{"$_.run"}) {
		print STDERR "ALREADY: $_\n";
	} else {
		my $rec = lookup_explanation($_);
		if (defined $rec) {
			print STDERR "MATCHED: $_ matched to ofile $rec->{scraped_from_ofile}\n";
			if (!$dryRun) {
				write_run_file("$_.run", $rec);
			}
		} else {
			print STDERR "UNMATCHED: $_\n";
		}
	}
}
