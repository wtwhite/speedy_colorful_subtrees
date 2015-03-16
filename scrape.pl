#!/usr/bin/env perl

use strict;
use warnings;
use autodie ':all';
use DBI;
use Digest::MD5;		# For calculating canon_tail_md5

#die "This is designed to run only on Unix." if exists $ENV{OS} && $ENV{OS} =~ /\AWin/;		#HACK

chomp(my $thisHost = `hostname`);

my $verbose = 0;		#HACK: Currently no way to turn this on

my $NUMREGEX = '[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?';		#HACK: Probably still wrong somehow, but it gets close to matching any valid FP number.

# Open DB connection and set up its closing
#my $dbh = DBI->connect('dbi:Pg:dbname=wtwhite;host=localhost;port=15433', undef, undef, { RaiseError => 1 });		# Works on my laptop (via SSH tunnel) and on darwin!  :)
my $dbPort = ($thisHost eq 'wallace' ? 5432 : 15433);			# On wallace we access the DB directly; on all other machines we need an SSH tunnel.
#my $dbh = DBI->connect("dbi:Pg:dbname=wtwhite;host=localhost;port=$dbPort", 'wtwhite', 'asdfasdf', { RaiseError => 1, ShowErrorStatement => 1, TraceLevel => 'SQL' });		#HACK: Disgusting that I include the password...
my $dbh = DBI->connect("dbi:Pg:dbname=wtwhite;host=localhost;port=$dbPort", 'wtwhite', 'asdfasdf', { RaiseError => 1, ShowErrorStatement => 1 });		#HACK: Disgusting that I include the password...
END {
	if (defined $dbh) {
		$dbh->disconnect;
		$dbh = undef;
	}
}

sub insert_into($$) {
	my ($tableName, $rec) = @_;
	
	my $sql = "INSERT INTO $tableName (" . join(", ", sort keys %$rec) . ") VALUES (" . join(", ", map { "?" } sort keys %$rec) . ")";
	my @paramValues = map { $rec->{$_} } sort keys %$rec;
	
	$dbh->do($sql, {}, @paramValues);
}

# INSERTs data into a temporary table corresponding to some other existing table.  Creates this table if necessary.
my $tempTableType = "";		# Set to " TEMP" to create all tmp_* tables as TEMP tables
my %alreadyCreatedTempTableFor;		# Caches the names of tables for which we have already created temporary tables whose names begin with "tmp_".
my $nTempTables = 0;		# How many temp tables have we created?  Also used for automatically ordering the updating of master tables (they'll be updated in the order of temp table creation).
sub insert_into_temp($$) {
	my ($tableName, $rec) = @_;
	
	my $tmpTableName = "tmp_$tableName";
	if (!exists $alreadyCreatedTempTableFor{$tableName}) {
		# It might already exist in the DB (e.g. if it was not created as a TEMP table at the DB level)
		$dbh->do("CREATE$tempTableType TABLE IF NOT EXISTS $tmpTableName (LIKE $tableName INCLUDING DEFAULTS)");		# Even nicer than SELECT * INTO, due to IF NOT EXISTS and INCLUDING DEFAULTS which lets us get sensible create_times :)
		$alreadyCreatedTempTableFor{$tableName} = ++$nTempTables;
	}
	
	return insert_into $tmpTableName, $rec;
}

#HACK: We can't use "*" as the SELECT list because we need to exclude the create_time field!
my %fieldsToIgnore = ('create_time' => 1);
sub get_relevant_fields_for_table($) {
	my ($table) = @_;
	return join(", ", sort grep { !exists $fieldsToIgnore{$_} } @{$dbh->selectcol_arrayref("select column_name from information_schema.columns where table_name = ?", {}, $table)});
}

# Any rows in tmp_XYZ not already present in XYZ will be copied across, and then deleted from tmp_XYZ.
# This whole shenanigans is because we don't have a nice way to UPSERT...
# Well: Actually we don't UPDATE any rows that already match (since they have to match EXACTLY anyway), so it's not quite an UPSERT.  A true UPSERT would require a 2nd temp table (I think).
# If there have been changes to some row that don't result in a new PK, then we'll get a referential integrity failure, which should be OK.
# DAMMIT, the create_time field means that EXCEPT ALL won't work!
sub update_table_from_temp($) {
	my ($tableName) = @_;
	my $relevantFieldsStr = get_relevant_fields_for_table $tableName;
	print STDERR "Updating table $tableName from tmp_$tableName...\n";
	#$dbh->do("INSERT INTO $tableName (SELECT * FROM tmp_$tableName EXCEPT ALL SELECT * FROM $tableName)");		# EXCEPT ALL is nice :)
	my $nRowsInserted = 0 + $dbh->do("INSERT INTO $tableName ($relevantFieldsStr) (SELECT $relevantFieldsStr FROM tmp_$tableName EXCEPT ALL SELECT $relevantFieldsStr FROM $tableName)");		# EXCEPT ALL is nice :)
	my $nRowsDeleted = 0 + $dbh->do("DELETE FROM tmp_$tableName");
	print STDERR "$nRowsInserted of $nRowsDeleted rows were completely new for table $tableName.\n";
}

# Updates master tables in the order their tmp_ versions were first INSERTed into.
sub update_all_tables_from_temp() {
	my @updateOrder = sort { $alreadyCreatedTempTableFor{$a} <=> $alreadyCreatedTempTableFor{$b} } keys %alreadyCreatedTempTableFor;
	foreach my $tn (@updateOrder) {
		update_table_from_temp $tn;
	}
}

# /usr/bin/time on some versions of Linux (e.g. those running on the cluster nodes) returns 4x the memory usage, and not on others!
#HACK: We guess whether to divide by 4 based on the hostname...
#HACK: The hostname might be undef (e.g. if no runfile exists), in which case we make the thoroughly unconservative guess that it was running on an exc-* node anyway and divide by 4...
sub calc_correct_mem($$) {
	my ($hostname, $kb) = @_;
	#if ($hostname =~ /\Aexc-\d+\z/) {
	if (!defined $hostname || $hostname =~ /\Aexc-\d+\z/) {
		return int($kb / 4);
	} else {
		return $kb;
	}
}

sub break_path($) {
	if ($_[0] =~ m!\A(.*)[\\/](.*)\z!) {
		return ($1, $2);
	} else {
		return (undef, $_);
	}
}

# Input: An hostname (needed for converting memory usage correctly -- yes, this is insane) and array of lines that should contain *only* timing data in an acceptable form (e.g. exactly 2 lines containing info from /usr/bin/time, or when
# I get around to implementing it, output lines from stopwatch.exe).
# Extract the timing and memory usage info.  If anything else (besides blank lines) appears, it is treated as an error, which is recorded by appending it to $rec{result_type}.
# Note that we don't actually extract much of the info there (like timings) because we already have a better breakdown from the output of the program itself.
sub extract_timing_info($@) {
	my $hostname = shift;		# This could be undef if no runfile was found!
	my %rec;
	my $line = 1;
	
	foreach my $x (@_) {
		if ($line == 1 && $x =~ /\A([\d\.]+)user\s+([\d\.]+)system\s+(?:(?:(\d+):)?(\d+):)?([\d\.]+)elapsed\s+([\d\.]+)\%CPU\s+\(\d+avgtext\+\d+avgdata\s+(\d+)maxresident\)k\s*\z/) {
			#@rec{qw(user_time system_time elapsed_time percent_cpu max_mem_kb)} = ($1, $2, ((defined($3) ? $3 : 0) * 60 + (defined($4) ? $4 : 0)) * 60 + $5, $6, $7);
			#$rec{max_mem_kb} = calc_correct_mem($rec{hostname}, $7);		# /usr/bin/time reports 4x the memory usage on some versions of Linux and not others!
			$rec{max_mem_kb} = calc_correct_mem($hostname, $7);		# /usr/bin/time reports 4x the memory usage on some versions of Linux and not others!
			#@rec{qw(user_secs system_secs percent_cpu max_mem_kb)} = ($1, $2, $6, $7);		# We now get elapsed_secs from the internally-reported action timing table
			@rec{qw(user_secs system_secs percent_cpu)} = ($1, $2, $6);		# We now get elapsed_secs from the internally-reported action timing table
			$line = 2;
		} elsif ($line == 1 && $x =~ /\Astopwatch: Terminated. Elapsed time: (\d+)ms\s*\z/) {
			#$rec{elapsed_time} = $1 / 1000;
			die "We don't handle this on Windows yet...";
		} elsif ($line == 2 && $x =~ /\A\d+inputs\+\d+outputs\s*\(\d+major\+\d+minor\)pagefaults\s+\d+swaps\s*\z/) {		#0inputs+5928outputs (0major+61121minor)pagefaults 0swaps
			$line = 3;		# Don't actually record anything
		} else {
			#warn "Second-to-last line does not contain the elapsed time etc.!\n";
			my $y = $x;		# Don't want to modify the caller's variable!
			$y =~ s/\A\s+//;
			$y =~ s/\s+\z//;
			if (length $y) {
				$rec{result_type} = '' if !exists $rec{result_type};
				$rec{result_type} .= "[Expected line $line of timing info, got: $y]";
			}
		}
	}
	
	return \%rec;
}

# Reductions and solves create a file with the same name as the output file but with .run appended (so *.reduced.txt.run, *.sol.run), to store important info like the SGE job ID, hostname and completion time.
# This is loaded into the following hash.
my %runFile;
my %importantRunFileKeys = map { ($_ => 1) } qw/job_id command_line hostname completed_time/;		# Don't bother keeping every last env var.  Note both job_id and JOB_ID will be there!  The former is "canonical".
sub load_runfile($) {
	my ($fn) = @_;
	print STDERR "Loading file '$fn' as a runfile...\n";
	open my $f, "<", $fn;
	$runFile{$fn} = {} if !exists $runFile{$fn};
	local $_;
	while (<$f>) {
		chomp;
		s/\r\z//;
		if (/\A([^=]+)=(.*)\z/ && exists $importantRunFileKeys{$1}) {
			$runFile{$fn}{$1} = $2;
		}
	}
}

# Load timing info for a orig_cplex, cut_cplex or grb run from a .stderr file.
# We store all this in the $runfile hash keyed by the ".run" extension, because it's convenient.
sub load_stderrfile($) {
	my ($fn) = @_;
	my $rffn = $fn;
	$rffn =~ s|\.stderr\z|.sol.run| or die "Could not massage '$fn' into the corresponding runfile name!";
	print STDERR "Loading file '$fn' as a timing-info-containing stderr file...\n";
	open my $f, "<", $fn;
	chomp(my @lines = <$f>);
	my $timingInfo = extract_timing_info $runFile{$rffn}{hostname}, @lines;			# Need to run load_runfile() first to get the hostname to pass to extract_timing_info()
	$runFile{$rffn} = {} if !exists $runFile{$rffn};
	$runFile{$rffn}{$_} = $timingInfo->{$_} foreach grep { $_ ne 'result_type' } keys %$timingInfo;		# Merge it in, but treat result_type specially -- we wanna concatenate not overwrite
	if (exists $timingInfo->{result_type}) {
		$runFile{$rffn}{result_type} = '' if !exists $runFile{$rffn}{result_type};
		$runFile{$rffn}{result_type} .= $timingInfo->{result_type};
	}
}

# CPLEX puts important info in the solutionStatusString field of the .sol file.
sub load_cplexsolfile($) {
	my ($fn) = @_;
	my $rffn = "$fn.run";
	print STDERR "Loading file '$fn' as a CPLEX raw solution file...\n";
	open my $f, "<", $fn;
	local $_;
	while (<$f>) {
		if (/solutionStatusString="([^"]*)"/) {
			if ($1 ne 'integer optimal solution' && $1 ne 'optimal') {			# 'optimal' seems to appear when the problem instance contains 0 edges -- but only *sometimes*.  In such cases, CPLEX doesn't list any variable values, but they're all 0, and canonicalise_cplex_solution.pl happens to do the right thing :)
				$runFile{$rffn}{result_type} = '' if !exists $runFile{$rffn}{result_type};
				$runFile{$rffn}{result_type} .= "[CPLEX solutionStatusString: $1]";
			}
			
			return;
		}
	}
	
	$runFile{$rffn}{result_type} = '' if !exists $runFile{$rffn}{result_type};
	$runFile{$rffn}{result_type} .= "[Couldn't find CPLEX solutionStatusString]";
}

# Extract the first line, and create an MD5 hash of the rest.
sub load_canonfile($) {
	my ($fn) = @_;
	my $rffn = $fn;
	$rffn =~ s/\.canon\z/.run/ or die "Couldn't change .canon file $fn into corresponding runfile!";
	print STDERR "Loading file '$fn' as a .canon solution file...\n";
	open my $f, "<", $fn;
	local $_;
	$_ = <$f>;
	chomp;
	s/\r\z//;
	$runFile{$rffn}{canon_first_line} = $_;
	my $md5 = Digest::MD5->new();
	$md5->addfile($f);		# Reads the rest of the file
	$runFile{$rffn}{canon_tail_md5} = $md5->hexdigest();
}

# Pass it an output path (filename) for either a reduction or a solution, and a hashref record.  It will add in the fields from the corresponding runfile to the hashref.
# Actually it may add info gathered from other places too (e.g. timings and memory usage taken from .stderr files).
# Because of a stupid design mistake, we don't produce .run files if there was an error during processing (e.g. out of memory).  Because I still want those records in the DB,
# I'll just create some dummy info instead of erroring out.
sub inject_common_runfile_info($$) {
	my ($fn, $rec) = @_;
	#die "No runfile loaded for output path '$fn'!" if !exists $runFile{"$fn.run"};
	if (!exists $runFile{"$fn.run"}) {
		#$rec->{job_id} = "UNKNOWN";
		#$rec->{command_line} = "UNKNOWN";
		#$rec->{hostname} = "UNKNOWN";
		#$rec->{completed_time} = "3000-01-01 00:00:00";		# We could actually guess this based on the logfile's mod time, but it's too much effort...
		# Actually, just let all these fields become NULL in the DB.  The ok_* and bad_* DB VIEWs can be used to find what's relevant.
	} else {
		#foreach my $v (qw(job_id command_line hostname completed_time)) {
		#foreach my $v (keys %{$runFile{"$fn.run"}}) {
		foreach my $v (qw(job_id command_line hostname completed_time user_secs system_secs percent_cpu max_mem_kb canon_first_line canon_tail_md5)) {
			#die "Could not find $v value in runfile for '$fn'!" if !exists $runFile{"$fn.run"}{$v};
			#$rec->{$v} = $runFile{"$fn.run"}{$v};
			$rec->{$v} = $runFile{"$fn.run"}{$v} if exists $runFile{"$fn.run"}{$v} && !exists $rec->{$v};		# Don't overwrite any field already there
		}
		
		# Treat result_type specially: concatenate, don't overwrite
		if (exists $runFile{"$fn.run"}{result_type}) {
			$rec->{result_type} = '' if !exists $rec->{result_type};
			$rec->{result_type} .= $runFile{"$fn.run"}{result_type};
		}
	}
}

# Scrape info from the .ftr.log files produced by saving ft_reduce's stderr.
# The corresponding instance row must already exist in instances.
# Table structure:
#create table reductions (
#	input_path varchar not null references instances (path),		-- Refers to the path of the instance that was input to the reduction.
#	output_path varchar references instances (path),	-- Identifies the reduced instance.
#	job_id varchar not null,		-- Generally integer, but doesn't have to be.  Could use my own bespoke string IDs if I want.  Mainly used to allow multiple reps.
#	command_line varchar not null,		-- Should contain the reduction "program"
#	hostname varchar not null,
#	completed_time timestamp not null,
#	elapsed_secs numeric not null,
#	elapsed_nonio_secs numeric not null,	-- Includes all time except that used by "read" and "write" actions
#	user_secs numeric,
#	system_secs numeric,
#	percent_cpu numeric,
#	max_mem_kb numeric,
#	create_time timestamp not null default current_timestamp,
#	primary key (output_path, job_id)
#);
sub scrape_ftr_log($) {
	my ($fn) = @_;
	local $_;
	
	print STDERR "Scraping file '$fn' as an ft_reduce logfile...\n";
	open my $f, "<", $fn;
	
	my %rec;
	$_ = $fn;
	s/\.txt\.ftr.log\z// or warn "Could not strip .txt.ftr.log suffix from ft_reduce logfile's name!\n";
	#@rec{'path', 'filename'} = break_path $_;
	$rec{output_path} = "$_.txt";
	if (s/(\.reduced\d*)\z//) {
		$rec{reduction} = $1;
	} else {
		die "Could not figure out the reduction used for '$_'!";		# There must always be a reduction for a .ftr.log file!
	}
	$rec{input_path} = "$_.txt";
	#@rec{'path', 'filename'} = break_path $fn;
	$rec{dataset} = $_;			# Whatever's left.
#	$rec{input_path} =~ s|\.reduced\d+\.|.| or die "Could not figure out the input path for the reduction that produced '$_' from its filename!";		# input_path must be in the instances table!
#	$rec{hostname} = $thisHost;		#HACK: We don't have this information, so we hackfully guess that we are being run on the same machine that the reduction was run on...
	# Certain info is extracted from the .run file instead, which must have already been loaded.
	inject_common_runfile_info $rec{output_path}, \%rec;
	
	my @final;		# Will hold last 2 lines
	while (<$f>) {
		#if (/\AAbout to parse the following program: <(.*)>.\s*\z/) {
		#	$rec{reduction_script} = $1;
		#}
		
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
				$rec{elapsed_secs} += $actionTotSecs;
				$rec{elapsed_nonio_secs} += $actionTotSecs if $action ne 'read' && $action ne 'write';
			}
		}
		
		push @final, $_;
		shift @final if @final > 2;
	}
	
	# Linux /usr/bin/time format example:
	#1122.50user 7.40system 18:51.17elapsed 99%CPU (0avgtext+0avgdata 10703040maxresident)k
	#159176inputs+188336outputs (0major+1777375minor)pagefaults 0swaps
	
	# Windows stopwatch example:
	#stopwatch: Terminated. Elapsed time: 4203ms
	#stopwatch: Process completed with exit code 0.
	#if ($final[0] =~ /\A([\d\.]+)user\s+([\d\.]+)system\s+(?:(?:(\d+):)?(\d+):)?([\d\.]+)elapsed\s+([\d\.]+)\%CPU\s+\(\d+avgtext\+\d+avgdata\s+(\d+)maxresident\)k\s*\z/) {
	#	#@rec{qw(user_time system_time elapsed_time percent_cpu max_mem_kb)} = ($1, $2, ((defined($3) ? $3 : 0) * 60 + (defined($4) ? $4 : 0)) * 60 + $5, $6, $7);
	#	$rec{max_mem_kb} = calc_correct_mem($rec{hostname}, $7);		# /usr/bin/time reports 4x the memory usage on some versions of Linux and not others!
	#	#@rec{qw(user_secs system_secs percent_cpu max_mem_kb)} = ($1, $2, $6, $7);		# We now get elapsed_secs from the internally-reported action timing table
	#	@rec{qw(user_secs system_secs percent_cpu)} = ($1, $2, $6);		# We now get elapsed_secs from the internally-reported action timing table
	#} elsif ($final[0] =~ /\Astopwatch: Terminated. Elapsed time: (\d+)ms\s*\z/) {
	#	#$rec{elapsed_time} = $1 / 1000;
	#	die "We don't handle this on Windows yet...";
	#} else {
	#	warn "Second-to-last line does not contain the elapsed time etc.!\n";
	#}
	my $timingInfo = extract_timing_info $rec{hostname}, @final;
	$rec{$_} = $timingInfo->{$_} foreach keys %$timingInfo;		# Merge in all the fields
	return \%rec;
}

sub scrape_stephan_cplex_stdout($) {
	my ($fn) = @_;
	local $_;
	
	print STDERR "Scraping file '$fn' as a Stephan CPLEX .stdout file...\n";
	open my $f, "<", $fn;
	
	my %rec;
	$_ = $fn;
	s/\.([^\.]+_cplex)\.stdout\z// or warn "Could not strip .grb.log suffix from Gurobi logfile's name!\n";		# So far we have orig_cplex and cut_cplex.
	$rec{solver} = $1;
	$rec{input_path} = "$_.txt";
	$rec{output_path} = "$_.$1.sol";
	$rec{reduction} = '';		# Assume an original instance at first...
	s/(\.reduced\d*)\z// and $rec{reduction} = $1;
	#@rec{'path', 'filename'} = break_path $fn;
	$rec{dataset} = $_;			# Whatever's left.
#	$rec{hostname} = $thisHost;		#HACK: We don't have this information, so we hackfully guess that we are being run on the same machine that the reduction was run on...
	# Certain info is extracted from the .run file instead, which must have already been loaded.
	inject_common_runfile_info $rec{output_path}, \%rec;		# This will also populate result_type if any unusual lines (usually indicating errors) were found in the .stderr file
	
	#$rec{found_opt} = 0;
	my $readSecs = 0;
	#$rec{result_type} = 'ERROR';		#HACK: Ugh, we can't be more specific...
	my $resultType = 'ERROR';		#HACK: Ugh, we can't be more specific...
	$rec{ncuts_generated} = 0;
	$rec{nthreads} = 1;
	while (<$f>) {
		chomp;
		s/\r\z//;
		
		if (/\ATime limit: (\d+)s/) {
			$rec{time_limit_secs} = $1;
		} elsif (/\AMemory limit: (\d+)MB/) {
			$rec{mem_limit_kb} = $1 * 1024;
		} elsif (/\A  -- initialization time: ([\d.]+) s/) {
			$readSecs = $1;			# Technically this includes setup of the solver model and constraints etc., but oh well
		} elsif (/\A  -- total time used to solve problem: ([\d.]+) seconds/) {
			$rec{elapsed_secs} = $1;
		} elsif (/\A  \*\* Cuts found (by MinCut): (\d+)/) {
			$rec{ncuts_generated} += $1;
		} elsif (/\A\t([\d.]+)\s*\z/) {
			$rec{opt_sol_val} = $1;
		} elsif (/\AWriting out solution to/) {
			#HACK: Assume that if we made it this far, we found the correct solution.  Should really somehow check the .stderr file too for strangeness (e.g. termination due to signal).
			$rec{opt_bound} = $rec{opt_sol_var};
			#$rec{result_type} = 'OPT';
			$resultType = 'OPT';
		}
	}
	
	# $rec{result_type} needs to reflect anything discovered in the .stderr file too.
	$rec{result_type} = '' if !exists $rec{result_type};		# May have been created by errors in the .stderr
	$rec{result_type} .= $resultType;
	$rec{orig_result_type} = $rec{result_type};		#HACK
	if (defined $rec{elapsed_secs} && defined $readSecs) {
		$rec{elapsed_nonio_secs} = $rec{elapsed_secs} - $readSecs;
	}
	return \%rec;
}

# We ignore large chunks of the Gurobi log output.
sub scrape_grb_log($) {
	my ($fn) = @_;
	local $_;
	
	print STDERR "Scraping file '$fn' as a Gurobi logfile...\n";
	open my $f, "<", $fn;
	
	my %rec;
	$_ = $fn;
	s/\.grb\.log\z// or warn "Could not strip .grb.log suffix from Gurobi logfile's name!\n";
	#@rec{'path', 'filename'} = break_path $fn;
	$rec{solver} = 'grb';
	$rec{input_path} = "$_.txt";
	$rec{output_path} = "$_.grb.sol";
	$rec{reduction} = '';		# Assume an original instance at first...
	s/(\.reduced\d*)\z// and $rec{reduction} = $1;
	#@rec{'path', 'filename'} = break_path $fn;
	$rec{dataset} = $_;			# Whatever's left.
#	$rec{hostname} = $thisHost;		#HACK: We don't have this information, so we hackfully guess that we are being run on the same machine that the reduction was run on...
	# Certain info is extracted from the .run file instead, which must have already been loaded.
	inject_common_runfile_info $rec{output_path}, \%rec;
	
	#$rec{found_opt} = 0;
	my $readSecs = 0;
	my $resultType = 'ERROR';
	while (<$f>) {
		chomp;
		s/\r\z//;
		
		if (/\A(Gurobi Optimizer version .*?)\s*\z/) {
			#$rec{solver} = $1;
			$rec{command_line} .= " # $1";		#HACK: Encode the Gurobi version as a shell comment!
		} elsif (/\ASet parameter Threads to value (-?\d+)\s*\z/) {
			$rec{nthreads} = $1;
		} elsif (/\ASet parameter TimeLimit to value (\d+)\s*\z/) {
			$rec{time_limit_secs} = $1;
		} elsif (/\ASet parameter NodefileStart to value ([\d.]+)\s*\z/) {
			$rec{mem_limit_kb} = $1 * 1024 * 1024;
		#} elsif (/\A(.*): (\d+) rows, (\d+) columns, (\d+) nonzeros\s*\z/ && !exists $rec{nrows_before_presolve}) {
		#	if ($1 eq 'Presolved') {
		#		@rec{qw(nrows_after_presolve  ncols_after_presolve  nnonzeros_after_presolve)}  = ($2, $3, $4);
		#	} else {
		#		@rec{qw(nrows_before_presolve ncols_before_presolve nnonzeros_before_presolve)} = ($2, $3, $4);
		#	}
		#} elsif (/\AExplored (\d+) nodes? \(\d+ simplex iterations?\) in ([\d\.]+) seconds?\s*\z/) {
		#	@rec{qw(nnodes_explored elapsed_time)} = ($1, $2);
		} elsif (/\AOptimal solution found \(tolerance ([-+\d\.eE]+)\)\s*\z/) {
			#$rec{tolerance} = $1 + 0;
			#$rec{found_opt} = 1;
			#$rec{result_type} = 'OPT';
			$resultType = 'OPT';
		} elsif (/\AERROR\b/ || /\ATime limit reached/) {
			#$rec{result_type} = $_;		# The only error I've seen so far is "ERROR 10001: Out of memory", but there could be others.  A (generally subopt) obj func score will still be reported for "Time limit reached" it seems.
			$resultType = $_;		# The only error I've seen so far is "ERROR 10001: Out of memory", but there could be others.  A (generally subopt) obj func score will still be reported for "Time limit reached" it seems.
		#} elsif (/\ABest objective ([-+\d\.eE]+), best bound ([-+\d\.eE]+), gap .*$/) {		# Need $ instead of \z to match newline!
		} elsif (/\ABest objective ($NUMREGEX), best bound ($NUMREGEX), gap .*$/) {		# Need $ instead of \z to match newline!  HACK: Need fancy $NUMREGEX because $1 can be just "-", if Gurobi runs out of memory!
			$rec{opt_sol_val} = -$1 + 0;		# Remember, Gurobi has been instructed to *minimise* the negative of the real objective function.
			$rec{opt_bound} = -$2 + 0;
		} elsif (/\AReading time = ([\d.]+) seconds/) {
			$readSecs = $1;
		} elsif (/\AExplored \d+ nodes .*in ([\d.]+) seconds/) {
			$rec{elapsed_nonio_secs} = $1;
		}
	}
	
	# $rec{result_type} needs to reflect anything discovered in the .stderr file too.
	$rec{result_type} = '' if !exists $rec{result_type};		# May have been created by errors in the .stderr
	$rec{result_type} .= $resultType;
	$rec{orig_result_type} = $rec{result_type};		#HACK
	if (defined $rec{elapsed_nonio_secs} && defined $readSecs) {
		$rec{elapsed_secs} = $rec{elapsed_nonio_secs} + $readSecs;
	}
	$rec{ncuts_generated} = 0;			# We know this because if cuts were used then we'd be using Stephan's program instead, which produces completely different output.
	return \%rec;
}

# Table structure:
#create table instances (
#	path varchar primary key,	-- So all instances must have distinct filenames.  Include directory paths if necessary.
#	parent_path varchar references instances (path),	-- Non-NULL only for reduced instances
#	nverts integer not null,
#	nedges integer not null,
#	ncols integer not null,
#	opt_sol_val numeric,		-- Contains the known optimal solution value (or the first number if 2 are separated by a '+') if this is specified as line 4 in the file
#	opt_sol_extra_val numeric,	-- Contains the number after the '+' (the weight of the edge to this subtree) if line 4 specifies a known optimal solution value in 2 parts
#	create_time timestamp not null default current_timestamp
#);

sub scrape_instance($) {
	my ($fn) = @_;
	local $_;
	
	print STDERR "Scraping file '$fn' as a problem instance...\n";
	open my $f, "<", $fn;
	
	my %rec;
	#@rec{'path', 'filename'} = break_path $fn;
	$_ = $fn;
	$rec{path} = $_;
	s/\.txt\z// or die "Could not strip .txt suffix from instance file '$fn'!";
	$rec{reduction} = '';		# Assume an original instance at first...
	s/(\.reduced\d*)\z// and $rec{reduction} = $1;
	$rec{dataset} = $_;			# Whatever's left.
	($rec{data_group}) = m%([^/]+)/b\.[0-9a-f][0-9a-f]/[^/]+$%;		# DB actually has a CHECK for this
	$rec{is_combined} = 0 + m%/__COMBINED__fixed_colours$%;		# DB actually has a CHECK for this
	
	chomp($_ = <$f>);
	$rec{nverts} = $_;
	chomp($_ = <$f>);
	$rec{nedges} = $_;
	chomp($_ = <$f>);
	$rec{ncols} = $_;
	
	# Next line could be either the first colour assignment or a line giving the known optimal score -- either as a single number, or two numbers separated by a '+'
	$_ = <$f>;
	if (/\A([\d\.]+)(?:\s*\+\s*([\d\.]+))?\s*\z/) {
		$rec{opt_sol_val} = $1;
		$rec{opt_sol_extra_val} = $2;		# This can be undef
	}
	
	# Recording which instances are reduced is now done by joining to the reductions table and using its input_path and output_path fields.
	## If this looks like a reduced instance, look for a parent instance
	#if ($fn =~ s/\.reduced\d+\././) {
	#	$rec{parent_path} = $fn;
	#}
	
	return \%rec;
}

# Try something, and if an error occurs, do something else.
sub try($$) {
	my ($f, $handler) = @_;
	my @retVal;		#HACK: Forces list context
	
	eval {
		@retVal = $f->();
	};
	
	if (length $@) {
		return $handler->($@);
	} else {
		return @retVal;
	}
}

# Returns a closure that can be called to record the problem.
my @filesWithProblems;
sub problem_with($) {
	my ($fn) = @_;
	return sub {
		chomp $_[0];
		push @filesWithProblems, "$fn: $_[0]";
	};
}



# Main program

#HACK: Process filenames containing ".reduced." after files that don't have this, so that DB lookups for parent rows will find them.
sub count_word($$) {
	my ($needle, $haystack) = @_;
	
	my $n = 0;
	++$n while $haystack =~ /$needle/g;
	return $n;
}

#my @files = sort { count_word("\\.reduced\\.", $a) <=> count_word("\\.reduced\\.", $b) } map { exists $ENV{OS} && $ENV{OS} =~ /\AWin/ ? glob $_ : $_ } @ARGV;
#
#foreach (@files) {
#	if (/\.txt\z/) {
#		insert_into "instances", scrape_instance $_;
#	} elsif (/\.ftr\.log\z/) {
#		insert_into "reductions", scrape_ftr_log $_;
#	} elsif (/\.grb\.log\z/) {
#		insert_into "solutions", scrape_grb_log $_;
#	} else {
#		die "Don't recognise filename '$_'.";
#	}
#}

my @files;
my %processed;
if (@ARGV && $ARGV[0] eq '-i') {
	# Read files from stdin instead of the command line.
	shift @ARGV;
	local $_;
	while (<>) {
		chomp;
		s/\r\z//;
		push @files, $_;
	}
} else {
	# Process any files listed on the command line in the right order so that DB INSERTs will have FKs satisfied
	@files = map { exists $ENV{OS} && $ENV{OS} =~ /\AWin/ ? glob $_ : $_ } @ARGV;
}

print STDERR scalar(@files), " files to process.\n";

#HACK: Discard non fixed-coloured versions of COMBINED instances.  That means anything that has "__COMBINED__" not immediately followed by "fixed_colours".
@files = grep { !/\b__COMBINED__\./ } @files;
print STDERR scalar(@files), " files to process after eliminating __COMBINED__ instances without fixed_colours.\n";

# The scrape_...() functions for reduction logs and solver logs require that .run files have already been loaded into memory.
foreach (grep { /\.run\z/ } @files) {
	try sub { load_runfile $_ }, problem_with $_;
	$processed{$_} = 1;
}

# The scrape_...() functions also need the .stderr files for CPLEX and Gurobi runs to be loaded into memory first...  What a hack.
# At least all 3 solver types (orig_cplex, cut_cplex and grb) store this info in the same way in files with the same naming system.
foreach (grep { /\.(?:[^\.]+_cplex|grb)\.stderr\z/ } @files) {
	try sub { load_stderrfile $_ }, problem_with $_;
	$processed{$_} = 1;
}

# CPLEX writes important info, like whether or not it timed out, to the .sol file.  Need to read this before the *_cplex.stdout files.
foreach (grep { /\.[^\.]+_cplex\.sol\z/ } @files) {
	try sub { load_cplexsolfile $_ }, problem_with $_;
	$processed{$_} = 1;
}

# Extract first line and MD5 hash of rest from .canon files
foreach (grep { /\.canon\z/ } @files) {
	try sub { load_canonfile $_ }, problem_with $_;
	$processed{$_} = 1;
}

# We now allow ".reduced" to be followed by one or more digits encoding the particular reduction script.
foreach (grep { !/\.reduced\d*\./ } grep { /\.txt\z/ } @files) {
	try sub { insert_into_temp "instances", scrape_instance $_; }, problem_with $_;
	$processed{$_} = 1;
}

foreach (grep { /\.reduced\d*\./ } grep { /\.txt\z/ } @files) {
	try sub { insert_into_temp "instances", scrape_instance $_; }, problem_with $_;
	$processed{$_} = 1;
}

foreach (grep { /\.ftr\.log\z/ } @files) {
	try sub { insert_into_temp "reductions", scrape_ftr_log $_; }, problem_with $_;
	$processed{$_} = 1;
}

foreach (grep { /\.[^\.]+_cplex\.stdout\z/ } @files) {
	try sub { insert_into_temp "solutions", scrape_stephan_cplex_stdout $_; }, problem_with $_;
	$processed{$_} = 1;
}

foreach (grep { /\.grb\.log\z/ } @files) {
	try sub { insert_into_temp "solutions", scrape_grb_log $_; }, problem_with $_;
	$processed{$_} = 1;
}

# Update the main DB tables from the temp ones
update_all_tables_from_temp();

# Did we miss any?
my $nUnprocessedFiles = 0;
foreach (@files) {
	if (!exists $processed{$_}) {
		if ($verbose) {
			warn "Did not process unrecognised file '$_'!\n" ;
		} else {
			++$nUnprocessedFiles;
		}
	}
}

print STDERR "$nUnprocessedFiles unprocessed files.\n";

# Were there problems with any?
foreach (@filesWithProblems) {
	warn "Problem processing file $_\n";
}
