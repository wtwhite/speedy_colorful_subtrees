#!/usr/bin/perl

use strict;
use warnings;
use autodie ':all';
use DBI;

#die "This is designed to run only on Unix." if exists $ENV{OS} && $ENV{OS} =~ /\AWin/;		#HACK

# Doesn't work because we still need the mpos* stuff, e.g. to distinguish COMBINED instances.
#		'fmm1/' || substring(dataset from '([^/]+)\$'),
my $newDatasetExpr = "'fmm1/' || substring(regexp_replace(src.dataset, '/b\\.[0-9a-f][0-9a-f]/', '/') from '/(.+)\$')";		#HACK: Horrifying.  Relies on the table always being called "src"...

chomp(my $thisHost = `hostname`);

# Open DB connection and set up its closing
my $dbPort = ($thisHost eq 'wallace' ? 5432 : 15433);			# On wallace we access the DB directly; on all other machines we need an SSH tunnel.
#my $dbh = DBI->connect("dbi:Pg:dbname=wtwhite;host=localhost;port=$dbPort", 'wtwhite', 'asdfasdf', { RaiseError => 1, ShowErrorStatement => 1, TraceLevel => 'SQL' });		#HACK: Disgusting that I include the password...
my $dbh = DBI->connect("dbi:Pg:dbname=wtwhite;host=localhost;port=$dbPort", 'wtwhite', 'asdfasdf', { RaiseError => 1, ShowErrorStatement => 1 });		#HACK: Disgusting that I include the password...
END {
	if (defined $dbh) {
		$dbh->disconnect;
		$dbh = undef;
	}
}

sub update_instances($) {
	my ($src) = @_;
	
	print "update_instances($src): Running SQL...\n";
	my $sql = <<THE_END;
insert into instances (
	dataset,
	reduction,
	data_group,
	is_combined,
	path,
	nverts,
	nedges,
	ncols,
	opt_sol_val,
	opt_sol_extra_val
) (select
		$newDatasetExpr,
		reduction,
		data_group,
		is_combined,
		path,
		nverts,
		nedges,
		ncols,
		opt_sol_val,
		opt_sol_extra_val
	from instances src
	where src.dataset like '$src/%'
	and not exists (select 1
		from instances already
		where already.dataset = $newDatasetExpr
		and already.reduction = src.reduction
	)
)
THE_END

	my $nRows = 0 + $dbh->do($sql);
	print "update_instances($src): $nRows rows inserted.\n";
}

sub update_reductions($) {
	my ($src) = @_;
	
	print "update_reductions($src): Running SQL...\n";
	my $sql = <<THE_END;
insert into reductions (
	dataset,
	reduction,
	input_path,
	output_path,
	job_id,
	command_line,
	hostname,
	completed_time,
	elapsed_secs,
	elapsed_nonio_secs,
	user_secs,
	system_secs,
	percent_cpu,
	max_mem_kb
) (select
		$newDatasetExpr,
		reduction,
		input_path,
		output_path,
		job_id,
		command_line,
		hostname,
		completed_time,
		elapsed_secs,
		elapsed_nonio_secs,
		user_secs,
		system_secs,
		percent_cpu,
		max_mem_kb
	from reductions src
	where src.dataset like '$src/%'
	and not exists (select 1
		from reductions already
		where already.dataset = $newDatasetExpr
		and already.reduction = src.reduction
	)
)
THE_END

	my $nRows = 0 + $dbh->do($sql);
	print "update_reductions($src): $nRows rows inserted.\n";
}

sub update_solutions($) {
	my ($src) = @_;
	
	print "update_solutions($src): Running SQL...\n";
	my $sql = <<THE_END;
insert into solutions (
	dataset,
	reduction,
	solver,
	input_path,
	output_path,
	job_id,
	command_line,
	hostname,
	completed_time,
	nthreads,
	result_type,
	orig_result_type,
	opt_sol_val,
	opt_bound,
	elapsed_secs,
	elapsed_nonio_secs,
	ncuts_generated,
	user_secs,
	system_secs,
	percent_cpu,
	max_mem_kb,
	time_limit_secs,
	mem_limit_kb,
	canon_first_line,
	canon_tail_md5
) (select
		$newDatasetExpr,
		reduction,
		solver,
		input_path,
		output_path,
		job_id,
		command_line,
		hostname,
		completed_time,
		nthreads,
		result_type,
		orig_result_type,
		opt_sol_val,
		opt_bound,
		elapsed_secs,
		elapsed_nonio_secs,
		ncuts_generated,
		user_secs,
		system_secs,
		percent_cpu,
		max_mem_kb,
		time_limit_secs,
		mem_limit_kb,
		canon_first_line,
		canon_tail_md5
	from solutions src
	where src.dataset like '$src/%'
	and not exists (select 1
		from solutions already
		where already.dataset = $newDatasetExpr
		and already.reduction = src.reduction
		and already.solver = src.solver
	)
)
THE_END

	my $nRows = 0 + $dbh->do($sql);
	print "update_solutions($src): $nRows rows inserted.\n";
}

sub update($) {
	my ($src) = @_;
	update_instances $src;
	update_reductions $src;
	update_solutions $src;
}

# Main program
update 'missing1';
update 'missing';
update 'fragmentationgraphs';
