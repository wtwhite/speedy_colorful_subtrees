#!/usr/bin/perl

use strict;
use warnings;
use autodie ':all';
use DBI;

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

my @cols = qw(
	.orig_cplex
	.reduced4.orig_cplex
	.reduced6.orig_cplex
	.cut_cplex
	.reduced4.cut_cplex
	.reduced6.cut_cplex
);

print STDERR "Running query...\n";
my @rows = @{$dbh->selectall_arrayref(<<THE_END, { Slice => {} })};
	select
		sq1.data_group,
		sq1.solver,
		sq1.reduction,
		comb_red_time + comb_sol_time as comb_time,
		sep_red_time + sep_sol_time as sep_time,
		ncands
	from (
		select
			data_group,
			solver,
			reduction,
			coalesce(r.elapsed_nonio_secs, 0) as comb_red_time,
			s.elapsed_nonio_secs as comb_sol_time
		from nice_instances i
		left join ok_reductions r using (dataset, reduction)
		join ok_solutions s using (dataset, reduction)
		where is_combined
	) sq1
	join (
		select
			data_group,
			solver,
			reduction,
			sum(coalesce(r.elapsed_nonio_secs, 0)) as sep_red_time,
			sum(s.elapsed_nonio_secs) as sep_sol_time,
			count(*) as ncands
		from nice_instances i
		left join ok_reductions r using (dataset, reduction)
		join ok_solutions s using (dataset, reduction)
		where not is_combined
		group by 1, 2, 3
	) sq2 on sq1.data_group = sq2.data_group and sq1.solver = sq2.solver and sq1.reduction = sq2.reduction
	order by 1, 2, 3
THE_END
print STDERR "Query produced " . scalar(@rows) . " rows.\n";

my %d;
foreach my $r (@rows) {
	#$d{$r->{data_group}}{"$r->{reduction}.$r->{solver}.sep"} = $r->{sep_time};
	#$d{$r->{data_group}}{"$r->{reduction}.$r->{solver}.comb"} = $r->{comb_time};
	$d{$r->{data_group}}{"$r->{reduction}.$r->{solver}.sep"} = sprintf '%.0f', $r->{sep_time};
	$d{$r->{data_group}}{"$r->{reduction}.$r->{solver}.comb"} = sprintf '%.0f', $r->{comb_time};
	$d{$r->{data_group}}{ncands} = sprintf '%.0f', $r->{ncands};
}

#print "\\begin{tabular}{|c" . ("|d{-1}" x 6) . "|}\\hline\n";
#print "\\begin{tabular}{|c" . ("|r" x 6) . "|}\\hline\n";		# "d{-1}" doesn't work if we try to bold with \textbf{}...
#print "\\begin{tabular}{c" . ("r" x 6) . "}\\hline\n";		# "d{-1}" doesn't work if we try to bold with \textbf{}...
print <<'THE_END'
\begin{tabular}{c r|rrrrrr|rrrrrr}
& & \multicolumn{6}{c|}{solve all candidate instances} &
  \multicolumn{6}{c}{solve global instances} \\
%
\emph{ID} & \#cand	&	\multicolumn{1}{c}{o}	&	oR2	&	oR3	&
\multicolumn{1}{c}{c}	&	cR2	&	cR3
&	\multicolumn{1}{c}{o}	&	oR2	&	oR3	&
\multicolumn{1}{c}{c}	&	cR2	&	cR3
\\\hline
THE_END

#print join("\t&\t", 'ID/\#cand', 'o\qquad', 'oR2', 'oR3', 'c\qquad', 'cR2', 'cR3'), "\\\\\\hline\n";
#print join("\t&\t", '\emph{ID}/(\#cand)', '\mbox{o\,\,\,\,}', 'oR2', 'oR3', '\mbox{c\,\,\,\,}', 'cR2', 'cR3'), "\\\\\\hline\n";		#HACK: the \,\,\,\, is disgusting

foreach my $dn (sort keys %d) {
	my $short = $dn;
	$short =~ s/\.ms$//;
	$short =~ s/^mpos/m/;		# Need all the horizontal space we can get!
	#print "$short";
	print "\\emph{$short}";
	print "\t&\t$d{$dn}{ncands}";		# Sneak the number of candidates in there!  :)
	
	# Find maximum over instance.
	#print STDERR join("|", @{$d{$dn}}{map { ("$_.sep", "$_.comb") } @cols}), "\n";		#DEBUG
	my $min = (sort { $a <=> $b } @{$d{$dn}}{map { ("$_.sep", "$_.comb") } @cols})[0];
	#print STDERR "min=<$min>\n";		#DEBUG
	
	foreach my $c (@cols) {
		my $key = "$c.sep";
		my $val = $d{$dn}{$key};
		$val = "\\textbf{$val}" if $val == $min;
		print "\t&\t$val";
	}
	
	#print "\\\\\n($d{$dn}{ncands})";		# Sneak the number of candidates in there!  :)
	foreach my $c (@cols) {
		my $key = "$c.comb";
		my $val = $d{$dn}{$key};
		$val = "\\textbf{$val}" if $val == $min;
		print "\t&\t$val";
	}
	#print "\\\\\\hline\n";
	print "\\\\\n";
}

print "\\hline\n";
print "\\end{tabular}\n";
