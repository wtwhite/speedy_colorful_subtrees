#!/usr/bin/perl

use strict;
use warnings;

$_ = join '', <>;
warn "Your SQL contains '--'...  If that's for a comment, then this will go horribly wrong!" if /--/;			#HACK: "--" comments cause havoc if we join everything into 1 line, but it's too hard to parse them out properly...
tr/\n\r\t/ /s;		# psql's backslash-commands can't handle multiline arguments (and tabs and carriage returns look ugly)
s/;\s*\z//;			# Can't have a ';' inside a backslash-command
s/'/'"'"'/g;		# Quote single quotes

# GLE interprets spaces as starting a new field unless the field is quoted.
# PostgreSQL's COPY TO doesn't quote fields with spaces unless you specify FORCE_QUOTE, in which case it quotes all fields
# (except NULL values), including numeric fields, which GLE will then try to interpret as strings...
#HACK: So we have to scrape COPY TO's output to strip quotes off anything that could be a number...
my $cmd = "psql -c '\\copy ($_) to stdout null as '\"'*'\"' csv header force quote *'";
#my $cmd = "psql -c '\\copy ($_) to stdout null as '\"'*'\"' csv header'";
print STDERR "$cmd\n";
#system $cmd and die "Failed with exit code $?!";
open my $f, "-|", $cmd or die;
while (<$f>) {
	chomp;
	my @f;
	while (s/\A("(?:[^"]*"")*[^"]*"|[^",]*)//) {
		my $x = $1;
		$x =~ s/\A"([^\s,"]*)"\z/$1/;			# Unquote if it's safe to do so
		if ($x =~ /\W/) {		# ... and quote if it's not safe not to do so.  This happens for header rows, which \\copy to doesn't quote.
			$x =~ s/"/""/g;
			$x = "\"$x\"";
		}
		push @f, $x;
		s/\A,// or last;
	}
	
	print join(",", @f), "\n";
}
