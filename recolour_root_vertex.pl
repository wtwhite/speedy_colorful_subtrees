#!/usr/bin/perl

# Not strictly necessary but without this the size can blow up quite a bit!
binmode STDIN;
binmode STDOUT;

print scalar <>;	# nVerts
print scalar <>;	# nEdges
my $nCols = 0 + <>;
print $nCols + 1, "\n";
$_ = <>;		# Could be known optimal objective function score, or first line of colours
if (/\A\s*\S+\s*\z/) {
	print;
	$_ = <>;	# It was just the score; read the next line
}
#if (/\A\s*0\s+0\s*\z/) {
if (/\A\s*0\s+\d+\s*\z/) {		# Recognise a root of any colour
	print "0 $nCols\n";	# Line recognised; print the new colour
} else {
	die "Did not recognise line $. -- expected to see '0 0', representing the root vertex coloured 0.";
}

while (<>) {
	print;
}
