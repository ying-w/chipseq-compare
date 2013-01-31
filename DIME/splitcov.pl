#!/usr/bin/perl
#Written by Ying Wu
#06/01/2012
#split coverage files into each chromosome
#the inefficient way, better way would be count from top the lines needed using
#head -498502 _.cov or head -(498502+486399) _.cov | tail +498502 w/cumulative sum

use strict;
use warnings;

my @cov = `ls ./data/*cov`;
my $winsize = 500;
my @chr = (1 .. 22);
push(@chr, "M", "X", "Y");

foreach my $covfile (@cov)
{
	chomp $covfile;
	my $prefix = "tmp";
	if($covfile =~ /\/(\w+)\.cov/) { $prefix = $1; }
	foreach my $chrnum (@chr)
	{
		my $filename = join("_", $prefix, $winsize, "chr$chrnum");
		print "grep \"chr$chrnum\t\" $covfile | cut -f 4 > ./data/$prefix/$filename\n";
		print `grep "chr$chrnum	" $covfile | cut -f 4 > ./data/$prefix/$filename`;
	#keep in mind that there is a tab after chrnum
	}
}
