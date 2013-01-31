#!/usr/bin/perl
#Written by Ying Wu
#06/01/2012
#generate coverage bed for DIME input

use strict;
use warnings;

die "Usage: perl covBEDforDIME.pl chrsizes.txt > 500bp_windows.bed
" unless($ARGV[0]);

open(FILE, $ARGV[0]) || die("Could not open chrsizes file \"$ARGV[0]\": $!");

my @chr = (1 .. 22);
push(@chr, "M", "X", "Y");
my %chrsize;
my $winsize = 500;

while(<FILE>)
{
	chomp;
	$chrsize{+(split("\t", $_))[0]} = +(split("\t", $_))[1];
	#print "KEY: ", +(split("\t", $_))[0], "\t VAL: ", +(split("\t", $_))[1],"\n";
}
#for my $key ( sort keys %chrsize ) { print "KEY: $key \t VAL: $chrsize{$key}\n"; }

for my $chrnum (@chr) {
	my $i = 0;
	while(exists $chrsize{"chr$chrnum"} && ($i+1)*$winsize < $chrsize{"chr$chrnum"})
	{ #not sure if +1 for winsize should be there
		print "chr$chrnum\t",$i*$winsize+1,"\t",($i+1)*$winsize,"\n";
		$i++;
	} #not quite sure how coveragedbed handles out of bounds so do things this way
	if(exists $chrsize{"chr$chrnum"}) { print "chr$chrnum\t",$i*$winsize+1,"\t",$chrsize{"chr$chrnum"},"\n"; }
}
