#! /usr/bin/perl

use warnings;
use strict;

if (@ARGV != 2){
	die("Syntax: $0 [coverage file] [window size]\n");
}

my $hash;

open(IN, $ARGV[0]);
while (<IN>){
	chomp;
	my @data = split /\t/;
	push(@{$hash->{$data[0]}}, $data[2]);
}
close(IN);

my $out = $ARGV[0] . "_in_" . $ARGV[1] . "_windows.txt";
open(OUT, ">$out");

foreach my $i (keys %{$hash}){
	my @data = @{$hash->{$i}};
	my $p = 1;
	for (my $a = 0; $a < @data-$ARGV[1]; $a++){
		my $t = 0;
		#print "Adding:";
		for (my $b = $a; $b < $a+$ARGV[1]; $b++){
		#	print " " . $data[$b];
			$t += $data[$b];
		}
		#print "\n";
		print OUT $i . "\t" . $p . "\t" .  ($t/$ARGV[1]) . "\n";
		$p++;
	}
}
close(OUT);

exit;
		