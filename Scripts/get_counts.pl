#!/usr/bin/perl
use strict;
use warnings;

#my $file = "counttimes.txt";
#my $times = &counttime($file);
#print $times . "\n";
#sub counttime{
#	my $infile = @_;
my $counttimes;
open IN, "counttimes.txt" || die;
open OUT, ">TEMPTIMES" || die;
while (my $line = <IN>) {
	$line=~/(\d+)/;
	$counttimes = $1 + 1;
	print OUT $counttimes;
}
close IN;
close OUT;
rename("TEMPTIMES", "counttimes.txt");
	#unlink("TEMP");
#	return $count;
#}
