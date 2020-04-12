#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;

use lib dirname(__FILE__) . "/lib";

use own;
use File::Spec;

system("mkdir -p Results/Toxins");
system("mkdir -p tmp");
my $input = shift;
my $output = shift;
my $seq_type = shift;
my @results;
my $download_dir = "./Results/Toxins/";

my $target_file = File::Spec->catfile("./tmp/", $output);
if(!system("cp $input $target_file"))  {
	@results = own::BTTCMP($output, $seq_type, $target_file);
	foreach (@results)  {
		system("mv $_ $download_dir");
	}
	system("rm -f tmp/*");
}
