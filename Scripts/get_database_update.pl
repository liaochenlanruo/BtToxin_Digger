#!/usr/bin/perl
use strict;
use warnings;

my %hash;
local $/ = ">";
open IN, "bt_toxin20220609.fas" || die;# BtToxin_Digger database
<IN>;
while (<IN>) {
	chomp;
	my ($id, $seq) = split (/\n/, $_, 2);
	$hash{$id}++;
}
close IN;


open INF, "all_app_cry_cyt_gpp_mcf_mpf_mpp_mtx_pra_prb_spp_tpp_vip_vpa_vpb_xpp_fasta_sequences.txt" || die;# The latest BPPRC database
open OUT, ">Update.txt" || die;# protein sequences in BPPRC database but not exist in BtToxin_Digger database

<INF>;
while (<INF>) {
	chomp;
	my ($id, $seq) = split (/\n/, $_, 2);
	if (not exists $hash{$id}) {
		print OUT ">$id\n$seq" . "\n";
	}
}

close INF;
close OUT;
