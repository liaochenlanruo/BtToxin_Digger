#!/usr/bin/env perl
use strict;
use warnings;


my %hash;
my %strain;
my %gene;
#my %lensvm;

my @files = glob("*.list");
foreach  (@files) {
	$_=~/(\S+).list/;
	$strain{$1}++;
	my $str = $1;
	my $toxin_type;
	open IN, "$_" || die;
	while (<IN>) {
		chomp;
		if (/Toxin\s+type:\s+(\S+)\s+protein/) {
			$toxin_type = $1;
		}
		if (/^[0-9]+/) {
			my @lines = split /\t/;
			#my $gene = $lines[6];
			#my $coverage = $lines[8];
			#my $identity = $lines[9];
			if ($lines[5] eq "YES") {#判断blast是否有结果
				if ($lines[8] >= 50) {#coverage 低于阈值的不考虑
					$gene{$lines[6]}++;
					#$hash{$str}{$lines[6]} = $lines[9];#若存在重复基因?
					$hash{$str}{$lines[6]} .= $lines[9] . ",";#若存在重复基因?2020/4/12
				}
			}elsif ($lines[-1] eq "YES") {
				#$lensvm{$lines[3]}++;
				my $predictor = "HMM_" . $toxin_type . "_len_" . $lines[3];
				$gene{$predictor}++;
				$hash{$str}{$predictor}++;
			}elsif ($lines[-2] eq "YES") {
				my $predictor = "SVM_" . $toxin_type . "_len_" . $lines[3];
				$gene{$predictor}++;
				$hash{$str}{$predictor}++;
			}
		}
	}
}

open OUT, ">Bt_all_genes.table" or die "open OUT file faild\n";
my @strains = sort keys %strain;
my @genes = sort keys %gene;
print OUT "-";
foreach  (@genes) {
	print OUT "\t" . $_;
}

print OUT "\n";
for (my $i = 0;$i < @strains;$i++) {
	print OUT $strains[$i];
	for (my $j = 0;$j < @genes;$j++) {
		if (exists $hash{$strains[$i]}{$genes[$j]}) {#2020/4/12
			if ($hash{$strains[$i]}{$genes[$j]}=~/(.+),/) {#2020/4/12
				print OUT "\t" . $1;#2020/4/12
			}else {#2020/4/12
				print OUT "\t" . "$hash{$strains[$i]}{$genes[$j]}";#2020/4/12
			}#2020/4/12
			#print OUT "\t" . "$hash{$strains[$i]}{$genes[$j]}";
		}else{
			print OUT "\t";
		}
	}
	print OUT "\n";
}
close IN;
close OUT;
