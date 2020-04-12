#!/usr/bin/env perl
use strict;
use warnings;

# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Pipeline for identification of Bacillus thuringiensis toxins from massive genomic data

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 FEEDBACK

=head2 Mailing Lists

=head2 Reporting Bugs

=head1 AUTHOR - Weixing YE & Hualin Liu

E-mail: yeweixing@yahoo.cn or liaochenlanruo@webmail.hzau.edu.cn

=head1 CONTRIBUTORS

=cut

#Let the code begin...

package own;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Bio::SearchIO;
use Bio::Seq;
use Bio::Tools::CodonTable;
use Cwd;
use Bio::SeqFeature::Generic;
use Bio::Tools::SeqStats;
use File::Basename;
use File::Spec;


=head2  BTTCMP
 Title   :
 Usage   :
 Function:
 Return  :
 Args    :
=cut

sub BTTCMP  {
	my ($output_filename, $seqtype, @seqfile) = @_;
	my @preprocess_out = Pre_Process($seqtype, @seqfile);
	my %step1_out = Step1($output_filename, @preprocess_out);
	my %step2_out = Step2($output_filename, %step1_out);
	my %domain = Step2x($output_filename);#2020/4/7
	#my @results = Step3($seqtype, $output_filename,\@seqfile, %step2_out);
	my @results = Step3($seqtype, $output_filename,\@seqfile, %step2_out, %domain);#2020/4/7
	return @results;
}


=head2  Pre_Process
 Title   : Pre_Process
 Usage   : @seq = &Pre_process($seq_type, @input);
 Function: based on different sequence input, the pre-process would handle everything
           so that the input would be suitable for the next step analysis
 Return  : An array contains preprocess file names and the program would be used
 Args    : An hash list contains all needed parameters 
=cut

sub Pre_Process  {
	my ($seqtype, @seqfile) = @_; 
	# For $seq_quality, 0 means bad sequence quality containing thousands of small fragments
	# while 1 represents good sequence quality containing sequences been assembled well
	my @preprocess_out;
	if($seqtype eq 'prot')   {
		@preprocess_out = &Length_Filter(@seqfile);
	}elsif($seqtype eq 'orfs')  {
		@preprocess_out = &Batch_translation(@seqfile);
	}elsif($seqtype eq 'nucl')  { 
		@preprocess_out = &six_translation(@seqfile);
	}else  {
		print "Error!! Wrong sequence type!! \n";
		exit;
	}
	return(@preprocess_out);
}


=head2 Step1 
 Title   : Step1 
 Usage   : %step1_out = &Step1($system, $param, @step1_in);
 Function: Blast against Bt toxin, and screen the proteins showing significant hits to Bt toxins
 Return  : An hash list contains the query id and blast results
 Args    : system, parameters for Blast analysis and pre-process inputs 
=cut

sub Step1  {
	my ($output_filename, @files) = @_;
	my %step1_out;
	my $dir;
	my $path = `which bttcmp`;
	if ($path=~/(.+)\/bttcmp/) {
		$dir = $1;
	}
	my @hmm_models = (
						$dir . "/BTTCMP_models/hmm_models/cry.hmm",
						$dir . "/BTTCMP_models/hmm_models/cyt.hmm",
						$dir . "/BTTCMP_models/hmm_models/vip.hmm",
					);

	#Combine all sequence files into seperate file
	my $step1_1_file = $output_filename . ".step1_in";
	my $seqout = Bio::SeqIO->new(-file => ">$step1_1_file", -format => "fasta");
	#file test to determine if the file size equal to zero
	foreach my $file (@files)  {
		my $in = Bio::SeqIO->new(-file => $file, -format => "fasta");
		while(my $seq = $in->next_seq())  {
			if($seq->id =~ /.*\|(.*)\|/)  {
				$seq->id("$1");
			}
			$seqout->write_seq($seq);
		}
	}

	#After that, delete sequence files
	system("rm @files");

	#Determine the file size
	my $size = -s "$step1_1_file";

	#If size is greater than 0, Implement BLAST, SVM and HMM prediction
	if($size)  {
		#print "Step1: Start Blast search.\n";
		my $Blast_file_out = &Blast_search('step1', $step1_1_file);

		#print "Step1: Blast completed.\n";
		#print "Step1: Now start Blast parsing.\n";
		my %blast_results = &Blast_parser($Blast_file_out);

		#Implement SVM prediction
		my %svm_results = &svm_prediction($step1_1_file);

		#Implement HMM prediction
		my %hmm_results = &hmm_prediction($step1_1_file, @hmm_models);

		%step1_out = (
						'blast' => \%blast_results,
						'svm'   => \%svm_results,
						'hmm'   => \%hmm_results,
					);


		#Extract IDs and sequences 
		my $db = &Index_maker(&Guess_header_format($step1_1_file), $step1_1_file);
		my $step1_out = $output_filename . ".step1";
		my $out = Bio::SeqIO->new(-file => ">$step1_out", -format => "fasta");
		my $idex1 = $step1_1_file . ".index";#modified
		#my $idex2 = $step1_1_file . ".index";#modified
		my %step1_ids;
		foreach  (sort keys %step1_out)  {
			my $temp = $step1_out{$_};
			foreach (sort keys %$temp)  {
				if(!(exists $step1_ids{$_}))  {
					$step1_ids{$_} = 1;
					my $seq = $db->get_Seq_by_acc($_);
					$out->write_seq($seq);
				}
			}
		}

		system("rm $step1_1_file $idex1");#modified
		return %step1_out;
	}else  {
		system("rm $step1_1_file @files");
		printnoresults($output_filename);
		#clear the temp files
	}
}


=head2 Step2
 Title   : Background elimination
 Usage   : 
 Function:
 Return  : An hash list contains the query id and the results received after background elimination
 Args    : 
=cut

sub Step2  {
	my ($output_filename, %step2_in) = @_;
	my $step1_in = $output_filename . ".step1";
	my $step2_out = $output_filename . ".step2";

	#If size is greater than 0, Implement BLAST, SVM and HMM prediction
	if(-s "$step1_in")  {
		my $Blast_result = &Blast_search('step2', $step1_in);
		my %step2_out = &Blast_parser1($Blast_result);
		my $seq_number = 0;
		my %seq_ids;
		foreach my $key(sort keys %step2_in)  {
			my $temp = $step2_in{$key};
			foreach (sort keys %$temp)  {
				if((exists $step2_out{$_}) && (!($$temp{$_} =~ /vip\.hmm|cyt\.hmm|cry\.hmm/i)))  {
					delete $$temp{$_};
				}else  {
					if(!(exists $seq_ids{$_}))  {
						$seq_ids{$_} = '';
						$seq_number++;
					}
				}
			}
		}
		if($seq_number == 0)  {
				#Clear the temp files
				system("rm $step1_in");
				printnoresults($output_filename);
		}else  {
			#Extract Sequence
			my $db = &Index_maker('normal', $step1_in);
			my $out = Bio::SeqIO->new(-file => ">$step2_out", -format => "fasta");
			foreach my $id(sort keys %seq_ids)  {
				my $seq = $db->get_Seq_by_id($id);
				$out->write_seq($seq);
			}
			my $index1 = $step1_in . ".index";#modified
			#my $index2 = $step1_in . ".index";#modified

			#Clear the temp file
			system("rm $index1 $step1_in");#modified
			return %step2_in;
		}

	}else  {
		#clear the temp files
		system("rm $step1_in");
		printnoresults($output_filename);
	}
}

sub hmm_prediction2  {
	my ($seq_file, @domain_models) = @_;
#	my %HMM_results;
#	my (%D1, %D2, %D3, %MTX2, %Toxin10) = ();
	my %d1=();
	my %d2=();
	my %d25=();
	my %d3=();
	my %mtx2=();
	my %toxin10=();
	my %domain1;#2020/4/8
	# Implement HMM searching and parsing
	foreach my $hmm(@domain_models)  {
		#For each hmm models, implement hmmsearch
		$hmm=~/.+\/(.+)/;
		my $hmm_file_out = $seq_file . ".$1";
		system("hmmsearch -E 1e-10 -o $hmm_file_out $hmm $seq_file");

		#HMM results parsing
		my $searchio = Bio::SearchIO->new(-file => $hmm_file_out, -format => "hmmer3");
		while(my $result = $searchio->next_result)  {
			if($result->num_hits)  {
				while(my $hit = $result->next_hit)  {
					my $hmm_name = $result->query_name;
					my $hmm_len = $result->query_length;
					my $hit_desc = $hit->description;
					my $hsp = $hit->next_hsp;
					my $evalue = $hsp->evalue;
					my $id = $hit->name;
					#print $hmm_name . "\n";
					#while($hsp = $hit->next_hsp)  {
					#	if($evalue > $hsp->evalue)  {
					#		$evalue = $hsp->evalue;
					#	}
					#}
					#print '$hmm_name: ' . $hmm_name . "\n";
					#print "sequence ID: " . $id . "\n";
					if ($hmm_name=~/EndotoxinN/i) {
						$d1{$id}++;
						#print $hmm_name . "\t". $id . "\n";
					}
					if ($hmm_name=~/EndotoxinM/i) {
						$d2{$id}++;
						#print $hmm_name . "\t". $id . "\n";
					}
					if ($hmm_name=~/EndotoxinC/i) {
						$d3{$id}++;
						#print $hmm_name . "\t". $id . "\n";
					}
					if ($hmm_name=~/Endotoxin_mid/i) {
						$d25{$id}++;
						#print $hmm_name . "\t". $id . "\n";
					}
					if ($hmm_name=~/ETXMTX2/i) {
						$mtx2{$id}++;
						#print $hmm_name . "\t". $id . "\n";
					}
					if ($hmm_name=~/Toxin10/i) {
						$toxin10{$id}++;
						#print $hmm_name . "\t". $id . "\n";
					}
					#$$hmm_name{$id}++;
				}
			}
		}
		system("rm $hmm_file_out");
	}
	%domain1 = (
				'en'   => \%d1,
				'em'   => \%d2,
				'ec'   => \%d3,
				'emid' => \%d25,
				'mtx'  => \%mtx2,
				'tox'  => \%toxin10,
	);
	return %domain1;
}

sub Step2x  {
	###my ($output_filename, %step2x_in) = @_;
	my ($output_filename) = @_;
	#my %step2x_out;
	my %domain;#2020/4/8
	my $dir;
	my $path = `which bttcmp`;
	if ($path=~/(.+)\/bttcmp/) {
		$dir = $1;
	}
	my @domain_models = (
						$dir . "/BTTCMP_models/cry_domains/EndotoxinN.hmm",
						$dir . "/BTTCMP_models/cry_domains/EndotoxinM.hmm",
						$dir . "/BTTCMP_models/cry_domains/EndotoxinC.hmm",
						$dir . "/BTTCMP_models/cry_domains/Endotoxin_mid.hmm",
						$dir . "/BTTCMP_models/cry_domains/ETXMTX2.hmm",
						$dir . "/BTTCMP_models/cry_domains/Toxin10.hmm",
					);
	my $input = $output_filename . ".step2";
	if (-s "$input") {
		my %domain = &hmm_prediction2($input, @domain_models);
		return %domain;#2020/4/8
	}else  {
		printnoresults($output_filename);
	}
}

=head2 Step3
 Title   : Step3
 Usage   : 
 Function: 
 Return  : none
 Args    : An hash list contains the informations about all of the 
=cut

sub Step3  {
	my ($seq_type, $output_file, $seqfile, %step3_in, %domain) = @_;
	my @results = Writer($seq_type, $output_file, $seqfile, %step3_in, %domain);
	return @results;
}

=head2 Length_Filter
 Title   : Length_Filter
 Usage   : @prot = &Length_Filter(@inputs)
 Function: Retain those proteins with sequences length greater than 115 amino acid residues
 Returns : An array of files contains multiple fasta format of protein sequences
 Args    : An array of files contains multiple fasta format of protein sequences
=cut

sub Length_Filter  {
	my @inputs = @_;
	my @outputs;

	#Set a sequence count
	my $seq_count = 0;

	foreach my $inputs (@inputs)  {
		my $output = $inputs . ".prot";
		my $in = Bio::SeqIO->new(-file => $inputs, -format => "fasta");
 		my $out = Bio::SeqIO->new(-file => ">$output", -format => "fasta");

		#Pick up those protein sequence with sequence length greater than 115 amino acid residues
		while(my $seq = $in->next_seq)  {
			$seq_count++;
			if(($seq->length >= 115) && !($seq->seq =~ /[^acdefghiklmnpqrstvwy]/i))  {
				$out->write_seq($seq);
			}
		}
		push @outputs, $output;
	}
	if($seq_count == 0)  {
		#Delete both the input file and the output file
		print "Empty Sequence file.\n";
		system ("rm @outputs");
		exit;
	}else  {
		return @outputs;
	}
}


=head2 Batch_translation
 Title   : Batch_translation
 Usage   : @prot_file = &prot_translate(@orf_file);
 Function: Translate open reading frames to proteins batchly
 Returns : An array of files contains multiple fasta format of proteins
 Args    : An array of files contains multiple fasta format of open reading frames 
=cut

sub Batch_translation  {

	my (@inputs) = @_;
	# Specify the eubacterial codon table during translation
	my $CodonTable  = Bio::Tools::CodonTable -> new ( -id => 11 );  
	my @outputs;

	my $seq_count = 0;
	foreach (@inputs)  {
		my $orf_fin = $_;
		my $prot_fout = $_ . ".prot";
		push @outputs, $prot_fout;
		my $orf_in = Bio::SeqIO->new(-file => "$orf_fin", -format => "fasta");
		my $prot_out = Bio::SeqIO->new(-file => ">>$prot_fout", -format => "fasta");

		while(my $seq = $orf_in->next_seq())  {
		$seq_count++;
		#Make sure that the sequences length is between 115 amino acid residues
		#Make sure that the sequences is a open reading frames, else return error
		if(($seq->length >= 348) && !($seq->seq =~ /[^atcg]/i))  {
			my $prot = $seq->translate(-complete => 1, -throw => 0, -codontable => $CodonTable);
			$prot_out->write_seq($prot);
		}
		}
	}

	if($seq_count == 0)  {
		#Delete both the input file and the output file
		print "Empty Sequence file.\n";
		system ("rm @outputs");
		exit;
	}else  {
		# Return translation results
		return(@outputs);
	}
}


=head2 Six_translation
 Title   : Six_translation
 Usage   : @Prot  = &six_translation(@input)
 Function: Implement six frame translation
 Return  : An array of files contains multiple fasta format of protein sequences
 Args    : An array of files contains multiple fasta format of normal nucleotide sequences
=cut

sub six_translation  {
	my @inputs = @_;
	my @outputs;
	my $seq_count = 0;

	foreach my $input(@inputs)  {
		my $output = $input . ".prot";
		push @outputs, $output;        

		my $in = Bio::SeqIO->new(-file => $input, -format => "fasta");
		my $out = Bio::SeqIO->new(-file => ">$output", -format => "fasta");

		#maximum orf number is up to 9999
		my $orf_width = 5;
		my @seqs;

		while(my $seq = $in->next_seq())  {
			$seq_count++;
			my $seq_id = $seq->id;
			my $seqstr = $seq->seq;
			$seqstr =~ s/[^atcg]//ig;
			my $desc = $seq->description;
			$seq = Bio::Seq->new(-id => $seq_id, -seq => $seqstr, -desc => $desc);
			my $seq_len = $seq->length;

			if(($seq->length >= 345) && !($seq->seq =~ /[^atcg]/i))  {
				my $frame_number = 0;
				my $id = '';
				my $pre_length = 0;
				my $frame_length = 0;

				my @strands = (-1, 1);
				my @frame_flags = (0, 1, 2);

				foreach my $strand (@strands)  {
					foreach my $frame_flag(@frame_flags)  {
						if($strand == 1)  {
							my $str = $seq->translate(-frame => $frame_flag);
							my $flag = $frame_flag + 1;
							my $subseq = $str->seq;
							my $start = 0;
							my $end = 0;
							my $repeat = 0;
							my $desc = '';
							my $frame_seq= '';
							my $frame_tag = '';
							my $total_len = 0;
							while($subseq =~ /(\w+)/i)  {
								$frame_length = length($&);
								$pre_length = length($`);
								$total_len += $pre_length; 
								$start = $total_len * 3 + $flag;
								$total_len += $frame_length;
								$end = $total_len * 3 + $flag - 1;
								if($frame_length >= 115)  {
									$frame_number++;
									$repeat = $orf_width - length($frame_number);
									$frame_tag = '0' x$repeat . $frame_number;
									$id = $seq_id . "_$frame_tag";
									if($end > $seq_len)  {
										$end = $end - 3;
										$desc = "+$flag " . "$start-$end" . " len=$frame_length";
										$frame_seq = Bio::Seq->new(-id => $id, -description => $desc, -seq => $&);
										my $temp = $frame_seq->subseq(1, $frame_length - 1);
										$frame_seq = Bio::Seq->new(-id => $id, -description => $desc, -seq => $temp);
										push @seqs, $frame_seq;
									}else  {
										$desc = "+$flag " . "$start-$end" . " len=$frame_length";
										$frame_seq = Bio::Seq->new(-id => $id, -description => $desc, -seq => $&);
										push @seqs, $frame_seq;
									}
								}
								$subseq = $';
							}
						}else  {
							my $str = $seq->revcom->translate(-frame => $frame_flag);
							my $flag = $frame_flag + 1;
							my $subseq = $str->seq;
							my $start = 0;
							my $end = 0;
							my $repeat = 0;
							my $desc = '';
							my $frame_seq= '';
							my $frame_tag = '';
							my $total_len = 0;

							while($subseq =~ /(\w+)/i)  {
								$frame_length = length($&);
								$pre_length = length($`);
								$total_len += $pre_length; 
								$start = $total_len * 3 + $flag;
								$total_len += $frame_length;
								$end = $total_len * 3 + $flag - 1 ;
								$start = $seq_len - $start + 1;
								$end = $seq_len - $end + 1;

							if($frame_length >= 115)  {
								$frame_number++;
								$repeat = $orf_width - length($frame_number);
								$frame_tag = '0' x$repeat . $frame_number;
								$id = $seq_id . "_$frame_tag";
								if($end <= 0)  {
									$end = $end + 3;
									$frame_length = $frame_length - 1;
									$desc = "-$flag " . "$end-$start" . " len=$frame_length";
									$frame_seq = Bio::Seq->new(-id => $id, -description => $desc, -seq => $&);
									my $temp = $frame_seq->subseq(1, $frame_length - 1);
									$frame_seq = Bio::Seq->new(-id => $id, -description => $desc, -seq => $temp);
									push @seqs, $frame_seq;
								}else  {
									$desc = "-$flag " . "$end-$start" . " len=$frame_length";
									$frame_seq = Bio::Seq->new(-id => $id, -description => $desc, -seq => $&);
									push @seqs, $frame_seq;
								}
							}
							$subseq = $';
							}
						}
					}
				}
			}
		}
		$out->write_seq(@seqs);
	}
	if($seq_count == 0)  {
		#Delete both the input file and the output file
		print "Empty Sequence file.\n";
		system ("rm @outputs");
		exit;
	}else  {
		# Return translation results
		return @outputs;
	}
}


=head2 Blast_search
 Title   : Blast analysis
 Usage   : @blast_result = &Blast($system, @files);
 Function: Invoke Blast program
 Returns : An array of files contains Blast results
 Args    :  
=cut

sub Blast_search  {
	my ($step, $input) = @_;
	my $output = $input . ".blast";
#	my $dir = getcwd();
	my $db;
	my $evalue;
	my $program;
	my $dir;
	my $path = `which bttcmp`;
	if ($path=~/(.+)\/bttcmp/) {
		$dir = $1;
	}
	if($step eq 'step1')  {
		$db = $dir . "/BTTCMP_db/bt_toxin/db/bt_toxin";
		$evalue = 1e-25;
		system("blastp -query $input -out $output -db $db -evalue $evalue -num_threads 4");
		return $output;
	}

	if($step eq 'step2')  {
		$db = $dir . "/BTTCMP_db/back/db/back";
		$evalue = 1e-30;
		system("blastp -query $input -out $output -db $db -evalue $evalue -num_threads 4");
		return $output;
	}
}



=head2 Blast parser
 Title   : Blast parser
 Usage   : @blast_result = &Blast();
 Function: Implement Blast analysis 
 Returns : An hash contains query id and detail Blast results
 Args    : An array of files contains multiple fasta format of sequences 
=cut

sub Blast_parser  {
	my ($in) = @_;
	my %step1_out;
	my $count = 0;

	my $fin = Bio::SearchIO->new(-format => 'blast', -file => $in);
	while(my $result = $fin->next_result)  {
		if($result->num_hits)  {
			my $query_id = $result->query_accession;
			my $query_desc = ($result->query_description eq "") ? "no" : $result->query_description;
			my $query_length = $result->query_length;
			my $hit = $result->next_hit;
			my $hit_desc = $hit->description;
			my $hit_id = $hit->accession;
			my $hit_length = $hit->length;
			my $hsp = $hit->next_hsp;
			my $evalue = $hsp->evalue;
			my $query_start = $hsp->start('query');
			my $query_end = $hsp->end('query');
			my $hit_start = $hsp->start('hit');
			my $hit_end = $hsp->end('hit');
			my $aln_length = $hsp->length;
			my $percent_identity = $hsp->percent_identity;
			my $desc = "Query_desc:$query_desc\tQuery_Length:$query_length\tQuery_Start-End:$query_start-$query_end\tHit_id:$hit_id\tHit_desc:$hit_desc\tHit_length:$hit_length\tHit_Start-End:$hit_start-$hit_end\tAln_length:$aln_length\tPercent_identity:$percent_identity\tE-value:$evalue";
			$step1_out{$query_id} = $desc;
		}
	}

	#delete blast file
	system("rm $in");
	#print "Step1 identifies about $count proteins.\n";
	return %step1_out;
}


=head2 Blast_parser1
 Title   :
 Usage   :
 Function:
 Returns :
 Args    :
=cut

sub Blast_parser1  {
	my ($blast_result) = @_;
	my %step2_out;

	my $fin = Bio::SearchIO->new(-file => $blast_result, -format => 'blast');
	while(my $result = $fin->next_result)  {
		my $query_id = $result->query_accession;
		if($result->num_hits)  {
			$step2_out{$query_id} = '';
		}
	}
	#Clear the temp file
	system("rm $blast_result");
	#print "$count proteins have been deleted.\n";
	return %step2_out;
}


=head2 svm_prediction
 Title    : svm_prediction
 usage    :
 Function :
 Return   :
 Arg      :
=cut

sub svm_prediction  {
	my ($input) = @_;
#	my $dir = getcwd;
	my %svm_results;
	my $dir;
	my $path = `which bttcmp`;
		if ($path=~/(.+)\/bttcmp/) {
		$dir = $1;
	}
	#set the parameter 
	my $svm_model = $dir . "/BTTCMP_models/svm_model/model";

	#set file holding sequence feature
	my $svm_input = $input . ".feat";
	my $svm_output = $input . ".svmout";
	my @ids;
	my @ids_desc;

	my $seqin = Bio::SeqIO->new(-file => $input, -format => "fasta");
	open(MYHAND, ">$svm_input") || die "could not open.\n";
	while(my $seq = $seqin->next_seq())  {
		my $len = $seq->length;
		my $id = $seq->id;
		my $desc = $seq->desc;
		my $i = 1;
		$desc = $desc . "\t" . $len;
		my $hash = count_dimers($seq);
		push @ids, $id;
		push @ids_desc, $desc;
		print MYHAND "0";
		foreach (@$hash)  {
			print MYHAND "\t$i:";
			printf MYHAND "%0.8f", $_;
			$i++;
		}
		print MYHAND "\n";
	}

	#Implement SVM prediction
	system("svm-predict $svm_input $svm_model $svm_output");

	#parse SVM prediction results
	open(MYHAND, $svm_output) || die "could not open.\n";
	my $line = 0;
	while(<MYHAND>)  {
		chomp;
		if($_ == 1)  {
			$svm_results{$ids[$line]} = $ids_desc[$line] . "\t" . "YES";
		}
		$line++;
	}
	close MYHAND;

	#delete SVM prediction input and output file
	system("rm $svm_input $svm_output");

	#return SVM results
	return %svm_results;
}


=head2 Guess header format
 Title   : Guess header format
 Usage   : $header_format = &Guess_header_format($file);
 Function: Guess the header format of the sequence, two kinds of header could be recoganized,
	         NCBI header format: >gi|30018279|ref|NP_829910.1|, or normal header format: Contig1_020, contig132
	         The major difference between the NCBI header and the normal header format is that the NCBI header contains
	         several keys with each sepearated with "|", while the normal header format do not contain such notation "|" 
 Returns : string, either ncbi or normal
 Args    : File contain the sample sequences
=cut

sub Guess_header_format  {
	my ($file) = @_;
	my $in = Bio::SeqIO->new(-file => $file, -format => "fasta");
	my $seq = $in->next_seq();
	my $format;
	if($seq->id =~ /^gi\|.+\|.+\|(.+)\|/)  {
		$format = 'ncbi';
	}else  {
		$format = 'normal';
	}
	$in->DESTROY();
	return $format;
}


=head2 Index_maker
 Title   : Index maker
 Usage   : $index = &Index_maker($format, $file);
 Function: Give the file and its format, the Idex_maker will give the index file
 Returns : The index file name
 Args    : The file and the correct format of its fasta header
=cut

sub Index_maker  {
	my ($format, $file) = @_;
	my $db;
	if($format eq 'normal')  {
		$db = Bio::DB::Fasta->new($file);
	}else  {
		$db = Bio::DB::Fasta->new($file, -makeid => \&make_ncbi_id);     
	}
	return $db;
}


=head2 Make_ncbi_id
 Title   : Make NCBI ID
 Usage   : $db = Bio::DB::Fasta->new($file, -makeid=>\&make_ncbi_id);
 Function: The make_ncbi_id option gives you a chance to modify sequence IDs during indexing
 Returns : A scalar result
 Args    : A code reference that will take a scalar argument
=cut

sub make_ncbi_id  {
	my $description_line = shift;
	# get a different id from the fasta header, e.g.
	$description_line =~ /^>gi\|.+\|.+\|(.+)\|/;  
	return $1;
}

=head2 hmm_prediction
 Title   : hmm_prediction
 Usage   : %hmm_result = hmm_prediction($seqfile, $hmm_models);
 Function: Return Hidden Markov model prediction results
 Args    : HMM model and sequence file
 Return  :
=cut

sub hmm_prediction  {
	my ($seq_file, @hmm_models) = @_;
	my %HMM_results;

	# Implement HMM searching and parsing
	foreach my $hmm(@hmm_models)  {
		#For each hmm models, implement hmmsearch
		my $hmm_file_out = $seq_file . ".hmmout";
		system("hmmsearch -E 1e-10 -o $hmm_file_out $hmm $seq_file");

		#HMM results parsing
		my $searchio = Bio::SearchIO->new(-file => $hmm_file_out, -format => "hmmer3");
		while(my $result = $searchio->next_result)  {
			if($result->num_hits)  {
				while(my $hit = $result->next_hit)  {
					my $hmm_name = $result->query_name;
					my $hmm_len = $result->query_length;
					my $hit_desc = $hit->description;
					my $hsp = $hit->next_hsp;
					my $evalue = $hsp->evalue;
					my $id = $hit->name;
					while($hsp = $hit->next_hsp)  {
						if($evalue > $hsp->evalue)  {
							$evalue = $hsp->evalue;
						}
					}
					if(exists $HMM_results{$id})  {
						my @hmm_results = split /\t/, $HMM_results{$id};
						my @temp = split /:/, $hmm_results[-1];
						my $evalue1 = $temp[1];
						if($evalue1 > $evalue)  {
							$HMM_results{$id} = "Hit_desc:$hit_desc\tHmm_name:$hmm_name\tHmm_len:$hmm_len\tE-value:$evalue";
						}
					}else  {
						$HMM_results{$id} = "Hit_desc:$hit_desc\tHmm_name:$hmm_name\tHmm_len:$hmm_len\tE-value:$evalue";
					}
				}
			}
		}
		system("rm $hmm_file_out");
	}

	#Return hmm prediction result
	return %HMM_results;
}


=head2


=cut

sub printnoresults  {
	my $str = @_;#
	open OUT, ">>Strains_without_toxins_found.txt" || die;#
	print OUT $str . "\n";#
	#print "No Bt toxin has been detected.\n";
	exit;
}


=head2


=cut

sub count_dimers  {
	my @dip_results;
	my %count = (
               'AA' => 0, 'AR' => 0, 'AN' => 0, 'AD' => 0, 'AC' => 0, 'AQ' => 0, 'AE' => 0, 'AG' => 0, 'AH' => 0, 'AI' => 0, 
               'AL' => 0, 'AK' => 0, 'AM' => 0, 'AF' => 0, 'AP' => 0, 'AS' => 0, 'AT' => 0, 'AW' => 0, 'AY' => 0, 'AV' => 0, 
               'RA' => 0, 'RR' => 0, 'RN' => 0, 'RD' => 0, 'RC' => 0, 'RQ' => 0, 'RE' => 0, 'RG' => 0, 'RH' => 0, 'RI' => 0,
               'RL' => 0, 'RK' => 0, 'RM' => 0, 'RF' => 0, 'RP' => 0, 'RS' => 0, 'RT' => 0, 'RW' => 0, 'RY' => 0, 'RV' => 0,
               'NA' => 0, 'NR' => 0, 'NN' => 0, 'ND' => 0, 'NC' => 0, 'NQ' => 0, 'NE' => 0, 'NG' => 0, 'NH' => 0, 'NI' => 0,
               'NL' => 0, 'NK' => 0, 'NM' => 0, 'NF' => 0, 'NP' => 0, 'NS' => 0, 'NT' => 0, 'NW' => 0, 'NY' => 0, 'NV' => 0,
               'DA' => 0, 'DR' => 0, 'DN' => 0, 'DD' => 0, 'DC' => 0, 'DQ' => 0, 'DE' => 0, 'DG' => 0, 'DH' => 0, 'DI' => 0,
               'DL' => 0, 'DK' => 0, 'DM' => 0, 'DF' => 0, 'DP' => 0, 'DS' => 0, 'DT' => 0, 'DW' => 0, 'DY' => 0, 'DV' => 0,
               'CA' => 0, 'CR' => 0, 'CN' => 0, 'CD' => 0, 'CC' => 0, 'CQ' => 0, 'CE' => 0, 'CG' => 0, 'CH' => 0, 'CI' => 0,
               'CL' => 0, 'CK' => 0, 'CM' => 0, 'CF' => 0, 'CP' => 0, 'CS' => 0, 'CT' => 0, 'CW' => 0, 'CY' => 0, 'CV' => 0,
               'QA' => 0, 'QR' => 0, 'QN' => 0, 'QD' => 0, 'QC' => 0, 'QQ' => 0, 'QE' => 0, 'QG' => 0, 'QH' => 0, 'QI' => 0,
               'QL' => 0, 'QK' => 0, 'QM' => 0, 'QF' => 0, 'QP' => 0, 'QS' => 0, 'QT' => 0, 'QW' => 0, 'QY' => 0, 'QV' => 0,
               'EA' => 0, 'ER' => 0, 'EN' => 0, 'ED' => 0, 'EC' => 0, 'EQ' => 0, 'EE' => 0, 'EG' => 0, 'EH' => 0, 'EI' => 0,
               'EL' => 0, 'EK' => 0, 'EM' => 0, 'EF' => 0, 'EP' => 0, 'ES' => 0, 'ET' => 0, 'EW' => 0, 'EY' => 0, 'EV' => 0,
               'GA' => 0, 'GR' => 0, 'GN' => 0, 'GD' => 0, 'GC' => 0, 'GQ' => 0, 'GE' => 0, 'GG' => 0, 'GH' => 0, 'GI' => 0,
               'GL' => 0, 'GK' => 0, 'GM' => 0, 'GF' => 0, 'GP' => 0, 'GS' => 0, 'GT' => 0, 'GW' => 0, 'GY' => 0, 'GV' => 0,
               'HA' => 0, 'HR' => 0, 'HN' => 0, 'HD' => 0, 'HC' => 0, 'HQ' => 0, 'HE' => 0, 'HG' => 0, 'HH' => 0, 'HI' => 0,
               'HL' => 0, 'HK' => 0, 'HM' => 0, 'HF' => 0, 'HP' => 0, 'HS' => 0, 'HT' => 0, 'HW' => 0, 'HY' => 0, 'HV' => 0,
               'IA' => 0, 'IR' => 0, 'IN' => 0, 'ID' => 0, 'IC' => 0, 'IQ' => 0, 'IE' => 0, 'IG' => 0, 'IH' => 0, 'II' => 0,
               'IL' => 0, 'IK' => 0, 'IM' => 0, 'IF' => 0, 'IP' => 0, 'IS' => 0, 'IT' => 0, 'IW' => 0, 'IY' => 0, 'IV' => 0,
               'LA' => 0, 'LR' => 0, 'LN' => 0, 'LD' => 0, 'LC' => 0, 'LQ' => 0, 'LE' => 0, 'LG' => 0, 'LH' => 0, 'LI' => 0,
               'LL' => 0, 'LK' => 0, 'LM' => 0, 'LF' => 0, 'LP' => 0, 'LS' => 0, 'LT' => 0, 'LW' => 0, 'LY' => 0, 'LV' => 0,
               'KA' => 0, 'KR' => 0, 'KN' => 0, 'KD' => 0, 'KC' => 0, 'KQ' => 0, 'KE' => 0, 'KG' => 0, 'KH' => 0, 'KI' => 0,
               'KL' => 0, 'KK' => 0, 'KM' => 0, 'KF' => 0, 'KP' => 0, 'KS' => 0, 'KT' => 0, 'KW' => 0, 'KY' => 0, 'KV' => 0,
               'MA' => 0, 'MR' => 0, 'MN' => 0, 'MD' => 0, 'MC' => 0, 'MQ' => 0, 'ME' => 0, 'MG' => 0, 'MH' => 0, 'MI' => 0,
               'ML' => 0, 'MK' => 0, 'MM' => 0, 'MF' => 0, 'MP' => 0, 'MS' => 0, 'MT' => 0, 'MW' => 0, 'MY' => 0, 'MV' => 0,
               'FA' => 0, 'FR' => 0, 'FN' => 0, 'FD' => 0, 'FC' => 0, 'FQ' => 0, 'FE' => 0, 'FG' => 0, 'FH' => 0, 'FI' => 0,
               'FL' => 0, 'FK' => 0, 'FM' => 0, 'FF' => 0, 'FP' => 0, 'FS' => 0, 'FT' => 0, 'FW' => 0, 'FY' => 0, 'FV' => 0,
               'PA' => 0, 'PR' => 0, 'PN' => 0, 'PD' => 0, 'PC' => 0, 'PQ' => 0, 'PE' => 0, 'PG' => 0, 'PH' => 0, 'PI' => 0,
               'PL' => 0, 'PK' => 0, 'PM' => 0, 'PF' => 0, 'PP' => 0, 'PS' => 0, 'PT' => 0, 'PW' => 0, 'PY' => 0, 'PV' => 0,
               'SA' => 0, 'SR' => 0, 'SN' => 0, 'SD' => 0, 'SC' => 0, 'SQ' => 0, 'SE' => 0, 'SG' => 0, 'SH' => 0, 'SI' => 0,
               'SL' => 0, 'SK' => 0, 'SM' => 0, 'SF' => 0, 'SP' => 0, 'SS' => 0, 'ST' => 0, 'SW' => 0, 'SY' => 0, 'SV' => 0,
               'TA' => 0, 'TR' => 0, 'TN' => 0, 'TD' => 0, 'TC' => 0, 'TQ' => 0, 'TE' => 0, 'TG' => 0, 'TH' => 0, 'TI' => 0,
               'TL' => 0, 'TK' => 0, 'TM' => 0, 'TF' => 0, 'TP' => 0, 'TS' => 0, 'TT' => 0, 'TW' => 0, 'TY' => 0, 'TV' => 0,
               'WA' => 0, 'WR' => 0, 'WN' => 0, 'WD' => 0, 'WC' => 0, 'WQ' => 0, 'WE' => 0, 'WG' => 0, 'WH' => 0, 'WI' => 0,
               'WL' => 0, 'WK' => 0, 'WM' => 0, 'WF' => 0, 'WP' => 0, 'WS' => 0, 'WT' => 0, 'WW' => 0, 'WY' => 0, 'WV' => 0,
               'YA' => 0, 'YR' => 0, 'YN' => 0, 'YD' => 0, 'YC' => 0, 'YQ' => 0, 'YE' => 0, 'YG' => 0, 'YH' => 0, 'YI' => 0,
               'YL' => 0, 'YK' => 0, 'YM' => 0, 'YF' => 0, 'YP' => 0, 'YS' => 0, 'YT' => 0, 'YW' => 0, 'YY' => 0, 'YV' => 0,
               'VA' => 0, 'VR' => 0, 'VN' => 0, 'VD' => 0, 'VC' => 0, 'VQ' => 0, 'VE' => 0, 'VG' => 0, 'VH' => 0, 'VI' => 0,
               'VL' => 0, 'VK' => 0, 'VM' => 0, 'VF' => 0, 'VP' => 0, 'VS' => 0, 'VT' => 0, 'VW' => 0, 'VY' => 0, 'VV' => 0,
	);
	my $seqobj = shift @_;
	my $seqstring = uc $seqobj->seq();
	my $seq_len = $seqobj->length;
	my $seq_stats = Bio::Tools::SeqStats->new(-seq => $seqobj);
	my $hash = $seq_stats->count_monomers();
	foreach my $aa (sort keys %$hash)  {
		if($aa =~ /[acdefghiklmnpqrstvwy]/i)  {
			push @dip_results, $$hash{$aa}/$seq_len;
		}
	}
	foreach my $element( keys %count)  {
		$count{$element} = ($seqstring =~ s/$element/$element/g);
		$count{$element} = $count{$element}/($seq_len - 1);
		push @dip_results, $count{$element};
	}
	return \@dip_results;
}


=head2 Writer
 Title   : Writer
 Usage   : 
 Function: 
 Return  :
 Args    :
=cut

sub Writer  {
	my ($seq_type, $final_output, $seqfile, %step3_in, %domain) = @_;#2020/4/7
	#make index for the original sequence file
	my $orign_seq = pop @$seqfile;
	my $db = &Index_maker(&Guess_header_format($orign_seq), $orign_seq);
	#*#my $index =  $orign_seq . ".index";
	
	#define final output files
	my $table_list = $final_output . ".list";
	my $seqout = $final_output . ".gbk";
	my @results = ($table_list, $seqout);
	my $seqin = $final_output . ".step2";
	my %seq_info;
	open(MYHAND, ">$table_list") || die "could not open.\n";
	
	my $in = Bio::SeqIO->new(-file=> "$seqin", -format=> "fasta");
	while(my $seq = $in->next_seq())  {
		$seq_info{$seq->id} = {
								"seq_length" => $seq->length,
								"seq_desc"   => $seq->desc,
							};
	}

	my %rank4;
	my %rank3;
	my %rank2;
	my %rank1;
	my %cry;
	my %cyt;
	my %vip;
	my %others;

	my %Cry_count =  (
						'Rank1' => 0,
						'Rank2' => 0,
						'Rank3' => 0,
						'Rank4' => 0,
					);

	my %Cyt_count =  (
						'Rank1' => 0,
						'Rank2' => 0,
						'Rank3' => 0,
						'Rank4' => 0,
					);

	my %Vip_count =  (
						'Rank1' => 0,
						'Rank2' => 0,
						'Rank3' => 0,
						'Rank4' => 0,
					);

	my %others_count =  (
						'Rank1' => 0,
						'Rank2' => 0,
						'Rank3' => 0,
						'Rank4' => 0,
					);
	my $Cyt_count = 0;
	my $Vip_count = 0;
	my $others_count = 0;
	my %Summary_count =  (
							'Cry'     => 0,
							'Cyt'     => 0,
							'Vip'     => 0,
							'Summary' => 0,
							'others'  => 0,
						);

	my $blast = $step3_in{'blast'};
	my $hmm = $step3_in{'hmm'};
	my $svm = $step3_in{'svm'};
	my $hmmem = $step3_in{'em'};#added
	my $hmmen = $step3_in{'en'};#added
	my $hmmec = $step3_in{'ec'};#added
	my $hmmemid = $step3_in{'emid'};#added
	my $hmmmtx = $step3_in{'mtx'};#added
	my $hmmtox = $step3_in{'tox'};#added
#	print "\nblast: $blast\nhmm: $hmm\nsvm: $svm\nhmmem: $hmmem\nhmmen: $hmmen\nhmmec: $hmmec\nhmmmtx: $hmmmtx\nhmmtox: $hmmtox\n";#2020/4/8

	#After that, Proteins in $blast variable contains only Cry proteins
	foreach (sort keys %$hmm)  {
		my @hmm_split = split /\t/, $$hmm{$_};
		my @protein_desc = split /:/, $hmm_split[0];
		#process cyt predicted by hmm method
		if($$hmm{$_} =~ /cyt.hmm/i)  {
			if(exists $$blast{$_})  {
				my %temp;
				my @blast_split = split /\t/, $$blast{$_};
				foreach (@blast_split)  {
					my ($key, $value) = split /:/, $_;
					$temp{$key} = ($value eq '') ? "" : $value;
				}
				$cyt{$_} = {
							'rank'               => 'ND',
							'hmm_prediction'     => 'YES',
							'protein_id'         => $_,
							'protein_desc'       => $temp{'Query_desc'},
							'protein_len'        => $temp{'Query_Length'},
							'blast_prediction'   => 'YES',
							'best_hit'           => $temp{'Hit_id'},
							'hit_length'         => $temp{'Hit_length'},
							'coverage'           => $temp{'Aln_length'}/$temp{'Hit_length'}*100,
							'Percent_identity'   => $temp{'Percent_identity'},
							'e-value'            => $temp{'E-value'},
							'blast_detail'       => $$blast{$_},
							'Endotoxin_N'        => 'NA',#2020/4/7
							'Endotoxin_M'        => 'NA',#2020/4/7
							'Endotoxin_C'        => 'NA',#2020/4/7
							'Endotoxin_mid'        => 'NA',#2020/4/8
							'ETX_MTX2'           => 'NA',#2020/4/7
							'Toxin_10'           => 'NA',#2020/4/7
							};
				if(exists $$svm{$_})  {
					$cyt{$_}->{'svm_prediction'} = 'YES';
					$cyt{$_}->{'svm_detail'} = $$svm{$_};
				}else  {
					$cyt{$_}->{'svm_prediction'} = 'NO';
				}
				if($temp{'Percent_identity'} >= 95)  {
					$Cyt_count{'Rank4'}++;
					$cyt{$_}->{'rank'} = "Rank4";
				}elsif($temp{'Percent_identity'} >= 78)  {
					$Cyt_count{'Rank3'}++;
					$cyt{$_}->{'rank'} = "Rank3";
				}elsif($temp{'Percent_identity'} >= 45)  {
					$Cyt_count{'Rank2'}++;
					$cyt{$_}->{'rank'} = "Rank2";
				}else  {
					$Cyt_count{'Rank1'}++;
					$cyt{$_}->{'rank'} = "Rank1";
				}
				delete $$blast{$_};
			}else {
				$Cyt_count++;
				$cyt{$_} = {
							  'protein_id'         => $_,
							  'protein_desc'       => $protein_desc[1],
							  'protein_len'        => $seq_info{$_}->{'seq_length'},
							  'rank'               => 'ND',
							  'hmm_prediction'     => 'YES',
							  'blast_prediction'   => 'ND',
							  'best_hit'           => 'ND',
							  'hit_length'         => 'ND',
							  'coverage'           => 'ND',
							  'Percent_identity'   => 'ND',
							  'e-value'            => 'ND',
							  'blast_detail'       => 'ND',
							  'hmm_detail'         => $$hmm{$_},
							  'Endotoxin_N'        => 'NA',#2020/4/7
							  'Endotoxin_M'        => 'NA',#2020/4/7
							  'Endotoxin_C'        => 'NA',#2020/4/7
							  'Endotoxin_mid'        => 'NA',#2020/4/8
							  'ETX_MTX2'           => 'NA',#2020/4/7
							  'Toxin_10'           => 'NA',#2020/4/7
							};
				if(exists $$svm{$_})  {
						$cyt{$_}->{'svm_prediction'} = 'YES';
						$cyt{$_}->{'svm_detail'} = $$svm{$_};
						delete  $$svm{$_};
				}else  {
					$cyt{$_}->{'svm_prediction'} = 'NO';
				}
			}
		}
		#process vip predicted by hmm method
		if($$hmm{$_} =~ /vip.hmm/i)  {
			if(!(exists $$blast{$_}))  {
				$Vip_count++;
				$vip{$_} = {
							'protein_id'         => $_,
							'protein_desc'       => $protein_desc[1],
							'protein_len'        => $seq_info{$_}->{'seq_length'},
							'rank'               => 'ND',
							'hmm_prediction'     => 'YES',
							'blast_prediction'   => 'ND',
							'best_hit'           => 'ND',
							'hit_length'         => 'ND',
							'coverage'           => 'ND',
							'Percent_identity'   => 'ND',
							'e-value'            => 'ND',
							'blast_detail'       => 'ND',
							'hmm_detail'         => $$hmm{$_},
							'Endotoxin_N'        => 'NA',#2020/4/7
							'Endotoxin_M'        => 'NA',#2020/4/7
							'Endotoxin_C'        => 'NA',#2020/4/7
							'Endotoxin_mid'        => 'NA',#2020/4/8
							'ETX_MTX2'           => 'NA',#2020/4/7
							'Toxin_10'           => 'NA',#2020/4/7
							};
				if(exists $$svm{$_})  {
					$vip{$_}->{'svm_prediction'} = 'YES';
					$vip{$_}->{'svm_detail'} = $$svm{$_};
				}else  {
					$vip{$_}->{'svm_prediction'} = 'NO';
				}
			}else  {
				my %temp;
				my @blast_split = split /\t/, $$blast{$_};
				foreach (@blast_split)  {
					my ($key, $value) = split /:/, $_;
					$temp{$key} = ($value eq '') ? "" : $value;
				}
				$vip{$_} = {
							'rank'               => 'ND',
							'hmm_prediction'     => 'YES',
							'protein_id'         => $_,
							'protein_desc'       => $temp{'Query_desc'},
							'protein_len'        => $temp{'Query_Length'},
							'blast_prediction'   => 'YES',
							'best_hit'           => $temp{'Hit_id'},
							'hit_length'         => $temp{'Hit_length'},
							'coverage'           => $temp{'Aln_length'}/$temp{'Hit_length'}*100,
							'Percent_identity'   => $temp{'Percent_identity'},
							'e-value'            => $temp{'E-value'},
							'blast_detail'       => $$blast{$_},
							'Endotoxin_N'        => 'NA',#2020/4/7
							'Endotoxin_M'        => 'NA',#2020/4/7
							'Endotoxin_C'        => 'NA',#2020/4/7
							'Endotoxin_mid'        => 'NA',#2020/4/8
							'ETX_MTX2'           => 'NA',#2020/4/7
							'Toxin_10'           => 'NA',#2020/4/7
							};
				if(exists $$svm{$_})  {
					$vip{$_}->{'svm_prediction'} = 'YES';
					$vip{$_}->{'svm_detail'} = $$svm{$_};
				}else  {
					$vip{$_}->{'svm_prediction'} = 'NO';
				}
				if($temp{'Percent_identity'} >= 95)  {
					$Vip_count{'Rank4'}++;
					$vip{$_}->{'rank'} = "Rank4";
				}elsif($temp{'Percent_identity'} >= 78)  {
					$Vip_count{'Rank3'}++;
					$vip{$_}->{'rank'} = "Rank3";
				}elsif($temp{'Percent_identity'} >= 45)  {
					$Vip_count{'Rank2'}++;
					$vip{$_}->{'rank'} = "Rank2";
				}else  {
					$Vip_count{'Rank1'}++;
					$vip{$_}->{'rank'} = "Rank1";
				}
				delete $$blast{$_};
			}
		}
	}

	#Obtain Cry protein ID
	foreach (sort keys %seq_info)  {
		if((!exists $vip{$_}) && !(exists $cyt{$_}))  {
			$cry{$_} = '';#sequences ids contain crys(blast and hmm and svm) and (vip and cyt predicted only by blast)
		}
	}

	foreach (sort keys %cry)  {
		my $EN = "NO";
		my $EM = "NO";
		my $EC = "NO";
		my $EMID = "NO";
		my $MTX = "NO";
		my $TOXIN = "NO";
		if (exists $$hmmem{$_}) {
			$EM = "YES";
		}
		if (exists $$hmmen{$_}) {
			$EN = "YES";
		}
		if (exists $$hmmec{$_}) {
			$EC = "YES";
		}
		if (exists $$hmmemid{$_}) {
			$EMID = "YES";
		}
		if (exists $$hmmmtx{$_}) {
			$MTX = "YES";
		}
		if (exists $$hmmtox{$_}) {
			$TOXIN = "YES";
		}
		if(exists $$blast{$_})  {
=pod
			if (exists $$hmmem{$_}) {
				$EM = "YES";
				#print "EM=YES\n";
			}
			
			if (exists $$hmmen{$_}) {
				$EN = "YES";
				#print "EN=YES\n";
			}
			if (exists $$hmmec{$_}) {
				$EC = "YES";
				#print "EC=YES\n";
			}
			if (exists $$hmmemid{$_}) {
				$EMID = "YES";
				#print "EC=YES\n";
			}
			if (exists $$hmmmtx{$_}) {
				$MTX = "YES";
				#print "MTX=YES\n";
			}
			if (exists $$hmmtox{$_}) {
				$TOXIN = "YES";
				#print "TOXIN=YES\n";
			}
=cut
			my %temp;
			my @blast_split = split /\t/, $$blast{$_};
				foreach (@blast_split)  {
					my ($key, $value) = split /:/, $_;
					$temp{$key} = ($value eq '') ? "" : $value;
				}
				#added by liu to excluding vip and cyt
				if ($temp{'Hit_id'} =~ /vip/i) {
					$vip{$_} = {
								'hmm_prediction'     => 'NO',
								'protein_id'         => $_,
								'protein_desc'       => $temp{'Query_desc'},
								'protein_len'        => $temp{'Query_Length'},
								'blast_prediction'   => 'YES',
								'best_hit'           => $temp{'Hit_id'},
								'hit_length'         => $temp{'Hit_length'},
								'coverage'           => $temp{'Aln_length'}/$temp{'Hit_length'}*100,
								'Percent_identity'   => $temp{'Percent_identity'},
								'e-value'            => $temp{'E-value'},
								'blast_detail'       => $$blast{$_},
								'Endotoxin_N'        => 'NA',#2020/4/7
								'Endotoxin_M'        => 'NA',#2020/4/7
								'Endotoxin_C'        => 'NA',#2020/4/7
								'Endotoxin_mid'        => 'NA',#2020/4/8
								'ETX_MTX2'           => 'NA',#2020/4/7
								'Toxin_10'           => 'NA',#2020/4/7
								};
					if(exists $$svm{$_})  {
						$vip{$_}->{'svm_prediction'} = 'YES';
						$vip{$_}->{'svm_detail'} = $$svm{$_};
					}else  {
						$vip{$_}->{'svm_prediction'} = 'NO';
					}
					if($temp{'Percent_identity'} >= 95)  {
						$Vip_count{'Rank4'}++;
						$vip{$_}->{'rank'} = "Rank4";
					}elsif($temp{'Percent_identity'} >= 78)  {
						$Vip_count{'Rank3'}++;
						$vip{$_}->{'rank'} = "Rank3";
					}elsif($temp{'Percent_identity'} >= 45)  {
						$Vip_count{'Rank2'}++;
						$vip{$_}->{'rank'} = "Rank2";
					}else  {
						$Vip_count{'Rank1'}++;
						$vip{$_}->{'rank'} = "Rank1";
					}
					delete $cry{$_};
				}elsif ($temp{'Hit_id'} =~ /cyt/i) {
					$cyt{$_} = {
								'hmm_prediction'     => 'NO',
								'protein_id'         => $_,
								'protein_desc'       => $temp{'Query_desc'},
								'protein_len'        => $temp{'Query_Length'},
								'blast_prediction'   => 'YES',
								'best_hit'           => $temp{'Hit_id'},
								'hit_length'         => $temp{'Hit_length'},
								'coverage'           => $temp{'Aln_length'}/$temp{'Hit_length'}*100,
								'Percent_identity'   => $temp{'Percent_identity'},
								'e-value'            => $temp{'E-value'},
								'blast_detail'       => $$blast{$_},
								'Endotoxin_N'        => 'NA',#2020/4/7
								'Endotoxin_M'        => 'NA',#2020/4/7
								'Endotoxin_C'        => 'NA',#2020/4/7
								'Endotoxin_mid'        => 'NA',#2020/4/8
								'ETX_MTX2'           => 'NA',#2020/4/7
								'Toxin_10'           => 'NA',#2020/4/7
								};
					if(exists $$svm{$_})  {
						$cyt{$_}->{'svm_prediction'} = 'YES';
						$cyt{$_}->{'svm_detail'} = $$svm{$_};
					}else  {
						$cyt{$_}->{'svm_prediction'} = 'NO';
					}
					if($temp{'Percent_identity'} >= 95)  {
						$Cyt_count{'Rank4'}++;
						$cyt{$_}->{'rank'} = "Rank4";
					}elsif($temp{'Percent_identity'} >= 78)  {
						$Cyt_count{'Rank3'}++;
						$cyt{$_}->{'rank'} = "Rank3";
					}elsif($temp{'Percent_identity'} >= 45)  {
						$Cyt_count{'Rank2'}++;
						$cyt{$_}->{'rank'} = "Rank2";
					}else  {
						$Cyt_count{'Rank1'}++;
						$cyt{$_}->{'rank'} = "Rank1";
					}
					delete $cry{$_};
				}elsif ($temp{'Hit_id'} =~ /other/i) {
					$others{$_} = {
								'hmm_prediction'     => 'NA',
								'protein_id'         => $_,
								'protein_desc'       => $temp{'Query_desc'},
								'protein_len'        => $temp{'Query_Length'},
								'blast_prediction'   => 'YES',
								'best_hit'           => $temp{'Hit_id'},
								'hit_length'         => $temp{'Hit_length'},
								'coverage'           => $temp{'Aln_length'}/$temp{'Hit_length'}*100,
								'Percent_identity'   => $temp{'Percent_identity'},
								'e-value'            => $temp{'E-value'},
								'blast_detail'       => $$blast{$_},
								'Endotoxin_N'        => $EN,#2020/4/7
								'Endotoxin_M'        => $EM,#2020/4/7
								'Endotoxin_C'        => $EC,#2020/4/7
								'Endotoxin_mid'        => $EMID,#2020/4/8
								'ETX_MTX2'           => $MTX,#2020/4/7
								'Toxin_10'           => $TOXIN,#2020/4/7
								};
					if(exists $$svm{$_})  {
						$others{$_}->{'svm_prediction'} = 'YES';
						$others{$_}->{'svm_detail'} = $$svm{$_};
					}else  {
						$others{$_}->{'svm_prediction'} = 'NO';
					}
					if($temp{'Percent_identity'} >= 95)  {
						$others_count{'Rank4'}++;
						$others{$_}->{'rank'} = "Rank4";
					}elsif($temp{'Percent_identity'} >= 78)  {
						$others_count{'Rank3'}++;
						$others{$_}->{'rank'} = "Rank3";
					}elsif($temp{'Percent_identity'} >= 45)  {
						$others_count{'Rank2'}++;
						$others{$_}->{'rank'} = "Rank2";
					}else  {
						$others_count{'Rank1'}++;
						$others{$_}->{'rank'} = "Rank1";
					}
					delete $cry{$_};
				}else {
					$cry{$_} = {
								'protein_id'         => $_,
								'protein_desc'       => $temp{'Query_desc'},
								'protein_len'        => $temp{'Query_Length'},
								'blast_prediction'   => 'YES',
								'best_hit'           => $temp{'Hit_id'},
								'hit_length'         => $temp{'Hit_length'},
								'coverage'           => $temp{'Aln_length'}/$temp{'Hit_length'}*100,
								'Percent_identity'   => $temp{'Percent_identity'},
								'e-value'            => $temp{'E-value'},
								'blast_detail'       => $$blast{$_},
								'Endotoxin_N'        => $EN,#2020/4/7
								'Endotoxin_M'        => $EM,#2020/4/7
								'Endotoxin_C'        => $EC,#2020/4/7
								'Endotoxin_mid'        => $EMID,#2020/4/8
								'ETX_MTX2'           => $MTX,#2020/4/7
								'Toxin_10'           => $TOXIN,#2020/4/7
								};
					if($temp{'Percent_identity'} >= 95)  {
						$Cry_count{'Rank4'}++;
						$cry{$_}->{'rank'} = "Rank4";
					}
					elsif($temp{'Percent_identity'} >= 78)  {
						$Cry_count{'Rank3'}++;
						$cry{$_}->{'rank'} = "Rank3";
					}
					elsif($temp{'Percent_identity'} >= 45)  {
						$Cry_count{'Rank2'}++;
						$cry{$_}->{'rank'} = "Rank2";
					}
					else  {
						$Cry_count{'Rank1'}++;
						$cry{$_}->{'rank'} = "Rank1";
					}
					if(exists $$hmm{$_})  {
						$cry{$_}->{'hmm_prediction'} = 'YES';
						$cry{$_}->{'hmm_detail'} = $$hmm{$_};
					}
					else  {
						$cry{$_}->{'hmm_prediction'} = 'NO';
					}
					if(exists $$svm{$_})  {
						$cry{$_}->{'svm_prediction'} = 'YES';
						$cry{$_}->{'svm_detail'} = $$svm{$_};
					}
					else  {
						$cry{$_}->{'svm_prediction'} = 'NO';
					}
				}
		}else  {
			$Cry_count{'Rank1'}++;
			$cry{$_} = {
						'protein_id'         => $_,
						'blast_prediction'   => 'NO',
						'best_hit'           => 'NO',
						'hit_length'         => 'NO',
						'coverage'           => 'NO',
						'Percent_identity'   => 'NO',
						'e-value'            => 'NO',
						'blast_detail'       => 'NO',
						'Endotoxin_N'        => $EN,#2020/4/7
						'Endotoxin_M'        => $EM,#2020/4/7
						'Endotoxin_C'        => $EC,#2020/4/7
						'Endotoxin_mid'        => $EMID,#2020/4/8
						'ETX_MTX2'           => $MTX,#2020/4/7
						'Toxin_10'           => $TOXIN,#2020/4/7
						'rank'               => 'Rank1',
						};
			if(exists $$hmm{$_})  {
				my @hmm_split = split /\t/, $$hmm{$_};
				my @hit_desc = split /:/, $hmm_split[0];
				my $hit_desc = $hit_desc[1];
				$cry{$_}->{'protein_desc'}   = $hit_desc;
				$cry{$_}->{'protein_len'}    = $seq_info{$_}->{'seq_length'};
				$cry{$_}->{'hmm_prediction'} = "YES";
				$cry{$_}->{'hmm_detail'}     = $$hmm{$_};
				if(exists $$svm{$_})  {
					$cry{$_}->{'svm_prediction'} = "YES";
				}
				else  {
					$cry{$_}->{'svm_prediction'} = 'NO';
				}
			}else  {
				my ($desc, $len, $mark) = split /\t/, $$svm{$_};
				$cry{$_}->{'hmm_prediction'} = 'NO';
				$cry{$_}->{'hmm_detail'} = 'NO';
				$cry{$_}->{'protein_desc'}   = $desc;
				$cry{$_}->{'protein_len'}    = $len;
				if(exists $$svm{$_})  {
					$cry{$_}->{svm_prediction} = "YES";
				}
			}
		}
	}

	$Summary_count{'Cry'}     = $Cry_count{'Rank1'} + $Cry_count{'Rank2'} + $Cry_count{'Rank3'} + $Cry_count{'Rank4'};
	$Summary_count{'Cyt'}     = $Cyt_count{'Rank1'} + $Cyt_count{'Rank2'} + $Cyt_count{'Rank3'} + $Cyt_count{'Rank4'} + $Cyt_count;
	$Summary_count{'Vip'}     = $Vip_count{'Rank1'} + $Vip_count{'Rank2'} + $Vip_count{'Rank3'} + $Vip_count{'Rank4'} + $Vip_count;
	$Summary_count{'others'}  = $others_count{'Rank1'} + $others_count{'Rank2'} + $others_count{'Rank3'} + $others_count{'Rank4'} + $others_count;
	my $Summary_count_R1      = $Cry_count{'Rank1'} + $Cyt_count{'Rank1'} + $Vip_count{'Rank1'} + $others_count{'Rank1'};
	my $Summary_count_R2      = $Cry_count{'Rank2'} + $Cyt_count{'Rank2'} + $Vip_count{'Rank2'} + $others_count{'Rank2'};
	my $Summary_count_R3      = $Cry_count{'Rank3'} + $Cyt_count{'Rank3'} + $Vip_count{'Rank3'} + $others_count{'Rank3'};
	my $Summary_count_R4      = $Cry_count{'Rank4'} + $Cyt_count{'Rank4'} + $Vip_count{'Rank4'} + $others_count{'Rank4'};
	my $Summary_count_HMMSVM  = $Cyt_count + $Vip_count;
	$Summary_count{'Summary'} = $Summary_count{'Cry'} + $Summary_count{'Cyt'} + $Summary_count{'Vip'} + $Summary_count{'others'};

	#Print Overview of prediction
	print MYHAND "Overview of prediction \t Sequence type: $seq_type\n";
	print MYHAND "Name\tCry\tCyt\tVip\tOthers\tSummary\n";
	print MYHAND "Rank1\t", $Cry_count{'Rank1'}, "\t$Cyt_count{'Rank1'}", "\t$Vip_count{'Rank1'}", "\t$others_count{'Rank1'}", "\t$Summary_count_R1", "\n";
	print MYHAND "Rank2\t", $Cry_count{'Rank2'}, "\t$Cyt_count{'Rank2'}", "\t$Vip_count{'Rank2'}", "\t$others_count{'Rank2'}", "\t$Summary_count_R2", "\n";
	print MYHAND "Rank3\t", $Cry_count{'Rank3'}, "\t$Cyt_count{'Rank3'}", "\t$Vip_count{'Rank3'}", "\t$others_count{'Rank3'}", "\t$Summary_count_R3", "\n";
	print MYHAND "Rank4\t", $Cry_count{'Rank4'}, "\t$Cyt_count{'Rank4'}", "\t$Vip_count{'Rank4'}", "\t$others_count{'Rank4'}", "\t$Summary_count_R4", "\n";
	print MYHAND "HMM_SVM\t", "0", "\t$Cyt_count", "\t$Vip_count\t", $Summary_count_HMMSVM, "\n";
	print MYHAND "Summary\t", $Summary_count{'Cry'}, "\t", $Summary_count{'Cyt'}, "\t", $Summary_count{'Vip'}, "\t", $Summary_count{'others'}, "\t", $Summary_count{'Summary'}, "\n";
	print MYHAND "//\n\n\n";

	my $ID = 1;
	print MYHAND "Toxin type: Cry protein\n";
	print MYHAND "ID\tProtein_ID\tProtein_description\tLength\tRank\tBLAST\tBest_hit\tHit length\tCoverage\tIdentity\tSVM\tHMM\n";
	
	foreach (sort keys %cry)  {
		print  MYHAND $ID, "\t";
		print  MYHAND $cry{$_}->{'protein_id'}, "\t";###Can't use string ("") as a HASH ref while "strict refs" in use
		print  MYHAND $cry{$_}->{'protein_desc'}, "\t";
		print  MYHAND $cry{$_}->{'protein_len'}, "\t";
		print  MYHAND $cry{$_}->{'rank'}, "\t";
		print  MYHAND $cry{$_}->{'blast_prediction'}, "\t";
		print  MYHAND $cry{$_}->{'best_hit'}, "\t";
		print  MYHAND $cry{$_}->{'hit_length'}, "\t";
		printf MYHAND ($cry{$_}->{'coverage'}=~ /\d+/)? ("%0.2f\t", $cry{$_}->{'coverage'}) : "$cry{$_}->{'coverage'}\t";
		printf MYHAND ($cry{$_}->{'Percent_identity'} =~ /\d+/)? ("%0.2f\t", $cry{$_}->{'Percent_identity'}):"$cry{$_}->{'Percent_identity'}\t";
		print  MYHAND $cry{$_}->{'svm_prediction'}, "\t";
		print  MYHAND $cry{$_}->{'hmm_prediction'}, "\n";
		$ID++;
	}
	print MYHAND "//\n\n\n";


	$ID = 1;
	print MYHAND "Toxin type: Cyt protein\n";
	print_toxin2($ID, %cyt); 
	
	$ID = 1;
	print MYHAND "Toxin type: Vip protein\n";
	print_toxin2($ID, %vip); #print_toxin

	$ID = 1;
	print MYHAND "Toxin type: Other toxins\n";
	print_toxin2($ID, %others); #print_toxin

	my %toxin_info;
	foreach (keys %cry)  {
		$toxin_info{$_} = $cry{$_};
	}
	foreach (keys %cyt)  {
		$toxin_info{$_} = $cyt{$_};
	}
	foreach (keys %vip)  {
		$toxin_info{$_} = $vip{$_};
	}
	foreach (keys %others)  {
		$toxin_info{$_} = $others{$_};
	}

	my $in1 = Bio::SeqIO->new(-file=> "$seqin", -format=> "fasta");
	my $out = Bio::SeqIO->new(-file=> ">$seqout", -format=> "genbank");
	while(my $seq = $in1->next_seq())  {
		if(exists $toxin_info{$seq->id})  {
			if($seq_type eq 'nucl')  {
				my $contig_id;
				if($seq->id =~ /(.*)_\d+/)  {
					$contig_id = $1;
				}
				my @temp =split /\s+/, $toxin_info{$seq->id}->{'protein_desc'};
				my $strand;
				my $frame;

				my $nucl_seq = $db->get_Seq_by_id($contig_id);
				my $seq1 = Bio::Seq->new(-id => $seq->id, -seq => $nucl_seq->seq, -accession_number => $seq->id);

				if($temp[0] > 0)  {
					$strand = 1;
				}else  {
					$strand = -1;
				}
				my ($start, $end) = split /-/, $temp[1];
				$contig_id = $contig_id . ":" . $start . "-" . $end;
				my $feat = Bio::SeqFeature::Generic->new(
															-start     => $start,
															-end       => $end,
															-primary   => 'Segment',
															-strand    => $strand,
															-tag =>  {
																		Segment_name      => $seq->id,
																		locus_tag         => $seq->id,
																		protein_id        => $toxin_info{$seq->id}->{'protein_id'},
																		protein_desc      => $toxin_info{$seq->id}->{'protein_desc'},
																		protein_len       => $toxin_info{$seq->id}->{'protein_len'},
																		Rank              => $toxin_info{$seq->id}->{'rank'},
																		svm_prediction    => $toxin_info{$seq->id}->{'svm_prediction'},
																		blast_prediction  => $toxin_info{$seq->id}->{'blast_prediction'},
																		blast_detail      => $toxin_info{$seq->id}->{'blast_detail'},
																		hmm_prediction    => $toxin_info{$seq->id}->{'hmm_prediction'},
																		hmm_detail        => $toxin_info{$seq->id}->{'hmm_detail'},
																		Endotoxin_N       => $toxin_info{$seq->id}->{'Endotoxin_N'},#2020/4/7
																		Endotoxin_M       => $toxin_info{$seq->id}->{'Endotoxin_M'},#2020/4/7
																		Endotoxin_C       => $toxin_info{$seq->id}->{'Endotoxin_C'},#2020/4/7
																		Endotoxin_mid     => $toxin_info{$seq->id}->{'Endotoxin_mid'},#2020/4/8
																		ETX_MTX2          => $toxin_info{$seq->id}->{'ETX_MTX2'},#2020/4/7
																		Toxin_10          => $toxin_info{$seq->id}->{'Toxin_10'},#2020/4/7
																		translation       => $seq->seq,
																	},
														);
				$seq1->add_SeqFeature($feat);
				$out->write_seq($seq1);                
			}
			elsif($seq_type eq 'orfs')  {
				my $nucl_seq = $db->get_Seq_by_id($seq->id);
				my $seq1 = Bio::Seq->new(-id => $seq->id, -seq => $nucl_seq->seq, -accession_number => $seq->id);
				my $feat = Bio::SeqFeature::Generic->new (
															-start      => 1,
															-end        => $nucl_seq->length,
															-primary    => 'CDS',
															-strand     => 1,
															-tag        => {
																			locus_tag         => $nucl_seq->id,
																			codon_start       => 1,
																			transl_table      => 11,
																			product           => 'Possible Bt toxin protein',
																			protein_id        => $seq->id,
																			translation       => $seq->seq,
																			protein_len       => $toxin_info{$seq->id}->{'protein_len'},
																			Rank              => $toxin_info{$seq->id}->{'rank'},
																			svm_prediction    => $toxin_info{$seq->id}->{'svm_prediction'},
																			blast_prediction  => $toxin_info{$seq->id}->{'blast_prediction'},
																			blast_detail      => $toxin_info{$seq->id}->{'blast_detail'},
																			hmm_prediction    => $toxin_info{$seq->id}->{'hmm_prediction'},
																			hmm_detail        => $toxin_info{$seq->id}->{'hmm_detail'},
																			Endotoxin_N       => $toxin_info{$seq->id}->{'Endotoxin_N'},#2020/4/7
																			Endotoxin_M       => $toxin_info{$seq->id}->{'Endotoxin_M'},#2020/4/7
																			Endotoxin_C       => $toxin_info{$seq->id}->{'Endotoxin_C'},#2020/4/7
																			Endotoxin_mid     => $toxin_info{$seq->id}->{'Endotoxin_mid'},#2020/4/8
																			ETX_MTX2          => $toxin_info{$seq->id}->{'ETX_MTX2'},#2020/4/7
																			Toxin_10          => $toxin_info{$seq->id}->{'Toxin_10'},#2020/4/7
																			},
														);
				$seq1->add_SeqFeature($feat);
				$out->write_seq($seq1);
			}
			else  {
				$seq->accession_number($seq->id);
				my $feat = Bio::SeqFeature::Generic->new (
															-start => 1,
															-end   => $seq->length, 
															-primary   => $seq->id,
															-tag =>  {
																		protein_id        => $toxin_info{$seq->id}->{'protein_id'},
																		protein_desc      => $toxin_info{$seq->id}->{'protein_desc'},
																		protein_len       => $toxin_info{$seq->id}->{'protein_len'},
																		Rank              => $toxin_info{$seq->id}->{'rank'},
																		svm_prediction    => $toxin_info{$seq->id}->{'svm_prediction'},
																		blast_prediction  => $toxin_info{$seq->id}->{'blast_prediction'},
																		blast_detail      => $toxin_info{$seq->id}->{'blast_detail'},
																		hmm_prediction    => $toxin_info{$seq->id}->{'hmm_prediction'},
																		hmm_detail        => $toxin_info{$seq->id}->{'hmm_detail'},
																		Endotoxin_N       => $toxin_info{$seq->id}->{'Endotoxin_N'},#2020/4/7
																		Endotoxin_M       => $toxin_info{$seq->id}->{'Endotoxin_M'},#2020/4/7
																		Endotoxin_C       => $toxin_info{$seq->id}->{'Endotoxin_C'},#2020/4/7
																		Endotoxin_mid     => $toxin_info{$seq->id}->{'Endotoxin_mid'},#2020/4/8
																		ETX_MTX2          => $toxin_info{$seq->id}->{'ETX_MTX2'},#2020/4/7
																		Toxin_10          => $toxin_info{$seq->id}->{'Toxin_10'},#2020/4/7
																	},
														);
				$seq->add_SeqFeature($feat);
				$out->write_seq($seq);
			}
		}
	}

	system("rm $seqin");
	return @results;
}


sub print_toxin  {
	my ($ID, %toxin) = @_;
	print MYHAND "ID\tProtein_ID\tProtein_desc\tLength\tRank\tBLAST\tSVM\tHMM\n";
	foreach (sort keys %toxin)  {
		print MYHAND $ID, "\t";
		print MYHAND $toxin{$_}->{'protein_id'}, "\t";
		print MYHAND $toxin{$_}->{'protein_desc'}, "\t";
		print MYHAND $toxin{$_}->{'protein_len'}, "\t";
		print MYHAND $toxin{$_}->{'rank'}, "\t";
		print MYHAND $toxin{$_}->{'blast_prediction'}, "\t";
		print MYHAND $toxin{$_}->{'svm_prediction'}, "\t";
		print MYHAND $toxin{$_}->{'hmm_prediction'}, "\n"; 
		$ID++;
	}
	print MYHAND "//\n\n\n";
}

sub print_toxin2  {
	my ($ID, %toxin) = @_;
	print MYHAND "ID\tProtein_ID\tProtein_description\tLength\tRank\tBLAST\tBest_hit\tHit length\tCoverage\tIdentity\tSVM\tHMM\n";
	foreach (sort keys %toxin)  {
		print MYHAND $ID, "\t";
		print MYHAND $toxin{$_}->{'protein_id'}, "\t";
		print MYHAND $toxin{$_}->{'protein_desc'}, "\t";
		print MYHAND $toxin{$_}->{'protein_len'}, "\t";
		print MYHAND $toxin{$_}->{'rank'}, "\t";
		print MYHAND $toxin{$_}->{'blast_prediction'}, "\t";
		print  MYHAND $toxin{$_}->{'best_hit'}, "\t";
		print  MYHAND $toxin{$_}->{'hit_length'}, "\t";
		printf MYHAND ($toxin{$_}->{'coverage'}=~ /\d+/)? ("%0.2f\t", $toxin{$_}->{'coverage'}) : "$toxin{$_}->{'coverage'}\t";
		printf MYHAND ($toxin{$_}->{'Percent_identity'} =~ /\d+/)? ("%0.2f\t", $toxin{$_}->{'Percent_identity'}):"$toxin{$_}->{'Percent_identity'}\t";
		print MYHAND $toxin{$_}->{'svm_prediction'}, "\t";###Use of uninitialized value in print
		print MYHAND $toxin{$_}->{'hmm_prediction'}, "\n"; #2020/4/12
		#print MYHAND $toxin{$_}->{'Endotoxin_N'}, "\t";#2020/4/7,2020/4/12
		#print MYHAND $toxin{$_}->{'Endotoxin_M'}, "\t";#2020/4/7,2020/4/12
		#print MYHAND $toxin{$_}->{'Endotoxin_C'}, "\t";#2020/4/7,2020/4/12
		#print MYHAND $toxin{$_}->{'ETX_MTX2'}, "\t";#2020/4/7,2020/4/12
		#print MYHAND $toxin{$_}->{'Toxin_10'}, "\n";#2020/4/7,2020/4/12
		$ID++;
	}
	print MYHAND "//\n\n\n";
}
1;
