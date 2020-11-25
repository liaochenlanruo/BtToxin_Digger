# High-throughput _Bacillus thuringiensis_ toxin mining pipeline

![Platform](https://badgen.net/badge/platform/WSL,Linux,macOS,Docker?list=|)
![License](https://badgen.net/github/license/liaochenlanruo/BtToxin_Digger)
[![GitHubversion](https://badge.fury.io/gh/liaochenlanruo%2FBtToxin_Digger.svg)](https://badge.fury.io/gh/liaochenlanruo%2FBtToxin_Digger)
![Downloads conda](https://img.shields.io/conda/dn/bioconda/bttoxin_digger.svg?style=flat)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/bttoxin_digger/README.html)


	
			______    _ _____         _                 ______ _                       
			| ___ \  | |_   _|       (_)                |  _  (_)                      
			| |_/ / _| |_| | _____  ___ _ __            | | | |_  __ _  __ _  ___ _ __ 
			| ___ \|__ __| |/ _ \ \/ / | '_ \           | | | | |/ _` |/ _` |/ _ \ '__|
			| |_/ /  | |_| | (_) >  <| | | | |  ______  | |/ /| | (_| | (_| |  __/ |   
			\____/   \___\_/\___/_/\_\_|_| |_| |______| |___/ |_|\__, |\__, |\___|_|   
									      __/ | __/ |          
									     |___/ |___/           


<font face="STCAIYUN" color=#0099ff size=4>A web server of BtToxin_Digger can be found at [http://bcam.hzau.edu.cn/BtToxin_Digger](http://bcam.hzau.edu.cn/BtToxin_Digger).</font>

## Contents

- [Introduction](#introduction)

- [Installation](#installation)

- [Usage](#usage)

- [Examples](#examples)

- [Outputs](#outputs)

- [License](#license)

- [Feedback](#feedback)

- [Citation](#citation)

- [FAQs](#faqs)

- [Updates](#updates)

## Introduction

- **What is BtToxin_Digger?**

    BtToxin_Digger is a high-throughput, automatic gene mining tool that can mine toxin genes, such as Cry, Cyt and Vip, etc, from _Bacillus thuringiensis_. The pipeline accepts multiple forms of input data including Reads, assembled genomes, CDSs, and protein sequences and can output rich and useful results. It is derived from the re-design of the tool [BtToxin_Digger](http://bcam.hzau.edu.cn/BtToxin_Digger/) we developed previously. Compared with BtToxin_Digger, BtToxin_Digger has many improvements, as follows:
  - Can be run in batches, suitable for large-scale genome analysis.
  - Added genome assembly functions, including second-generation short-reads assembly, third-generation long-reads assembly, and hybrid assembly of short-reads and long-reads, to realize the full-automatic mining of genes from Reads to virulence factors. The previous three input files (assembled genomes, ORFs and protein sequences) are still supported, and genome assembly can be used independently.
  - Fixed a bug where BtToxin_Digger often reported errors when processing assembled genomes.
  - Added support for CDSs and not limited to ORFs.
  - The database has been updated, adding support for App, Gpp, Mcf, Mpf, Mpp, Mtx, Pra, Prb, Spp, Tpp, Cyt, Vip, Vpa, Vpb, Xpp and other virulence factors.
  - BtToxin_Digger generates comprehensive and readable outputs including toxin list and sequence for each input; a matrix of all strains and the virulence factors it contains (behavior strain names, listed as virulence factor names), which can be used as virulence factors contained in the strain Database; and a file writes the information and sequences of all toxins (Table 1) to facilitate centralized processing and downstream analysis and experiment designs.
  - Added multi-thread support, greatly improving the running speed of the pipeline.


## Installation

- Required dependencies
  - [BioPerl](http://metacpan.org/pod/BioPerl)
  - [HMMER](https://www.ebi.ac.uk/Tools/hmmer/)
  - [libsvm](https://github.com/cjlin1/libsvm)
  - [NCBI-blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
  - [Perl](http://www.perl.org/get.html)
  - [PGCGAP](https://liaochenlanruo.hzaubmb.org/pgcgap/)

- Install with Bioconda - OSX/Linux/WSL
```
conda create -n toxin python=3
conda activate toxin
conda install bttoxin_digger
```

- Install with the docker container
```
docker pull quay.io/biocontainers/bttoxin_digger:<tag>
```
  (See [bttoxin_digger/tags](https://quay.io/repository/biocontainers/bttoxin_digger?tab=tags) for valid values for \<tag\>)

## Usage
```
BtToxin_Digger [Options]
```

Options:

    [--help]                      Print the help message and exit

    [--version]                   Show version number of BtToxin_Digger and exit

    [--threads (INT)]             Number of threads to be used ( Default 4 )

    [--SeqPath (PATH)]            [Required] The path of input sequences ( Default "the current directory" )

    [--SequenceType (STRING)]     [Required] Sequence type for inputs. "reads", "nucl", "orfs", and "prot" avaliable ( Default nucl )

    [--platform (STRING)]         [Required] Sequencing Platform, "illumina", "pacbio", "oxford" and "hybrid" available ( Default illumina )

    [--assemble_only (STRING)]    Only perform genome assembly without predicting toxins.

    [--reads1 (STRING)]           [Required by "reads"] The suffix name of reads 1 ( for example: if the name of reads 1 is "YBT-1520_L1_I050.R1.clean.fastq.gz", "YBT-1520" is the strain same, so the suffix name should be ".R1.clean.fastq.gz" )

    [--reads2 (STRING)]           [Required by "reads"] The suffix name of reads 2( not required by "oxford" and "pacbio". For example: if the name of reads 2 is "YBT-1520_2.fq", the suffix name should be _2.fq" )

    [--suffix_len (INT)]          [Required by "reads"] (Strongly recommended) The suffix length of the reads file, that is the length of the reads name minus the length of the strain name. For example the --suffix_len of "YBT-1520_L1_I050.R1.clean.fastq.gz" is 26 ( "YBT-1520" is the strain name ) ( Default 0 )

    [--short1 (STRING)]           [Required] FASTQ file of first short reads in each pair. Needed by hybrid assembly ( Default Unset )

    [--short2 (STRING)]           [Required] FASTQ file of second short reads in each pair. Needed by hybrid assembly ( Default Unset )

    [--long (STRING)]             [Required] FASTQ or FASTA file of long reads. Needed by hybrid assembly ( Default Unset )

    [--hout (STRING)]             [Required] Output directory for hybrid assembly ( Default "../../Results/Assembles/Hybrid" )

    [--genomeSize (STRING)]       [Required] An estimate of the size of the genome. Common suffixes are allowed, for example, 3.7m or 2.8g. Needed by PacBio data and Oxford data ( Default 6.07m )

    [--Scaf_suffix (STRING)]      The suffix of scaffolds or genomes ( Default ".filtered.fas" )

    [--orfs_suffix (STRING)]      The suffix of orfs files ( Default ".ffn" )

    [--prot_suffix (STRING)]      The suffix of protein files ( Default ".faa" )

## Examples

- Processing Illumina paired-end Reads
```
BtToxin_Digger --SeqPath <Illumina Reads PATH> --SequenceType reads --platform illumina --reads1 <suffix name of reads 1> -reads2 <suffix name of reads 2> --threads <INT> --suffix_len <INT>
```

- Processing PacBio long Reads
```
BtToxin_Digger --SeqPath <PacBio Reads PATH> --SequenceType reads --platform pacbio --reads1 <suffix name of PacBio reads> --threads <INT> --suffix_len <INT>
```

- Processing Oxford long Reads
```
BtToxin_Digger --SeqPath <Oxford Reads PATH> --SequenceType reads --platform oxford --reads1 <suffix name of Oxford reads> --threads <INT> --suffix_len <INT>
```

- Processing Hybrid Reads (Long reads + illumina short reads)
```
BtToxin_Digger --SeqPath <Reads PATH> --SequenceType reads --platform hybrid --short1 <short reads 1> --short2 <short reads 2> --long <long reads> --threads <INT>
```

- Processing assembled genomes
```
BtToxin_Digger --SeqPath <Assembled genome PATH> --SequenceType nucl --Scaf_suffix <suffix of genomes> --threads <INT>
```

- Processing protein sequences
```
BtToxin_Digger --SeqPath <Protein file PATH> --SequenceType prot --prot_suffix <suffix of protein files> --threads <INT>
```

- Processing orfs sequences
```
BtToxin_Digger --SeqPath <orfs file PATH> --SequenceType orfs --orfs_suffix <suffix of orfs files> --threads <INT>
```

## Outputs

- __Results/Assembles/*:__ Genome assembly results;
- __Results/Toxins/*.list:__ Toxin list of each strain;
- __Results/Toxins/*.gbk:__ Toxin sequences in Genbank format of each strain;
- __Results/Toxins/Bt_all_genes.table:__ A matrix describes Strains vs. Toxins. See Table 3 for details;
- __Results/Toxins/All_Toxins.txt:__ A table containing all information and sequences of all toxin genes. See Table 1 and Table 2 for details.


__Contents of \'*.list\':__
```
Overview of prediction 	 Sequence type: nucl
Name	Cry	Cyt	Vip	Others	Summary
Rank1	1	0	0	0	1
Rank2	0	0	0	0	0
Rank3	0	0	0	0	0
Rank4	10	0	1	1	12
HMM_SVM	0	0	0	0
Summary	11	0	1	1	13
//


Toxin type: Cry protein
ID	Protein_ID		Protein_description	Length	Rank	BLAST	Best_hit	Hit length	Coverage	Identity	SVM	HMM
1	NHPK02000085.1_00018	+2 19034-20206 len=391	391	Rank1	YES	Cry78Aa1	359		97.21		32.66		NO	NO
2	NHPK02000099.1_00001	-1 5374-7407 len=678	678	Rank4	YES	Cry2Ab35	633		100.00		100.00		NO	YES
3	NHPK02000115.1_00002	+2 2-1342 len=447	447	Rank4	YES	Cry1Da2		1165		38.37		100.00		NO	YES
4	NHPK02000115.1_00006	+3 4275-7898 len=1208	1208	Rank4	YES	Cry1Ca15	1189		100.00		100.00		NO	YES
5	NHPK02000153.1_00001	-1 4176-6365 len=730	730	Rank4	YES	Cry1Ia40	719		100.00		100.00		YES	YES
6	NHPK02000168.1_00002	-3 1673-5146 len=1158	1158	Rank4	YES	Cry9Ea9		1150		100.00		100.00		YES	YES
7	NHPK02000196.1_00001	+3 723-2942 len=740	740	Rank4	YES	Cry1Da2		1165		63.43		100.00		NO	YES
8	NHPK02000240.1_00001	+2 2-1813 len=604	604	Rank4	YES	Cry1Ac16	1178		51.27		100.00		NO	YES
9	NHPK02000243.1_00001	+2 2-1738 len=579	579	Rank4	YES	Cry1Aa18	1225		47.27		100.00		NO	YES
10	NHPK02000256.1_00001	+2 2-1294 len=431	431	Rank4	YES	Cry1Ac30	1178		36.59		100.00		NO	YES
11	NHPK02000386.1_00001	-2 3-530 len=176	176	Rank4	YES	Cry1Ab11	695		25.32		99.43		NO	YES
//


Toxin type: Cyt protein
ID	Protein_ID	Protein_description	Length	Rank	BLAST	Best_hit	Hit length	Coverage	Identity	SVM	HMM
//


Toxin type: Vip protein
ID	Protein_ID		Protein_description	Length	Rank	BLAST	Best_hit	Hit length	Coverage	Identity	SVM	HMM
1	NHPK02000099.1_00004	-3 10013-12397 len=795	795	Rank4	YES	Vip3Aa12	789		100.00		100.00		NO	YES
//


Toxin type: Other toxins
ID	Protein_ID		Protein_description	Length	Rank	BLAST	Best_hit	Hit length	Coverage	Identity	SVM	HMM
1	NHPK02000017.1_00027	-3 15272-16258 len=329	329	Rank4	YES	Zwa5B-other	322		100.00		99.07		NO	NA
//

```

__Contents of \'*.gbk\' (Partial):__
```
LOCUS       NHPK02000017.1_00027        77664 bp    dna     linear   UNK
ACCESSION   NHPK02000017.1_00027
FEATURES             Location/Qualifiers
     Segment         complement(15272..16258)
                     /ETX_MTX2="NO"
                     /Endotoxin_C="NO"
                     /Endotoxin_M="NO"
                     /Endotoxin_N="NO"
                     /Endotoxin_mid="NO"
                     /Rank="Rank4"
                     /Segment_name="NHPK02000017.1_00027"
                     /Toxin_10="NO"
                     /blast_detail="Query_desc:-3 15272-16258 len=329
                     Query_Length:329	Query_Start-End:8-329	Hit_id:Zwa5B-other
                     Hit_desc:AAZ67352.1	Hit_length:322	Hit_Start-End:1-322
                     Aln_length:322	Percent_identity:99.0683229813665
                     E-value:0.0"
                     /blast_prediction="YES"
                     /hmm_detail=""
                     /hmm_prediction="NA"
                     /locus_tag="NHPK02000017.1_00027"
                     /protein_desc="-3 15272-16258 len=329"
                     /protein_id="NHPK02000017.1_00027"
                     /protein_len="329"
                     /svm_prediction="NO"
                     /translation="YRNEEAIMMYLNTKHINEMGVNWEETINVISKAVQSLDAEDFSQ
                     PIKPYLRFDDPANRIIAMPAYIGGEFKVSGIKWIASFPKNIEKGIQRAHSVTILNDAM
                     TGKPFATLNTAMVSVIRTASVTGLMIREFAKLRDLNNVKVGIIGFGPIGQMHLKMVTA
                     LLGDKIESVYLYDINGIKDELIPEDIYSKTQKVNAYEEAYNDADIFITCTVSAEGYID
                     KKPKDGALLLNVSLRDFKPDILEYTKSLVVDNWEEVCREKTDVERMHLERGLQKEDTV
                     SIADVVIRGALQNFPYDKAILFNPMGMAIFDVAIAAYYYQRAREKEIGVLLED"
ORIGIN      
        1 cttctaaccg caaggactgc ggggttagcc taatcaatta gagttcgata gaactcttta
       61 cttaggaatc cctcacttct aaacgaagtg aaagtggggt tagttcaaaa ctattaacta
      121 taatataccc tttcaagaaa ttttaaaaag gttgaagtag ccaaaaaata ttttccggaa
      181 taattagatt aatttctctt ttttgtatat ttatgttaaa ttattgttat cactacaaat
      241 ttattgaata attttaatac tctccccttt atactatggt aatatgtttt tcacaataca
      301 tattaccact ataattgcaa acatatataa acccatattt agaattttta agattctttc
      361 atagcattaa gatatttttt accttttagt ttgtttattc ttaattttta aactaaaata
      421 atatatattg gtaataatta aataaaattc caataattat aggaaggaga aaattatgaa
```

__Table 1: Description of \'All_Toxins.txt\'__

|Header|Description|
|:-----|:----------|
|Strain|The name of your input|
|Protein_id|The protein ID|
|Protein_len|The length of protein sequence|
|Strand|Positive or negative strand where the gene comes from|
|Gene location on scaffold|Gene coordinates on the genome|
|SVM|Is the protein predicted by SVM|
|BLAST|Is the protein predicted by BLAST|
|HMM|Is the protein predicted by HMM|
|Hit_id|The subject sequence ID|
|Hit_length|The length of subject sequence|
|Aln_length|alignment length (sequence overlap)|
|Query start-end|Start and end of alignment in query|
|Hit stard-end|Start and end of alignment in subject|
|Identity|Percentage of identical matches|
|Evalue of blast|Expect value of BLAST|
|Hmm hit|The subject model ID|
|Hmm hit length|The length of subject model sequence|
|Evalue of Hmm|Expect value of HMM|
|Nomenclature|[Bt nomenclature](http://www.lifesci.sussex.ac.uk/home/Neil_Crickmore/Bt/) containing 4 Ranks|
|Endotoxin_N|Whether the Cry protein contain Endotoxin_N domain|
|Endotoxin_M|Whether the Cry protein contain Endotoxin_M domain|
|Endotoxin_C|Whether the Cry protein contain Endotoxin_C domain|
|Endotoxin_mid|Whether the Cry protein contain Endotoxin_mid domain|
|Toxin_10|Whether the Cry protein contain Toxin_10 domain|
|ETX_MTX2|Whether the Cry protein contain ETX_MTX2 domain|
|Gene sequence|The nucleotide sequence of the toxin|
|Protein sequence|Amino acid sequence of the toxin|
|Scaffold sequence|The scaffold sequence where the toxin gene is located|

__Table 2. The contents of \'All_Toxins.txt\' (Partial)__

|Strain|Protein_id|Protein_len|Strand|Gene location on scaffold|SVM|BLAST|HMM|Hit_id|Hit_length|Aln_length|Query start-end|Hit stard-end|Identity|Evalue of blast|Hmm hit|Hmm hit length|Evalue of Hmm|Nomenclature|Endotoxin_N|Endotoxin_M|Endotoxin_C|Endotoxin_mid|Toxin_10|ETX_MTX2|Gene sequence|Protein sequence|Scaffold sequence|
|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
|1126_1|NHPK02000017.1_00027|329|-3|10001-10987|NO|YES|NA|Zwa5B-other|322|322|8-329|1-322|99.0683229813665|0.0|NA|NA|NA|Rank4|NO|NO|NO|NO|NO|NO|GTCTTCTAATAAAACACCTATTTCCTTTTCCCTAGCCCTTTGATAGTAATAAGCAGCAATGGCTACATCAAAAATGGCCATACCCATTGGATTAAATAATATTGCTTTGTCATAAGGGAAGTTTTGCAGTGCTCCTCGGATTACAACGTCAGCAATTGATACTGTATCTTCTTTTTGTAAACCTCTTTCAAGGTGCATCCTTTCAACATCAGTTTTTTCACGACAGACTTCTTCCCAATTATCAACAACTAGGGATTTTGTATATTCTAAGATATCAGGCTTAAAATCACGAAGAGAAACATTGAGTAGTAGCGCTCCATCTTTAGGTTTCTTATCAATATAACCTTCTGCAGAAACAGTACAAGTAATAAAAATATCAGCATCATTATAGGCTTCTTCATAAGCATTAACCTTTTGAGTTTTAGAATAAATATCCTCAGGAATTAATTCGTCTTTAATTCCATTGATATCATATAGGTAAACACTTTCGATTTTATCGCCTAATAGGGCAGTTACCATCTTTAAATGCATCTGACCAATAGGTCCAAATCCAATAATTCCAACTTTAACATTATTTAAATCTCTTAACTTGGCAAACTCTCGAATCATTAACCCTGTTACAGAAGCTGTTCGTATTACACTGACCATTGCTGTATTAAGCGTCGCAAATGGTTTTCCAGTCATAGCATCATTTAATATAGTAACTGAATGAGCACGTTGAATACCTTTCTCAATATTCTTCGGAAAACTAGCTATCCATTTGATTCCTGAAACCTTAAATTCTCCTCCAATATAAGCTGGCATAGCAATAATTCTATTAGCCGGATCATCAAAACGTAAGTATGGCTTAATTGGCTGAGAAAAATCTTCTGCATCAAGACTTTGTACTGCCTTTGAAATGACATTTATTGTTTCTTCCCAATTAACACCCATTTCATTGATATGTTTAGTATTTAAATACATCATAATTGCTTCCTCATTTCTATA|YRNEEAIMMYLNTKHINEMGVNWEETINVISKAVQSLDAEDFSQPIKPYLRFDDPANRIIAMPAYIGGEFKVSGIKWIASFPKNIEKGIQRAHSVTILNDAMTGKPFATLNTAMVSVIRTASVTGLMIREFAKLRDLNNVKVGIIGFGPIGQMHLKMVTALLGDKIESVYLYDINGIKDELIPEDIYSKTQKVNAYEEAYNDADIFITCTVSAEGYIDKKPKDGALLLNVSLRDFKPDILEYTKSLVVDNWEEVCREKTDVERMHLERGLQKEDTVSIADVVIRGALQNFPYDKAILFNPMGMAIFDVAIAAYYYQRAREKEIGVLLED|TCTCAATACGAACCAGTAAAAATGGTTCAGGGCGTTGAAGATACTGGACGTAATTCTTTTAGCCAAGCTAGTTTTTCAAACGTACTTCTAAGAACAGATACATCTGGAGCTACCTACAAAAGTTGGAACGGTTCAATTTCTTCCTCATTTTCTAAACACGGTACTTATGGTAACAGTTTCACAGTTTATTCTCAAGTGCCATTAAGCGCAAGTTTACGCTAAAATTGTTCTTGTATACGAATAATTTTGTGAAAAAGGTAGCTGTACTTTTGTACAGCTACCTTTTTCATATGGTGAATTTTTAATAGATATTGTATCAGACATAAGTTACATTCTTCATATTTTAACTTCTTTGAATTCAATAATACTTAATCAGCTTGTCTTATTACTTCACTTAACGTATCGAGATACAGCTTATAGCAACTACAACACGCTTTATTATTAATAATACCACTCTTCTTAATTTAACTTTAATTCCATTGCATTGATTTTTTAATAATAAATATACAACTTAGTAATGATAATAGATATACAATAAGAATAATATACAACTGTTTTGGTTCAATGCTATCGATATAAAGACTCATAATTGGATATGACCAAGGAAACAGTATCCAATAGTTTGTATGCGAGATAAAAAAGCTACCCAACCCGCCTATAACCCCTACTAATAATGAAAGTAACGGGATTTTATAAAATATTATATTTATTAAAAAATGTAAATTCACAATAGTTAAAGAGCTAATCACCATAAAAAGTATAGAAAAAATTAATCTTTGATAAACAGTCCCAAGATTATCAACCGTTAATACTATAGGTATACTAGTTATGACTAAAGCACATATTAAATAACCAAAAGTAATTATTAAGCAAAACAAAAACTTATTACTTAATATAGTTTTTAAAGTGATTCCGTGATAAATAGCAGACATATATCCTTTTCTTCTGAACTCGTTTACAAAAGAAAATGACACTAATAGACCTATTGCAATAGGTAAAAATAATCTACAAAAAAGACTTAGAATCGCTATATGCATGAAAACGCCATCTTTTGATTGAGATGATACAGTTAGAAAGATTTGAAAACTAATTGGTATCAACAGAAATACTAGTAACAGGATATACAATCCATTATGAAGACTTTTTTTATATTCACCTAATATACCCACTTATTTCACCCTTTATTATTCATATGTATAATATTTAAAACTCACAGCTAATAACATAATTCCTAAAATATAGATTGAAGAACTTGCAAGATATATTAAAATAGTTGTTTCTTTCATCATTACATCACTTACTATATTCAATGCATAAAATATTGGAAATGTCTTTGTGAATACTACAGGAAATAATCCATAAGTTACGCTAAAAAAAATCCAAATAACTGACAAGGACGCCGCAAGAATTCCCTTTTTTATAATCAAATGAAAAAGTAATTGAGAATATGGTATTGCTAAAGATAGAATTGTACAAAGCATCAACATTTTCAAAAAATAAAGTATAGATATTTCATATCCATTTATAAAAATATAGATACCCATAACTGAAAAATTTAAAATATTAAAGATTAAATAAAACAGTATAAACAATAAATTTTTAGCTACCAATAATTTAAATTTTGAAATTGGCTTTGTATAGATAATGGCCCAAGTACCATATTTATTTTCAGAAGAAATGATAAAATACCACATACAAGCAGAAAAAATGCTTAATGTTAGTGTATTAAAGTTTACTGTAATTGCTATAAGATTGGTAGATTGAAACATATCGGGATTATTCTTGTTAAAAATAAAATAACTTATACACAAAAAGCTTGATAATAATGAAAAAAACAAAGAAAATAACAATAATTTTCTCCCACCAATTTTCCTCATCTCACTATTCACACTTTTAAAGATAAATTGCATTATGATAGTAACTCCTTTTGTTTAGTTAAAGATAAAAATACTTCCTCTAATGTTAGCTTTTTAGGGCTAACTTCAAAAATATCTATATCATTTTCAATTAATATTCTGGTTAAATTAGGAACACTATTTTGATCTATCATTAGTTTGAATACCCCTTCTTTTTCTATGCTATAAGGTATATCCTTTGTATTTAAAATCCTTTCTGTGGCAGTAGTGTCATTAACATGAACAATATATTCTCTAGAAGTAATAATACTGTTTAAATTCCCTTCATAGAGTAAGTTACCCTTATGCAGTAATGCTATGTAATCTGCTATCATTTCGATTTCAGCTAAATTATGACTAGAAATAAATATTGTCTTCCCCTCATTTTTTACTAAATGAAGTAATAAGCGCCTTATCTCGTGTACTCCCTCTGGGTCTAATCCATTTGTCGGCTCATCAAGAATTAGTATTTTAGGATCATTAAGAATAGCAAAGGAAATACCTAGCCTTTGTTTCATTCCTAAAGAAAAATCTTTCACTTTTTTATCTTTTGCTTTATATATTCCTGTAATTTTTAATACTCTGTCTATTTCATCTAAAGGCTTATTAATCATTTTTTGTACATAGGCTAAATTTTCATATGCTGTCAAATTAGGATAGTAGGTAGGTGAATCAATAAAACAACCCACATTTTCAAATAATTCTTCTTTCCAATCTTTGATATCTTTATTATTAAATAATATCGTTCCTTTATCTAATGGAATTAATCCTAATAGCATTCGCATTGTTGTTGTTTTTCCTGCACCATTTGGTCCAATAAAAGCATATAAACTTCCTTCAGGAACAGATAAGTTTAAATCCCTTACAATGGATTGATTGTCAAAAGACTTTGTTAAATTTTTAGTTTGGATTGCTAATTTCATTTTGTTCCTCACTTTCTTATTCTCAAATAAAAGTAAGTATCTATTAAATTCTTGATCTATTAGATATTCTATTTACCCTTTACCTAATGTTTGCAATTTGGATTTAGAAAATATTTCTCGTTGTTTTAATTTTTTATCAGTAAGTTCATAAACTCTATCTACACGCTCTAAAATATTTTCATCATGAGTGATAATAATCTTAGTGGATGGTATGGCAAATATATTTTCAATAACAGCATTTGCCGTTTCCCTATCTAAATGGTTAGTTGGTTCATCCAAAATAAGAATCGGATATTCTTTCATGAACAATCGTAACATAGCCAACCTTTGTCTTTGGCCACCGGAGAAATTTGTTCCTCCCTCTGAGATTAAAATATCATTTATACTTTGAAATTGAATATAATTACTCAGATTCAATTTATCAATAGACTCCTGGATTTTTGAAGGATCATAATTTATATCAAAGTAAGATAAATTCTCAGATAAACTTCCTCTAAATAAAACTCCATCCTGAGTCATTAGTCCTACATTTTTTCTTAGTTCCTCACTGTTCCAATCTTTAATATCTCTATCGTCTATTAATATATTACCCTTATCGGGAGTATGCAATTTTAATAATAAGTTAGCTAAAGTGGACTTACCACTGCCACTTTCTCCAACAATAGCAATACTTTCGTTATTTTTAATGTTCAAACTTAAGTTTTCTATTATTTTTTTGTCATCAAAACTAAACGAAATATTATTAAATTGAATCGATCCTCTAAGAATGTTTACTTTAGTATCCGATCCATCTTTTTCTAGTTCCTCTTCCAAGATATCTAAAGTTCTTAATAAAAGTGGTTTAACTGAATTCCAATTTAAAATTGATCCTATAATTGTTGAAACCGGAGAAAATATAATACTAGATATAGCAATAAAGAAGAAAATATTACTTAAATCCCCTTCTCCTTTACTAAAATTATAAAATCCTACAGCAAGCATTATGATAATAAATATTGTTGATAGGGAACTATTAAATGATTCTAAAGCATTTTGTTTTAGCGCTCGTTTATGTACTGTATCTAAGTATTTATTATAATCTTTTTGCCATACAGATTTAATTCTAGAAAAAATCCCTGTAGCCTTTACATAATACATCGCTTTTAATGATTCATTAATTTTCCCTCTATGTTCTCCTAGGTAATACGTTTCGACTAAATTAGCATTATATATGAATTTCCCTAAGGAAATATTGATAAATAGGGAACTAGCAATAAAAAACATCAATATTAAAGTATAAAAAGTCTCAGTATATAATAGATACATCAATAATGGAATGATAACTACTATTGATATAATCATTTCAATTAGATCCTCTAAAATATATCCTCGTATTGATTCTAGGCTCATGATTCTTACATTTAAATCCCCTGCATGACGTGTAGTAATTTTTTCTATGGGGGTATTAAACATATGCTCAAAAACTTTTGTTGTACTTTTTCTATCAATTTCTTTCTGAAAGTTTACTTTTATAATAGAAGTTATAGATTTTATCCCTATAACAACGATTAAACTTCCCAAAATAATTAAAAGGTATCTTACATCCACTTCCTTTGAAGCATTCAAAATATCCATAAACTTTTTCAAAATAAAGGGGAAAGCAAGAATTGAGATTTGTGCTAACAATGCTAATAGTGCTAATAAAACAATCTGGCTAGGAGATACAAAATGATCAATATATCTAAGTAATATATGTTTTTCACGCCTTTTAGTATTTGGCTGTATATCTTTTTTGTTTTCTATTACTAGGATATAGTCACAGAACATTTCCTCAAAACTTCTTATATCAATTATCATCGGCCCTGCCTTAGGATCTACAATATAAAATTTATTTGCTTTAACCCTTTCAATAATTACAAAGTGATTAAAGCCCCAAAATGCTATCATCGGCTTTTGTATCCTTTGTAAATCAGTTATATTTTCTATCCTATATGCAGTTGCAGAAAGTCCATACGATTTACTAACTTTTTTAATATCCAACAAACTCCATGCCGATTTTTTATTGATAATTTTACTATTTCTAATTTCCGAAGCTTGTACTGCTGAACCATAATAATGCATAATCATTGCTAGACAACTAGGCCCACAATCAAAAGATGTCATTTGCGGTATAAATTTAATTTTCTTCCTTAACATATTTATCACCTATTCTAGGCATTAAAGTATAAAATGCTATTTAGTTTATTATCCTCTAACCTTTGTAAAAAATACAAAATACTTCCTAGCCCTGTAAATAATGCGATTTTATGAGGATTAAAATCTTCTATATCAAAATTATTCAAAAAGTTATTTGTTACTTTCATAACAATTTGTTTTTCATTGCGCTCGAGCATACCTCTATTTTCTAGCTCCAATAAAAAGTCTAGTTTCCCCATTAAACCATGACATAAACAATAATCACTATTTTCATGTAATGTACTTAAGATTTGTTTTTTAAAAAAAATTAAATCTGATTGTATATCCTCATTTTGATACGATTCTAAAATTTTAATTCTCGAGATTCCCATTCCACTAAATCCTTTACACCAGCTATCTGTATACTCATATTCACTCATACTCATTTTCTCACTTTTATGTATACCATGTATTACTCTATCGTACTTATTTGTTTTTAATATGCCTTTTAACTTATAAAATGCGAGAAGTTGTCCTGATAATCCATGAACAAACCCCTTTAATGAATCCTCTAAATTTAGAGTATCGTGAGCTTTCCAATATACAAATTCATTATGTGTTTGCATATTATTGACAATGTAGTCACCGAGATTTTCAGCCACTTTAAGCAATATTTTATAATTAGGATTAATTTGATAGAACTTAACTGCCGTTAAAATTAGACCAGCTGAACCCCCTAATATATCAAAATGCTTATCGCTATATTTACCTTTTTCAATATCAGAGTTGACATCCATAAATATTTGCAATATCAATCTTTCATCGATAGATGTTATATGTGGTAAATTTGTTAAATAGTAAGAAATAATTGCATTAACACCATTTACAAAACCATAATTATTTTTATTTTTGTTTAATATTTCTTGCAGGGCTTTTTGTTCAAGTTTATTTATAATTTCAATCTGATTAGGATCAGCATTATTTCTCAATAATGAACTTATTCCCAATAATCCATCATATATTCCACTATTCAAATGGTTTATCTCAAGCTCCCCTTGCCAATTTGTTTTTAATCCAATAAAATCTACAGCACCGTCTTTTGATATATAACTGCCATCAAGAATTTTCCCTGTAATTGCTTCCTTTATTTTATTTGTAGGTAAATTAGTTCGAAAATTAATTGTTTCCTTTACTGTATTCAGCTTATATAGATTGTCATTTTGATTAGAAAGGGATATTTCTATAATTTTTCTTTGTAATAAAAAGTCAGCATATGAAATTTTTTTAAATTTTGTTTCTAGATCTATTTCATTATTTGATAAACTTTCGCAAGACTCCCATTCATACTTATTAAAGAAATATGGAACATCATAGGATAATAGTGCCTCTATCTCTTTTTGAATAATATCTTTCCTAGTAGCCAAGTATATGTTTTTTTCAAGAATATTCCGCAAATATTTTATAGTCTTTAATTTATTCTCTAATAAATCCGGAACATTTAGTCTATCTAATAAAAGGGTATAGACTTTTGTATTTCGATATACTTTTCGATAATTCCATTTAGAAACACTATTTTTTAATATACACATATACTCCTCTTTGTTTTTTAAAACTGAGTTATATGCATTCTCAAAACCTTTAATAATATCATCTATATAATTATCATACTCGTATATATCATCGTTAAGAAACGGTAAATTCTGAACCGTATCTTCTACTATTGTCTCTTCTCTTATTTCAATTGGATTCGATGTAAATTCATTTTTTATTTTAATTGTTTTTACCAAAGACTTCATTACTCTCCCAATTGCACTGTAATCTCTTATTACCTCCCTACCTCCACTAAACCTTATAGGTAGCATTCTAGAATTTAATACTGAATTAGATAATATAAAACTACTGGTTTCTTCAGTAACACTTGGATTTAATGTACTTCCTAATGTCTCAAAATCAATTATGACTGGTGAGGATTGTTTAGCTATGATATTTTCATTATGTATGTCACTTCCATTAATCAAATACATTACCGAAAGCACTACACCTAAATTATAATAATAATCTTTAATTTGCTCATTAGATTCACATGATTCATGATTAATGTATTCCATATATCCATAATTTCTTTTATTAACTACTTGCGGAATAAAGATATCCGTATTTAAATCTTTATTCAGTTGATTAATTAAACTATGATATAGTATAGAATTTTTCAAATCGACTGGTTTTAAAATTACCTTATTATTAAATTTATCTTCTAAAAGAATTGTTCTCTTTCCATTATTATGGTAGTCACCTAGGTTAATAACTTTTTCTACATGCATAATATTAAAACCCGACATTTTATGAATTTCTTTTTGATCTTTTTCAAGAGCACCTTGCAACCAAATCATATAATTTATTTCATTTTCAATGATTGTAAAAATTTTTTCGTTTAATACGGGTAATATTTTAAAAAAATTCATATAAAATTTTTTAAAATGAAAGAAGTTCTTTGTAAAATACTTATATTCTTCATTGGTATCCATTCCTCTTAAAGCATTTTCTACCCTAAGTTTATTTATGTATACAATTAAGGATGGCCTTAATACTTGGTATAGTCTTTTTTCTAAATCAAACATTATCGATTCAGAAACCTTCTCAAAATTACCAAAAAAATCAATTTTATTTATTATCATATTTTTAAAGAATACTACATATGTTTGTATTATTTTTTGAAAGTTATAGTTTTCGTCATCCTTAATGACATCATTACATGAAAAACCCTCTAAAAAGATAACAAAATTTTTTTCTAGATATTCTTTTGTTACTAAAGAATGTGATCCATTATCTATTGTATTCAAAATTCCATGTGAATTATTTTTCATAATAACCTCTTCTTTTTATAGGTAGGAGAAACAGTTAGAAATTCATTCTCATTAAAGAGAGTAAGGCCTTTCGCCTTACTCTCTTTAACTATTAGTAGCATTTTTTACAAATTCTTTGTCCCATAATCGTGACACATGCTACAGGCACATTCCAAGGACAATCTTTAGTTAAAGTTTCAACCCAACCAGGTCCTGCTCCACCAACTAATTGTTCAATCTCTTCATCTTCTACGACTTGTAAATATTTTTCAGTTTCCATTTCAATACCTACCTTTTATTTTTTTGAAAATGATGTTATAAATTCAATAATTTACCACAACATTATCTAAAAACAATATACCATTAAATGTTAATAAAAAAAAGATATATTTATAGAAAAATAAAATAATAGGATATTTAACTATATTAATATAATTAAAGGAGACTCACTTACTATATCATTAAAATTTAAAGGCTATTTTCCAAGTTAATATTTTTAAAATTTTATAAATTTATTAATAATCAATGCTGTTTAATAATTTAAATAATTCTTCTAATACCTGAATAACGATTAAGAATATTAATTAAACTTGTAAAAAATAGATTTTCAGAAGATATAGCTCAAGAATATTTTACCATAAAATAGGATAATCAATTCTTTTTTTAATAGTATAATGCTTCGTCAGAAAATTTTAATCATGAATTGCAAAAAGCATGTTAAAAATAAAAAGATGAAATTGCAAGTATATTCACTTCACAATCCCATCTCTTCTAAAAATTGTATATCATCAACTCTCTTGATCCTCGCAATTCTACCAAATTAATAGCAACACATTATGTAGATATAGGTTTATACTATTCTACATTTAAATACATGTAATATTGCCTATACTTTGTTATTTTCATGTGAAACTAATTCACTATATATCTCTTTCTCTATATTTTTCAAGAGAATTTCCAACCTATTTTTTGTAGCTAGATTAAAATGTGGATTTAATTCTTTTATTATAAATTCAAGTATCTCAAAGCTTATGTAATCTCCTAATTCAAGTGAATATTGATCACCCAAGAAAGTTTTTAATTTTTTAGAGATCTCCACTTTATCTTTTTTGGTAATATTTACCCCATCCGATAACAAACCTCTCTTAGATTAGTCTTAAATACTAATTCAAAAAGTTTTTATAACATAATAAAAATATTTTGGAATTAATATAAATCATATAATGAAATTTTTACAAATAAATAAAAAATATACATCATAATCATTAAACAATATATTAAGCATCGTAAACATATATTGGCTCGTCATAAAGGCGTAAACATTTTGACCAATTTCCACATAACGTTTCCCATTCACCAGCAGGATATACTTCCGTACGCAATATATGAAGTGGTATACGTAAATTCAATAACCGACCACTAAGCGTCGTACGAGTAACTCCCGCTGGCATTAAAAACTTTGCTTCAATTACCTTTTTCAGTTCCATCAACGAGTATTCTGGATAACTCATCAATATTTCATTTGAATGAGGCCAGACCTCAGTATCATTTGAATGACGGCGTACCGTATACTCATGTTGATATAAGCTGACAATACGATGCCATATCTCTAAATGATTAAATAAATCGCTATTTAAATCCCTTGCATATAAAAAATGTACTTGTTTATTTGCTTCTGTAATTTTTGCAAGCAAACGTTTATTAGCAATAATCTTACCTGACCAAATTATTGTAGGTTCTTTCTTTAAATTATTTACCCATTCCCCCATTGGTAATAAATGATTCCAAGCACGAATTTTTATTTTTTCATCAGGTACTACCTGTACAGCAACGCGTTTACAACCTAATGCCATCATAGCCTTCATACGATGGGCACCATCAAGTAACAGATAATTACCATCACTAAGCTCAGATACTAAGGGAGGATTCCGTAAAACACCTTCGCTTTCCATTACACGACAAATTTGTTCCAATCTTTTGTGTTCATATGATTCATGAAAGCGAATTTGCTTGGGATGTAAGAGATCTAAGGCCGAAATAATTTCACTCATAAAAACAGCTCCTCAAATAATTACATAAATTCTAAGCTTATCAAATAACTAATTTTCTATAATGATTAGTCGTGTTTACCGATATTTAGCAGCTATTAGACAATTATTTCAAACTACGTTTAATAATTTTACTAGTCTTCTAATAAAACACCTATTTCCTTTTCCCTAGCCCTTTGATAGTAATAAGCAGCAATGGCTACATCAAAAATGGCCATACCCATTGGATTAAATAATATTGCTTTGTCATAAGGGAAGTTTTGCAGTGCTCCTCGGATTACAACGTCAGCAATTGATACTGTATCTTCTTTTTGTAAACCTCTTTCAAGGTGCATCCTTTCAACATCAGTTTTTTCACGACAGACTTCTTCCCAATTATCAACAACTAGGGATTTTGTATATTCTAAGATATCAGGCTTAAAATCACGAAGAGAAACATTGAGTAGTAGCGCTCCATCTTTAGGTTTCTTATCAATATAACCTTCTGCAGAAACAGTACAAGTAATAAAAATATCAGCATCATTATAGGCTTCTTCATAAGCATTAACCTTTTGAGTTTTAGAATAAATATCCTCAGGAATTAATTCGTCTTTAATTCCATTGATATCATATAGGTAAACACTTTCGATTTTATCGCCTAATAGGGCAGTTACCATCTTTAAATGCATCTGACCAATAGGTCCAAATCCAATAATTCCAACTTTAACATTATTTAAATCTCTTAACTTGGCAAACTCTCGAATCATTAACCCTGTTACAGAAGCTGTTCGTATTACACTGACCATTGCTGTATTAAGCGTCGCAAATGGTTTTCCAGTCATAGCATCATTTAATATAGTAACTGAATGAGCACGTTGAATACCTTTCTCAATATTCTTCGGAAAACTAGCTATCCATTTGATTCCTGAAACCTTAAATTCTCCTCCAATATAAGCTGGCATAGCAATAATTCTATTAGCCGGATCATCAAAACGTAAGTATGGCTTAATTGGCTGAGAAAAATCTTCTGCATCAAGACTTTGTACTGCCTTTGAAATGACATTTATTGTTTCTTCCCAATTAACACCCATTTCATTGATATGTTTAGTATTTAAATACATCATAATTGCTTCCTCATTTCTATATCATTTTTATAAAGAAACTAATTGGTCTTCTACTGATTTTTTTGTATTAAGCCATTCTACCCATTCTACGTTGTAAATAGTTGATGTGTAAGCTTGTCCATTATCAGGACACAAGAAAACAACATTAGGTGTATTTTGAACATCCCTATTTTCAAAATACTTTTGGATTGCATAATATGAAGTCCCTGATGATCCACCAGCAAATATTGCATGTTTATTAAATAGCTCATAACAACCTGCAACAGTATGGACCTCCGGTACAATCATTACATCATCAATCAATGCTTTTTTAACCATACCAGGTATCATACTGGCACCAATTCCTGGTATGTACCGTTTACGAGGCTTATCACCAAAAATAATGGACCCTTGACTATCTACCGCAATAATTTTAATATTAGGGAATTTTTCTTTTAATCGTGTAGATACCCCAGCGATAGTTCCTCCAGTACTCACTCCAATGAAGGCATAATCTAACTGTTTAAAATCATTAGATATTTCTTCACCAATTCCTTGATAATGGGCTTCAAAGTTATCAGCATTATTATATTGATTTGTCCAGTATGCATTAGGAATAGTATTTAAAAGTTCTTCCACCTTATTTAAACGCGTTAATAAATAGCCACCTGTTTCATCTCGTTCATCCACTTTAGCTACTTGATAGGAAGTTGCTCTCAAAAAATTCTCATAACTATCATTAATATTTGGATCAATAACTGGTATGAATTTCAGGCCGATATACCTACAAAGAGTAGCAAGAGCAACCGCAAAGTTACCAGATGAAGATTCGATAATTGTAGAATTTTCAGTCACTTCACCGCGACTTATTGCTGATTTTAAAATATGGTGAGCAGCACGCACTTTAACACTATTCATTAAATTGTGATATTCTAGTTTTGCATACAGATTAATCTTTTCATGCTCTAATTTAATCATAGGTGTGTTGCCGATTACCCTTTCTAAACTCTCTAATTTCTTGAGCATAATTTGCTTCCTTTCTTGGTTAAAATCTAAAAAAATAATTAAAGATAAATAGTTGTAAATGTTCTTTCGAATGTATTTCAAATAGAACTTATACCTGAAAGACAATAAATGCTAGATTCTTTAGTATCGATAAGCTAAGCACCATTCAAGTACTGCCATAGCACTATTTAGCTTGTTTTCGGCCTGGGTAAATACAATACTTTTCGGTCCATCAATAACTGTTGCTGCTACTTCTTCCCCCCGGTGTGCTGGGAGGTCATGCATGAAAGTAGCTTCCGGATATAAAGTCAAAATTCTCTCTGTCACTGCAAATGGATTAAAGGCATCTTTCCATAAAGGATCTATTTTGCTTGTCCCTGTTGTTTGCCACCTTGTTGTATATACAACATCCACTTCACTAGGTAACATTTTTATATCATGACATTCAATAACTTTAGCACCAGATAATCGTGCATATTCTTTTGATTTTTCAAGCACACTTGGAGACACTCCATAACCAGCAGGTGTAAAAAGATATAACTCCGTACCTGGAAAACGAGATAAAGAAAGAGCTAGAGCTGCAGCACTGTTATTACCTTCTCCCATATATAAAACCCGTAATCCAGCTATTTCACCAAATTTTTGTTTCATTGTTGTTAAGTCTGTTAGCGCTTGTGTTGGATGTTCTTCCGTACTCATTGCATTAATGACTGCCATACGACTTTGAGATGCCAGCTTTTCCATTTCCTTTTGACTATCTGCAGTTCTTGCTACTAGTGCATCCAGCATTCGCGAAAGTACTTGAGTAGTGTCTTCTATAGACTCCCCTGTATTTTCTTGCAAATCATTTGGACCATACGTTACAATTTGTGCACCCATCTTTAATGCTGCAACAGAAAATGCCGATCGAGTTCTAGTCGATGTTTTACGGAAGTAAATCCCAATTATATTTCCCGCTAATATTTGATTAGGTTGTGATTTACCTGTGGCAAACTCTACCCCACGAGTTACAATTTGGTTAATATCACTATCAGTAAGATCCTCAAGAGTAATTAAATGTTTTCGAATAGCCATTATTATATTACCTCCCTTTTGAATAATTGCATGTAACAAAAATGTATTTACATATGTATGTAAATAAAAACTACTAAACAGCTAAAGTAAAAATAATAGATTGGAAATATAATTCATGTGGATCCATTACCACACTTTTCTATAATTTTTTTTAAACTAATGTATAGCAGTTGATATATAGCCTATATGAACTCGAATACTATTTTTACTCAGAAATACTTACATATAATTCTATTGTTAGTCTATATTCCTATTTTTTATGATATAGCTCTGCTGAATATATTACACACTACAACACTAAAAATATACTCAAAATACAGAATTTTCAGCTTTCTCAATTCTTTCATCAATTTTCTGTCAATTCATACTCAAATTATCTTTTTTACTAATAATATACATACAATGTTGGCTGTTTCTTTTATAGCCCACTTCTAAAAAATCCACCCATTCAAAAATAACGTTTTCCGGATTGTCAGTTAAGATACAAATTACTGCTGCAGTCATTAAAGAGTCCTCTTTAAAAGAGATATACTCCTTCTACTGGCTCGCTATTTTTAGTTATGATAACCTCTTAATTGTTACAGCTATAACGCTAACAATTGTTATTGCAACCATTGATCTAAGTAAACAAACAATGAATGACATACTAATATATACACAAAGTAATCTATTAAGCGTAAAATGATTGCATCTCCAGCAAGTTTAATTACTTATTCAGATCTTTGCTTTTATAATTTTCAAAATTAAATTTCCGAAAATTTTCTTTCCTTAAGGAAATATCGATCATTTTAAACTATGTAAGAATGCTTATCTTCTATAAAGTAATTACACAAATCAACAAAGTCAATTTCAATTATTGATTCTGGAAAACTATCTAAAGAAGCACAACAAGAAAGCTTATAATCTTTACCAATAGTAAACTTTTTATAATAAACCTCATTAGATAGCATGTCACCATAAAATATTATATTTTTTTTAATAACAAAAGAGTTTAAAGGTATAGTTAGTCCTTTCCCTATCCATTTAACAAAGCTTTCTTTTATAGTCCATAATTCATAAAACACATCTAATTGTTGTTCTACATTTAATGAATTTAAATAATTAAATTCTTCCTCAGTAAATATATTTTTAATCATTGTAGTATCCATGGGCCGGACTTTTTCTACATCTATTCCTACTTCTTCCTCATGAATAGCTCCAACAACCCAACTTTCGGAATGCGATATATTAAAATGAAAATTATTTAACTCATCCACATAAGGCTTTCCGTATTCATTGTATTTATATTTAATATCTTTATTTTCTTTTGAGAAATTTGTAATAATTAAATATCTTATTATAATGTCTCCAATTAAGGCGCTATATTGGTCGGGCTTTCTTTTATATTTCTGTATTTTCATCTGCTTTTCTTTTGAGACCAAACTCATAAGTTTCTGCATAATGTTATGTTCTATATTTCCTGGCACACGTACCTTATATATATTCATTATTTCCTTCCCACCTAACTTCTTATATAAATATAAAAAACTCTCAATGAATAATTCCTATTTATTTTTGATTAGTAAAAAACATTTTCAGAGGATAATATTCTATTTTACATTCATATTTTCGACTTATTTTCTAACCTCTTTATGTCTATTATTGAATGGGCTCAAAGAAAAAAGTCTCATAATATGTTATACACGAATTTATGAGACTTTTAGAAAACTAATAGTGTACAATTTCAAATAATATTTCACCTATTTACATGTGCTAATAAGAGCATTTGTTTAAAACTCATTTTTATACTCTTTACCAACTCATGAATATGCGGTTTCTCTAACATACTCAGATGGGAACCTGGCACTTGCAAAGCGTATACTTCCCCTGATGTATATTCATTCCATCTGTTATAATCTACTAGTGGGTGAATGTCATTAATAGATGCATTAAATAGAAAGATATCTGCTTTTATTTTTTGTTTACAATTATATTTTAGGTATGCATATCTATTAGCTATCATTACTTTTAACTTATTCATCATTGGATCTTCAAAATTTTGCTGACAACTATTTTTATTCAATGTAAATTTTTTCAGTAAAGATTCTATTAATTGTTCTTCACTCATTTGCTCAAAACTAATTTTCTCAATCCCTAACTGATCATTGAATTTTTCCAATTCTTCAAAGGCATTCTTAATATTTAAACTAAGAATCTCCTTACCTTGTTCAATAGGATGAACATCAAGTAAACCTAGAAAACTAACTTTATCACCTAGTTCTTCTAATTTTCTAGCCATTTCAAATGCAACTATTCCTCCAAAAGACCATCCTAATAAGGCATATGGCCCTTCTTTTTTTACTTGTTTAATCTCTTCTATATATCTTACCGCCATTTCCTCTACAGATAAGTTTGAAAATCTATTATCATCATATCCTATAGATTGTAACCCATATACAGTCTTATCCTCTCCTAATTCTCTAGCCAAATCATAATAGTTTAATATGCCCCCTCCTTGTCCATGAACAATGAACCATTGACTATCCTTGTTTGTCCCATTTTGAATTGGAATTAAACACTCACTGTCTATTCCTTTATTCCTACTTATTACATCACTTAATTGTTCAATAGTAGCATTTTGAAATAATAAACTTAATGGCAGTTGCACATTAAACATTCTCTTGATATTTTCGAATAACTTTAATCCCTTAAGAGAATGACCACCTAATTCAAAGAAATTATCATTTATACCAATATTATTTACGCCTAATATACTACTCCAAATATCAATTAAACTACTATCTATGTCATTTCTTGGTGGTACATAATTACTATTGTCCAGTGTATTTAGCTTCGGTAACTTACTTCTATCTATTTTCCCATTTTGTGTTAATGGTATATTGTGAATTGGTATGAGTTGTTGAGGTATCATATAATGTGGTAACTTTGTTGCCAAATATGCTCTTACCTCTGGAATAGGAATATCTTTTTCTGTAACTACGTATGCACATAAATACTTTTCCCCAGCTTCATCTTCTTGATCTATTACTACGGCCGTTTTGATTGTCTCATATTTTAACAAACTCGCTTCTATCTCACCTAGTTCTATTCGATATCCCCTTATTTTTACTTGATGATCGACTCGGCCCAAGTATTCAATGTTACCATCAGGTAGCCACCTTGCTATATCCCCTGTTTTATATAGTTTCTCACCTCGTTCAAATGGATGATCAATAAATTTATCTGCTGTTAATTCTGTTCGATTGATGTATCCTCTAGCTAACCCTATGCCAGAGATACATATTTCACCTGGAACACCTATAGGTTGTATCCTATAAAATGAATCCAATATATAAATTTTAGTATTTAACAATGGTGAACCTATCGGTACGTTTTGTGTTGTAATTTCTTTATTACTCTCATATCGATAAGATGTACAACAAACTGTAGCCTCTGTAGGTCCGTATCCATTTAATATTTGGAGGTTCCCCCTAAATAGATGATCGTATTTTGCGAGTAATTCTGTTTTAATGGGTTCTACTCCCACAAGTAGCTTATTTAACACTATCTTCTGGTTATCTCTAACAAAATAATCATATATTTCATTTAATAAAGTAGGTGGAATATATGATAATGTAACCTGTTCTTCAAGAATAACTTGTACAAGTTTTGATACATCAAACTTCTCACCTTGATAAATAGTCATTCTTGCGCCATATATCAATGGGACAAATATCTCAAAAATAGTAACATCAAAAGAAATACTACTTGAGAATAGAACGTTATCAGTTATCCCTATATCTTGAGAGAAATCTTCATACATTGCACACAAAAAATTAGTCAAGGATCGATGTTCAATCATTACTCCTTTAGGTTGTCCTGTAGAACCTGATGTATAAATAACATACGCTAAGTTGTGAGGTTCTATCATCATTTGCATGTCTTCTCCTGGTTCTTCTTCAAAAGACATATCCATTAGATCTATTACATTACCTTGGAATTCTATCCCCTTTATAATAGAGTTTTGATGTACTAATACATGTGAACACCCACTGTCTGTCAGCATATATTCCACTCTTTGTTTTGGTAAGGCGGTATCAATTGGTAAATACGCCCCACCTGCTTTTAAAACACCTAATATACCTATAACCATCTCGATGGATCGTTCCATCATCACACCAACAATGGATTCTCTTTTAACTCCTTGGTCTAATAACCTTCTCGCCAACTGATTAGCTTTTATATTTAATTCATTATATGTAATTCCTTTTTCATTACATACAACGGCTATTTGATTAGGGTTCCGTTTTACTTGTTCTTCAAACATTTTATGCACTAATAAATGATTAGAATTTGAATTCTCTTTTTTGTTAAATTCATTCATAATACAATGTTCTTCTTCTATAGATAACATATTAATATTACGCAATCGTACTCTAGGGTTATTAGTTACTTCCTCAACTATATTTGTAAAATGTACCATTAACCTTTCTATTGTCTCCGCTTTAAATAACTTAGTACTATATTCTACTTTTAAATGAATGTTATTATCTATTTCCGTTGCTACTAATGACAAATCAAATTTTGAAACTGACTGCTTAAACGGATACGGTGTAAATTCTAATTCACCAATAGATATTGGATTCATATCCATGTTTTGAAAAACAAACATGGTATCAAATAATGGATTTCTACTTGTATCCCTATGCAAGTCTAAACCTTCTAATAGTTCTTCAAAAGGATAGTCTTGATTTTCATAAGCTTCTAATGTATTGAGTTTTAATCTACTCAAAAACTCAATAAACTCATCGTCATTTTCTAGATAATTTCTCATTACCAAAGTATTAATAAACATACCAATCATATGATTAGTATCAGAATGAGACCTTCCAGCAATAGGCGAACCTACAATAATGTCTTCTTGACCTGTGTATCTAGATAAAAGTATGTTATAAATGGCCAATAAAATCATATATGGCGTAGTACCAGTTTCAGTTGCTAGCTTATTTACTTTAAAAGTTAAATCTCTTCCCAAATTAAAAGAACAAACATTACCTTTAAAGCTTTGTATAGTCGGTCTTTGAAAATCGGTTGGAAAATTTAAAACCGGGAGTTCTCCTTTTAAAGTTGTTAACCAATAATTCTTTTGTTCACTAATTAGATTCTTATAGTAAGGTCCATTTTGCCACATCACATAGTCTTTGTACTGCACTCTCAATTTTGGAAGCTCATTTCCTTTATACAATTCTACAAACTCTTTTATTAATATCCCCATTGATAAACCATCAGATATTATATGATGCATATCTACTACAAGGATATGTCTTTCTTCTGCTATCCTTAAAAGCAACACCCTTAATAATGGAGGTTTTGATAAATCAAATGGACTTATAAACTCATGTATTAAATAATCCGCATCTTTTTCATTTACATGAACGTATTCAATATTGAAATCTACATTAGGCTCAATTTTCTGCACTAATTCCCCATCTAAAATTTGAAAAGAAGTTCTTAATATTTCATGTCTCTCAATTAAAGATTGAAATATATTTTCAAACTTATCTTTACAAATATCCCCTTCTACTTTAAGTATTGTGGGCATATTATAAGTTGTATTTGTACCATCTTCGAATTGATCTACTATAAACATTCTCTTTTGTGACGTAGAAGCTAAATAATACTCTTGTTGTTTTACAGGTTCTATGGAAATATAATTACTTTTTTCCATTTCTAATATACATTTTGAAAAATCAACCAAAATAGGAAACTTAAATAGGGATTTAATAGATAATTGCACATTAAATTCTTTATTAACGATAGAAATTAGCCTAGCAGCCTTTAATGAATGACCACCTATCTCAAAAAAATTATCCCGTATTCCTACCCTTTGTATTCCTAATACATCTTTCCAAATCTCAACTAATTTTCTTTCTGTAGAATTTGTCGGCTCTAGATGACTAGATTTCAAATTATTTATAGGTTGAGGTAATTTTTTTCTATCTATTTTTCCGTTTTGTGTTAACGGTATGTTTTGAATGGATATGATTTGTTGAGGTATCATATAATATGGTAACTTGGTTGCTAAATATGCTCTTACCTCTGGAATAGGAATATCTTTTTCGGTAACTACATACGCACATAAATACTTTTCTCCACTTTCATCTTCTCGCTGTATTACTACGGCGGTTTTGATTGTCTCATATTTTAACAAACTTGCTTCTATCTCACCTAGTTCTATTCGATATCCCCTTATTTTCACTTGATGATCAACTCGTCCCAAGTATTCTATGTTACCATCAGGTAGCCATCTCGCTATATCTCCTGTTTTATATAATTTCTCACCATGTTCAAATGGATGATCAATAAATTTATCTGCTGTTAATTCTTTTCGATTGATGTATCCCCTAGCTAACCCTATGCCAGAAATACATATTTCTCCTGGAACACCTATAGGTTGTAGCCTGTGAAATGAATCCAATATATAAATTTTAGTATTTAACAATGGCGAACCTATTGGTACGTTTTGTATTGTAATTTCTTTATCCCTCTCATATTGATAAGATGTACAACAGACTGTAGCCTCTGTAGGTCCGTATAAATTTAATATTTGGAGGTTCCCCCTAAATAGATGATCGTATTTTGCAAGTAATTCTGTTTTAATTGGTTCTACTCCCACGAAAAGTTTATTTAACGATATCTTCTGATTCGCCCTTACAAAATAATCATATATTTCATTTAATAAGGTAGGTGGAATATATGCTAACGTGACCTGTTCCTCAAGAATGACTTGTACAAGTTTTGGTACATCAAACTTCTCACCTTGATAAATAGTCATTCTTGCTCCATATACCAATGGGACAAATATTTCAAATATAGTAACATCAAAAGAAATACTACTTGAGAATAGTACATTATCACTTATCCCTATATCTTGAGAGAAGTCCTCATACATTGCACACAAAAAATTAGTTAAGGAACGATGTTCAATCATTACCCCTTTGGGCTGGCCTGTAGAACCTGATGTATAAATAACATATGCCAAATTTTGAGGTTCCATCGTTATCTGCAAATCTTCTACTTGTTCTTCTTCAAAAGGAATATCCATTAGATTTATTACACTACCTTGAAATGCTACTCCCTTTATAATAGAGTTTTGATATGTTAGTACATGTGAACACCCGCTGTCTGTCAGCATATATTCCACTCTTTGTTTTGGTAACTCGGTATCAATCGGTAAATACGCCCCACCCGCTTTTAAGATACCTAATATGCCTACAATCATCTCAATAGAGCGTTCCATCATCACACCAACAATGGATTCTCTTTTAACTCCCTGGTCTAATAACCTTCTCGCCAACTGATTAGCTTTTATATTTAATTGTTTATATGTAATTTCTTTTCCATTACATACAATCGCTATTTGATTAGGATTCTGTTTTACCTGCTCTTCAAATAATTGAGGAGCTGTTACAGTCTCACAAAGTAATGTTTTATGAACTCTAGTCGTATGATTGAAATCAAATAATATTTGATTCTTTTCTGTTTTAGGCATTACATCTAGATCCATTGCAGATTTATTTGGATCTTTCATTAATATATCTAAAATATTATATAAATGATTAACTATTCTACTAACTAATCCTTCACTATATAAAATAGAATTGTAATCCACTTGAACCTTTAACTGTTCTTCATTTTTCATAAACCTAATAACCATATCAGAGTTTATTTTATCTGTACTCTCATAACAATGAATATCATCTAACATTACTATTGTATTTAATAAGGGAAGATTATTACTCTCTCCATCTAGACTTAATAATTGAGTAAGTTTATTAAAAGGAAAATGACAATGTTCATTTGACTCTAAAATTGTTTCTTTAATTTTATATATTATTTCTTTAAAATTATCTTCTTGATTTATTTGAGTACGTAATAATAAAAAGTTGTTTTGAAAAACCGTCTCTTCTTGACCTTGTTTAAA|

__Table 3. The contents of \'Bt_all_genes.table\'__

|-|Cry1Ac16|Cry1Ca15|Cry1Da2|Cry1Ia40|Cry2Aa10|Cry2Ab35|Cry78Aa1|Cry9Ea9|HMM_Cry_len_492|HMM_Cyt_len_531|HMM_Cyt_len_533|HMM_Cyt_len_615|Sip1A2-other|Vip1Aa2|Vip3Aa12|Zwa5B-other|
|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
|1126_1|100.00|100.00|100.00|100.00||100.00|32.66|100.00|||||||100.00|99.07|
|AFS094730|||||||36.14||1|1|||92.41|29.46|||
|AFS095482||||||||||1||1|||||
|AK47||||100.00|100.00|100.00|||||1||||100.00|100.00|

__Footnote:__ The float number represent the identity of blast search, and the integer number represent the number of toxins predicted by HMM and SVM.

## License

BtToxin_Digger is free software, licensed under [GPLv3](https://github.com/liaochenlanruo/BtToxin_Digger/blob/master/LICENSE).

## Feedback/Issues

Please report any issues about usage of the software to the [issues page](https://github.com/liaochenlanruo/BtToxin_Digger/issues).

## Citation

- If you use this software please cite: Liu H, Zheng J, Yu Y, Ye W, Peng D, Sun M. BtToxin_Digger: a comprehensive and high-throughput pipeline for mining toxin protein genes from _Bacillus thuringiensis_. _bioRxiv_, 2020. DOI: [10.1101/2020.05.26.114520](https://doi.org/10.1101/2020.05.26.114520).

- If you used the genome assembly function, please also cite: Liu H, Xin B, Zheng J, Zhong H, Yu Y, Peng D, Sun M. Build a bioinformatics analysis platform and apply it to routine analysis of microbial genomics and comparative genomics. _Protocol exchange_, 2020. DOI: [10.21203/rs.2.21224/v3+](https://dx.doi.org/10.21203/rs.2.21224/v3).

## FAQs

## Updates

- v1.0.2
  - Fixed a "Can not find path" error.

- v1.0.3
  - Fixed a bug of "get_all_info_nucl.pl", which can not get the gene location and strand information of some toxins.

- v1.0.4
  - Updated the database and models to support [the latest clasiffication of Bt toxins](https://www.bpprc.org).

- v1.0.5
  - The name of strains with no toxin found will be outputed into the file "Strains_without_toxins_found.txt".