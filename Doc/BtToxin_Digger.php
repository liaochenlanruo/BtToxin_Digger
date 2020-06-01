<html xmlns="http://www.w3.org/1999/xhtml">
<!--
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<link href="commen.css" type="text/css" rel="stylesheet" />
<title>BtToxin_Digger</title>
</head>
-->
<head>
	<meta http-equiv="content-type" content="text/html;charset=utf-8"/>
	<link href="BtToxin_scanner/css/pageStyle.css" rel="stylesheet" type="text/css" />
	<link href="BtToxin_scanner/css/Bttoxin2.css" rel="stylesheet" type="text/css" />
</head>

<title>BtToxin_Digger</title>

<body class="body">

<!--
              <!--Start top-->             
              <div class="top_head">
                    <div class="top_img"><img src="BtToxin_scanner/img/head.jpg" width="900px" height="80px"></div>
                    <div class="top_menu">
                          <ul>
                                <li><a href="./index.php">Home</a></li>
								<li><a href="./Database.php">Database</a></li>
								<li><a href="./Analysis.php">Analysis</a></li>
								<li><a href="./BioTools.php">BioTools</a></li>
								<li><a href="./Publication.php">Publication</a></li>
								<li><a href="./About.php" class="last">About</a></li>
                          </ul>
                    </div>
              </div>
              <!--End top-->
-->

<div class="page_content">
   <div class="main_content">
       <b class="b1"></b><b class="b2"></b><b class="b3"></b> 
	   <div class="main_top"><a href="BtToxin_Digger.php">Advanced</a></div>

    <div class="Introduction">
         <div class="Introduction">

			<h2 id="highthroughput_bacillusthuringiensis_toxinminingpipeline">BtToxin_Digger: High-throughput <em>Bacillus thuringiensis</em> toxin mining pipeline</h2>

			<div class="top2_menu">
				<ul>
					<li><p><a href="#introduction">Introduction</a></p></li>
					<li><p><a href="#installation">Installation</a></p></li>
					<li><p><a href="#usage">Usage</a></p></li>
					<li><p><a href="#examples">Examples</a></p></li>
					<li><p><a href="#outputs">Outputs</a></p></li>
					<!--<li><p><a href="#license">License</a></p></li>-->
					<li><p><a href="#feedback">Feedback</a></p></li>
					<li><p><a href="#citation">Citation</a></p></li>
					<li><p><a href="#faqs">FAQs</a></p></li>
					<li><p><a href="#updates">Updates</a></p></li>
				</ul>
			</div>

			<br/><br/><br/><h3 id="introduction">Introduction</h3>

			<ul>
			<li><p><strong>What is BtToxin_Digger?</strong></p>

			<p>BtToxin_Digger is a high-throughput, automatic gene mining tool that can mine toxin genes, such as Cry, Cyt and Vip, etc, from <em>Bacillus thuringiensis</em>. The pipeline accepts multiple forms of input data including Reads, assembled genomes, CDSs, and protein sequences and can output rich and useful results. It is derived from the re-design of the tool <a href="http://bcam.hzau.edu.cn/BtToxin_scanner/" target="blank">BtToxin_Scanner</a> we developed previously. Compared with BtToxin_Scanner, BtToxin_Digger has many improvements, as follows:</p>

			<ul>
			<li>Can be run in batches, suitable for large-scale genome analysis.</li>

			<li>Added genome assembly functions, including second-generation short-reads assembly, third-generation long-reads assembly, and hybrid assembly of short-reads and long-reads, to realize the full-automatic mining of genes from Reads to virulence factors. The previous three input files (assembled genomes, ORFs and protein sequences) are still supported, and genome assembly can be used independently.</li>

			<li>Fixed a bug where BtToxin_Scanner often reported errors when processing assembled genomes.</li>

			<li>Added support for CDSs and not limited to ORFs.</li>

			<li>The database has been updated, adding support for Cyt, Vip and other virulence factors.</li>

			<li>BtToxin_Digger generates comprehensive and readable outputs including toxin list and sequence for each input; a matrix of all strains and the virulence factors it contains (behavior strain names, listed as virulence factor names), which can be used as virulence factors contained in the strain Database; and a file writes the information and sequences of all toxins (Table 1) to facilitate centralized processing and downstream analysis and experiment designs.</li>

			<li>Added multi-thread support, greatly improving the running speed of the pipeline.</li></ul></li>

			</ul>

			<h3 id="installation">Installation</h3>

			<ul>
			<li><p>Required dependencies</p>

			<ul>
			<li><a href="http://metacpan.org/pod/BioPerl">BioPerl</a></li>

			<li><a href="https://www.ebi.ac.uk/Tools/hmmer/">HMMER</a></li>

			<li><a href="https://github.com/cjlin1/libsvm">libsvm</a></li>

			<li><a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&amp;PAGE_TYPE=BlastDocs&amp;DOC_TYPE=Download">NCBI-blast+</a></li>

			<li><a href="http://www.perl.org/get.html">Perl</a></li>

			<li><a href="https://liaochenlanruo.hzaubmb.org/pgcgap/">PGCGAP</a></li></ul></li>

			<li><p>Source codes</p>

			<p>The BtToxin_Digger codes can be downloaded <a href="http://bcam.hzau.edu.cn/BtToxin_Digger-latest.tar.gz">Here</a> or from <a href="https://github.com/liaochenlanruo/BtToxin_Digger" target="blank">GitHub</a>.</p></li>

			<li><p>Install with Bioconda - OSX/Linux/WSL</p>

			<pre style="color:black;background-color:#eeeeee;line-height: 22px;font-size: 14px;overflow-x:scroll;padding-left: 5px;padding-right: 5px;"><code>conda create -n toxin python=3<br/>conda activate toxin<br/>conda install bttoxin_digger</code></pre>
			</li>
			</ul>

			<ul>
			<li>Install with the docker container

			<pre style="color:black;background-color:#eeeeee;line-height: 22px;font-size: 14px;overflow-x:scroll;padding-left: 5px;padding-right: 5px;"><code>docker pull quay.io/biocontainers/bttoxin_digger:&lt;tag&gt;</code></pre>

			<p>(See <a href="https://quay.io/repository/biocontainers/bttoxin_digger?tab=tags">bttoxin_digger/tags</a> for valid values for &lt;tag&gt;)</p>
			</li>
			</ul>

			<h3 id="usage">Usage</h3>

			<ul>
			<pre style="color:black;background-color:#eeeeee;line-height: 22px;font-size: 14px;overflow-x:scroll;padding-left: 5px;padding-right: 5px;"><code>BtToxin_Digger [Options]</code></pre>


			<p>Options:</p>

			<pre style="color:black;background-color:#eeeeee;line-height: 22px;font-size: 14px;overflow-x:scroll;padding-left: 5px;padding-right: 5px;"><code>[--help]                      Print the help message and exit<br/>[--version]                   Show version number of BtToxin_Digger and exit<br/>[--threads (INT)]             Number of threads to be used ( Default 4 )<br/>[--SeqPath (PATH)]            [Required] The path of input sequences ( Default "the current directory" )<br/>[--SequenceType (STRING)]     [Required] Sequence type for inputs. "reads", "nucl", "orfs", and "prot" avaliable ( Default nucl )<br/>[--platform (STRING)]         [Required] Sequencing Platform, "illumina", "pacbio", "oxford" and "hybrid" available ( Default illumina )<br/>[--assemble_only (STRING)]    Only perform genome assembly without predicting toxins.<br/>[--reads1 (STRING)]           [Required by "reads"] The suffix name of reads 1 ( for example: if the name of reads 1 is "YBT-1520_L1_I050.R1.clean.fastq.gz", "YBT-1520" is the strain same, so the suffix name should be ".R1.clean.fastq.gz" )<br/>[--reads2 (STRING)]           [Required by "reads"] The suffix name of reads 2( not required by "oxford" and "pacbio". For example: if the name of reads 2 is "YBT-1520_2.fq", the suffix name should be _2.fq" )<br/>[--suffix_len (INT)]          [Required by "reads"] (Strongly recommended) The suffix length of the reads file, that is the length of the reads name minus the length of the strain name. For example the --suffix_len of "YBT-1520_L1_I050.R1.clean.fastq.gz" is 26 ( "YBT-1520" is the strain name ) ( Default 0 )<br/>[--short1 (STRING)]           [Required] FASTQ file of first short reads in each pair. Needed by hybrid assembly ( Default Unset )<br/>[--short2 (STRING)]           [Required] FASTQ file of second short reads in each pair. Needed by hybrid assembly ( Default Unset )<br/>[--long (STRING)]             [Required] FASTQ or FASTA file of long reads. Needed by hybrid assembly ( Default Unset )<br/>[--hout (STRING)]             [Required] Output directory for hybrid assembly ( Default "../../Results/Assembles/Hybrid" )<br/>[--genomeSize (STRING)]       [Required] An estimate of the size of the genome. Common suffixes are allowed, for example, 3.7m or 2.8g. Needed by PacBio data and Oxford data ( Default 6.07m )<br/>[--Scaf_suffix (STRING)]      The suffix of scaffolds or genomes ( Default ".filtered.fas" )<br/>[--orfs_suffix (STRING)]      The suffix of orfs files ( Default ".ffn" )<br/>[--prot_suffix (STRING)]      The suffix of protein files ( Default ".faa" )</code></pre>
			</ul>

			<h3 id="examples">Examples</h3>

			<ul>
			<li>Processing Illumina paired-end Reads

			<pre style="color:black;background-color:#eeeeee;line-height: 22px;font-size: 14px;overflow-x:scroll;padding-left: 5px;padding-right: 5px;"><code>BtToxin_Digger --SeqPath &lt;Illumina Reads PATH&gt; --SequenceType reads --platform illumina --reads1 &lt;suffix name of reads 1&gt; -reads2 &lt;suffix name of reads 2&gt; --threads &lt;INT&gt; --suffix_len &lt;INT&gt;</code></pre>
			</li>
			</ul>

			<ul>
			<li>Processing PacBio long Reads

			<pre style="color:black;background-color:#eeeeee;line-height: 22px;font-size: 14px;overflow-x:scroll;padding-left: 5px;padding-right: 5px;"><code>BtToxin_Digger --SeqPath &lt;PacBio Reads PATH&gt; --SequenceType reads --platform pacbio --reads1 &lt;suffix name of PacBio reads&gt; --threads &lt;INT&gt; --suffix_len &lt;INT&gt;</code></pre>
			</li>
			</ul>

			<ul>
			<li>Processing Oxford long Reads

			<pre style="color:black;background-color:#eeeeee;line-height: 22px;font-size: 14px;overflow-x:scroll;padding-left: 5px;padding-right: 5px;"><code>BtToxin_Digger --SeqPath &lt;Oxford Reads PATH&gt; --SequenceType reads --platform oxford --reads1 &lt;suffix name of Oxford reads&gt; --threads &lt;INT&gt; --suffix_len &lt;INT&gt;</code></pre>
			</li>
			</ul>

			<ul>
			<li>Processing Hybrid Reads (Long reads + illumina short reads)

			<pre style="color:black;background-color:#eeeeee;line-height: 22px;font-size: 14px;overflow-x:scroll;padding-left: 5px;padding-right: 5px;"><code>BtToxin_Digger --SeqPath &lt;Reads PATH&gt; --SequenceType reads --platform hybrid --short1 &lt;short reads 1&gt; --short2 &lt;short reads 2&gt; --long &lt;long reads&gt; --threads &lt;INT&gt;</code></pre>
			</li>
			</ul>

			<ul>
			<li>Processing Assembled genomes

			<pre style="color:black;background-color:#eeeeee;line-height: 22px;font-size: 14px;overflow-x:scroll;padding-left: 5px;padding-right: 5px;"><code>BtToxin_Digger --SeqPath &lt;Assembled genome PATH&gt; --SequenceType nucl --Scaf_suffix &lt;suffix of genomes&gt; --threads &lt;INT&gt;</code></pre>
			</li>
			</ul>

			<ul>
			<li>Processing Protein sequences

			<pre style="color:black;background-color:#eeeeee;line-height: 22px;font-size: 14px;overflow-x:scroll;padding-left: 5px;padding-right: 5px;"><code>BtToxin_Digger --SeqPath &lt;Protein file PATH&gt; --SequenceType prot --prot_suffix &lt;suffix of protein files&gt; --threads &lt;INT&gt;</code></pre>
			</li>
			</ul>

			<ul>
			<li>Processing Coding sequences

			<pre style="color:black;background-color:#eeeeee;line-height: 22px;font-size: 14px;overflow-x:scroll;padding-left: 5px;padding-right: 5px;"><code>BtToxin_Digger --SeqPath &lt;CDSs file PATH&gt; --SequenceType orfs --orfs_suffix &lt;suffix of orfs files&gt; --threads &lt;INT&gt;</code></pre>
			</li>
			</ul>

			<h3 id="outputs">Outputs</h3>

			<ul>
			<li><strong>Results/Assembles/*:</strong> Genome assembly results;</li>

			<li><strong>Results/Toxins/*.list:</strong> Toxin list of each strain;</li>

			<li><strong>Results/Toxins/*.gbk:</strong> Toxin sequences in Genbank format of each strain;</li>

			<li><strong>Results/Toxins/Bt<em>all</em>genes.table:</strong> A matrix describes Strains vs. Toxins;</li>

			<li><strong>Results/Toxins/All_Toxins.txt:</strong> A table containing all information and sequences of all toxin genes. See table 1 for details.</li>
			</ul>

			<p><br><strong>Contents of *.list:</strong><br><pre style="color:black;background-color:white;line-height: 22px;font-size: 16px;overflow-x:scroll;padding-left: 5px;padding-right: 5px;"><code><figure class="highlight plain"><table><tr><td class="code"><pre><span class="line">Overview of prediction 	 Sequence type: nucl</span><br><span class="line">Name	Cry	Cyt	Vip	Others	Summary</span><br><span class="line">Rank1	1	0	0	0	1</span><br><span class="line">Rank2	0	0	0	0	0</span><br><span class="line">Rank3	0	0	0	0	0</span><br><span class="line">Rank4	10	0	1	1	12</span><br><span class="line">HMM_SVM	0	0	0	0</span><br><span class="line">Summary	11	0	1	1	13</span><br><span class="line">//</span><br><span class="line"></span><br><span class="line"></span><br><span class="line">Toxin type: Cry protein</span><br><span class="line">ID	Protein_ID		Protein_description	Length	Rank	BLAST	Best_hit	Hit length	Coverage	Identity	SVM	HMM</span><br><span class="line">1	NHPK02000085.1_00018	+2 19034-20206 len=391	391	Rank1	YES	Cry78Aa1	359		97.21		32.66		NO	NO</span><br><span class="line">2	NHPK02000099.1_00001	-1 5374-7407 len=678	678	Rank4	YES	Cry2Ab35	633		100.00		100.00		NO	YES</span><br><span class="line">3	NHPK02000115.1_00002	+2 2-1342 len=447	447	Rank4	YES	Cry1Da2		1165		38.37		100.00		NO	YES</span><br><span class="line">4	NHPK02000115.1_00006	+3 4275-7898 len=1208	1208	Rank4	YES	Cry1Ca15	1189		100.00		100.00		NO	YES</span><br><span class="line">5	NHPK02000153.1_00001	-1 4176-6365 len=730	730	Rank4	YES	Cry1Ia40	719		100.00		100.00		YES	YES</span><br><span class="line">6	NHPK02000168.1_00002	-3 1673-5146 len=1158	1158	Rank4	YES	Cry9Ea9		1150		100.00		100.00		YES	YES</span><br><span class="line">7	NHPK02000196.1_00001	+3 723-2942 len=740	740	Rank4	YES	Cry1Da2		1165		63.43		100.00		NO	YES</span><br><span class="line">8	NHPK02000240.1_00001	+2 2-1813 len=604	604	Rank4	YES	Cry1Ac16	1178		51.27		100.00		NO	YES</span><br><span class="line">9	NHPK02000243.1_00001	+2 2-1738 len=579	579	Rank4	YES	Cry1Aa18	1225		47.27		100.00		NO	YES</span><br><span class="line">10	NHPK02000256.1_00001	+2 2-1294 len=431	431	Rank4	YES	Cry1Ac30	1178		36.59		100.00		NO	YES</span><br><span class="line">11	NHPK02000386.1_00001	-2 3-530 len=176	176	Rank4	YES	Cry1Ab11	695		25.32		99.43		NO	YES</span><br><span class="line">//</span><br><span class="line"></span><br><span class="line"></span><br><span class="line">Toxin type: Cyt protein</span><br><span class="line">ID	Protein_ID	Protein_description	Length	Rank	BLAST	Best_hit	Hit length	Coverage	Identity	SVM	HMM</span><br><span class="line">//</span><br><span class="line"></span><br><span class="line"></span><br><span class="line">Toxin type: Vip protein</span><br><span class="line">ID	Protein_ID		Protein_description	Length	Rank	BLAST	Best_hit	Hit length	Coverage	Identity	SVM	HMM</span><br><span class="line">1	NHPK02000099.1_00004	-3 10013-12397 len=795	795	Rank4	YES	Vip3Aa12	789		100.00		100.00		NO	YES</span><br><span class="line">//</span><br><span class="line"></span><br><span class="line"></span><br><span class="line">Toxin type: Other toxins</span><br><span class="line">ID	Protein_ID		Protein_description	Length	Rank	BLAST	Best_hit	Hit length	Coverage	Identity	SVM	HMM</span><br><span class="line">1	NHPK02000017.1_00027	-3 15272-16258 len=329	329	Rank4	YES	Zwa5B-other	322		100.00		99.07		NO	NA</span><br><span class="line">//</span><br></pre></td></tr></table></figure></code></pre></p>

			<p><br><strong>Contents of *.gbk (Partial):</strong><br><figure class="highlight plain"><table><tr><td class="code"><pre><span class="line">LOCUS       NHPK02000017.1_00027        77664 bp    dna     linear   UNK</span><br><span class="line">ACCESSION   NHPK02000017.1_00027</span><br><span class="line">FEATURES             Location/Qualifiers</span><br><span class="line">     Segment         complement(15272..16258)</span><br><span class="line">                     /ETX_MTX2=&quot;NO&quot;</span><br><span class="line">                     /Endotoxin_C=&quot;NO&quot;</span><br><span class="line">                     /Endotoxin_M=&quot;NO&quot;</span><br><span class="line">                     /Endotoxin_N=&quot;NO&quot;</span><br><span class="line">                     /Endotoxin_mid=&quot;NO&quot;</span><br><span class="line">                     /Rank=&quot;Rank4&quot;</span><br><span class="line">                     /Segment_name=&quot;NHPK02000017.1_00027&quot;</span><br><span class="line">                     /Toxin_10=&quot;NO&quot;</span><br><span class="line">                     /blast_detail=&quot;Query_desc:-3 15272-16258 len=329</span><br><span class="line">                     Query_Length:329	Query_Start-End:8-329	Hit_id:Zwa5B-other</span><br><span class="line">                     Hit_desc:AAZ67352.1	Hit_length:322	Hit_Start-End:1-322</span><br><span class="line">                     Aln_length:322	Percent_identity:99.0683229813665</span><br><span class="line">                     E-value:0.0&quot;</span><br><span class="line">                     /blast_prediction=&quot;YES&quot;</span><br><span class="line">                     /hmm_detail=&quot;&quot;</span><br><span class="line">                     /hmm_prediction=&quot;NA&quot;</span><br><span class="line">                     /locus_tag=&quot;NHPK02000017.1_00027&quot;</span><br><span class="line">                     /protein_desc=&quot;-3 15272-16258 len=329&quot;</span><br><span class="line">                     /protein_id=&quot;NHPK02000017.1_00027&quot;</span><br><span class="line">                     /protein_len=&quot;329&quot;</span><br><span class="line">                     /svm_prediction=&quot;NO&quot;</span><br><span class="line">                     /translation=&quot;YRNEEAIMMYLNTKHINEMGVNWEETINVISKAVQSLDAEDFSQ</span><br><span class="line">                     PIKPYLRFDDPANRIIAMPAYIGGEFKVSGIKWIASFPKNIEKGIQRAHSVTILNDAM</span><br><span class="line">                     TGKPFATLNTAMVSVIRTASVTGLMIREFAKLRDLNNVKVGIIGFGPIGQMHLKMVTA</span><br><span class="line">                     LLGDKIESVYLYDINGIKDELIPEDIYSKTQKVNAYEEAYNDADIFITCTVSAEGYID</span><br><span class="line">                     KKPKDGALLLNVSLRDFKPDILEYTKSLVVDNWEEVCREKTDVERMHLERGLQKEDTV</span><br><span class="line">                     SIADVVIRGALQNFPYDKAILFNPMGMAIFDVAIAAYYYQRAREKEIGVLLED&quot;</span><br><span class="line">ORIGIN      </span><br><span class="line">        1 cttctaaccg caaggactgc ggggttagcc taatcaatta gagttcgata gaactcttta</span><br><span class="line">       61 cttaggaatc cctcacttct aaacgaagtg aaagtggggt tagttcaaaa ctattaacta</span><br><span class="line">      121 taatataccc tttcaagaaa ttttaaaaag gttgaagtag ccaaaaaata ttttccggaa</span><br><span class="line">      181 taattagatt aatttctctt ttttgtatat ttatgttaaa ttattgttat cactacaaat</span><br><span class="line">      241 ttattgaata attttaatac tctccccttt atactatggt aatatgtttt tcacaataca</span><br><span class="line">      301 tattaccact ataattgcaa acatatataa acccatattt agaattttta agattctttc</span><br><span class="line">      361 atagcattaa gatatttttt accttttagt ttgtttattc ttaattttta aactaaaata</span><br><span class="line">      421 atatatattg gtaataatta aataaaattc caataattat aggaaggaga aaattatgaa</span><br><span class="line"></span><br></pre></td></tr></table></figure></p>


			<p><br><strong>Table 1. Description of "All_Toxins.txt"</strong></p>
			<pre style="color:black;background-color:white;line-height: 22px;font-size: 16px;overflow-x:scroll;padding-left: 5px;padding-right: 5px;">
			<table border="2" bordercolor="green" width="600px">
			<thead>
			<tr>
			<th style="text-align:left">Header</th>
			<th style="text-align:left">Description</th>
			</tr>
			</thead>
			<tbody>
			<tr>
			<td style="text-align:left">Strain</td>
			<td style="text-align:left">The name of your input</td>
			</tr>
			<tr>
			<td style="text-align:left">Protein_id</td>
			<td style="text-align:left">The protein ID</td>
			</tr>
			<tr>
			<td style="text-align:left">Protein_len</td>
			<td style="text-align:left">The length of protein sequence</td>
			</tr>
			<tr>
			<td style="text-align:left">Strand</td>
			<td style="text-align:left">Positive or negative strand where the gene comes from</td>
			</tr>
			<tr>
			<td style="text-align:left">Gene location on scaffold</td>
			<td style="text-align:left">Gene coordinates on the genome</td>
			</tr>
			<tr>
			<td style="text-align:left">SVM</td>
			<td style="text-align:left">Is the protein predicted by SVM</td>
			</tr>
			<tr>
			<td style="text-align:left">BLAST</td>
			<td style="text-align:left">Is the protein predicted by BLAST</td>
			</tr>
			<tr>
			<td style="text-align:left">HMM</td>
			<td style="text-align:left">Is the protein predicted by HMM</td>
			</tr>
			<tr>
			<td style="text-align:left">Hit_id</td>
			<td style="text-align:left">The subject sequence ID</td>
			</tr>
			<tr>
			<td style="text-align:left">Hit_length</td>
			<td style="text-align:left">The length of subject sequence</td>
			</tr>
			<tr>
			<td style="text-align:left">Aln_length</td>
			<td style="text-align:left">alignment length (sequence overlap)</td>
			</tr>
			<tr>
			<td style="text-align:left">Query start-end</td>
			<td style="text-align:left">Start and end of alignment in query</td>
			</tr>
			<tr>
			<td style="text-align:left">Hit stard-end</td>
			<td style="text-align:left">Start and end of alignment in subject</td>
			</tr>
			<tr>
			<td style="text-align:left">Identity</td>
			<td style="text-align:left">Percentage of identical matches</td>
			</tr>
			<tr>
			<td style="text-align:left">Evalue of blast</td>
			<td style="text-align:left">Expect value of BLAST</td>
			</tr>
			<tr>
			<td style="text-align:left">Hmm hit</td>
			<td style="text-align:left">The subject model ID</td>
			</tr>
			<tr>
			<td style="text-align:left">Hmm hit length</td>
			<td style="text-align:left">The length of subject model sequence</td>
			</tr>
			<tr>
			<td style="text-align:left">Evalue of Hmm</td>
			<td style="text-align:left">Expect value of HMM</td>
			</tr>
			<tr>
			<td style="text-align:left">Nomenclature</td>
			<td style="text-align:left"><a href="http://www.lifesci.sussex.ac.uk/home/Neil_Crickmore/Bt/" target="_blank" rel="noopener">Bt nomenclature</a> containing 4 Ranks</td>
			</tr>
			<tr>
			<td style="text-align:left">Endotoxin_N</td>
			<td style="text-align:left">Whether the Cry protein contain Endotoxin_N domain</td>
			</tr>
			<tr>
			<td style="text-align:left">Endotoxin_M</td>
			<td style="text-align:left">Whether the Cry protein contain Endotoxin_M domain</td>
			</tr>
			<tr>
			<td style="text-align:left">Endotoxin_C</td>
			<td style="text-align:left">Whether the Cry protein contain Endotoxin_C domain</td>
			</tr>
			<tr>
			<td style="text-align:left">Endotoxin_mid</td>
			<td style="text-align:left">Whether the Cry protein contain Endotoxin_mid domain</td>
			</tr>
			<tr>
			<td style="text-align:left">Toxin_10</td>
			<td style="text-align:left">Whether the Cry protein contain Toxin_10 domain</td>
			</tr>
			<tr>
			<td style="text-align:left">ETX_MTX2</td>
			<td style="text-align:left">Whether the Cry protein contain ETX_MTX2 domain</td>
			</tr>
			<tr>
			<td style="text-align:left">Gene sequence</td>
			<td style="text-align:left">The nucleotide sequence of the toxin</td>
			</tr>
			<tr>
			<td style="text-align:left">Protein sequence</td>
			<td style="text-align:left">Amino acid sequence of the toxin</td>
			</tr>
			<tr>
			<td style="text-align:left">Scaffold sequence</td>
			<td style="text-align:left">The scaffold sequence where the toxin gene is located</td>
			</tr>
			</tbody>
			</table>
			</pre>

			<p><br><strong>Table 2. The contents of "All_Toxins.txt"  (Partial)</strong></p>
			<pre style="color:black;background-color:white;line-height: 22px;font-size: 14px;overflow-x:scroll;padding-left: 5px;padding-right: 5px;"><code>
			<table border="2" bordercolor="green">
			<thead>
			<tr>
			<th style="text-align:center">Strain</th>
			<th style="text-align:center">Protein_id</th>
			<th style="text-align:center">Protein_len</th>
			<th style="text-align:center">Strand</th>
			<th style="text-align:center">Gene location on scaffold</th>
			<th style="text-align:center">SVM</th>
			<th style="text-align:center">BLAST</th>
			<th style="text-align:center">HMM</th>
			<th style="text-align:center">Hit_id</th>
			<th style="text-align:center">Hit_length</th>
			<th style="text-align:center">Aln_length</th>
			<th style="text-align:center">Query start-end</th>
			<th style="text-align:center">Hit stard-end</th>
			<th style="text-align:center">Identity</th>
			<th style="text-align:center">Evalue of blast</th>
			<th style="text-align:center">Hmm hit</th>
			<th style="text-align:center">Hmm hit length</th>
			<th style="text-align:center">Evalue of Hmm</th>
			<th style="text-align:center">Nomenclature</th>
			<th style="text-align:center">Endotoxin_N</th>
			<th style="text-align:center">Endotoxin_M</th>
			<th style="text-align:center">Endotoxin_C</th>
			<th style="text-align:center">Endotoxin_mid</th>
			<th style="text-align:center">Toxin_10</th>
			<th style="text-align:center">ETX_MTX2</th>
			<th style="text-align:center">Gene sequence</th>
			<th style="text-align:center">Protein sequence</th>
			<th style="text-align:center">Scaffold sequence</th>
			</tr>
			</thead>
			<tbody>
			<tr>
			<td style="text-align:center">1126_1</td>
			<td style="text-align:center">NHPK02000017.1_00027</td>
			<td style="text-align:center">329</td>
			<td style="text-align:center">-3</td>
			<td style="text-align:center">10001-10987</td>
			<td style="text-align:center">NO</td>
			<td style="text-align:center">YES</td>
			<td style="text-align:center">NA</td>
			<td style="text-align:center">Zwa5B-other</td>
			<td style="text-align:center">322</td>
			<td style="text-align:center">322</td>
			<td style="text-align:center">8-329</td>
			<td style="text-align:center">1-322</td>
			<td style="text-align:center">99.0683229813665</td>
			<td style="text-align:center">0.0</td>
			<td style="text-align:center">NA</td>
			<td style="text-align:center">NA</td>
			<td style="text-align:center">NA</td>
			<td style="text-align:center">Rank4</td>
			<td style="text-align:center">NO</td>
			<td style="text-align:center">NO</td>
			<td style="text-align:center">NO</td>
			<td style="text-align:center">NO</td>
			<td style="text-align:center">NO</td>
			<td style="text-align:center">NO</td>
			<td style="text-align:center">GTCTTCTAATAAAACACCTATTTCCTTTTCCCTAGCCCTTTGATAGTAATAAGCAGCAATGGCTACATCAAAAATGGCCATACCCATTGGATTAAATAATATTGCTTTGTCATAAGGGAAGTTTTGCAGTGCTCCTCGGATTACAACGTCAGCAATTGATACTGTATCTTCTTTTTGTAAACCTCTTTCAAGGTGCATCCTTTCAACATCAGTTTTTTCACGACAGACTTCTTCCCAATTATCAACAACTAGGGATTTTGTATATTCTAAGATATCAGGCTTAAAATCACGAAGAGAAACATTGAGTAGTAGCGCTCCATCTTTAGGTTTCTTATCAATATAACCTTCTGCAGAAACAGTACAAGTAATAAAAATATCAGCATCATTATAGGCTTCTTCATAAGCATTAACCTTTTGAGTTTTAGAATAAATATCCTCAGGAATTAATTCGTCTTTAATTCCATTGATATCATATAGGTAAACACTTTCGATTTTATCGCCTAATAGGGCAGTTACCATCTTTAAATGCATCTGACCAATAGGTCCAAATCCAATAATTCCAACTTTAACATTATTTAAATCTCTTAACTTGGCAAACTCTCGAATCATTAACCCTGTTACAGAAGCTGTTCGTATTACACTGACCATTGCTGTATTAAGCGTCGCAAATGGTTTTCCAGTCATAGCATCATTTAATATAGTAACTGAATGAGCACGTTGAATACCTTTCTCAATATTCTTCGGAAAACTAGCTATCCATTTGATTCCTGAAACCTTAAATTCTCCTCCAATATAAGCTGGCATAGCAATAATTCTATTAGCCGGATCATCAAAACGTAAGTATGGCTTAATTGGCTGAGAAAAATCTTCTGCATCAAGACTTTGTACTGCCTTTGAAATGACATTTATTGTTTCTTCCCAATTAACACCCATTTCATTGATATGTTTAGTATTTAAATACATCATAATTGCTTCCTCATTTCTATA</td>
			<td style="text-align:center">YRNEEAIMMYLNTKHINEMGVNWEETINVISKAVQSLDAEDFSQPIKPYLRFDDPANRIIAMPAYIGGEFKVSGIKWIASFPKNIEKGIQRAHSVTILNDAMTGKPFATLNTAMVSVIRTASVTGLMIREFAKLRDLNNVKVGIIGFGPIGQMHLKMVTALLGDKIESVYLYDINGIKDELIPEDIYSKTQKVNAYEEAYNDADIFITCTVSAEGYIDKKPKDGALLLNVSLRDFKPDILEYTKSLVVDNWEEVCREKTDVERMHLERGLQKEDTVSIADVVIRGALQNFPYDKAILFNPMGMAIFDVAIAAYYYQRAREKEIGVLLED</td>
			<td style="text-align:center">TCTCAATACGAACCAGTAAAAATGGTTCAGGGCGTTGAAGATACTGGACGTAATTCTTTTAGCCAAGCTAGTTTTTCAAACGTACTTCTAAGAACAGATACATCTGGAGCTACCTACAAAAGTTGGAACGGTTCAATTTCTTCCTCATTTTCTAAACACGGTACTTATGGTAACAGTTTCACAGTTTATTCTCAAGTGCCATTAAGCGCAAGTTTACGCTAAAATTGTTCTTGTATACGAATAATTTTGTGAAAAAGGTAGCTGTACTTTTGTACAGCTACCTTTTTCATATGGTGAATTTTTAATAGATATTGTATCAGACATAAGTTACATTCTTCATATTTTAACTTCTTTGAATTCAATAATACTTAATCAGCTTGTCTTATTACTTCACTTAACGTATCGAGATACAGCTTATAGCAACTACAACACGCTTTATTATTAATAATACCACTCTTCTTAATTTAACTTTAATTCCATTGCATTGATTTTTTAATAATAAATATACAACTTAGTAATGATAATAGATATACAATAAGAATAATATACAACTGTTTTGGTTCAATGCTATCGATATAAAGACTCATAATTGGATATGACCAAGGAAACAGTATCCAATAGTTTGTATGCGAGATAAAAAAGCTACCCAACCCGCCTATAACCCCTACTAATAATGAAAGTAACGGGATTTTATAAAATATTATATTTATTAAAAAATGTAAATTCACAATAGTTAAAGAGCTAATCACCATAAAAAGTATAGAAAAAATTAATCTTTGATAAACAGTCCCAAGATTATCAACCGTTAATACTATAGGTATACTAGTTATGACTAAAGCACATATTAAATAACCAAAAGTAATTATTAAGCAAAACAAAAACTTATTACTTAATATAGTTTTTAAAGTGATTCCGTGATAAATAGCAGACATATATCCTTTTCTTCTGAACTCGTTTACAAAAGAAAATGACACTAATAGACCTATTGCAATAGGTAAAAATAATCTACAAAAAAGACTTAGAATCGCTATATGCATGAAAACGCCATCTTTTGATTGAGATGATACAGTTAGAAAGATTTGAAAACTAATTGGTATCAACAGAAATACTAGTAACAGGATATACAATCCATTATGAAGACTTTTTTTATATTCACCTAATATACCCACTTATTTCACCCTTTATTATTCATATGTATAATATTTAAAACTCACAGCTAATAACATAATTCCTAAAATATAGATTGAAGAACTTGCAAGATATATTAAAATAGTTGTTTCTTTCATCATTACATCACTTACTATATTCAATGCATAAAATATTGGAAATGTCTTTGTGAATACTACAGGAAATAATCCATAAGTTACGCTAAAAAAAATCCAAATAACTGACAAGGACGCCGCAAGAATTCCCTTTTTTATAATCAAATGAAAAAGTAATTGAGAATATGGTATTGCTAAAGATAGAATTGTACAAAGCATCAACATTTTCAAAAAATAAAGTATAGATATTTCATATCCATTTATAAAAATATAGATACCCATAACTGAAAAATTTAAAATATTAAAGATTAAATAAAACAGTATAAACAATAAATTTTTAGCTACCAATAATTTAAATTTTGAAATTGGCTTTGTATAGATAATGGCCCAAGTACCATATTTATTTTCAGAAGAAATGATAAAATACCACATACAAGCAGAAAAAATGCTTAATGTTAGTGTATTAAAGTTTACTGTAATTGCTATAAGATTGGTAGATTGAAACATATCGGGATTATTCTTGTTAAAAATAAAATAACTTATACACAAAAAGCTTGATAATAATGAAAAAAACAAAGAAAATAACAATAATTTTCTCCCACCAATTTTCCTCATCTCACTATTCACACTTTTAAAGATAAATTGCATTATGATAGTAACTCCTTTTGTTTAGTTAAAGATAAAAATACTTCCTCTAATGTTAGCTTTTTAGGGCTAACTTCAAAAATATCTATATCATTTTCAATTAATATTCTGGTTAAATTAGGAACACTATTTTGATCTATCATTAGTTTGAATACCCCTTCTTTTTCTATGCTATAAGGTATATCCTTTGTATTTAAAATCCTTTCTGTGGCAGTAGTGTCATTAACATGAACAATATATTCTCTAGAAGTAATAATACTGTTTAAATTCCCTTCATAGAGTAAGTTACCCTTATGCAGTAATGCTATGTAATCTGCTATCATTTCGATTTCAGCTAAATTATGACTAGAAATAAATATTGTCTTCCCCTCATTTTTTACTAAATGAAGTAATAAGCGCCTTATCTCGTGTACTCCCTCTGGGTCTAATCCATTTGTCGGCTCATCAAGAATTAGTATTTTAGGATCATTAAGAATAGCAAAGGAAATACCTAGCCTTTGTTTCATTCCTAAAGAAAAATCTTTCACTTTTTTATCTTTTGCTTTATATATTCCTGTAATTTTTAATACTCTGTCTATTTCATCTAAAGGCTTATTAATCATTTTTTGTACATAGGCTAAATTTTCATATGCTGTCAAATTAGGATAGTAGGTAGGTGAATCAATAAAACAACCCACATTTTCAAATAATTCTTCTTTCCAATCTTTGATATCTTTATTATTAAATAATATCGTTCCTTTATCTAATGGAATTAATCCTAATAGCATTCGCATTGTTGTTGTTTTTCCTGCACCATTTGGTCCAATAAAAGCATATAAACTTCCTTCAGGAACAGATAAGTTTAAATCCCTTACAATGGATTGATTGTCAAAAGACTTTGTTAAATTTTTAGTTTGGATTGCTAATTTCATTTTGTTCCTCACTTTCTTATTCTCAAATAAAAGTAAGTATCTATTAAATTCTTGATCTATTAGATATTCTATTTACCCTTTACCTAATGTTTGCAATTTGGATTTAGAAAATATTTCTCGTTGTTTTAATTTTTTATCAGTAAGTTCATAAACTCTATCTACACGCTCTAAAATATTTTCATCATGAGTGATAATAATCTTAGTGGATGGTATGGCAAATATATTTTCAATAACAGCATTTGCCGTTTCCCTATCTAAATGGTTAGTTGGTTCATCCAAAATAAGAATCGGATATTCTTTCATGAACAATCGTAACATAGCCAACCTTTGTCTTTGGCCACCGGAGAAATTTGTTCCTCCCTCTGAGATTAAAATATCATTTATACTTTGAAATTGAATATAATTACTCAGATTCAATTTATCAATAGACTCCTGGATTTTTGAAGGATCATAATTTATATCAAAGTAAGATAAATTCTCAGATAAACTTCCTCTAAATAAAACTCCATCCTGAGTCATTAGTCCTACATTTTTTCTTAGTTCCTCACTGTTCCAATCTTTAATATCTCTATCGTCTATTAATATATTACCCTTATCGGGAGTATGCAATTTTAATAATAAGTTAGCTAAAGTGGACTTACCACTGCCACTTTCTCCAACAATAGCAATACTTTCGTTATTTTTAATGTTCAAACTTAAGTTTTCTATTATTTTTTTGTCATCAAAACTAAACGAAATATTATTAAATTGAATCGATCCTCTAAGAATGTTTACTTTAGTATCCGATCCATCTTTTTCTAGTTCCTCTTCCAAGATATCTAAAGTTCTTAATAAAAGTGGTTTAACTGAATTCCAATTTAAAATTGATCCTATAATTGTTGAAACCGGAGAAAATATAATACTAGATATAGCAATAAAGAAGAAAATATTACTTAAATCCCCTTCTCCTTTACTAAAATTATAAAATCCTACAGCAAGCATTATGATAATAAATATTGTTGATAGGGAACTATTAAATGATTCTAAAGCATTTTGTTTTAGCGCTCGTTTATGTACTGTATCTAAGTATTTATTATAATCTTTTTGCCATACAGATTTAATTCTAGAAAAAATCCCTGTAGCCTTTACATAATACATCGCTTTTAATGATTCATTAATTTTCCCTCTATGTTCTCCTAGGTAATACGTTTCGACTAAATTAGCATTATATATGAATTTCCCTAAGGAAATATTGATAAATAGGGAACTAGCAATAAAAAACATCAATATTAAAGTATAAAAAGTCTCAGTATATAATAGATACATCAATAATGGAATGATAACTACTATTGATATAATCATTTCAATTAGATCCTCTAAAATATATCCTCGTATTGATTCTAGGCTCATGATTCTTACATTTAAATCCCCTGCATGACGTGTAGTAATTTTTTCTATGGGGGTATTAAACATATGCTCAAAAACTTTTGTTGTACTTTTTCTATCAATTTCTTTCTGAAAGTTTACTTTTATAATAGAAGTTATAGATTTTATCCCTATAACAACGATTAAACTTCCCAAAATAATTAAAAGGTATCTTACATCCACTTCCTTTGAAGCATTCAAAATATCCATAAACTTTTTCAAAATAAAGGGGAAAGCAAGAATTGAGATTTGTGCTAACAATGCTAATAGTGCTAATAAAACAATCTGGCTAGGAGATACAAAATGATCAATATATCTAAGTAATATATGTTTTTCACGCCTTTTAGTATTTGGCTGTATATCTTTTTTGTTTTCTATTACTAGGATATAGTCACAGAACATTTCCTCAAAACTTCTTATATCAATTATCATCGGCCCTGCCTTAGGATCTACAATATAAAATTTATTTGCTTTAACCCTTTCAATAATTACAAAGTGATTAAAGCCCCAAAATGCTATCATCGGCTTTTGTATCCTTTGTAAATCAGTTATATTTTCTATCCTATATGCAGTTGCAGAAAGTCCATACGATTTACTAACTTTTTTAATATCCAACAAACTCCATGCCGATTTTTTATTGATAATTTTACTATTTCTAATTTCCGAAGCTTGTACTGCTGAACCATAATAATGCATAATCATTGCTAGACAACTAGGCCCACAATCAAAAGATGTCATTTGCGGTATAAATTTAATTTTCTTCCTTAACATATTTATCACCTATTCTAGGCATTAAAGTATAAAATGCTATTTAGTTTATTATCCTCTAACCTTTGTAAAAAATACAAAATACTTCCTAGCCCTGTAAATAATGCGATTTTATGAGGATTAAAATCTTCTATATCAAAATTATTCAAAAAGTTATTTGTTACTTTCATAACAATTTGTTTTTCATTGCGCTCGAGCATACCTCTATTTTCTAGCTCCAATAAAAAGTCTAGTTTCCCCATTAAACCATGACATAAACAATAATCACTATTTTCATGTAATGTACTTAAGATTTGTTTTTTAAAAAAAATTAAATCTGATTGTATATCCTCATTTTGATACGATTCTAAAATTTTAATTCTCGAGATTCCCATTCCACTAAATCCTTTACACCAGCTATCTGTATACTCATATTCACTCATACTCATTTTCTCACTTTTATGTATACCATGTATTACTCTATCGTACTTATTTGTTTTTAATATGCCTTTTAACTTATAAAATGCGAGAAGTTGTCCTGATAATCCATGAACAAACCCCTTTAATGAATCCTCTAAATTTAGAGTATCGTGAGCTTTCCAATATACAAATTCATTATGTGTTTGCATATTATTGACAATGTAGTCACCGAGATTTTCAGCCACTTTAAGCAATATTTTATAATTAGGATTAATTTGATAGAACTTAACTGCCGTTAAAATTAGACCAGCTGAACCCCCTAATATATCAAAATGCTTATCGCTATATTTACCTTTTTCAATATCAGAGTTGACATCCATAAATATTTGCAATATCAATCTTTCATCGATAGATGTTATATGTGGTAAATTTGTTAAATAGTAAGAAATAATTGCATTAACACCATTTACAAAACCATAATTATTTTTATTTTTGTTTAATATTTCTTGCAGGGCTTTTTGTTCAAGTTTATTTATAATTTCAATCTGATTAGGATCAGCATTATTTCTCAATAATGAACTTATTCCCAATAATCCATCATATATTCCACTATTCAAATGGTTTATCTCAAGCTCCCCTTGCCAATTTGTTTTTAATCCAATAAAATCTACAGCACCGTCTTTTGATATATAACTGCCATCAAGAATTTTCCCTGTAATTGCTTCCTTTATTTTATTTGTAGGTAAATTAGTTCGAAAATTAATTGTTTCCTTTACTGTATTCAGCTTATATAGATTGTCATTTTGATTAGAAAGGGATATTTCTATAATTTTTCTTTGTAATAAAAAGTCAGCATATGAAATTTTTTTAAATTTTGTTTCTAGATCTATTTCATTATTTGATAAACTTTCGCAAGACTCCCATTCATACTTATTAAAGAAATATGGAACATCATAGGATAATAGTGCCTCTATCTCTTTTTGAATAATATCTTTCCTAGTAGCCAAGTATATGTTTTTTTCAAGAATATTCCGCAAATATTTTATAGTCTTTAATTTATTCTCTAATAAATCCGGAACATTTAGTCTATCTAATAAAAGGGTATAGACTTTTGTATTTCGATATACTTTTCGATAATTCCATTTAGAAACACTATTTTTTAATATACACATATACTCCTCTTTGTTTTTTAAAACTGAGTTATATGCATTCTCAAAACCTTTAATAATATCATCTATATAATTATCATACTCGTATATATCATCGTTAAGAAACGGTAAATTCTGAACCGTATCTTCTACTATTGTCTCTTCTCTTATTTCAATTGGATTCGATGTAAATTCATTTTTTATTTTAATTGTTTTTACCAAAGACTTCATTACTCTCCCAATTGCACTGTAATCTCTTATTACCTCCCTACCTCCACTAAACCTTATAGGTAGCATTCTAGAATTTAATACTGAATTAGATAATATAAAACTACTGGTTTCTTCAGTAACACTTGGATTTAATGTACTTCCTAATGTCTCAAAATCAATTATGACTGGTGAGGATTGTTTAGCTATGATATTTTCATTATGTATGTCACTTCCATTAATCAAATACATTACCGAAAGCACTACACCTAAATTATAATAATAATCTTTAATTTGCTCATTAGATTCACATGATTCATGATTAATGTATTCCATATATCCATAATTTCTTTTATTAACTACTTGCGGAATAAAGATATCCGTATTTAAATCTTTATTCAGTTGATTAATTAAACTATGATATAGTATAGAATTTTTCAAATCGACTGGTTTTAAAATTACCTTATTATTAAATTTATCTTCTAAAAGAATTGTTCTCTTTCCATTATTATGGTAGTCACCTAGGTTAATAACTTTTTCTACATGCATAATATTAAAACCCGACATTTTATGAATTTCTTTTTGATCTTTTTCAAGAGCACCTTGCAACCAAATCATATAATTTATTTCATTTTCAATGATTGTAAAAATTTTTTCGTTTAATACGGGTAATATTTTAAAAAAATTCATATAAAATTTTTTAAAATGAAAGAAGTTCTTTGTAAAATACTTATATTCTTCATTGGTATCCATTCCTCTTAAAGCATTTTCTACCCTAAGTTTATTTATGTATACAATTAAGGATGGCCTTAATACTTGGTATAGTCTTTTTTCTAAATCAAACATTATCGATTCAGAAACCTTCTCAAAATTACCAAAAAAATCAATTTTATTTATTATCATATTTTTAAAGAATACTACATATGTTTGTATTATTTTTTGAAAGTTATAGTTTTCGTCATCCTTAATGACATCATTACATGAAAAACCCTCTAAAAAGATAACAAAATTTTTTTCTAGATATTCTTTTGTTACTAAAGAATGTGATCCATTATCTATTGTATTCAAAATTCCATGTGAATTATTTTTCATAATAACCTCTTCTTTTTATAGGTAGGAGAAACAGTTAGAAATTCATTCTCATTAAAGAGAGTAAGGCCTTTCGCCTTACTCTCTTTAACTATTAGTAGCATTTTTTACAAATTCTTTGTCCCATAATCGTGACACATGCTACAGGCACATTCCAAGGACAATCTTTAGTTAAAGTTTCAACCCAACCAGGTCCTGCTCCACCAACTAATTGTTCAATCTCTTCATCTTCTACGACTTGTAAATATTTTTCAGTTTCCATTTCAATACCTACCTTTTATTTTTTTGAAAATGATGTTATAAATTCAATAATTTACCACAACATTATCTAAAAACAATATACCATTAAATGTTAATAAAAAAAAGATATATTTATAGAAAAATAAAATAATAGGATATTTAACTATATTAATATAATTAAAGGAGACTCACTTACTATATCATTAAAATTTAAAGGCTATTTTCCAAGTTAATATTTTTAAAATTTTATAAATTTATTAATAATCAATGCTGTTTAATAATTTAAATAATTCTTCTAATACCTGAATAACGATTAAGAATATTAATTAAACTTGTAAAAAATAGATTTTCAGAAGATATAGCTCAAGAATATTTTACCATAAAATAGGATAATCAATTCTTTTTTTAATAGTATAATGCTTCGTCAGAAAATTTTAATCATGAATTGCAAAAAGCATGTTAAAAATAAAAAGATGAAATTGCAAGTATATTCACTTCACAATCCCATCTCTTCTAAAAATTGTATATCATCAACTCTCTTGATCCTCGCAATTCTACCAAATTAATAGCAACACATTATGTAGATATAGGTTTATACTATTCTACATTTAAATACATGTAATATTGCCTATACTTTGTTATTTTCATGTGAAACTAATTCACTATATATCTCTTTCTCTATATTTTTCAAGAGAATTTCCAACCTATTTTTTGTAGCTAGATTAAAATGTGGATTTAATTCTTTTATTATAAATTCAAGTATCTCAAAGCTTATGTAATCTCCTAATTCAAGTGAATATTGATCACCCAAGAAAGTTTTTAATTTTTTAGAGATCTCCACTTTATCTTTTTTGGTAATATTTACCCCATCCGATAACAAACCTCTCTTAGATTAGTCTTAAATACTAATTCAAAAAGTTTTTATAACATAATAAAAATATTTTGGAATTAATATAAATCATATAATGAAATTTTTACAAATAAATAAAAAATATACATCATAATCATTAAACAATATATTAAGCATCGTAAACATATATTGGCTCGTCATAAAGGCGTAAACATTTTGACCAATTTCCACATAACGTTTCCCATTCACCAGCAGGATATACTTCCGTACGCAATATATGAAGTGGTATACGTAAATTCAATAACCGACCACTAAGCGTCGTACGAGTAACTCCCGCTGGCATTAAAAACTTTGCTTCAATTACCTTTTTCAGTTCCATCAACGAGTATTCTGGATAACTCATCAATATTTCATTTGAATGAGGCCAGACCTCAGTATCATTTGAATGACGGCGTACCGTATACTCATGTTGATATAAGCTGACAATACGATGCCATATCTCTAAATGATTAAATAAATCGCTATTTAAATCCCTTGCATATAAAAAATGTACTTGTTTATTTGCTTCTGTAATTTTTGCAAGCAAACGTTTATTAGCAATAATCTTACCTGACCAAATTATTGTAGGTTCTTTCTTTAAATTATTTACCCATTCCCCCATTGGTAATAAATGATTCCAAGCACGAATTTTTATTTTTTCATCAGGTACTACCTGTACAGCAACGCGTTTACAACCTAATGCCATCATAGCCTTCATACGATGGGCACCATCAAGTAACAGATAATTACCATCACTAAGCTCAGATACTAAGGGAGGATTCCGTAAAACACCTTCGCTTTCCATTACACGACAAATTTGTTCCAATCTTTTGTGTTCATATGATTCATGAAAGCGAATTTGCTTGGGATGTAAGAGATCTAAGGCCGAAATAATTTCACTCATAAAAACAGCTCCTCAAATAATTACATAAATTCTAAGCTTATCAAATAACTAATTTTCTATAATGATTAGTCGTGTTTACCGATATTTAGCAGCTATTAGACAATTATTTCAAACTACGTTTAATAATTTTACTAGTCTTCTAATAAAACACCTATTTCCTTTTCCCTAGCCCTTTGATAGTAATAAGCAGCAATGGCTACATCAAAAATGGCCATACCCATTGGATTAAATAATATTGCTTTGTCATAAGGGAAGTTTTGCAGTGCTCCTCGGATTACAACGTCAGCAATTGATACTGTATCTTCTTTTTGTAAACCTCTTTCAAGGTGCATCCTTTCAACATCAGTTTTTTCACGACAGACTTCTTCCCAATTATCAACAACTAGGGATTTTGTATATTCTAAGATATCAGGCTTAAAATCACGAAGAGAAACATTGAGTAGTAGCGCTCCATCTTTAGGTTTCTTATCAATATAACCTTCTGCAGAAACAGTACAAGTAATAAAAATATCAGCATCATTATAGGCTTCTTCATAAGCATTAACCTTTTGAGTTTTAGAATAAATATCCTCAGGAATTAATTCGTCTTTAATTCCATTGATATCATATAGGTAAACACTTTCGATTTTATCGCCTAATAGGGCAGTTACCATCTTTAAATGCATCTGACCAATAGGTCCAAATCCAATAATTCCAACTTTAACATTATTTAAATCTCTTAACTTGGCAAACTCTCGAATCATTAACCCTGTTACAGAAGCTGTTCGTATTACACTGACCATTGCTGTATTAAGCGTCGCAAATGGTTTTCCAGTCATAGCATCATTTAATATAGTAACTGAATGAGCACGTTGAATACCTTTCTCAATATTCTTCGGAAAACTAGCTATCCATTTGATTCCTGAAACCTTAAATTCTCCTCCAATATAAGCTGGCATAGCAATAATTCTATTAGCCGGATCATCAAAACGTAAGTATGGCTTAATTGGCTGAGAAAAATCTTCTGCATCAAGACTTTGTACTGCCTTTGAAATGACATTTATTGTTTCTTCCCAATTAACACCCATTTCATTGATATGTTTAGTATTTAAATACATCATAATTGCTTCCTCATTTCTATATCATTTTTATAAAGAAACTAATTGGTCTTCTACTGATTTTTTTGTATTAAGCCATTCTACCCATTCTACGTTGTAAATAGTTGATGTGTAAGCTTGTCCATTATCAGGACACAAGAAAACAACATTAGGTGTATTTTGAACATCCCTATTTTCAAAATACTTTTGGATTGCATAATATGAAGTCCCTGATGATCCACCAGCAAATATTGCATGTTTATTAAATAGCTCATAACAACCTGCAACAGTATGGACCTCCGGTACAATCATTACATCATCAATCAATGCTTTTTTAACCATACCAGGTATCATACTGGCACCAATTCCTGGTATGTACCGTTTACGAGGCTTATCACCAAAAATAATGGACCCTTGACTATCTACCGCAATAATTTTAATATTAGGGAATTTTTCTTTTAATCGTGTAGATACCCCAGCGATAGTTCCTCCAGTACTCACTCCAATGAAGGCATAATCTAACTGTTTAAAATCATTAGATATTTCTTCACCAATTCCTTGATAATGGGCTTCAAAGTTATCAGCATTATTATATTGATTTGTCCAGTATGCATTAGGAATAGTATTTAAAAGTTCTTCCACCTTATTTAAACGCGTTAATAAATAGCCACCTGTTTCATCTCGTTCATCCACTTTAGCTACTTGATAGGAAGTTGCTCTCAAAAAATTCTCATAACTATCATTAATATTTGGATCAATAACTGGTATGAATTTCAGGCCGATATACCTACAAAGAGTAGCAAGAGCAACCGCAAAGTTACCAGATGAAGATTCGATAATTGTAGAATTTTCAGTCACTTCACCGCGACTTATTGCTGATTTTAAAATATGGTGAGCAGCACGCACTTTAACACTATTCATTAAATTGTGATATTCTAGTTTTGCATACAGATTAATCTTTTCATGCTCTAATTTAATCATAGGTGTGTTGCCGATTACCCTTTCTAAACTCTCTAATTTCTTGAGCATAATTTGCTTCCTTTCTTGGTTAAAATCTAAAAAAATAATTAAAGATAAATAGTTGTAAATGTTCTTTCGAATGTATTTCAAATAGAACTTATACCTGAAAGACAATAAATGCTAGATTCTTTAGTATCGATAAGCTAAGCACCATTCAAGTACTGCCATAGCACTATTTAGCTTGTTTTCGGCCTGGGTAAATACAATACTTTTCGGTCCATCAATAACTGTTGCTGCTACTTCTTCCCCCCGGTGTGCTGGGAGGTCATGCATGAAAGTAGCTTCCGGATATAAAGTCAAAATTCTCTCTGTCACTGCAAATGGATTAAAGGCATCTTTCCATAAAGGATCTATTTTGCTTGTCCCTGTTGTTTGCCACCTTGTTGTATATACAACATCCACTTCACTAGGTAACATTTTTATATCATGACATTCAATAACTTTAGCACCAGATAATCGTGCATATTCTTTTGATTTTTCAAGCACACTTGGAGACACTCCATAACCAGCAGGTGTAAAAAGATATAACTCCGTACCTGGAAAACGAGATAAAGAAAGAGCTAGAGCTGCAGCACTGTTATTACCTTCTCCCATATATAAAACCCGTAATCCAGCTATTTCACCAAATTTTTGTTTCATTGTTGTTAAGTCTGTTAGCGCTTGTGTTGGATGTTCTTCCGTACTCATTGCATTAATGACTGCCATACGACTTTGAGATGCCAGCTTTTCCATTTCCTTTTGACTATCTGCAGTTCTTGCTACTAGTGCATCCAGCATTCGCGAAAGTACTTGAGTAGTGTCTTCTATAGACTCCCCTGTATTTTCTTGCAAATCATTTGGACCATACGTTACAATTTGTGCACCCATCTTTAATGCTGCAACAGAAAATGCCGATCGAGTTCTAGTCGATGTTTTACGGAAGTAAATCCCAATTATATTTCCCGCTAATATTTGATTAGGTTGTGATTTACCTGTGGCAAACTCTACCCCACGAGTTACAATTTGGTTAATATCACTATCAGTAAGATCCTCAAGAGTAATTAAATGTTTTCGAATAGCCATTATTATATTACCTCCCTTTTGAATAATTGCATGTAACAAAAATGTATTTACATATGTATGTAAATAAAAACTACTAAACAGCTAAAGTAAAAATAATAGATTGGAAATATAATTCATGTGGATCCATTACCACACTTTTCTATAATTTTTTTTAAACTAATGTATAGCAGTTGATATATAGCCTATATGAACTCGAATACTATTTTTACTCAGAAATACTTACATATAATTCTATTGTTAGTCTATATTCCTATTTTTTATGATATAGCTCTGCTGAATATATTACACACTACAACACTAAAAATATACTCAAAATACAGAATTTTCAGCTTTCTCAATTCTTTCATCAATTTTCTGTCAATTCATACTCAAATTATCTTTTTTACTAATAATATACATACAATGTTGGCTGTTTCTTTTATAGCCCACTTCTAAAAAATCCACCCATTCAAAAATAACGTTTTCCGGATTGTCAGTTAAGATACAAATTACTGCTGCAGTCATTAAAGAGTCCTCTTTAAAAGAGATATACTCCTTCTACTGGCTCGCTATTTTTAGTTATGATAACCTCTTAATTGTTACAGCTATAACGCTAACAATTGTTATTGCAACCATTGATCTAAGTAAACAAACAATGAATGACATACTAATATATACACAAAGTAATCTATTAAGCGTAAAATGATTGCATCTCCAGCAAGTTTAATTACTTATTCAGATCTTTGCTTTTATAATTTTCAAAATTAAATTTCCGAAAATTTTCTTTCCTTAAGGAAATATCGATCATTTTAAACTATGTAAGAATGCTTATCTTCTATAAAGTAATTACACAAATCAACAAAGTCAATTTCAATTATTGATTCTGGAAAACTATCTAAAGAAGCACAACAAGAAAGCTTATAATCTTTACCAATAGTAAACTTTTTATAATAAACCTCATTAGATAGCATGTCACCATAAAATATTATATTTTTTTTAATAACAAAAGAGTTTAAAGGTATAGTTAGTCCTTTCCCTATCCATTTAACAAAGCTTTCTTTTATAGTCCATAATTCATAAAACACATCTAATTGTTGTTCTACATTTAATGAATTTAAATAATTAAATTCTTCCTCAGTAAATATATTTTTAATCATTGTAGTATCCATGGGCCGGACTTTTTCTACATCTATTCCTACTTCTTCCTCATGAATAGCTCCAACAACCCAACTTTCGGAATGCGATATATTAAAATGAAAATTATTTAACTCATCCACATAAGGCTTTCCGTATTCATTGTATTTATATTTAATATCTTTATTTTCTTTTGAGAAATTTGTAATAATTAAATATCTTATTATAATGTCTCCAATTAAGGCGCTATATTGGTCGGGCTTTCTTTTATATTTCTGTATTTTCATCTGCTTTTCTTTTGAGACCAAACTCATAAGTTTCTGCATAATGTTATGTTCTATATTTCCTGGCACACGTACCTTATATATATTCATTATTTCCTTCCCACCTAACTTCTTATATAAATATAAAAAACTCTCAATGAATAATTCCTATTTATTTTTGATTAGTAAAAAACATTTTCAGAGGATAATATTCTATTTTACATTCATATTTTCGACTTATTTTCTAACCTCTTTATGTCTATTATTGAATGGGCTCAAAGAAAAAAGTCTCATAATATGTTATACACGAATTTATGAGACTTTTAGAAAACTAATAGTGTACAATTTCAAATAATATTTCACCTATTTACATGTGCTAATAAGAGCATTTGTTTAAAACTCATTTTTATACTCTTTACCAACTCATGAATATGCGGTTTCTCTAACATACTCAGATGGGAACCTGGCACTTGCAAAGCGTATACTTCCCCTGATGTATATTCATTCCATCTGTTATAATCTACTAGTGGGTGAATGTCATTAATAGATGCATTAAATAGAAAGATATCTGCTTTTATTTTTTGTTTACAATTATATTTTAGGTATGCATATCTATTAGCTATCATTACTTTTAACTTATTCATCATTGGATCTTCAAAATTTTGCTGACAACTATTTTTATTCAATGTAAATTTTTTCAGTAAAGATTCTATTAATTGTTCTTCACTCATTTGCTCAAAACTAATTTTCTCAATCCCTAACTGATCATTGAATTTTTCCAATTCTTCAAAGGCATTCTTAATATTTAAACTAAGAATCTCCTTACCTTGTTCAATAGGATGAACATCAAGTAAACCTAGAAAACTAACTTTATCACCTAGTTCTTCTAATTTTCTAGCCATTTCAAATGCAACTATTCCTCCAAAAGACCATCCTAATAAGGCATATGGCCCTTCTTTTTTTACTTGTTTAATCTCTTCTATATATCTTACCGCCATTTCCTCTACAGATAAGTTTGAAAATCTATTATCATCATATCCTATAGATTGTAACCCATATACAGTCTTATCCTCTCCTAATTCTCTAGCCAAATCATAATAGTTTAATATGCCCCCTCCTTGTCCATGAACAATGAACCATTGACTATCCTTGTTTGTCCCATTTTGAATTGGAATTAAACACTCACTGTCTATTCCTTTATTCCTACTTATTACATCACTTAATTGTTCAATAGTAGCATTTTGAAATAATAAACTTAATGGCAGTTGCACATTAAACATTCTCTTGATATTTTCGAATAACTTTAATCCCTTAAGAGAATGACCACCTAATTCAAAGAAATTATCATTTATACCAATATTATTTACGCCTAATATACTACTCCAAATATCAATTAAACTACTATCTATGTCATTTCTTGGTGGTACATAATTACTATTGTCCAGTGTATTTAGCTTCGGTAACTTACTTCTATCTATTTTCCCATTTTGTGTTAATGGTATATTGTGAATTGGTATGAGTTGTTGAGGTATCATATAATGTGGTAACTTTGTTGCCAAATATGCTCTTACCTCTGGAATAGGAATATCTTTTTCTGTAACTACGTATGCACATAAATACTTTTCCCCAGCTTCATCTTCTTGATCTATTACTACGGCCGTTTTGATTGTCTCATATTTTAACAAACTCGCTTCTATCTCACCTAGTTCTATTCGATATCCCCTTATTTTTACTTGATGATCGACTCGGCCCAAGTATTCAATGTTACCATCAGGTAGCCACCTTGCTATATCCCCTGTTTTATATAGTTTCTCACCTCGTTCAAATGGATGATCAATAAATTTATCTGCTGTTAATTCTGTTCGATTGATGTATCCTCTAGCTAACCCTATGCCAGAGATACATATTTCACCTGGAACACCTATAGGTTGTATCCTATAAAATGAATCCAATATATAAATTTTAGTATTTAACAATGGTGAACCTATCGGTACGTTTTGTGTTGTAATTTCTTTATTACTCTCATATCGATAAGATGTACAACAAACTGTAGCCTCTGTAGGTCCGTATCCATTTAATATTTGGAGGTTCCCCCTAAATAGATGATCGTATTTTGCGAGTAATTCTGTTTTAATGGGTTCTACTCCCACAAGTAGCTTATTTAACACTATCTTCTGGTTATCTCTAACAAAATAATCATATATTTCATTTAATAAAGTAGGTGGAATATATGATAATGTAACCTGTTCTTCAAGAATAACTTGTACAAGTTTTGATACATCAAACTTCTCACCTTGATAAATAGTCATTCTTGCGCCATATATCAATGGGACAAATATCTCAAAAATAGTAACATCAAAAGAAATACTACTTGAGAATAGAACGTTATCAGTTATCCCTATATCTTGAGAGAAATCTTCATACATTGCACACAAAAAATTAGTCAAGGATCGATGTTCAATCATTACTCCTTTAGGTTGTCCTGTAGAACCTGATGTATAAATAACATACGCTAAGTTGTGAGGTTCTATCATCATTTGCATGTCTTCTCCTGGTTCTTCTTCAAAAGACATATCCATTAGATCTATTACATTACCTTGGAATTCTATCCCCTTTATAATAGAGTTTTGATGTACTAATACATGTGAACACCCACTGTCTGTCAGCATATATTCCACTCTTTGTTTTGGTAAGGCGGTATCAATTGGTAAATACGCCCCACCTGCTTTTAAAACACCTAATATACCTATAACCATCTCGATGGATCGTTCCATCATCACACCAACAATGGATTCTCTTTTAACTCCTTGGTCTAATAACCTTCTCGCCAACTGATTAGCTTTTATATTTAATTCATTATATGTAATTCCTTTTTCATTACATACAACGGCTATTTGATTAGGGTTCCGTTTTACTTGTTCTTCAAACATTTTATGCACTAATAAATGATTAGAATTTGAATTCTCTTTTTTGTTAAATTCATTCATAATACAATGTTCTTCTTCTATAGATAACATATTAATATTACGCAATCGTACTCTAGGGTTATTAGTTACTTCCTCAACTATATTTGTAAAATGTACCATTAACCTTTCTATTGTCTCCGCTTTAAATAACTTAGTACTATATTCTACTTTTAAATGAATGTTATTATCTATTTCCGTTGCTACTAATGACAAATCAAATTTTGAAACTGACTGCTTAAACGGATACGGTGTAAATTCTAATTCACCAATAGATATTGGATTCATATCCATGTTTTGAAAAACAAACATGGTATCAAATAATGGATTTCTACTTGTATCCCTATGCAAGTCTAAACCTTCTAATAGTTCTTCAAAAGGATAGTCTTGATTTTCATAAGCTTCTAATGTATTGAGTTTTAATCTACTCAAAAACTCAATAAACTCATCGTCATTTTCTAGATAATTTCTCATTACCAAAGTATTAATAAACATACCAATCATATGATTAGTATCAGAATGAGACCTTCCAGCAATAGGCGAACCTACAATAATGTCTTCTTGACCTGTGTATCTAGATAAAAGTATGTTATAAATGGCCAATAAAATCATATATGGCGTAGTACCAGTTTCAGTTGCTAGCTTATTTACTTTAAAAGTTAAATCTCTTCCCAAATTAAAAGAACAAACATTACCTTTAAAGCTTTGTATAGTCGGTCTTTGAAAATCGGTTGGAAAATTTAAAACCGGGAGTTCTCCTTTTAAAGTTGTTAACCAATAATTCTTTTGTTCACTAATTAGATTCTTATAGTAAGGTCCATTTTGCCACATCACATAGTCTTTGTACTGCACTCTCAATTTTGGAAGCTCATTTCCTTTATACAATTCTACAAACTCTTTTATTAATATCCCCATTGATAAACCATCAGATATTATATGATGCATATCTACTACAAGGATATGTCTTTCTTCTGCTATCCTTAAAAGCAACACCCTTAATAATGGAGGTTTTGATAAATCAAATGGACTTATAAACTCATGTATTAAATAATCCGCATCTTTTTCATTTACATGAACGTATTCAATATTGAAATCTACATTAGGCTCAATTTTCTGCACTAATTCCCCATCTAAAATTTGAAAAGAAGTTCTTAATATTTCATGTCTCTCAATTAAAGATTGAAATATATTTTCAAACTTATCTTTACAAATATCCCCTTCTACTTTAAGTATTGTGGGCATATTATAAGTTGTATTTGTACCATCTTCGAATTGATCTACTATAAACATTCTCTTTTGTGACGTAGAAGCTAAATAATACTCTTGTTGTTTTACAGGTTCTATGGAAATATAATTACTTTTTTCCATTTCTAATATACATTTTGAAAAATCAACCAAAATAGGAAACTTAAATAGGGATTTAATAGATAATTGCACATTAAATTCTTTATTAACGATAGAAATTAGCCTAGCAGCCTTTAATGAATGACCACCTATCTCAAAAAAATTATCCCGTATTCCTACCCTTTGTATTCCTAATACATCTTTCCAAATCTCAACTAATTTTCTTTCTGTAGAATTTGTCGGCTCTAGATGACTAGATTTCAAATTATTTATAGGTTGAGGTAATTTTTTTCTATCTATTTTTCCGTTTTGTGTTAACGGTATGTTTTGAATGGATATGATTTGTTGAGGTATCATATAATATGGTAACTTGGTTGCTAAATATGCTCTTACCTCTGGAATAGGAATATCTTTTTCGGTAACTACATACGCACATAAATACTTTTCTCCACTTTCATCTTCTCGCTGTATTACTACGGCGGTTTTGATTGTCTCATATTTTAACAAACTTGCTTCTATCTCACCTAGTTCTATTCGATATCCCCTTATTTTCACTTGATGATCAACTCGTCCCAAGTATTCTATGTTACCATCAGGTAGCCATCTCGCTATATCTCCTGTTTTATATAATTTCTCACCATGTTCAAATGGATGATCAATAAATTTATCTGCTGTTAATTCTTTTCGATTGATGTATCCCCTAGCTAACCCTATGCCAGAAATACATATTTCTCCTGGAACACCTATAGGTTGTAGCCTGTGAAATGAATCCAATATATAAATTTTAGTATTTAACAATGGCGAACCTATTGGTACGTTTTGTATTGTAATTTCTTTATCCCTCTCATATTGATAAGATGTACAACAGACTGTAGCCTCTGTAGGTCCGTATAAATTTAATATTTGGAGGTTCCCCCTAAATAGATGATCGTATTTTGCAAGTAATTCTGTTTTAATTGGTTCTACTCCCACGAAAAGTTTATTTAACGATATCTTCTGATTCGCCCTTACAAAATAATCATATATTTCATTTAATAAGGTAGGTGGAATATATGCTAACGTGACCTGTTCCTCAAGAATGACTTGTACAAGTTTTGGTACATCAAACTTCTCACCTTGATAAATAGTCATTCTTGCTCCATATACCAATGGGACAAATATTTCAAATATAGTAACATCAAAAGAAATACTACTTGAGAATAGTACATTATCACTTATCCCTATATCTTGAGAGAAGTCCTCATACATTGCACACAAAAAATTAGTTAAGGAACGATGTTCAATCATTACCCCTTTGGGCTGGCCTGTAGAACCTGATGTATAAATAACATATGCCAAATTTTGAGGTTCCATCGTTATCTGCAAATCTTCTACTTGTTCTTCTTCAAAAGGAATATCCATTAGATTTATTACACTACCTTGAAATGCTACTCCCTTTATAATAGAGTTTTGATATGTTAGTACATGTGAACACCCGCTGTCTGTCAGCATATATTCCACTCTTTGTTTTGGTAACTCGGTATCAATCGGTAAATACGCCCCACCCGCTTTTAAGATACCTAATATGCCTACAATCATCTCAATAGAGCGTTCCATCATCACACCAACAATGGATTCTCTTTTAACTCCCTGGTCTAATAACCTTCTCGCCAACTGATTAGCTTTTATATTTAATTGTTTATATGTAATTTCTTTTCCATTACATACAATCGCTATTTGATTAGGATTCTGTTTTACCTGCTCTTCAAATAATTGAGGAGCTGTTACAGTCTCACAAAGTAATGTTTTATGAACTCTAGTCGTATGATTGAAATCAAATAATATTTGATTCTTTTCTGTTTTAGGCATTACATCTAGATCCATTGCAGATTTATTTGGATCTTTCATTAATATATCTAAAATATTATATAAATGATTAACTATTCTACTAACTAATCCTTCACTATATAAAATAGAATTGTAATCCACTTGAACCTTTAACTGTTCTTCATTTTTCATAAACCTAATAACCATATCAGAGTTTATTTTATCTGTACTCTCATAACAATGAATATCATCTAACATTACTATTGTATTTAATAAGGGAAGATTATTACTCTCTCCATCTAGACTTAATAATTGAGTAAGTTTATTAAAAGGAAAATGACAATGTTCATTTGACTCTAAAATTGTTTCTTTAATTTTATATATTATTTCTTTAAAATTATCTTCTTGATTTATTTGAGTACGTAATAATAAAAAGTTGTTTTGAAAAACCGTCTCTTCTTGACCTTGTTTAAA</td>
			</tr>
			</tbody>
			</table></code></pre>


			<p><br><strong>Table 3. The contents of "Bt_all_genes.table"</strong></p>
			<pre style="color:black;background-color:white;line-height: 22px;font-size: 14px;overflow-x:scroll;padding-left: 5px;padding-right: 5px;"><code>
			<table border="2" bordercolor="green">
			<thead>
			<tr>
			<th style="text-align:center">-</th>
			<th style="text-align:center">Cry1Ac16</th>
			<th style="text-align:center">Cry1Ca15</th>
			<th style="text-align:center">Cry1Da2</th>
			<th style="text-align:center">Cry1Ia40</th>
			<th style="text-align:center">Cry2Aa10</th>
			<th style="text-align:center">Cry2Ab35</th>
			<th style="text-align:center">Cry78Aa1</th>
			<th style="text-align:center">Cry9Ea9</th>
			<th style="text-align:center">HMM_Cry_len_492</th>
			<th style="text-align:center">HMM_Cyt_len_531</th>
			<th style="text-align:center">HMM_Cyt_len_533</th>
			<th style="text-align:center">HMM_Cyt_len_615</th>
			<th style="text-align:center">Sip1A2-other</th>
			<th style="text-align:center">Vip1Aa2</th>
			<th style="text-align:center">Vip3Aa12</th>
			<th style="text-align:center">Zwa5B-other</th>
			</tr>
			</thead>
			<tbody>
			<tr>
			<td style="text-align:center">1126_1</td>
			<td style="text-align:center">100.00</td>
			<td style="text-align:center">100.00</td>
			<td style="text-align:center">100.00</td>
			<td style="text-align:center">100.00</td>
			<td style="text-align:center"></td>
			<td style="text-align:center">100.00</td>
			<td style="text-align:center">32.66</td>
			<td style="text-align:center">100.00</td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center">100.00</td>
			<td style="text-align:center">99.07</td>
			</tr>
			<tr>
			<td style="text-align:center">AFS094730</td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center">36.14</td>
			<td style="text-align:center"></td>
			<td style="text-align:center">1</td>
			<td style="text-align:center">1</td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center">92.41</td>
			<td style="text-align:center">29.46</td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			</tr>
			<tr>
			<td style="text-align:center">AFS095482</td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center">1</td>
			<td style="text-align:center"></td>
			<td style="text-align:center">1</td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			</tr>
			<tr>
			<td style="text-align:center">AK47</td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center">100.00</td>
			<td style="text-align:center">100.00</td>
			<td style="text-align:center">100.00</td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center">1</td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center"></td>
			<td style="text-align:center">100.00</td>
			<td style="text-align:center">100.00</td>
			</tr>
			</tbody>
			</table>
			<p><strong>Footnote:</strong> The float number represent the identity of blast search, and the integer number represent the number of toxins predicted by HMM and SVM.</p>
			</code></pre>

			<h3 id="license">License</h3>

			<p>BtToxin_Digger is free software, licensed under <a href="https://github.com/liaochenlanruo/BtToxin_Digger/blob/master/LICENSE">GPLv3</a>.</p>

			<h3 id="feedback">Feedback/Issues</h3>

			<p>Please report any issues about usage of the software to the <a href="https://github.com/liaochenlanruo/BtToxin_Digger/issues">issues page</a>.</p>

			<h3 id="citation">Citation</h3>

			<ul>
			<li><p>If you use this software please cite: Liu H, Zheng J, Yu Y, Ye W, Peng D, Sun M. BtToxin_Digger: a comprehensive and high-throughput pipeline for mining toxin protein genes from <em>Bacillus thuringiensis</em>. <em>bioRxiv</em>, 2020. <a href="https://doi.org/10.1101/2020.05.26.114520">10.1101/2020.05.26.114520</a>.</p></li>

			<li><p>If you used the genome assembly function, please also cite: Liu H, Xin B, Zheng J, Zhong H, Yu Y, Peng D, Sun M. Build a bioinformatics analysis platform and apply it to routine analysis of microbial genomics and comparative genomics. <em>Protocol exchange</em>, 2020. DOI: <a href="https://dx.doi.org/10.21203/rs.2.21224/v2">10.21203/rs.2.21224/v2</a>.</p></li>
			</ul>

			<h3 id="faqs">FAQs</h3>

			<h3 id="updates">Updates</h3>

			<ul>
			<li><p>v1.0.2
			-Fixed a "Can not find path" error.</p></li>

			<li><p>v1.0.3
			-Fixed a bug of "get_all_info_nucl.pl", which can not get the gene location and strand information of some toxins.</p></li>
			</ul>

         </div>
       </div>
        <b class="b3"></b><b class="b2"></b><b class="b1"></b> 
   </div>

<div align="center">
<?php
	include_once 'bottom2.html';
?>
</div>
</body>
</html>
