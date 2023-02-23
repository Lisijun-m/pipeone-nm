## Parameter Set in PipeOne-NM

Parameters that are not mentioned specifically in this file are the default parameters of each software/tool.

#### Ribodetector
* -e norrna
  * norrna: output only high confident non-rRNAs, the rest are clasified as rRNAs;
* --chunk_size 256:
  * chunk_size * 1024 reads to load each time.

#### fastp
* -f 10
  * trimming how many bases in front for read1, default is 0 (int [=0])
* -t 2 
  * trimming how many bases in tail for read1, default is 0 (int [=0])
* -q 15
  * the quality value that a base is qualified. Default 15 means phred quality >=Q15
* -u 20
  *  how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
* -n 0
  * if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
* -l 50
  * reads shorter than length_required will be discarded, default is 15. (int [=15]
* -w 4
  * worker thread number, default is 2 (int [=2])

#### HISAT2
* --dta -x
  * reports alignments tailored for transcript assemblers

#### samtools
* -m 4G
  * Set maximum memory per thread; suffix K/M/G recognized [768M]
* -@ 8
  * Number of additional threads to use [0]

#### gffread
* -M
  * cluster the input transcripts into loci, discarding "redundant" transcripts (those with the same exact introns and fully contained or equal boundaries)
* -K
  * for -M option: also discard as redundant the shorter, fully contained transcripts (intron chains matching a part of the container)
* --table @geneid
  * output a simple tab delimited format instead of GFF, with columns having the values of GFF attributes given in <attrlist>; special pseudo-attributes (prefixed by @) are recognized

#### trinityrnaseq
* --seqType fq
  * sequence type, fastq or fasta
* --est_method salmon
  * abundance estimation method. alignment_based:RSEM; alignment_free: kallisto|salmon
* --thread_count 32
  * number of threads to use (default = 4)
* --prep_reference
  * prep reference (builds target index)

#### blastp
* -max_target_seqs 1 
  * Maximum number of aligned sequences to keep (value of 5 or more is recommended) Default = `500'
* -outfmt 6 
  * alignment view options: 6 = Tabular
* -evalue 1e-5 
  * Expectation value (E) threshold for saving hits, Default = `10'

#### Trinotate: extract_GO_assignments_from_Trinotate_xls.pl
* -G
  * gene-mode

#### BWA MEM
* -T 19 
* -t 10

#### CIRI_AS
* -D yes
  * if output all processing info (Choose 'yes' would require more disk space. default: no)
