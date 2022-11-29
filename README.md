# PipeOne-NM

**PipeOne-NM** is a tool for genome functional annotation, non-coding RNA identification, transcripts alternative splicing analysis and gene differential expression in non-model organisms. 

PipeOne-NM is applicable to species with reference genome and can provide comprehensive insights into transcriptome of non-model organisms. It integrates twenty-one different tools to extract information from RNA-seq data and utilizes five R packages to carry out downstream analysis and data visualization. 

<p align="center">
<img src="https://github.com/Lisijun-m/pipeone-nm/blob/main/Figures/workflow.png" width="900px">
</p>

# Index
- [Prerequisites](#prerequisites)
- [Installation](#installation)
  - [1.Download PipeOne-NM](#1-download-pipeone-nm)
  - [2.Setup](#2-setup)
  - [3.Configuration](#3-configuration)
- [Run the Pipeline](#run-the-pipeline)
- [Results](#results) 
- [Rstudio scripts](#rstudio-scripts) 
- [Contact](#contact) 


## Prerequisites

1. [Docker](https://www.docker.com/) or [conda](https://docs.conda.io/en/latest/miniconda.html)
2. Java (version >= 1.7)
3. [Nextflow](https://www.nextflow.io/) (version >= 20.07.1.5413)

## Installation

### __1. Download PipeOne-NM__

```
git clone https://github.com/nongbaoting/pipeone-nm.git
```

### __2. Setup__

##### 2.1 pull docker image or install conda environments

**Option 1: pull docker images run by dock**

```bash
docker pull nongbaoting/pipeone_nm:latest
```

**Option 2: install conda environments** 

```
cd pipeone-nm/INSTALL
bash ./install.sh
```
To ensure a successful installation, we recommend that you run each command step by step in the shell script install.sh

##### __2.2 softwares need to be registeration__

* Trinotate
* tmhmm v2 (free academic download) [https://services.healthtech.dtu.dk/service.php?TMHMM-2.0](https://services.healthtech.dtu.dk/service.php?TMHMM-2.0)

* RNAMMER (free academic download) [https://services.healthtech.dtu.dk/service.php?RNAmmer-1.2](https://services.healthtech.dtu.dk/service.php?RNAmmer-1.2)
* signalP v4 (free academic download) [https://services.healthtech.dtu.dk/service.php?SignalP-5.0](https://services.healthtech.dtu.dk/service.php?SignalP-5.0)

The above four softwares are not included in the docker container and need to be installed by users.

### __3. Configuration__

#### basic configuration
Modify the program configuration file `pipeone-nm/conf/genomes.config`,  change the line below and store preparation files for *PipeOne-NM*:

```
params {
  genomes {
  "fish_test" {
	refgenome = "/path/to/genome/genomic.fna"
	sample = "/path/to/genome/samples.txt"
	ref_miRNA = '/path/to/genome/public5/l/FishmiRNA-June2021-dre-mature.fasta'			  // for miranda
	uniprot = "/path/to/genome/uniprot_sprot.pep"                                                     // for Trinotate
	pfam = "/path/to/genome/Pfam-A.hmm"                                                               // for Trinotate
	sqlite = "/path/to/genome/Trinotate_20210616.sqlite"                                         	  // for Trinotate 
	min_length = 150
    }
  }
}
```
* Replace *fish_test* with the name of your project
* *refgenome* is the file path of reference genome (fasta format), replace it with your path of reference genome
* *ref_miRNA* is the miRNA reference for non-model species. Reference from genetic similar model organisms can be used as alternative. Replace it with your path of miRNA reference. In the test data of grass carp, the miRNA reference was downloaded from [fishmirna](http://fishmirna.org/data/download/FishmiRNA-June2021-dre-mature.fasta)
* *uniprot* is a input database, downloading of *uniprot* can be found in [Trinotate](https://github.com/Trinotate/Trinotate/wiki/Software-installation-and-data-required)
* *pfam* is a input database, downloading of *pfam* can be found in [Trinotate](https://github.com/Trinotate/Trinotate/wiki/Software-installation-and-data-required)
* *sqlite* is input sqlite file, downloading of *sqlite* can be found in [Trinotate](https://github.com/Trinotate/Trinotate/wiki/Software-installation-and-data-required)
* *sample* is the sample information (including biological conditions and repetitions) of your project in txt format, for example:

```
Condition_A    sample1    sample1_left.fq    sample1_right.fq
Condition_A    sample2    sample2_left.fq    sample2_right.fq
Condition_B    sample3    sample3_left.fq    sample3_right.fq
Condition_B    sample4    sample4_left.fq    sample4_right.fq
```
*Condition_A* and *Condition_B* stand for biological conditions (treatments), *sample1* stands for sample name and *sample1_left.fq* stands for sequencing file in fastq format.

#### Configuraton for alternative splicing analysis
For alternative splicing analysis, provide information file for two biological conditions in txt format and name the two file as *AS_con1.txt* and *AS_con2.txt*. 
For example, *AS_con1.txt* is file that contains sample alignment file path for the first biological condition. Each bam file path is separated by comma.

```
/result/hisat2/SRR16151827.bam,/result/hisat2/SRR16151828.bam,/result/hisat2/SRR16151829.bam
```

Open the `pipeone-nm/main.nf` and replace line 25 with corresponding file path

```
as_files = tuple('/path/to/AS_con1.txt', '/path/to/AS_con2.txt')
```


## __Run the Pipeline__

```
nextflow run /path/to/pipeOne-nm  -profile docker --sra './sra/*.sra' --genome fish_test
```

* Download RNA-seq data in sra format for test from NCBI
* Required Input Files
  * RNA-seq data
    * --sra, RNA-seq data in sra format downloaded from GEO, store in directory *sra*, e.g. *SRR16151823.sra*
    * or --fastq, RNA-seq data in fastq format
  * --genome, genome defined in `conf/genomes.config`


## Results
The results generated by operating **PipeOne-NM** will be organized like the following.
Use grass carp RNA-seq data GSE185170 as analysis example.

* directory *fasterq-dump*
  * fastq format file that has been converted by fastq-dump from sratoolkit. For pair end sequencing, each SRR corresponds to two files. 
    * e.g. *SRR16151823_1.fastq* and *SRR16151823_2.fastq*
* directory *fastp*
  * Quality control result file. A directory named fastp is created to store data filtering and trimming results. 
  * Fastp will trim 10 bp from the beginning and 2 bp from the tail of each raw read. PipeOne-NM applys to both single-end and pair-end sequencing data.
  * Note that for pair-end (PE) reads, the forward and reversed reads are trimmed separately. 
  * Fastp also removes reads with phred score less than 15, unqualified percent higher than 20% (percents of bases are allowed to be unqualified (0~100)), N base number greater than 0 or length less than 50 bp. 
  * Data filtering and trimming results are stored in fq.gz format.
    * e.g. *SRR16151823_fastp.html* 
    * e.g. *SRR16151823_fastp.json* 
    * e.g. *SRR16151823_fastp.log*
  * fq.gz format file
    * e.g. *SRR16151823_1.fq.gz* and *SRR16151823_2.fq.gz*
* directory *hisat2*
  * Reads mapping results. 
  * directory *hisat2_index*
    * index file for genome alignment with HISAT2
      * e.g. *genomic.1.ht2*
  * bam format alignment file by hisat2
    * e.g. *SRR16151823.bam*
  * directory *un_conz_fastq*
    * e.g. *SRR16151823_um_1.fq.gz* and *SRR16151823_um_2.fq.gz*
* directory *stringtie*
  * gtf format file generated by stringtie
    * e.g. *SRR16151823.gtf*
* directory *TACO*
  * For more information about taco results, visit [https://tacorna.github.io/](https://tacorna.github.io/).
  * *assembly.gtf*
  * *assembly.bed*
  * *args.pickle*
  * *change_points.gtf*
  * *expr.neg.bedgraph*
  * *expr.none.bedgraph*
  * *expr.pos.bedgraph*
  * *path_graph_stats.txt*
  * *loci.txt*
  * *samples.txt*
  * *splice_graph.gtf*
  * *splice_junctions.bed*
  * *status.json*
  * *transfrags.bed*
  * *transfrags.filtered.bed*
  * *TACO.log*
* directory *gffread*
  * *Transcripts.fasta*
  * *Transctipts.fasta.gene_trans_map*
  * *GTFs.txt*
* directory *salmon*
  * For more information about salmon results, visit [https://combine-lab.github.io/salmon/](https://combine-lab.github.io/salmon/).
  * directory *salmon_index*
  * *matrix.log*
  * directory *samples*
    * *align_and_estimate_abundance.log*
    * quantify results for each sample (separately direcotry)
      * e.g. SRR16151823
    * *Transcripts.gene.counts.matrix*
    * *Transcripts.gene.TMM.EXPR.matrix*
    * *Transcripts.gene.TPM.not_cross_norm*
    * *Transcripts.gene.TPM.not_cross_norm.runTMM.R*
    * *Transcripts.gene.TPM.not_cross_norm.TMM_info.txt*
    * *Transcripts.isoform.counts.matrix*
    * *Transcripts.isoform.TMM.EXPR.matrix*
    * *Transcripts.isoform.TPM.not_cross_norm*
    * *Transcripts.isoform.TPM.not_cross_norm.runTMM.R*
    * *Transcripts.isoform.TPM.not_cross_norm.TMM_info.txt*
* directory *transdecoder*
  * For more information about transdecoder, visit [https://github.com/TransDecoder/TransDecoder/wiki](https://github.com/TransDecoder/TransDecoder/wiki).
  * *blastp.outfmt6*
  * *pfam.domtblout*
  * *Transcripts.fasta.transdecoder.bed*
  * *Transcripts.fasta.transdecoder.cds*
  * directory *Transcripts.fasta.transdecoder_dir*
  * directory *Transcripts.fasta.transdecoder_dir.__checkpoints*
  * directory *Transcripts.fasta.transdecoder_dir.__checkpoints_longorfs*
  * *TransDecoder.LongOrfs.log*
  * *TransDecoder.Predict.log*
* directory *trinotate*
  * For more information about trinotate results, visit [https://rnabio.org/module-07-trinotate/0007/02/01/Trinotate/](https://rnabio.org/module-07-trinotate/0007/02/01/Trinotate/).
  * *blastp.outfmt6*
  * *blastx.outfmt6*
  * *go_annotations.txt*
  * *signalp.out*
  * *tmhmm.out*
  * *Transcripts.fasta.rnammer.gff*
  * *transcriptsSuperScaffold.bed*
  * *transcriptsSuperScaffold.fasta*
  * *Trinotate_lncRNA.fa*
  * *Trinotate_lncRNA.fa_n50_stat*
  * *Trinotate_lncRNA.txt*
  * *Trinotate_mRNA.fa*
  * *Trinotate_mRNA.fa_n50_stat*
  * *Trinotate_mRNA.txt*
  * *TrinotatePFAM.out*
  * *Trinotate_rRNA.fa*
  * *Trinotate_rRNA.fa_n50_stat*
  * *Trinotate_rRNA.txt*
  * *Trinotate.xls*
* directory *bwaindex*
* directory *bwamem*
  * sam format alignment file by bwa-mem
    * e.g. *SRR16151823.sam*
* directory *ciri*
  * For more detailed description, please visit [https://ciri-cookbook.readthedocs.io/en/latest/CIRI2.html](https://ciri-cookbook.readthedocs.io/en/latest/CIRI2.html).
  * circRNA detection result based on sam format alignment file
    * e.g. *SRR16151823.ciri*
  * circRNA alternative splicing events detection
    * e.g. *SRR16151823.as_AS.list*
    * e.g. *SRR16151823.as_coverage.list*
    * e.g. *SRR16151823.as_jav.list*
    * e.g. *SRR16151823.as_library_length.list*
    * e.g. *SRR16151823.as.list*
    * e.g. *SRR16151823.as_splice.list*
  * Visualization of circRNA alternative splicing detection in pdf format
    * e.g. directory *SRR16151823*
  * directory *ciriquant*
    * directory *diff*
    * directory of each SRR, containing quantification result for each sample
* directory *bedtools*
  * fasta format circRNA extracted according to circRNA detection results
    * e.g. *SRR16151823.fasta*
* directory *miranda*
  * For more details about miRanda result, visit [https://cbio.mskcc.org/miRNA2003/miranda.html](https://cbio.mskcc.org/miRNA2003/miranda.html) for more information.
  * circRNA-miRNA interaction prediction file
    * e.g. *SRR16151823_circ_miRNA.txt*
* directory *rmats*
  * For details about rmats result, visit RMATS website [http://rnaseq-mats.sourceforge.net/](http://rnaseq-mats.sourceforge.net/) for more information.
  * A3SS.MATS.JCEC.txt*
  * *A3SS.MATS.JC.txt*
  * *A5SS.MATS.JCEC.txt*
  * *A5SS.MATS.JC.txt*
  * *fromGTF.A3SS.txt*
  * *fromGTF.A5SS.txt*
  * *fromGTF.MXE.txt*
  * *fromGTF.novelEvents.A3SS.txt*
  * *fromGTF.novelEvents.A5SS.txt*
  * *fromGTF.novelEvents.MXE.txt*
  * *fromGTF.novelEvents.RI.txt*
  * *fromGTF.RI.txt*
  * *fromGTF.SE.txt*
  * *JCEC.raw.input.A3SS.txt*
  * *JCEC.raw.input.A5SS.txt*
  * *JCEC.raw.input.MXE.txt*
  * *JCEC.raw.input.RI.txt*
  * *JCEC.raw.input.SE.txt*
  * *JC.raw.input.A3SS.txt*
  * *JC.raw.input.A5SS.txt*
  * *JC.raw.input.MXE.txt*
  * *JC.raw.input.RI.txt*
  * *JC.raw.input.SE.txt*
  * *MXE.MATS.JCEC.txt*
  * *MXE.MATS.JC.txt*
  * *RI.MATS.JCEC.txt*
  * *RI.MATS.JC.txt*
  * *SE.MATS.JCEC.txt*
  * *SE.MATS.JC.txt*

## Supplementary scripts
Compared with above process, differential expression analysis needs to be tailored to different datasets. It may involve removing of batch effect, differential analysis between multiple biological conditions and retrieving differential expression genes according to individualized threshold. 

Therefore, the differential analysis pipeline is not included in our PipeOne-NM, but here we present an example of performing differential analysis on grass carp circRNA host genes in folder Supplementary-example-scripts.

For the analysis of co-experssion of lncRNA-mRNA pair, we provide the matlab scripts in folder Supplementary-example-scripts.

## Contact
For any problem, write to lisj35@mail2.sysu.edu.cn


