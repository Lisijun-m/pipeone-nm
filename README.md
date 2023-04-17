# PipeOne-NM

**PipeOne-NM** is a tool for genome functional annotation, non-coding RNA identification, transcripts alternative splicing analysis and gene differential expression in non-model organisms. 

PipeOne-NM is applicable to species with reference genome and can provide comprehensive insights into transcriptome of non-model organisms. It integrates twenty-one different tools to extract information from RNA-seq data and utilizes five R packages to carry out downstream analysis and data visualization. 

<p align="center">
<img src="https://github.com/Lisijun-m/pipeone-nm/blob/main/Figures/workflow.png" width="900px">
</p>

# Index
- [Prerequisites and Configuration](#prerequisites-and-configuration)
- [Installation](#installation)
  - [1.Download PipeOne-NM](#1-download-pipeone-nm)
  - [2.Setup](#2-setup)
  - [3.Input Preparation](#3-input-preparation)
- [Run the Pipeline](#run-the-pipeline)
- [Results](#results) 
- [Supplementary scripts](#supplementary-scripts)


## Prerequisites and Configuration

1. [Docker](https://www.docker.com/) or [conda](https://docs.conda.io/en/latest/miniconda.html)
2. Java (version >= 1.7)
3. [Nextflow](https://www.nextflow.io/) (version >= 20.07.1.5413) (can be installed via conda)

* Operation System: Linux
* CPU Core Number: 33
* RAM: 55G

## Installation

### __1. Download PipeOne-NM__

```
git clone https://github.com/Lisijun-m/pipeone-nm.git
```

### __2. Setup__

##### 2.1 pull docker image or install conda environments

**Option 1: pull docker images run by docker**

```bash
docker pull sijunli/pipeone_nm:latest
```

**Option 2: install conda environments** 

```
cd pipeone-nm/INSTALL
bash ./install.sh
```
To ensure a successful installation, we recommend that you run each command step by step in the shell script install.sh

After conda enviroments are settled, install Trinotate, trinityrnaseq, miRanda and CIRI2 in pipeone directory.

An example may be like this:

```
## Trinotate
cd pipeone
git clone https://github.com/Trinotate/Trinotate.git

## set up trinotate databases
/path/to/pipeone/Trinotate/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate

## trinityrnaseq
cd pipeone
wget https://github.com/trinityrnaseq/trinityrnaseq/archive/refs/tags/Trinity-v2.13.2.tar.gz
tar -xzvf Trinity-v2.13.2.tar.gz

## CIRI2
cd pipeone
wget https://sourceforge.net/projects/ciri/files/CIRI-full/CIRI-full_v2.0.zip
unzip CIRI-full_v2.0.zip

## miRanda
cd pipeone
wget http://cbio.mskcc.org/microrna_data/miRanda-aug2010.tar.gz
tar zxvf miRanda-aug2010.tar.gz
mv miRanda-3.3a  miRanda
./configure --prefix=$HOME/software/miRanda
make
make install
echo 'export PATH=/path/to/pipeone/miRanda-3.3a/src:$PATH' >> ~/.bashrc
source ~/.bashrc

```
The docker images contains Trinotate, trinityseq, miRanda and CIRI2, so the installation is needed only in option 2.

##### __2.2 softwares need to be registerated__

* tmhmm v2 (free academic download) [http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm)
* RNAMMER (free academic download) [http://www.cbs.dtu.dk/cgi-bin/sw_request?rnammer](http://www.cbs.dtu.dk/cgi-bin/sw_request?rnammer)
* signalP v4 (free academic download) [http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp)

The above four softwares are not included in the docker container nor the conda install.sh script and need to be installed by users.

**In Docker Images**

Run the docker container in interaction mode, and install tmhmm, RNAMMER and signalP in directory _/pipeone_. The following codes can be used as a guide.

```
docker run -it sijunli/pipeone_nm
cd /pipeone
tar -zxvf signalp-4.1g.Linux.tar.gz
tar -zxvf rnammer-1.2.src.tar.gz
tar -zxvf tmhmm-2.0c.Linux.tar.gz
echo 'export PATH=/pipeone/signalp-4.1:$PATH' >> ~/.bashrc
echo 'export PATH=/pipeone/rnammer-1.2:$PATH' >> ~/.bashrc
echo 'export PATH=/pipeone:$PATH' >> ~/.bashrc
echo 'export PATH=/pipeone/tmhmm-2.0c/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
exit
```

**In conda environments**

Install tmhmm, RNAMMER and signalP in _/path/to/pipeone_. The following codes can be used as a guide.

```
cd /path/to/pipeone
tar -zxvf signalp-4.1g.Linux.tar.gz
tar -zxvf rnammer-1.2.src.tar.gz
tar -zxvf tmhmm-2.0c.Linux.tar.gz
echo 'export PATH=/pipeone/signalp-4.1:$PATH' >> ~/.bashrc
echo 'export PATH=/pipeone/rnammer-1.2:$PATH' >> ~/.bashrc
echo 'export PATH=/pipeone:$PATH' >> ~/.bashrc
echo 'export PATH=/pipeone/tmhmm-2.0c/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

### __3. Input Preparation__

#### sample information file

sample information (including biological conditions and repetitions) of your project in txt format, for example:

```
Condition_A    sample1    sample1_left.fq    sample1_right.fq
Condition_A    sample2    sample2_left.fq    sample2_right.fq
Condition_B    sample3    sample3_left.fq    sample3_right.fq
Condition_B    sample4    sample4_left.fq    sample4_right.fq
```
*Condition_A* and *Condition_B* stand for biological conditions (treatments), *sample1* stands for sample name and *sample1_left.fq* stands for sequencing file in fastq format.

#### Alternative Splicing file

For alternative splicing analysis, provide information file for two biological conditions in txt format and name the two file as *AS_con1.txt* and *AS_con2.txt*. 

For example, *AS_con1.txt* is file that contains sample alignment file path for the first biological condition. Each bam file path is separated by comma.

```
/result/hisat2/sample1.bam,/result/hisat2/SRR16151828.bam,/result/hisat2/sample2.bam
```

## __Run the Pipeline__

**Option 1: run with nextflow (Docker image required)**

#### basic configuration for running with nextflow
Modify the program configuration file `pipeone-nm/conf/genomes.config`,  change the line below and store preparation files for *PipeOne-NM*:

```
params {
  genomes {
  "test" {
	refgenome = "/path/to/genome/genomic.fna"
	sample = "/path/to/genome/samples.txt"
	ref_miRNA = '/path/to/genome/miRNA_ref.fasta'	// for miranda
	uniprot = "/path/to/genome/uniprot_sprot.pep"                                                     // for Trinotate
	pfam = "/path/to/genome/Pfam-A.hmm"                                                               // for Trinotate
	sqlite = "/path/to/genome/Trinotate.sqlite"                                         	  // for Trinotate 
	min_length = 150
    }
  }
}
```

* Replace *test* with the name of your project
* *refgenome* is the file path of reference genome (fasta format), replace it with your path of reference genome
* *ref_miRNA* is the miRNA reference for non-model species. Reference from genetic similar model organisms can be used as alternative. Replace it with your path of miRNA reference. In the test data of grass carp, the miRNA reference was downloaded from [fishmirna](http://fishmirna.org/data/download/FishmiRNA-June2021-dre-mature.fasta)
* *uniprot* is a input database, downloading of *uniprot* can be found in [Trinotate](https://github.com/Trinotate/Trinotate/wiki/Software-installation-and-data-required)
* *pfam* is a input database, downloading of *pfam* can be found in [Trinotate](https://github.com/Trinotate/Trinotate/wiki/Software-installation-and-data-required)
* *sqlite* is input sqlite file, downloading of *sqlite* can be found in [Trinotate](https://github.com/Trinotate/Trinotate/wiki/Software-installation-and-data-required)
* *sample* is the file mentioned in Input Preparation part

Open the `pipeone-nm/main.nf` and replace line 25 with corresponding file path

```
as_files = tuple('/path/to/AS_con1.txt', '/path/to/AS_con2.txt')
```
Then run the pipeline with nextflow:
```
docker pull nongbaoting/pipeone_nm:latest
nextflow run /path/to/pipeOne-nm  -profile docker --sra './sra/*.sra' --genome test
```

* Required Input Files
  * RNA-seq data
    * --sra, RNA-seq data in sra format downloaded from GEO, store in directory *sra*, e.g. *SRR16151823.sra*
    * or --fastq, RNA-seq data in fastq format
  * --genome, genome defined in `conf/genomes.config`

**Option 2: run with python script (Docker image required)**

```
docker run -it -v /Users/project:/work_dir sijunli/pipeone_nm /bin/bash 
python basicPipeline.py -i /work_dir/path/to/sra -o /work_dir/path/to/output
exit
```
* directory /Users/project is the directory contains all the raw data and will store the output, it will be projected in docker images as directory /work_dir
* After the execution is done, exit the docker images and the analysis result will be stored in directory /Users/project/

**Option 2: run with python script (Conda enviroments required)**

The python script of pipeline has multiple parameters, you can see a full description of the parameters typing "python basicPipeline.py -h".

An example might be this:
```
python /path/to/pipeone/basicPipeline.py -i /path/to/sra -o /path/to/output
```

## Parameters
The parameters set in the pipeline can be viewed in *Parameters.md*. Other parameters are default parameters for each tool.


## Supplementary scripts
Compared with above process, differential expression analysis needs to be tailored to different datasets. It may involve removing of batch effect, differential analysis between multiple biological conditions and retrieving differential expression genes according to individualized threshold. 

Therefore, the differential analysis pipeline is not included in our PipeOne-NM, but here we present an example of performing differential analysis on grass carp circRNA host genes in folder Supplementary-example-scripts.
