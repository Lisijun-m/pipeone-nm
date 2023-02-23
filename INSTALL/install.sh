conda create --name pipeone_nm python=3 -y
conda activate pipeone_nm
conda install fastp=0.22.0 -y
conda install hisat2 -y
conda install samtools=1.6 -y
conda install stringtie=2.1.6 -y
conda install gffread=0.12.7 -y
conda install blast=2.10.1 -y
conda install transdecoder=5.5.0 -y
conda install pfam_scan=1.6 -y
conda install salmon -y
#conda install bowtie=1.3.1
conda install bowtie2 -y
conda install hmmer=3.3.2
conda install bbmap bedtools -y
#conda install taco
conda install bwa -y
conda install -y pilon
conda install -y canu
conda install -c bioconda seqkit -y
conda install -y zlib
conda install -y R
pip install openpyxl
pip install xlrd
conda deactivate

conda create --name pipeone_nm_2 python=2 -y
conda activate pipeone_nm_2
#pip install CIRIquant
conda install -c bioconda taco=0.7.3  samtools -y
conda deactivate

conda create --name pipeone_nm_4 python=3 -y
conda activate pipeone_nm_4
#conda install -c bioconda trinotate=3.2.2 -y
conda install -c bioconda perl-dbi -y  #install dbi modules
conda deactivate

conda create -n ribodetector python=3.8 -y
conda activate ribodetector
conda install -c bioconda ribodetector -y
conda deactivate

conda create -n rmats -y
conda activate rmats
conda install -c bioconda rmats -y

conda create -n asgal -y
conda activate asgal
conda install python=3.6 -y
conda install biopython -y
pip install pysam
pip install pandas
conda install samtools -y
conda install cmake -y
conda install gffutils -y
conda install zlib
conda deactivate

conda env create -f environment.yml
# this activates the conda environment
conda activate CIRI
