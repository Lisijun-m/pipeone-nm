conda create --name pipeone_nm python=3
conda activate pipeone_nm
conda install -c bioconda sra-tools=2.11.0 fastp=0.22.0 \
    hisat2 samtools=1.6  stringtie=2.1.6  gffread=0.12.7 \
     transdecoder=5.5.0  blast=2.10.1  pfam_scan=1.6  salmon=1.9.0 \
      bowtie=1.3.1  bowtie2=2.4.5  hmmer=3.3.2  trinotate=3.2.2  bbmap   bedtools  rmats
conda deactivate

conda create --name pipeone_nm_2 python=2
conda activate pipeone_nm_2
pip install CIRIquant
conda install -c bioconda taco=0.7.3  samtools
conda deactivate


# conda create --name pipeone_nm_3 python=3
# conda activate pipeone_nm_3
# conda install -c predector tmhmm=2.0c
# conda install -c predector signalp4 #但是因为没注册装了也等于没装。。
# conda install bedtools
# conda deactivate

conda create --name pipeone_nm_4 python=3
conda activate pipeone_nm_4
conda install -c bioconda trinotate=3.2.2
conda install -c bioconda perl-dbi   #install dbi modules
conda deactivate

#cpan -i install  XML::Simple