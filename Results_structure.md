## Results
The results generated by operating **PipeOne-NM** will be organized like the following.

* directory *fasterq-dump*
  * fastq format file that has been converted by fastq-dump from sratoolkit. For pair end sequencing, each SRR corresponds to two files. 
    * e.g. *SRRXXXXX_1.fastq* and *SRRXXXXX_2.fastq*
* directory *fastp*
  * Quality control result file. A directory named fastp is created to store data filtering and trimming results. 
  * Fastp will trim 10 bp from the beginning and 2 bp from the tail of each raw read. PipeOne-NM applys to both single-end and pair-end sequencing data.
  * Note that for pair-end (PE) reads, the forward and reversed reads are trimmed separately. 
  * Fastp also removes reads with phred score less than 15, unqualified percent higher than 20% (percents of bases are allowed to be unqualified (0~100)), N base number greater than 0 or length less than 50 bp. 
  * Data filtering and trimming results are stored in fq.gz format.
    * e.g. *SRRXXXXX_fastp.html* 
    * e.g. *SRRXXXXX_fastp.json* 
    * e.g. *SRRXXXXX_fastp.log*
  * fq.gz format file
    * e.g. *SRRXXXXX_1.fq.gz* and *SRRXXXXX_2.fq.gz*
* directory *hisat2*
  * Reads mapping results. 
  * directory *hisat2_index*
    * index file for genome alignment with HISAT2
      * e.g. *genomic.1.ht2*
  * bam format alignment file by hisat2
    * e.g. *SRRXXXXX.bam*
  * directory *un_conz_fastq*
    * e.g. *SRRXXXXX_um_1.fq.gz* and *SRRXXXXX_um_2.fq.gz*
* directory *stringtie*
  * gtf format file generated by stringtie
    * e.g. *SRRXXXXX.gtf*
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
    * e.g. *SRRXXXXX.sam*
* directory *ciri*
  * For more detailed description, please visit [https://ciri-cookbook.readthedocs.io/en/latest/CIRI2.html](https://ciri-cookbook.readthedocs.io/en/latest/CIRI2.html).
  * circRNA detection result based on sam format alignment file
    * e.g. *SRRXXXXX.ciri*
  * circRNA alternative splicing events detection
    * e.g. *SRRXXXXX.as_AS.list*
    * e.g. *SRRXXXXX.as_coverage.list*
    * e.g. *SRRXXXXX.as_jav.list*
    * e.g. *SRRXXXXX.as_library_length.list*
    * e.g. *SRRXXXXX.as.list*
    * e.g. *SRRXXXXX.as_splice.list*
  * Visualization of circRNA alternative splicing detection in pdf format
    * e.g. directory *SRRXXXXX*
  * directory *ciriquant*
    * directory *diff*
    * directory of each SRR, containing quantification result for each sample
* directory *bedtools*
  * fasta format circRNA extracted according to circRNA detection results
    * e.g. *SRRXXXXX.fasta*
* directory *miranda*
  * For more details about miRanda result, visit [https://cbio.mskcc.org/miRNA2003/miranda.html](https://cbio.mskcc.org/miRNA2003/miranda.html) for more information.
  * circRNA-miRNA interaction prediction file
    * e.g. *SRRXXXXX_circ_miRNA.txt*
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
