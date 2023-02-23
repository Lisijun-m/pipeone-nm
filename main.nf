#!/usr/bin/env nextflow
/*  @Date         : 2022/08/09 10:13:34
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting
    @Version      : 
    @Description  : 
*/
nextflow.enable.dsl=2


// defined params
params.sra = ''
params.fastq = ''
params.genome
params.cleaned  = false
params.singleEnd =  false
params.read_length = 100
params.rm_rRNA = true

params.refgenome = params.genome ? params.genomes[ params.genome ].refgenome ?: false : false
params.min_length = params.genome ? params.genomes[ params.genome ].min_length ?: false : false
params.ref_miRNA = params.genome ? params.genomes[ params.genome ].ref_miRNA ?: false : false
params.sqlite = params.genome ? params.genomes[ params.genome ].sqlite ?: false : false
params.uniprot = params.genome ? params.genomes[ params.genome ].uniprot ?: false : false
params.pfam = params.genome ? params.genomes[ params.genome ].pfam ?: false : false

// config_file = '/public5/lisj/Software/CIRI-quant/CIRI-quantContig.yaml'
as_files = tuple('/public5/lisj/18_AS/AS_con1.txt', '/public5/lisj/18_AS/AS_con2.txt')
// lst = tuple("/public5/lisj/14_exp/SRR16151823.lst", "/public5/lisj/14_exp/SRR16151824.lst")

// include functions must be placed after params !
include {Header; check_file; input_reads; input_sra; get_base_index; get_dir_files; input_fna} from "${baseDir}/modules/functions.nf"

include { Fasterq_dump; } from  "${baseDir}/modules/process/fasterq_dump.nf"
include { Fastp; } from  "${baseDir}/modules/process/fastp.nf"
include { Hisat2; Hisat2_build; } from  "${baseDir}/modules/process/hisat2.nf"
include { Stringtie; } from  "${baseDir}/modules/process/stringtie.nf"
include { Taco;  } from "${baseDir}/modules/process/taco.nf"
include { Gffread;  } from "${baseDir}/modules/process/gffread.nf"
include { Align_and_estimate_abundance; Abundance_estimates_to_matrix } from "${baseDir}/modules/process/trinity.nf"
include { Salmon_idx_cal  } from "${baseDir}/modules/process/salmon.nf"
include { Transdecoder_LongOrfs; Blastp_LongOrfs;Hmmpress; Hmmscan_LongOrfs; Transdecoder_Predict  } from "${baseDir}/modules/process/transdecoder.nf"
include { Search_Transcripts; Blastp_predict; Rnammer; Hmmscan_Predict; Signalp; Tmhmm; Trinotate; Extract_GO; RNAs; Extract_seqs_by_IDs } from "${baseDir}/modules/process/trinotate.nf"
include { Bwaindex; Bwamem } from "${baseDir}/modules/process/bwa.nf"
include { Ciri; CiriAs; CIRIvis; CIRIquant;PrepDE_circRNA_Count_gene_reads } from "${baseDir}/modules/process/ciri.nf"
include { Reformat } from "${baseDir}/modules/process/reformat.nf"
include { Bedtools } from "${baseDir}/modules/process/bedtools.nf"
include { Miranda } from "${baseDir}/modules/process/miranda.nf"
include { Rmats } from "${baseDir}/modules/process/rmats.nf"
include {Asgal} from "${baseDir}/modules/process/asgal.nf"
include {Ribodetector} from "${baseDir}/modules/process/ribodetector.nf"

// include { PrepDE } from "${baseDir}/modules/process/prepDE.nf"


log.info Header()
def refgenome = check_file(params.refgenome,'refgenome')
def uniprot = check_file(params.uniprot, 'uniprot')
def pfam = check_file(params.pfam,'pfam')
def ref_miRNA = check_file(params.ref_miRNA,'ref_miRNA')
def sqlite = check_file(params.sqlite, 'sqlite')

workflow {

    if(params.sra && ! params.fastq){
        sraFiles = input_sra(params.sra)
        Fasterq_dump(sraFiles)
        reads = Fasterq_dump.out.reads.map{it -> [it[0], [it[1], it[2]]] }
    }else if(!params.sra &&  params.fastq){
        reads = input_reads(params.fastq)
    }else{
        log.info "\033[0;35m" + "Please provide sra/fasq files with parameter: --sra or --fastq !" + "\033[0m"
        System.exit(1)
    }
    
    
    Fastp(reads)
    cleaned_reads_fp = Fastp.out.cleaned_reads.map{it -> [it[0], [it[1], it[2]]] }
    
    if(params.rm_rRNA){

        Ribodetector(cleaned_reads_fp)
        cleaned_reads = Ribodetector.out
    }else{
        cleaned_reads =cleaned_reads_fp
    }


    Hisat2_build(refgenome)
    hisat2_idx = Hisat2_build.out.ht2.collect()

    Hisat2(cleaned_reads, "genomic", hisat2_idx)
   

    Stringtie(Hisat2.out.bam)
    gtfs = Stringtie.out.gtf.collect()
    
    Taco(gtfs)
    TACO = Taco.out.TACO
    
    Gffread(refgenome,TACO)
    exon = Gffread.out.exon
    
   Salmon_idx_cal(reads, refgenome, TACO)
    // Align_and_estimate_abundance(exon, reads)
    // aea_outdir = Align_and_estimate_abundance.out.outdir
    gene_trans_map = Salmon_idx_cal.out.gene_trans_map
    
    // Abundance_estimates_to_matrix(aea_outdir, gene_trans_map)

    Transdecoder_LongOrfs(exon)
    dir = Transdecoder_LongOrfs.out.dir
   

    Blastp_LongOrfs( dir, uniprot)
    blastp_LongOrfs = Blastp_LongOrfs.out.blastp_LongOrfs
    makeblastu = Blastp_LongOrfs.out.makeblastu.collect()
    
    Hmmpress(pfam)
    press = Hmmpress.out.press

    Hmmscan_LongOrfs( dir, press)
    domtblout = Hmmscan_LongOrfs.out.domtblout
    // press = Hmmscan_LongOrfs.out.press

    Transdecoder_Predict(exon, domtblout, blastp_LongOrfs,  dir)
    predict = Transdecoder_Predict.out.predict
    
    Search_Transcripts(exon, makeblastu)
    blastx = Search_Transcripts.out.blastx
    
    Rnammer(exon)
    rnammer = Rnammer.out.rnammer
    
    Hmmscan_Predict(predict, press)
    TPFAM = Hmmscan_Predict.out.TPFAM
    
    Blastp_predict(predict, makeblastu)
    blastp_predict = Blastp_predict.out.blastp_predict
    
    Signalp(predict)
    signalp = Signalp.out.signalp
    
    Tmhmm(predict)
    tmhmm = Tmhmm.out.tmhmm
    
    Trinotate(sqlite, gene_trans_map, exon, predict, blastp_predict, blastx, TPFAM, tmhmm, signalp, rnammer)
    trinotatefile = Trinotate.out.trinotatefile
    
    Extract_GO(trinotatefile)
    go_annotations =  Extract_GO.out.go

    RNAs(trinotatefile)
    mRNA = RNAs.out.mRNA
    rRNA = RNAs.out.rRNA
    lncRNA = RNAs.out.lncRNA

    Extract_seqs_by_IDs(exon, rRNA, lncRNA, mRNA)

    Bwaindex(refgenome)
    bwaindex = Bwaindex.out.bwaindex
    
    Reformat(reads, params.min_length)
    trimmed_reads = Reformat.out.trimmed_reads.map{it -> [it[0], [it[1], it[2]]] }
    no_id = Reformat.out.no_id.collect()

    Bwamem(bwaindex, trimmed_reads)
    bwamem = Bwamem.out.bwamem
    
    Ciri(refgenome,  TACO, bwamem,)
    ciri = Ciri.out.ciri

    bwabam_ciri =  bwamem.combine(ciri, by:0)
    CiriAs(refgenome,  TACO, bwabam_ciri)
    cirias = CiriAs.out.cirias

    CIRIvis(cirias, refgenome)

    id_bam_reads_ciri = Hisat2.out.bam
		.combine(trimmed_reads, by: 0)
        .combine(ciri, by:0)
   
    // id_bam_reads_ciri.view()

    CIRIquant(id_bam_reads_ciri, refgenome, TACO, bwaindex, hisat2_idx)
    ciriquant = CIRIquant.out.ciriquant
    PrepDE_circRNA_Count_gene_reads(CIRIquant.out.ciriquant_res.collect() )

    Bedtools(refgenome, ciriquant)
    getfasta = Bedtools.out.getfasta

    Miranda(getfasta, ref_miRNA)

    // Rmats(TACO, as_files)

    Asgal(refgenome,TACO, exon, cleaned_reads)
}   
