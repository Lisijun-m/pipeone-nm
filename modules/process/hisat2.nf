params.unstranded = true // defualt unstranded
params.reverse_stranded = true // dUTP if stranded
params.forward_stranded =false
params.singleEnd = false //
params.refgenome = ''


process Hisat2 {
		label "bigCPU"
        tag {id}
        publishDir "${params.outdir}/hisat2/", mode: 'copy',
            saveAs: {filename -> 
                if(filename =~ /bam/ )  "bam/$filename"
                else if (filename =~/log/) "logs/${filename}"
				else if (filename =~/.fq.gz/) "un_conz_fastq/${filename}"
                else null
                }
        
        input:
        tuple val(id), path(reads)
        val hisat2_base
        path "hisat2_index/*"
        
        output:
        tuple val(id), path("${id}.bam") , emit: bam
        tuple val(id), path("${id}_*.fq.gz") , emit: un_conz_fastq
        
        
        script:
        def rnastrandness = ''
        if (params.forward_stranded && !params.unstranded){
            rnastrandness =params.singleEnd ?  '--rna-strandness F' :  '--rna-strandness FR'
        } else if (params.reverse_stranded && !params.unstranded){
            rnastrandness = params.singleEnd ? '--rna-strandness R'  :  '--rna-strandness RF'
        }
        
    
        """
        set +u; source activate pipeone_nm; set -u
        hisat2 -p $task.cpus --dta $rnastrandness -x hisat2_index/$hisat2_base -1 ${reads[0]} -2 ${reads[1]} \\
        --un-conc-gz ${id}_%.fq.gz --new-summary --summary-file ${id}.align.log \\
        | samtools sort -@ 2 -o ${id}.bam
        
        """
       
 
}


process Hisat2_build {
    label "bigCPU"
    publishDir "${params.outdir}/hisat2/hisat2_index/", mode: 'copy'
      

    input:
    path("genome.fna")

    output:
    path("genomic.*.ht2") , emit: ht2

    script:
    """
    set +u; source activate pipeone_nm; set -u
    hisat2-build -p ${task.cpus } genome.fna genomic
    """

}