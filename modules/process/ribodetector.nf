process Ribodetector {
    publishDir "${params.outdir}/Ribodetecto/", mode: 'copy'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_rmrRNA_*.fastq")

    script:
     if(! params.singleEnd){
        """
        conda run -n asgal ribodetector_cpu -t 20 -l ${params.read_length} -i ${reads} -e norrna --chunk_size 256 -o ${id}_rmrRNA_1.fastq ${id}_rmrRNA_2.fastq
        """
        }else{
        """
        conda run -n asgal ribodetector_cpu -t 20 -l ${params.read_length} -i ${reads} -e norrna --chunk_size 256 -o ${id}_rmrRNA_1.fastq
        """
    }

}