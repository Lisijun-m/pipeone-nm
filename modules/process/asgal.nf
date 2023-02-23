process Asgal {
    publishDir "${params.outdir}/AS_asgal", mode: 'copy'

    input:
    path "genome.fasta"
    path gtf
    path transcripts
    tuple val(id), path(reads)

    output:
    path "${id}.asgal"

    script:
   
    if(! params.singleEnd){
        """
        mkdir -p  ${id}.asgal
        python ${baseDir}/bin/gtf.py to_gene_gtf assembly.gtf gene.gtf
        cat assembly.gtf >> gene.gtf
        conda run -n asgal  asgal --multi -g genome.fasta -a  gene.gtf -s ${reads[0]} -s2 ${reads[1]}  -o ${id}.asgal -t $transcripts
        """
        }else{
        """
        mkdir -p  ${id}.asgal
        python ${baseDir}/bin/gtf.py to_gene_gtf assembly.gtf gene.gtf
        cat assembly.gtf >> gene.gtf
        conda run -n asgal  asgal --multi -g genome.fasta -a  gene.gtf -s ${reads[0]}  -o ${id}.asgal -t $transcripts
        """
        }
 
}