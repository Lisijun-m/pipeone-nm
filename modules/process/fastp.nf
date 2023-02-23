
process Fastp {
    tag {id}

    publishDir "${params.outdir}/fastp/", mode: 'copy',
        saveAs: {filename -> 
            if(filename =~ /fastp.fq.gz/) "fastq/${filename}"
            else "report/${id}.${filename}" 
            }

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_1.fastp.fq"), path("${id}_2.fastp.fq"), emit: cleaned_reads
    tuple val(id), path("${id}_fastp.json"), path("${id}_fastp.html"), path("${id}_fastp.log"), emit: log
    

    script:
    """
    set +u; source activate pipeone_nm; set -u
    fastp -i ${id}_1.fastq -I ${id}_2.fastq -o ${id}_1.fastp.fq -O ${id}_2.fastp.fq \\
    -f 10 -t 2 -q 15 -u 20 -n 0 -l 50 -w 4 \\
    -j ${id}_fastp.json -h ${id}_fastp.html -R "${id} fastp report" 2> ${id}_fastp.log

    
    """
}
