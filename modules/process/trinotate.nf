// signalp_path = params.signalp_path
// tmhmm_path = params.tmhmm_path
// rnammer_path = params.rnammer_path

process Search_Transcripts {
        label "bigCPU"
		
        publishDir "${params.outdir}/trinotate/", mode: 'link',
            saveAs: {filename -> 
                if (filename =~/blastx/) "blastx/${filename}"
                
                else null
                }
        
        input:
        path("Transcripts.fasta")
        path(makeblastu)
        
        
        output:
        path("blastx.outfmt6"), emit: blastx
        
        
        
        script:

        """
        set +u; source activate pipeone_nm; set -u
        blastx -query Transcripts.fasta -db blast_uniprot_sprot -max_target_seqs 1 -outfmt 6 \\
        -evalue 1e-3 -out blastx.outfmt6 -num_threads ${task.cpus }

        
        """
}

process Blastp_predict {
	label "bigCPU"	
        publishDir "${params.outdir}/trinotate/", mode: 'link',
            saveAs: {filename -> 
                if (filename =~/blastp/) "blastp/${filename}"
                else if (filename =~/log/) "logs/${filename}"
                else null
                }
        
        input:
        path("Transcripts.fasta.transdecoder.pep")
        path(makeblastu)
        
        
        output:
        path("blastp_predict.outfmt6"), emit: blastp_predict
        
        
        
        script:

        """
        set +u; source activate pipeone_nm; set -u
        blastp -query Transcripts.fasta.transdecoder.pep -db blast_uniprot_sprot \\
        -max_target_seqs 1 -evalue 1e-3 -outfmt 6 -out blastp_predict.outfmt6 -num_threads ${task.cpus }

        
        """
}

process Hmmscan_Predict {
		
        publishDir "${params.outdir}/trinotate/", mode: 'link',
            saveAs: {filename -> 
                if (filename =~/out/) "hmmscan/${filename}"
                
                else null
                }
        
        input:
        path("Transcripts.fasta.transdecoder.pep")
        path(press)
        
        
        output:
        path("TrinotatePFAM.out"), emit: TPFAM
        
        
        
        script:

        """
        set +u; source activate pipeone_nm; set -u
        hmmscan --cpu ${task.cpus } --domtblout TrinotatePFAM.out Pfam-A.hmm \\
        Transcripts.fasta.transdecoder.pep > /dev/null

        
        """
}

process Signalp {
	label "local"
        publishDir "${params.outdir}/trinotate/", mode: 'link',
            saveAs: {filename -> 
                if (filename =~/signalp/) "signalp/${filename}"
                else if (filename =~/log/) "logs/${filename}"
                else null
                }
        
        input:
        path("Transcripts.fasta.transdecoder.pep")
        
        
        
        output:
        path("signalp.out"), emit: signalp
        path("signalp.log")
        
        
        
        script:

        """
       
        signalp -f short -n signalp.out -T . Transcripts.fasta.transdecoder.pep > signalp.log 2>&1
        
        
        """
}

process Tmhmm {
	label "local"
        publishDir "${params.outdir}/trinotate/", mode: 'link',
            saveAs: {filename -> 
                if (filename =~/tmhmm/) "tmhmm/${filename}"
                
                else null
                }
        
        input:
        path("Transcripts.fasta.transdecoder.pep")
        
        
        
        output:
        path("tmhmm.out"), emit:tmhmm
        
        
        
        script:

        """
       
        cat Transcripts.fasta.transdecoder.pep | tmhmm --short > tmhmm.out
        
        
        """
}

process Rnammer {
	label "local"
        label "bigCPU"
        label "bigMEM"

        publishDir "${params.outdir}/trinotate/", mode: 'link',
            saveAs: {filename -> 
                if (filename =~/gff/) "rnammer/${filename}"
                else if (filename =~/log/) "logs/${filename}"
                else null
                }
        
        input:
        path("Transcripts.fasta")
        
        
        
        output:
        path("Transcripts.fasta.rnammer.gff"), emit: rnammer
        path("rnammer.log"), emit:log
        
        
        
        script:

        """
        set +e
        rnammer_path=\$(which rnammer)
        RnammerTranscriptome.pl --transcriptome Transcripts.fasta --path_to_rnammer \$rnammer_path > rnammer.log 2>&1
        set -e
        if [ ! -f "Transcripts.fasta.rnammer.gff" ];then
                touch Transcripts.fasta.rnammer.gff
        fi
        """
}

process Trinotate {

	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'terminate' }		
        publishDir "${params.outdir}/trinotate/", mode: 'copy',
            saveAs: {filename -> 
                if (filename =~/log/) "logs/${filename}"
                else "${filename}"}
        
        input:
        path("Trinotate_20210616.sqlite")
        path("Transcripts.fasta.gene_trans_map")
        path("Transcripts.fasta")
        path("Transcripts.fasta.transdecoder.pep")
        path("blastp_predict.outfmt6")
        path("blastx.outfmt6")
        path("TrinotatePFAM.out")
        path("tmhmm.out")
        path("signalp.out")
        path("Transcripts.fasta.rnammer.gff")

        output:
        path("Trinotate.xls"), emit: trinotatefile
        tuple path("init.log"), path("load_blastp.log"), path("load_blastx.log"), path("load_pfam.log"), path("load_tmhmm.log"), path("load_signalp.log"), path("load_rnammer.log"), path("report.log")
        path("go_annotations.txt")
        
        script:

        """
        set +u; source activate pipeone_nm_4
        
        Trinotate Trinotate_20210616.sqlite init --gene_trans_map Transcripts.fasta.gene_trans_map --transcript_fasta Transcripts.fasta --transdecoder_pep Transcripts.fasta.transdecoder.pep > init.log 2>&1
        
        Trinotate Trinotate_20210616.sqlite LOAD_swissprot_blastp blastp_predict.outfmt6 > load_blastp.log 2>&1

        Trinotate Trinotate_20210616.sqlite LOAD_swissprot_blastx blastx.outfmt6 > load_blastx.log 2>&1

        Trinotate Trinotate_20210616.sqlite LOAD_pfam TrinotatePFAM.out > load_pfam.log 2>&1

        Trinotate Trinotate_20210616.sqlite LOAD_tmhmm tmhmm.out > load_tmhmm.log 2>&1

        Trinotate Trinotate_20210616.sqlite LOAD_signalp signalp.out > load_signalp.log 2>&1

        if [ ! -s "Transcripts.fasta.rnammer.gff" ];then
                
                Trinotate Trinotate_20210616.sqlite LOAD_rnammer Transcripts.fasta.rnammer.gff > load_rnammer.log 2>&1
        else
                echo "no rRNA detected"> load_rnammer.log
        fi
       

        Trinotate Trinotate_20210616.sqlite report --incl_pep --incl_trans > Trinotate.xls 2> report.log
        extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls Trinotate.xls -G --include_ancestral_terms > go_annotations.txt
        source activate pipeone_nm
        """
}



process Extract_GO {
		
        publishDir "${params.outdir}/trinotate/", mode: 'link',
            saveAs: {filename -> 
                if (filename =~/go/) "go/${filename}"
                else null}
        
        input:
        path("Trinotate.xls")

        output:
        path("go_annotations.txt"), emit: go
        path("go_forStats.txt")
        
        
        script:

        """
        set +u; source activate pipeone_nm; set -u
        extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls Trinotate.xls -G \\
        --include_ancestral_terms > go_annotations.txt

        python3 ${baseDir}/bin/create_go_stats.py

        
        """
}

process RNAs {
		
        publishDir "${params.outdir}/trinotate/", mode: 'link',
            saveAs: {filename -> 
                if (filename =~/RNA/) "RNAs/${filename}"
                else null}
        
        input:
        path("Trinotate.xls")

        output:
        path("Trinotate_rRNA.txt"), emit: rRNA
        path("Trinotate_lncRNA.txt"), emit: lncRNA
        path("Trinotate_mRNA.txt"), emit: mRNA
        
        
        script:

        """
        set +u; source activate pipeone_nm; set -u
        python3 ${baseDir}/bin/trinotateRNAs.py
        """
}

process Extract_seqs_by_IDs {
		
        publishDir "${params.outdir}/trinotate/", mode: 'link',
            saveAs: {filename -> 
                if (filename =~/RNA/) "RNAs/${filename}"
                else null}
        
        input:
        path("Transcripts.fasta")
        path("Trinotate_rRNA.txt")
        path("Trinotate_lncRNA.txt")
        path("Trinotate_mRNA.txt")

        output:
        path("Trinotate_lncRNA.fa")
        path("Trinotate_mRNA.fa")
        path("Trinotate_rRNA.fa")
        path("Trinotate_lncRNA.fa_n50_stat")
        path("Trinotate_mRNA.fa_n50_stat")
        path("Trinotate_rRNA.fa_n50_stat")
        
        
        script:

        """
        
        
        perl ${baseDir}/bin/extract_seqs_by_IDs.pl -f Transcripts.fasta -i Trinotate_lncRNA.txt -o Trinotate_lncRNA.fa

        perl ${baseDir}/bin/N50Stat.pl -i Trinotate_lncRNA.fa

        perl ${baseDir}/bin/extract_seqs_by_IDs.pl -f Transcripts.fasta -i Trinotate_mRNA.txt -o Trinotate_mRNA.fa

        perl ${baseDir}/bin/N50Stat.pl -i Trinotate_mRNA.fa

        perl ${baseDir}/bin/extract_seqs_by_IDs.pl -f Transcripts.fasta -i Trinotate_rRNA.txt -o Trinotate_rRNA.fa

        perl ${baseDir}/bin/N50Stat.pl -i Trinotate_rRNA.fa

        
        """
}
