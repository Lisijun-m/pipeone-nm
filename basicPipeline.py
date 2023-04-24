
import argparse

parser = argparse.ArgumentParser(description="RNA-seq analysis pipeline for non-model organisms.")
parser.add_argument("-i", "--raw", required = True, help = "input raw RNA-seq data directory (raw data in fastq or sra format) e.g. /public5/lisj/00_raw. Attention: Please name the raw data as *_1.fastq and *_2.fastq" )
parser.add_argument("-f", "--reference", type=str, help = "reference genome fasta file e.g. /public5/lisj/Genome/fasta.fna")
parser.add_argument("-o", "--output", required = True, help = "output directory e.g. /public5/lisj/ ")
parser.add_argument("--samples", required = True, help = "Basic sample information of raw data in txt format")
parser.add_argument("--pair_end", required = False, default = True, type = bool, help = "Whether raw data is pair-end sequencing data, default= True")
parser.add_argument("--rm_rRNA", required = False, default = False, type = bool, help = "Whether raw data needs to remove rRNA sequences, default: False")
parser.add_argument("--read_length", required = False, default = 100, type= int, help = "Mean read length of raw data. Default:100 ")
parser.add_argument("--total_RNA", required = True, default = True, type= bool, help = "Whether raw data contains total RNA. If raw data contains total RNA, the detection of circRNA and lncRNA-mRNA interaction will be performed. Default = True")
parser.add_argument("--AS_transcript", required = False, type=str, default = None, nargs='+',help = "The transcript to be analyzed alternative splicing by ASGAL and the order of the analyzed sample in samples information txt. e.g. /public5/lisj/transcripts.fasta 1")
parser.add_argument("--miRNA", required = False, type= str, help = "Fasta format reference miRNA file for circRNA and miRNA interaction prediction. ")
parser.add_argument("--conda_path", required = False, type= str, help = "Conda path. e.g. /dat1/lisijun/miniconda3/ ")
parser.add_argument("--species_name", required = False, type= str, help = "Species name. e.g. grassCarp ")

args = parser.parse_args()

print("checking modules")
import glob,os,sys,time,pandas,numpy
print("imports done")

### Dump sra to fastq
wk_dir = args.raw
os.chdir(wk_dir)
for filename in glob.glob('*.sra'):
	id = os.path.splitext(filename)[0]
	if os.path.isfile(id+'_1.fastq') == False:
		print(time.asctime(time.localtime(time.time()))+' Dump sra to fastq: '+filename)
		cmd = "fasterq-dump "+filename+' -e 8 2> '+id+'_dump.log'
		print('CMD: '+cmd)
		os.system(cmd)
	else: print(filename+' already dumped to '+id+'_1.fastq and '+id+'_2.fastq')
os.chdir(args.output)

## RiboDetector
in_dir = wk_dir
if args.rm_rRNA == True and args.pair_end == True:
	for filename in glob.glob(in_dir + '/*_1.fastq'):
		basename = os.path.basename(filename)
		id = basename.split('_1.fastq')[0]
		if os.path.isfile(id+'_rmrRNA_1.fastq')==False:
			cmd = 'conda run -n ribodetector ribodetector_cpu -t 20 -l '+str(args.read_length)+' -i '+in_dir+'/'+id+'_1.fastq '+ in_dir+'/'+id+'_2.fastq -e norrna --chunk_size 256 -o '+in_dir+'/'+id+'_rmrRNA_1.fastq '+in_dir+'/'+id+'_rmrRNA_2.fastq'
			print('remove rRNA sequences from raw data in '+id+'_1.fastq and '+id+'_2.fastq')
			print(cmd)
			os.system(cmd)
		else: print(filename+' rRNA has already been removed')
if args.rm_rRNA == True and args.pair_end == False:
	for filename in glob.glob(in_dir + '/*_1.fastq'):
		basename = os.path.basename(filename)
		id = basename.split('_1.fastq')[0]
		if os.path.isfile(id+'_rmrRNA_1.fastq')==False:
			cmd = 'conda run -n ribodetector ribodetector_cpu -t 20 -l '+str(args.read_length)+' -i '+in_dir+'/'+id+'_1.fastq -e norrna --chunk_size 256 -o '+in_dir+'/'+id+'_rmrRNA_1.fastq'
			print('remove rRNA sequences from raw data in '+id+'_1.fastq and '+id+'_2.fastq')
			print(cmd)
			os.system(cmd)
		else:print(filename + ' rRNA has already been removed')

### fastp
in_dir = wk_dir
out_dir = os.path.join(args.output+'01_fastp')
if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
os.chdir(out_dir)
if args.rm_rRNA == False:
	for filename in glob.glob(in_dir+'/*_1.fastq'):
		basename = os.path.basename(filename)
		id = basename.split('_1.fastq')[0]
		if os.path.isfile(id+'_1.fq.gz') == False and args.pair_end == True:
			print(time.asctime(time.localtime(time.time()))+' Trim: '+id+'.fastq')
			cmd = 'conda run -n pipeone_nm fastp -i '+in_dir+'/'+id+'_1.fastq -I '+in_dir+'/'+id+'_2.fastq -o '+id+'_1.fq.gz -O '+id+'_2.fq.gz -f 10 -t 2 -q 15 -u 20 -n 0 -l 50 -w 4 -j '+id+'_fastp.json -h '+id+'_fastp.html -R "'+id+' fastp report" 2> '+id+'_fastp.log'
			print('CMD: '+cmd)
			os.system(cmd)
		elif os.path.isfile(id+'_1.fq.gz') == False and args.pair_end == False:
			print(time.asctime(time.localtime(time.time())) + ' Trim: ' + id + '.fastq')
			cmd = 'conda run -n pipeone_nm fastp -i ' + in_dir + '/' + id + '_1.fastq -o ' + id + '_1.fq.gz -f 10 -t 2 -q 15 -u 20 -n 0 -l 50 -w 4 -j ' + id + '_fastp.json -h ' + id + '_fastp.html -R "' + id + ' fastp report" 2> ' + id + '_fastp.log'
			print('CMD: ' + cmd)
			os.system(cmd)
		else: print('FASTQ already trimmed before: '+id+'.fastq.gz')
if args.rm_rRNA == True:
	for filename in glob.glob(in_dir+'/*_rmrRNA_1.fastq'):
		basename = os.path.basename(filename)
		id = basename.split('_rmrRNA_1.fastq')[0]
		if os.path.isfile(id+'_1.fq.gz') == False and args.pair_end == True:
			print(time.asctime(time.localtime(time.time()))+' Trim: '+id+'.fastq')
			cmd = 'conda run -n pipeone_nm fastp -i '+in_dir+'/'+id+'_rmrRNA_1.fastq -I '+in_dir+'/'+id+'_rmrRNA_2.fastq -o '+id+'_1.fq.gz -O '+id+'_2.fq.gz -f 10 -t 2 -q 15 -u 20 -n 0 -l 50 -w 4 -j '+id+'_fastp.json -h '+id+'_fastp.html -R "'+id+' fastp report" 2> '+id+'_fastp.log'
			print('CMD: '+cmd)
			os.system(cmd)
		elif os.path.isfile(id+'_1.fq.gz') == False and args.pair_end == False:
			print(time.asctime(time.localtime(time.time())) + ' Trim: ' + id + '.fastq')
			cmd = 'conda run -n pipeone_nm fastp -i ' + in_dir + '/' + id + '_rmrRNA_1.fastq -o ' + id + '_1.fq.gz -f 10 -t 2 -q 15 -u 20 -n 0 -l 50 -w 4 -j ' + id + '_fastp.json -h ' + id + '_fastp.html -R "' + id + ' fastp report" 2> ' + id + '_fastp.log'
			print('CMD: ' + cmd)
			os.system(cmd)
		else: print('FASTQ already trimmed before: '+id+'.fastq.gz')
os.chdir(args.output)

## HISAT2 INDEX
genome_dir = os.path.dirname(str(args.reference))
os.chdir(genome_dir)
if os.path.isfile('genomic.8.ht2') == False:
	cmd = 'conda run -n pipeone_nm hisat2-build -p 8 '+args.reference+' '+'genomic'
	print('Building HISAT2 Index!')
	print(cmd)
	os.system(cmd)
else:print('HISAT2 index has already been built')
os.chdir(args.output)

### HISAT2
in_dir = os.path.join(args.output+'01_fastp')
out_dir = os.path.join(args.output+'02_HISAT2')
out_dir2 = os.path.join(args.output+'B02_Unmapped')
out_dir3 = os.path.join(args.output+'03_BAM')
if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
if os.path.isdir(out_dir2) == False: os.mkdir(out_dir2)
if os.path.isdir(out_dir3) == False: os.mkdir(out_dir3)
os.chdir(out_dir)
for filename in glob.glob(in_dir+'/*_1.fq.gz'):
	basename = os.path.basename(filename)
	id = basename.split('_1.fq.gz')[0]
	if args.pair_end ==True:
		if os.path.isfile(id+'.align.log') == False:
			print(time.asctime(time.localtime(time.time()))+' Align to genome by HISAT2: '+id+'_1.fq.gz and '+id+'_2.fq.gz')
			cmd = 'conda run -n pipeone_nm hisat2 -p 8 --dta -x '+genome_dir+'/genomic -1 '+in_dir+'/'+id+'_1.fq.gz -2 '+in_dir+'/'+id+'_2.fq.gz --un-conc-gz '+out_dir2+'/'+id+'_um_%.fq.gz --new-summary --summary-file '+id+'.align.log | samtools sort -@ 8 -m 4G -o '+out_dir3+'/'+id+'.bam'
			print('CMD: '+cmd)
			os.system(cmd)
		else: print('Aligning already done before: '+id+'_1.fq.gz and '+id+'_2.fq.gz')
	if args.pair_end ==False:
		if os.path.isfile(id+'.align.log') == False:
			print(time.asctime(time.localtime(time.time()))+' Align to genome by HISAT2: '+id+'_1.fq.gz')
			cmd = 'conda run -n pipeone_nm hisat2 -p 8 --dta -x '+genome_dir+'/genomic -U '+in_dir+'/'+id+'_1.fq.gz --un-conc-gz '+out_dir2+'/'+id+'_um_%.fq.gz --new-summary --summary-file '+id+'.align.log | samtools sort -@ 8 -m 4G -o '+out_dir3+'/'+id+'.bam'
			print('CMD: '+cmd)
			os.system(cmd)
os.chdir(args.output)

### StringTie
in_dir = os.path.join(args.output+'03_BAM')
out_dir = os.path.join(args.output+'04_StringTie')
if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
for filename in glob.glob(in_dir+'/*.bam'):
	basename = os.path.basename(filename)
	id = os.path.splitext(basename)[0]
	if os.path.isfile(out_dir+'/'+id+'.gtf') == False:
		print(time.asctime(time.localtime(time.time()))+' Reconstructing transcriptome: '+id+'.bam')
		cmd = 'conda run -n pipeone_nm stringtie '+filename+' -o '+out_dir+'/'+id+'.gtf'
		print('CMD: '+cmd)
		os.system(cmd)
	else: print('Transcriptome already reconstructed before: '+filename)
os.chdir(args.output)

### TACO
in_dir = os.path.join(args.output+'04_StringTie')
out_dir = os.path.join(args.output+'05_TACO')
if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
os.chdir(out_dir)
if os.path.isfile('GTFs.txt') == False:
	print(time.asctime(time.localtime(time.time()))+' Create GTF list: GTFs.txt')
	gtf_list = open('GTFs.txt', 'w')
	for filename in glob.glob(in_dir+'/*.gtf'): gtf_list.write(filename+'\n')
else: print('GTF list already created: GTFs.txt')
if os.path.isfile('TACO/assembly.gtf') == False:
	print(time.asctime(time.localtime(time.time()))+'Merge GTFs by TACO: ')
	cmd = 'conda run -n pipeone_nm_2 taco_run -p 10 -o TACO GTFs.txt 2> TACO.log'
	print('CMD: '+cmd)
	os.system(cmd)
else: print('TACO result already exists!')
if os.path.isfile('Transcripts.fasta') == False:
	print(time.asctime(time.localtime(time.time()))+' Extract transcripts sequences:')
	cmd = 'conda run -n pipeone_nm gffread -g '+args.reference+' -w Transcripts.fasta -M -K --table @geneid TACO/assembly.gtf'
	print('CMD: '+cmd)
	os.system(cmd)
else: print('Transcripts sequences already extracted!')
if os.path.isfile('Transcripts.fasta.gene_trans_map') == False:
	print(time.asctime(time.localtime(time.time()))+' Create gene-trans map:')
	fa = open ('Transcripts.fasta')
	gtm = open ('Transcripts.fasta.gene_trans_map', 'w')
	for line in fa.readlines():
		if line[0] == '>':
			line = line.lstrip('>')
			line = line.rstrip('\n')
			line = line.split('\t')
			tran = line[0]
			gene = line[1]
			gtm.write(gene+'\t'+tran+'\n')
else: print('Gene-trans map already extracted!')
os.chdir(args.output)

## Salmon
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

in_dir = os.path.join(args.output+'05_TACO')
out_dir = os.path.join(args.output+'06_Quantify')
trinity_dir = os.path.join(script_dir+'/trinityrnaseq-Trinity-v2.13.2/util/')

if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
os.chdir(out_dir)
if os.path.isdir('Salmon') == False: os.mkdir('Salmon')
os.chdir('Salmon')
if os.path.isfile('align_and_estimate_abundance.log') == False:
	print(time.asctime(time.localtime(time.time()))+' Run Salmon: ')
	cmd = trinity_dir+'align_and_estimate_abundance.pl --transcripts '+in_dir+'/Transcripts.fasta --seqType fq --gene_trans_map '+in_dir+'/Transcripts.fasta.gene_trans_map --samples_file '+args.samples+' --est_method salmon --thread_count 32 --prep_reference > align_and_estimate_abundance.log 2>&1'
	print('CMD: '+cmd)
	os.system(cmd)
else: print('Please check the Salmon output folder!')
os.chdir(out_dir)
if os.path.isfile('Transcripts.gene.TMM.EXPR.matrix') == False:
	print(time.asctime(time.localtime(time.time()))+' Get matrix: ')
	cmd = trinity_dir+'abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map '+in_dir+'/Transcripts.fasta.gene_trans_map --out_prefix Transcripts --name_sample_by_basedir Salmon/*/quant.sf 2> matrix.log'
	print('CMD: '+cmd)
	os.system(cmd)
else: print('Matrix already created from abundance estimates!')
os.chdir(out_dir)


### TransDecoder
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
db_dir = os.path.join(script_dir+'/Trinotate_Database/')

in_dir = os.path.join(args.output+'05_TACO')
out_dir = os.path.join(args.output+'07_TransDecoder')
if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
os.chdir(out_dir)

if os.path.isfile('Transcripts.fasta.transdecoder_dir/longest_orfs.pep') == False:
	print(time.asctime(time.localtime(time.time()))+' Extract the long open reading frames:')
	cmd = 'conda run -n pipeone_nm TransDecoder.LongOrfs -t '+in_dir+'/Transcripts.fasta > TransDecoder.LongOrfs.log 2>&1'
	print('CMD: '+cmd)
	os.system(cmd)
else: print('Long ORFs already extracted before!')

if os.path.isfile('blastp.outfmt6') == False:
	print(time.asctime(time.localtime(time.time()))+' Identify ORFs with homology to known proteins via blast:')
	cmd = 'conda run -n pipeone_nm blastp -query Transcripts.fasta.transdecoder_dir/longest_orfs.pep -db '+db_dir+'uniprot_sprot.pep -out blastp.outfmt6 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 32'
	print('CMD: '+cmd)
	os.system(cmd)
else: print('Long ORFs already identified via blast before!')

if os.path.isfile('pfam.domtblout') == False:
	print(time.asctime(time.localtime(time.time()))+' Identify ORFs with homology to known proteins via pfam:')
	cmd = 'hmmscan --cpu 8 --domtblout pfam.domtblout '+db_dir+'Pfam-A.hmm Transcripts.fasta.transdecoder_dir/longest_orfs.pep > /dev/null'
	print('CMD: '+cmd)
	os.system(cmd)
else: print('Long ORFs already identified via pfam before!')

if os.path.isfile('Transcripts.fasta.transdecoder.pep') == False:
	print(time.asctime(time.localtime(time.time()))+' Predict the long open reading frames:')
	cmd = 'conda run -n pipeone_nm TransDecoder.Predict -t '+in_dir+'/Transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 1> TransDecoder.Predict.log 2>&1'
	print('CMD: '+cmd)
	os.system(cmd)
else: print('Long ORFs already predicted before!')
os.chdir(args.output)

### Trinotate
trinotate_soft_path = os.path.join(script_dir+'/Trinotate/Trinotate')
trinotate_soft_dir = os.path.join(script_dir+'/Trinotate')
in_dir = os.path.join(args.output+'05_TACO')
in_dir2 = os.path.join(args.output+'07_TransDecoder')
out_dir = os.path.join(args.output+'08_Trinotate')
if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
os.chdir(out_dir)

if os.path.isfile('blastx.outfmt6') == False:
	print(time.asctime(time.localtime(time.time()))+' Search transcripts:')
	cmd = 'conda run -n pipeone_nm blastx -query '+in_dir+'/Transcripts.fasta -db '+db_dir+'uniprot_sprot.pep -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -out blastx.outfmt6 -num_threads 32'
	print('CMD: '+cmd)
	os.system(cmd)
else: print('Transcripts already searched before!')

if os.path.isfile('blastp.outfmt6') == False:
	print(time.asctime(time.localtime(time.time()))+' Search Transdecoder-predicted proteins:')
	cmd = 'conda run -n pipeone_nm blastp -query '+in_dir2+'/Transcripts.fasta.transdecoder.pep -db '+db_dir+'uniprot_sprot.pep -max_target_seqs 1 -evalue 1e-3 -outfmt 6 -out blastp.outfmt6 -num_threads 32'
	print('CMD: '+cmd)
	os.system(cmd)
else: print('Transdecoder-predicted proteins already searched before!')

if os.path.isfile('TrinotatePFAM.out') == False:
	print(time.asctime(time.localtime(time.time()))+' Identify protein domains by HMMER:')
	cmd = 'hmmscan --cpu 8 --domtblout TrinotatePFAM.out '+db_dir+'Pfam-A.hmm '+in_dir2+'/Transcripts.fasta.transdecoder.pep > /dev/null'
	print('CMD: '+cmd)
	os.system(cmd)
else: print('Protein domain already identified before!')

if os.path.isfile('signalp.out') == False:
	print(time.asctime(time.localtime(time.time()))+' Predict signal peptides by SingalP:')
	cmd = 'signalp -f short -n signalp.out -T . '+in_dir2+'/Transcripts.fasta.transdecoder.pep > signalp.log 2>&1'
	print('CMD: '+cmd)
	os.system(cmd)
else: print('Signal peptides already predicted before!')

if os.path.isfile('tmhmm.out') == False:
	print(time.asctime(time.localtime(time.time()))+' Predict transmembrane regions by tmHMM:')
	cmd = 'tmhmm --short < '+in_dir2+'/Transcripts.fasta.transdecoder.pep > tmhmm.out'
	print('CMD: '+cmd)
	os.system(cmd)
else: print('Transmembrane regions already predicted before!')

rnammer = os.path.join(script_dir+'/rnammer-1.2/rnammer')
if os.path.isfile('Transcripts.fasta.rnammer.gff') == False:
	print(time.asctime(time.localtime(time.time()))+' Predict rRNA transcripts by RNAMMER:')
	cmd = trinotate_soft_dir+'/util/rnammer_support/RnammerTranscriptome.pl --transcriptome '+in_dir+'/Transcripts.fasta --path_to_rnammer '+rnammer+' > rnammer.log 2>&1'
	print('CMD: '+cmd)
	os.system(cmd)
else: print('rRNA transcripts already predicted before!')

if os.path.isfile('Trinotate.xls') == False:
	print(time.asctime(time.localtime(time.time()))+' Copy original Trinotate_20210616.sqlite:')
	cmd = 'cp '+db_dir+'Trinotate.sqlite .'
	print('CMD: '+cmd)
	os.system(cmd)

	print(time.asctime(time.localtime(time.time()))+' Load transcripts and coding regions:')
	cmd = trinotate_soft_path+' Trinotate.sqlite init --gene_trans_map '+in_dir+'/Transcripts.fasta.gene_trans_map --transcript_fasta '+in_dir+'/Transcripts.fasta --transdecoder_pep '+in_dir2+'/Transcripts.fasta.transdecoder.pep > init.log 2>&1'
	print('CMD: '+cmd)
	os.system(cmd)

	print(time.asctime(time.localtime(time.time()))+' Load protein hits:')
	cmd = trinotate_soft_path+' Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6 > load_blastp.log 2>&1'
	print('CMD: '+cmd)
	os.system(cmd)

	print(time.asctime(time.localtime(time.time()))+' Load transcript hits:')
	cmd = trinotate_soft_path+' Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6 > load_blastx.log 2>&1'
	print('CMD: '+cmd)
	os.system(cmd)

	print(time.asctime(time.localtime(time.time()))+' Load Pfam domain entries:')
	cmd = trinotate_soft_path+' Trinotate.sqlite LOAD_pfam TrinotatePFAM.out > load_pfam.log 2>&1'
	print('CMD: '+cmd)
	os.system(cmd)

	print(time.asctime(time.localtime(time.time()))+' Load transmembrane domains:')
	cmd = trinotate_soft_path+' Trinotate.sqlite LOAD_tmhmm tmhmm.out > load_tmhmm.log 2>&1'
	print('CMD: '+cmd)
	os.system(cmd)

	print(time.asctime(time.localtime(time.time()))+' Load signal peptide predictions:')
	cmd = trinotate_soft_path+' Trinotate.sqlite LOAD_signalp signalp.out > load_signalp.log 2>&1'
	print('CMD: '+cmd)
	os.system(cmd)

	print(time.asctime(time.localtime(time.time()))+' Load rRNA predictions:')
	cmd = trinotate_soft_path+' Trinotate.sqlite LOAD_rnammer Transcripts.fasta.rnammer.gff > load_rnammer.log 2>&1'
	print('CMD: '+cmd)
	os.system(cmd)

	print(time.asctime(time.localtime(time.time()))+' Generate report:')
	cmd = trinotate_soft_path+' Trinotate.sqlite report --incl_pep --incl_trans > Trinotate.xls 2> report.log'
	print('CMD: '+cmd)
	os.system(cmd)
else: print('Trinotate.xls already generated before!')

if os.path.isfile('go_annotations.txt') == False:
	cmd = trinotate_soft_dir+'/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls Trinotate.xls -G --include_ancestral_terms > go_annotations.txt'
	print('CMD: '+cmd)
	os.system(cmd)

if os.path.isfile('go_forStats.txt')==False:
	f1 = open('go_annotations.txt', 'r')
	f2 = open('go_forStats.txt', 'w')
	for i in f1.readlines():
		j = i.split('\t')
		for k in j[1].split(','):
			m = j[0] + '\t' + k
			if (m[-1] != "\n"):
				m = m + "\n"
			print(m)
			f2.write(m)
	f1.close()
	f2.close()

if os.path.isfile('go_anno.txt')==False:
	file = pandas.read_table('Trinotate.xls')
	out_file = open('go_anno.txt','w')
	for i in range(len(file)):
		if file['gene_ontology_BLASTX'][i] != '.' :
			out_file.write(file['gene_ontology_BLASTX'][i]+"\n")
		if file['gene_ontology_BLASTP'][i] !='.':
			out_file.write(file['gene_ontology_BLASTP'][i]+"\n")
		if file['gene_ontology_Pfam'][i] !='.':
			out_file.write(file['gene_ontology_Pfam'][i]+"\n")
	out_file.close()

if os.path.isfile('go_anno_2.txt')==False:
	f1 = open('go_anno.txt', 'r')
	f2 = open('go_anno_2.txt', 'w')
	for i in f1.readlines():
		for j in i.split('`'):
			k = j.split('^')
			m = k[0] + '\t' + k[1] + '\t' + k[2]
			if (m[-1] != "\n"):
				m = m + "\n"
			print(m)
			f2.write(m)
	f1.close()
	f2.close()

if os.path.isfile('go_anno_3.txt') == False:
	fi = open('go_anno_2.txt', 'r')
	txt = fi.readlines()
	with open('go_anno_3.txt', 'a') as f:
		f.close()
	for w in txt:
		fi2 = open('go_anno_3.txt', 'r')
		txt2 = fi2.readlines()
		with open('go_anno_3.txt', 'a') as f:
			if w not in txt2:
				f.write(w)
			f.close()
	fi.close()

if os.path.isfile('Trinotate_lncRNA.txt') == False:
	mRNA_list = []
	tn = open('Trinotate.xls')
	rRNA = open('Trinotate_rRNA.txt', 'w')
	lncRNA = open('Trinotate_lncRNA.txt', 'w')
	tn.readline()
	for line in tn.readlines():
		unit = line.split('\t')
		if unit[3] != '.': rRNA.write(unit[1]+'\n')
		elif unit[2] == '.' and unit[16] == '.\n': lncRNA.write(unit[1]+'\n')
		else: mRNA_list.append(unit[1])
	tn.close()
	rRNA.close()
	lncRNA.close()
	mRNA_list = list(set(mRNA_list))
	mRNA = open('Trinotate_mRNA.txt', 'w')
	for G in mRNA_list: mRNA.write(G+'\n')
	mRNA.close()

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
pl_script_dir = os.path.join(script_dir+'/scripts/')

out_file = 'Trinotate_lncRNA.fa'
if os.path.isfile(out_file) == False:
	cmd = pl_script_dir+'extract_seqs_by_IDs.pl -f '+in_dir+'/Transcripts.fasta -i Trinotate_lncRNA.txt -o '+out_file
	print('CMD: '+cmd)
	os.system(cmd)
if os.path.isfile(out_file+'_n50_stat') == False:
	cmd = pl_script_dir+'N50Stat.pl -i '+out_file
	print('CMD: '+cmd)
	os.system(cmd)
	
out_file = 'Trinotate_mRNA.fa'
if os.path.isfile(out_file) == False:
	cmd = pl_script_dir+'extract_seqs_by_IDs.pl -f '+in_dir+'/Transcripts.fasta -i Trinotate_mRNA.txt -o '+out_file
	print('CMD: '+cmd)
	os.system(cmd)
if os.path.isfile(out_file+'_n50_stat') == False:
	cmd = pl_script_dir+'N50Stat.pl -i '+out_file
	print('CMD: '+cmd)
	os.system(cmd)

out_file = 'Trinotate_rRNA.fa'
if os.path.isfile(out_file) == False:
	cmd = pl_script_dir+'extract_seqs_by_IDs.pl -f '+in_dir+'/Transcripts.fasta -i Trinotate_rRNA.txt -o '+out_file
	print('CMD: '+cmd)
	os.system(cmd)
if os.path.isfile(out_file+'_n50_stat') == False:
	cmd = pl_script_dir+'N50Stat.pl -i '+out_file
	print('CMD: '+cmd)
	os.system(cmd)
os.chdir(args.output)

### BWA
if args.total_RNA == True:
	os.chdir(genome_dir)
	cmd = 'bwa index '+args.reference
	os.system('CMD: '+cmd)
	print('Generating bwa index for reference genome')
	in_dir = os.path.join(args.output+'01_fastp')
	out_dir = os.path.join(args.output+'SAM')
	if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
	os.chdir(out_dir)
	if args.rm_rRNA == True and args.pair_end == True:
		for filename in glob.glob(in_dir + '/*_rmrRNA_1.fastq'):
			basename = os.path.basename(filename)
			id = basename.split('_rmrRNA_1.fastq')[0]
			if os.path.isfile(id+'.sam')==False:
				cmd = 'conda run -n pipeone_nm bwa mem -T 19 -t 10 '+args.reference+' '+in_dir+'/'+id+'_rmrRNA_1.fastq '+in_dir+'/'+id+'_rmrRNA_2.fastq > '+id+'.sam'
				print(cmd)
				os.system(cmd)
			else:print(id+'.sam has already been created! ')
	if args.rm_rRNA == True and args.pair_end == False:
		for filename in glob.glob(in_dir + '/*_rmrRNA_1.fastq'):
			basename = os.path.basename(filename)
			id = basename.split('_rmrRNA_1.fastq')[0]
			if os.path.isfile(id+'.sam')==False:
				cmd = 'conda run -n pipeone_nm bwa mem -T 19 -t 10 '+args.reference+' '+in_dir+'/'+id+'_rmrRNA_1.fastq > '+id+'.sam'
				print(cmd)
				os.system(cmd)
			else:print(id+'.sam has already been created! ')
	if args.rm_rRNA == False and args.pair_end == True:
		for filename in glob.glob(in_dir + '/*_1.fastq'):
			basename = os.path.basename(filename)
			id = basename.split('_1.fastq')[0]
			if os.path.isfile(id+'.sam')==False:
				cmd = 'conda run -n pipeone_nm bwa mem -T 19 -t 10 '+args.reference+' '+in_dir+'/'+id+'_1.fastq '+in_dir+'/'+id+'_2.fastq > '+id+'.sam'
				print(cmd)
				os.system(cmd)
			else:print(id + '.sam has already been created! ')
	if args.rm_rRNA == False and args.pair_end == False:
		for filename in glob.glob(in_dir + '/*_1.fastq'):
			basename = os.path.basename(filename)
			id = basename.split('_1.fastq')[0]
			if os.path.isfile(id + '.sam') == False:
				cmd = 'conda run -n pipeone_nm bwa mem -T 19 -t 10 '+args.reference+' '+in_dir+'/'+id+'_1.fastq > '+id+'.sam'
				print(cmd)
				os.system(cmd)
			else:print(id + '.sam has already been created! ')

## circRNA detection
if args.total_RNA == True:
	in_dir = os.path.join(args.output+'SAM')
	out_dir = os.path.join(args.output+'CIRI')
	gtf_dir = os.path.join(args.output+'05_TACO/TACO')
	if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
	cmd= 'chmod 777 '+out_dir
	print('CMD: '+cmd)
	os.system(cmd)
	os.chdir(out_dir)
	for filename in glob.glob(in_dir+'/*.sam'):
		basename = os.path.basename(filename)
		id = os.path.splitext(basename)[0]
		if os.path.isfile(id+'.ciri') == False:
			print(time.asctime(time.localtime(time.time()))+' Detecting circRNA using CIRI2: '+id+'.sam')
			cmd = 'conda run -n CIRI perl '+script_dir+'/CIRI-full_v2.0/bin/CIRI_v2.0.6/CIRI2.pl -T 16 -I '+in_dir+'/'+id+'.sam -O '+id+'.ciri -F '+args.reference+' -A '+gtf_dir+'/assembly.gtf'
			print('CMD: '+cmd)
			os.system(cmd)
		else: print('circRNA already detected before: '+filename)
	for filename in glob.glob(in_dir+'/*.sam'):
		basename = os.path.basename(filename)
		id = os.path.splitext(basename)[0]
		if os.path.isfile(id+'.as_AS.list') == False:
			print(time.asctime(time.localtime(time.time()))+' Detecting circRNA alternative splicing using CIRI2: '+id+'.sam')
			cmd = 'conda run -n CIRI perl '+script_dir+'/CIRI-full_v2.0/bin/CIRI_AS_v1.2/CIRI_AS_v1.2.pl -D yes -T 16 -S '+in_dir+'/'+id+'.sam -O '+id+'.as'+' -C '+id+'.ciri -F '+args.reference+' -A '+gtf_dir+'/assembly.gtf'
			print('CMD: '+cmd)
			os.system(cmd)
		else: print('circRNA alternative splicing already detected before: '+filename)
	for filename in glob.glob(in_dir+'/*.sam'):
		basename = os.path.basename(filename)
		id = os.path.splitext(basename)[0]
		if os.path.isfile(id+'/'+id+'.list') == False:
			print(time.asctime(time.localtime(time.time()))+' Visualizing circRNA detection using CIRI2_vis: '+id+'.sam')
			cmd = 'conda run -n CIRI java -jar '+script_dir+'/CIRI-full_v2.0/CIRI-vis.jar -i '+id+'.as_jav.list -l '+id+'.as_library_length.list -d '+id+' -r '+args.reference+' -type pdf -o '+id
			print('CMD: '+cmd)
			os.system(cmd)
		else: print('Visulization of circRNA done before: '+filename)
	os.chdir(args.output)

## circRNA Quantification
if args.total_RNA ==True:
	yaml_path = os.path.join(args.output+'CIRI-quantContig.yaml')
	yaml_file = open(yaml_path,'a')
	yaml_file.write('names: '+args.species_name+"\n")
	yaml_file.write('tools:'+"\n")
	yaml_file.write(' bwa: '+args.conda_path+'envs/pipeone_nm/bin/bwa'+"\n")
	yaml_file.write(' hisat2: ' + args.conda_path + 'envs/pipeone_nm/bin/hisat2' + "\n")
	yaml_file.write(' srtingtie: ' + args.conda_path + 'envs/pipeone_nm/bin/stringtie' + "\n")
	yaml_file.write(' samtools: ' + args.conda_path + 'envs/pipeone_nm/bin/samtools' + "\n")
	yaml_file.write("\n")
	yaml_file.write('reference: ' + "\n")
	yaml_file.write(' fasta: ' +args.reference + "\n")
	yaml_file.write(' gtf: ' + os.path.join(args.output+'05_TACO/TACO/assembly.gtf') + "\n")
	yaml_file.write(' bwa_index: ' +args.reference + "\n")
	yaml_file.write(' hisat_index: ' +os.path.join(genome_dir+'/genomic')+ " \n")
	yaml_file.close()
	in_dir = os.path.join(args.output+'01_fastp')
	circ_dir = os.path.join(args.output+'CIRI')
	out_dir = os.path.join(args.output+'CIRIquantify')
	bam_dir = os.path.join(args.output+'03_BAM')
	if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
	os.chdir(out_dir)
	if args.pair_end == True:
		for filename in glob.glob(in_dir + '/*_1.fq.gz'):
			basename = os.path.basename(filename)
			id = basename.split('_1.fq.gz')[0]
			if os.path.isfile(id + '/' + id + '.log') == False:
				print(time.asctime(time.localtime(
					time.time())) + ' quantify CIRI detected circRNA: ' + id + '_1.fq.gz and ' + id + '_2.fq.gz')
				cmd = 'conda run -n CIRI CIRIquant --config ' + yaml_path + ' -t 16 -1 ' + in_dir + '/' + id + '_1.fq.gz -2 ' + in_dir + '/' + id + '_2.fq.gz -o ' + id + ' -p ' + id + ' --circ ' + circ_dir + '/' + id + '.ciri --bam ' + bam_dir + '/' + id + '.bam --tool CIRI2'
				print('CMD: ' + cmd)
				os.system(cmd)
			else: print('CIRI already quanted before: ' + id + '_1.fq.gz and ' + id + '_2.fq.gz')
		os.chdir(args.output)

## circRNA-miRNA interaction prediction
if args.total_RNA ==True:
	file_dir = os.path.join(args.output+'SAM')
	out_dir = os.path.join(args.output+'FASTA')
	bed_dir = os.path.join(args.output+'CIRIquantify')
	if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
	os.chdir(out_dir)
	for filename in glob.glob(file_dir + '/*.sam'):
		basename = os.path.basename(filename)
		id = os.path.splitext(basename)[0]
		if os.path.isfile(id + '.fasta') == False:
			print(time.asctime(time.localtime(time.time())) + ' Exstracting FASTA from BED: ' + id + '.bed')
			cmd = 'conda run -n pipeone_nm bedtools getfasta -fi ' + args.reference + ' -bed ' + bed_dir + '/' + id + '/' + id + '.bed -fo '+out_dir+'/'+id + '.fasta'
			print('CMD: ' + cmd)
			os.system(cmd)
		else:
			print('FASTA already extracted before: ' + filename)
	os.chdir(args.output)
	in_dir = os.path.join(args.output+'FASTA')
	mi_dir = os.path.join(args.output+'miRanda')
	os.chdir(mi_dir)
	for filename in glob.glob(in_dir + '/*.fasta'):
		basename = os.path.basename(filename)
		id = os.path.splitext(basename)[0]
		if os.path.isfile(id + '_circ_miRNA.txt') == False:
			print(time.asctime(time.localtime(time.time())) + ' Analyzing circRNA targeted miRNA: ' + id + '.fasta')
			cmd = 'miranda ' + args.miRNA +' '+ in_dir+'/' + id + '.fasta > '+mi_dir+'/' + id + '_circ_miRNA.txt'
			print('CMD: ' + cmd)
			os.system(cmd)
		else:
			print('Interaction of circRNA and miRNA have already analyzed before: ' + id + '.fasta')
	os.chdir(args.output)


## Alternative Splicing
def getline(the_file_path, line_number):
	if line_number < 1:
		return ''
	for cur_line_number, line in enumerate(open(the_file_path, 'rU')):
		if cur_line_number == line_number-1:
			return line
	return ''

if args.AS ==True:
	out_dir = os.path.join(args.output+'AS')
	GTF_file = os.path.join(args.output+'05_TACO/TACO/assembly.gtf')
	if os.path.isdir(out_dir) == False: os.mkdir(out_dir)
	os.chdir(out_dir)
	ASGAL_dir = os.path.join(script_dir+'galig/')
	transcript_fasta = args.AS_transcript[0]
	object_AS = args.AS_transcript[1:len(args.AS_transcript)]
	sample_path = args.samples
	for i in object_AS:
		line = getline(sample_path,int(i)).strip("\n").split('\t')
		if args.pair_end == True:
			cmd ='conda run -n asgal python3 ' +ASGAL_dir +'asgal --multi -g '+args.reference+' -a '+GTF_file+' -s '+line[2]+' -s2 '+line[3]+' -t '+args.AS_transcript[0]+' -o '+out_dir
		if args.pair_end == False:
			cmd ='conda run -n asgal python3 ' +ASGAL_dir +'asgal --multi -g '+args.reference+' -a '+GTF_file+' -s '+line[2]+' -t '+args.AS_transcript[0]+' -o '+out_dir

