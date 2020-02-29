<b>RNACocktail: A comprehensive framework for accurate and efficient RNA-Seq analysis</b>

See http://bioinform.github.io/rnacocktail/ for help and downloads. 
RNACocktail
A comprehensive framework for accurate and efficient RNA-Seq analysis
The RNACocktail pipeline is composed of a high-accuracy tools for different steps of RNA-Seq analysis. It performs a broad spectrum RNA-Seq analysis on both short- and long-read technologies to enable meaningful insights from transcriptomic data. It was developed after analyzing a variety of RNA-Seq samples (ranging from germline, cancer to stem cell datasets) and technologies using a multitude of tool combinations to determine a pipeline which is comprehensive, fast and accurate.
RNACocktail supports: 
short-read
long-read
alignment
error correction
transcriptome reconstruction
alignment
denovo transcriptome assembly
transcriptome reconstruction
alignment-free quantification
fusion prediction
differential expression analysis

fusion prediction

variant calling

RNA editing prediction


For more information contact us at bioinformatics.red@roche.com 
Publication
If you use RNACocktail in your work, please cite the following:
Sayed Mohammad Ebrahim Sahraeian, Marghoob Mohiyuddin, Robert Sebra, Hagen Tilgner, Pegah T. Afshar, Kin Fai Au, Narges Bani Asadi, Mark B. Gerstein, Wing Hung Wong, Michael P. Snyder, Eric Schadt, and Hugo Y. K. Lam
Gaining comprehensive biological insight into the transcriptome by performing a broad-spectrum RNA-seq analysis
Nature Communications 8, Article number: 59 (2017). doi:10.1038/s41467-017-00050-4 
Download RNACocktail
Latest version: https://github.com/bioinform/RNACocktail/archive/v0.3.0.tar.gz
For other versions, see "releases". https://github.com/bioinform/RNACocktail/releases
RNACocktail Docker Image
The docker image with all the packages installed can be found at https://hub.docker.com/repository/docker/rssbred/rnacocktail/
Older docker images (versions < 0.2.2) can be found here.
The dockerfile is also available at docker/Dockerfile for local build.
System Requirements
The current implementation of RNACocktail is tested with Python 2.7 and the following Python packages: 
Tool
Version tested
Pipeline modes used in
pybedtools
0.8.0
editing
pysam
0.15.0
editing
numpy
1.16.5
editing
scipy
1.2.2
long_fusion
biopython
1.74
fusion
openpyxl
2.6.4
fusion
pandas
0.24.2
fusion
xlrd
1.1.0
fusion
Note that bedtools (v2.29.0) has to be installed separately in order for pybedtools to work.
In addition, paths to the following tools must be provided as RNACocktail arguments. Alternatively, the executables can be on PATH environmental variable or defined on defaults.py:
Tool
Version tested
Pipeline modes used in
SAMtools
1.2
align, reconstruct, long_align, long_reconstruct, and editing
HISAT2
2.1.0
align
StringTie
2.0.4
reconstruct and diff
Salmon
0.11.0
quantify
Oases
0.2.09
assembly
Velvet
1.2.10
assembly
R with DESeq2, readr, and tximport libraries
3.6.1
diff, editing
featureCounts
2.0.0
diff
LoRDEC
0.9
long_correct
STAR
2.7.0f
long_align
IDP
0.1.9
long_reconstruct
IDP-fusion
1.1.1
long_fusion
GATK
4.1.4.0
variant and editing
Picard
2.19.0
variant
GIREMI
0.2.1
editing
HTSlib
1.3
editing
FusionCatcher
1.10
fusion
bowtie
1.2.2
fusion
bowtie2
2.2.9
fusion, long_fusion
bwa
0.7.17
fusion
sra toolkit
2.9.6
fusion
coreutils
8.27
fusion
pigz
2.3.1
fusion
blat
0.35
fusion
faToTwoBit

fusion
liftOver

fusion
SeqTK
1.2-r101c
fusion
gmap
2019-09-12
long_fusion

Installing RNACocktail
RNACocktail is a python package and can be installed using pip. To install type pip install https://github.com/bioinform/RNACocktail/archive/v0.3.0.tar.gz. The current version of RNACocktail is v0.3.0. In general, the install source would be https://github.com/bioinform/RNACocktail/archive/version.tar.gz
Running RNACocktail
Type run_rnacocktail.py -h for help.
Type run_rnacocktail.py align -h for short-read alignment help.
Type run_rnacocktail.py reconstruct -h for short-read transcriptome reconstruction help.
Type run_rnacocktail.py quantify -h for short-read quantification help.
Type run_rnacocktail.py diff -h for short-read differential expression help.
Type run_rnacocktail.py denovo -h for short-read de novo assembly help.
Type run_rnacocktail.py long_correct -h for long-read error correction help.
Type run_rnacocktail.py long_align -h for long-read alignment help.
Type run_rnacocktail.py long_reconstruct -h for long-read transcriptome reconstruction help.
Type run_rnacocktail.py long_fusion -h for long-read fusion detection help.
Type run_rnacocktail.py variant -h for variant calling help.
Type run_rnacocktail.py editing -h for RNA editing detection help.
Type run_rnacocktail.py fusion -h for RNA fusion detection help.
Type run_rnacocktail.py all -h for running all RNACocktail pipeline steps help.
The all mode for RNACocktail will automatically perform the most comprehensive analysis possible given the input data, which includes steps from alignment to differential expression analysis.
Testing RNACoktail
Small test
cd test
./test_run.sh
Extensive test of all modes on Docker image
cd test
./docker_test.sh
Analysis scripts
Several IPython Notebook and .py scripts to analyze the predictions in different tasks can be found at analaysis_scripts folder
Output files
The table below summarizes the output files generated by each mode of RNACocktail.
Task
Command
Default Tool
Output Files
Short-read alignment
align
HISAT2
alignments: alignments.sorted.bam
junctions: splicesites.tab
Short-read transcriptome reconstruction
reconstruct
StringTie
trasncripts: transcripts.gtf
expressions: gene_abund.tab
Short-read quantification
quantify
Salmon-SMEM
expressions: quant.sf
Short-read differential expression
diff
DESeq2
differential expressions: deseq2_res.tab
Short-read de novo assembly
denovo
Oases
trasncripts: transcripts.fa
Long-read error correction
long_correct
LoRDEC
corrected reads long_corrected.fa
Long-read alignment
long_align
STARlong
alignments Aligned.out.psl
Long-read transcriptome reconstruction
long_reconstruct
IDP
trasncripts: isoform.gtf
expressions: isoform.exp
Long-read fusion detection
long_fusion
IDP-fusion
fusions: fusion_report.tsv
Variant calling
variant
GATK
variants: variants_filtered.vcf
RNA editing detection
editing
GIREMI
edits: giremi_out.txt.res
RNA Fusion detection
fusion
FusionCatcher
fusions: final-list_candidate-fusion-genes.txt
Running all steps
all
whole pipeline
all outputs of the successful steps.
Examples
Some example command-lines for running RNACocktail with various modes and data type (short- and long-reads) are shown below. In particular, examples 17 and 18 show how to use the all mode for the most comprehensive analysis. Note that RNACocktail requires pre-built indexes for the genomic and transcriptomic references.
Example 1 (align):
Run of RNACocktail for alignment of paired-end short-read sequences (HISAT2). 
run_rnacocktail.py align --align_idx hisat2-idx --outdir out --workdir work --ref_gtf genes.GRCh37.gtf --1 seq_1.fq.gz --2 seq_2.fq.gz --hisat2 /path/to/hisat2 --hisat2_sps /path/to/hisat2_extract_splice_sites.py --samtools /path/to/samtools --threads 10 --sample A 
Example 2 (align):
Run of RNACocktail for alignment of single-end short-read sequences (HISAT2). 
run_rnacocktail.py align --align_idx hisat2-idx --outdir out --workdir work --ref_gtf genes.GRCh37.gtf --U seq.fq.gz --hisat2 /path/to/hisat2 --hisat2_sps /path/to/hisat2_extract_splice_sites.py --samtools /path/to/samtools --threads 10 --sample A 
Example 3 (reconstruct):
Run of RNACocktail for short-read transcriptome reconstruction (StringTie). 
run_rnacocktail.py reconstruct --alignment_bam work/hisat2/A/alignments.sorted.bam --outdir out --workdir work --ref_gtf genes.GRCh37.gtf --stringtie /path/to/stringtie --threads 10 --sample A 
Example 4 (quantify):
Run of RNACocktail for (alignment-free) quantification of paired-end short-read sequences (Salmon-SMEM). 
run_rnacocktail.py quantify --quantifier_idx salmon_fmd_idx --1 seq_1.fq.gz --2 seq_2.fq.gz --libtype IU --salmon_k 19 --outdir out --workdir work --salmon /path/to/salmon --threads 10 --sample A --unzip 
Example 5 (quantify):
Run of RNACocktail for (alignment-free) quantification of single-end short-read sequences (Salmon-SMEM). 
run_rnacocktail.py quantify --quantifier_idx salmon_fmd_idx --U seq.fq.gz --libtype U --salmon_k 19 --outdir out --workdir work --salmon /path/to/salmon --threads 10 --sample A --unzip 
Example 6 (diff):
Run of RNACocktail for differential expression analysis of quantifications computed using Salmon-SMEM (DESeq2). 
run_rnacocktail.py diff --quant_files work/salmon_smem/A1/quant.sf,work/salmon_smem/A2/quant.sf work/salmon_smem/B1/quant.sf,work/salmon_smem/B2/quant.sf --sample A1,A2 B1,B2 --ref_gtf genes.GRCh37.gtf --outdir out --workdir work 
Example 7 (diff):
Run of RNACocktail for differential expression analysis of reads aligned using HISAT2 on reference transcriptome (DESeq2). 
run_rnacocktail.py diff --alignments work/hisat2/A1/alignments.sorted.bam,work/hisat2/A2/alignments.sorted.bam work/hisat2/B1/alignments.sorted.bam,work/hisat2/B2/alignments.sorted.bam --sample A1,A2 B1,B2 --ref_gtf genes.GRCh37.gtf --outdir out --workdir work --featureCounts /path/to/featureCounts 
Example 8 (diff):
Run of RNACocktail for differential expression analysis of reads aligned using HISAT2 on StringTie computed transcriptome (DESeq2). 
run_rnacocktail.py diff --alignments work/hisat2/A1/alignments.sorted.bam,work/hisat2/A2/alignments.sorted.bam work/hisat2/B1/alignments.sorted.bam,work/hisat2/B2/alignments.sorted.bam --transcripts_gtfs work/stringtie/A1/transcripts.gtf,work/stringtie/A2/transcripts.gtf work/stringtie/B1/transcripts.gtf,work/stringtie/B2/transcripts.gtf --sample A1,A2 B1,B2 --ref_gtf genes.GRCh37.gtf --outdir out --workdir work --featureCounts /path/to/featureCounts 
Example 9 (denovo):
Run of RNACocktail for de novo assembly (Oases). 
run_rnacocktail.py denovo --1 seq_1.fq.gz --2 seq_2.fq.gz --outdir out --workdir work --oases /path/to/oases --velveth /path/to/velveth --velvetg /path/to/velvetg --threads 4 --sample A --file_format fastq.gz 
Example 10 (long_correct):
Run of RNACocktail for long-read error correction (LoRDEC). 
run_rnacocktail.py long_correct --kmer 23 --solid 3 --short seq.fq.gz --long seq_long.fa --outdir out --workdir work --lordec /path/to/lordec-correct --threads 4 --sample A 
Example 11 (long_align):
Run of RNACocktail for long-read alignment (STARlong). 
run_rnacocktail.py long_align --long work/lordec/A/long_corrected.fa --outdir out --workdir work --starlong /path/to/STARlong --threads 4 --sample A --sam2psl /path/to/sam2psl.py --samtools /path/to/samtools --genome_dir /path/to/STAR/genome_idx 
Example 12 (long_reconstruct):
Run of RNACocktail for long-read transcriptome reconstruction (IDP). 
run_rnacocktail.py long_reconstruct --alignment work/hisat2/A/alignments.sorted.bam --short_junction work/hisat2/A/splicesites.bed --long_alignment work/starlong/A/Aligned.out.psl --outdir out --workdir work --idp /path/to/runIDP.py --threads 4 --sample A --read_length 100 --ref_genome genome.GRCh37.fa --ref_all_gpd hg19.all.refSeq_gencode_ensemble_EST_known.gpd --ref_gpd genes.GRCh37.refFlat.txt --samtools /path/to/samtools --idp_cfg idp.cfg 
Example 13 (long_fusion):
Run of RNACocktail for long-read fusion detection (IDP-fusion). 
run_rnacocktail.py long_fusion --alignment work/hisat2/A/alignments.sorted.bam --short_junction work/hisat2/A/splicesites.bed --short_fasta seq.fa--long_fasta work/lordec/A/long_corrected.fa --outdir out --workdir work --threads 4 --sample A --ref_genome genome.GRCh37.fa --ref_all_gpd hg19.all.refSeq_gencode_ensemble_EST_known.gpd --ref_gpd genes.GRCh37.refFlat.txt --read_length 100 --genome_bowtie2_idx genome.bt2_idx --transcriptome_bowtie2_idx genes.bt2_idx --uniqueness_bedgraph uniqueness.bedGraph --gmap_idx gmap_idx --idpfusion /path/to/runIDP.py --samtools /path/to/samtools --idpfusion_cfg idpfusion.cfg 
Example 14 (variant):
Run of RNACocktail for RNA-Seq variant calling (GATK). 
run_rnacocktail.py variant --alignment work/hisat2/A/alignments.sorted.bam --outdir out --workdir work --picard /path/to/picard.jar --gatk /path/to/gatk.jar --threads 10 --sample A --ref_genome genome.GRCh37.fa --knownsites dbsnp_138.b37.vcf 
Example 15 (editing):
Run of RNACocktail for RNA editing detection (GIREMI) 
run_rnacocktail.py editing --alignment work/gatk/A/bsqr.bam --variant work/gatk/A/variants_filtered.vcf --strand_pos test/GRCh37_strand_pos.bed --genes_pos test/GRCh37_genes_pos.bed --outdir out --workdir work --giremi_dir /path/to/giremi/directory/ --gatk /path/to/gatk.jar --samtools /path/to/samtools --htslib_dir /path/to/htslib/directory/ --threads 10 --sample A --ref_genome genome.GRCh37.fa --knownsites dbsnp_138.b37.vcf 
Example 16 (fusion):
Run of RNACocktail for RNA fusion detection (FusionCatcher) 
run_rnacocktail.py fusion --data_dir /path/to/fusioncatcher/ensembl/data/directory/ --input seq_1.fq.gz,seq_2.fq.gz --outdir out --workdir work --fusioncatcher /path/to/fusioncatcher --threads 4 --sample A 
Example 17 (all):
Run all pipeline steps (Short-read example) 
run_rnacocktail.py all --outdir out --workdir work --threads 10 --1 A1_1.fq.gz,A2_1.fq.gz B1_1.fq.gz,B2_1.fq.gz --2 A1_2.fq.gz,A2_2.fq.gz B1_2.fq.gz,B2_2.fq.gz --sample all_A1,all_A2 all_B1,all_B2 --ref_gtf genes.GRCh37.gtf --ref_genome genome.GRCh37.fa --align_idx hisat2-idx --quantifier_idx salmon_fmd_idx --unzip --file_format fastq.gz --CleanSam --knownsites dbsnp_138.b37.vcf --strand_pos test/GRCh37_strand_pos.bed --genes_pos test/GRCh37_genes_pos.bed --data_dir /path/to/fusioncatcher/ensembl/data/directory/ --giremi_dir /path/to/giremi/directory/ --gatk /path/to/gatk.jar --htslib_dir /path/to/htslib/directory/ --picard /path/to/picard.jar --samtools /path/to/samtools --hisat2 /path/to/hisat2 --hisat2_sps /path/to/hisat2_extract_splice_sites.py --stringtie /path/to/stringtie --salmon /path/to/salmon --featureCounts /path/to/featureCounts --oases /path/to/oases --velveth /path/to/velveth --velvetg /path/to/velvetg --lordec /path/to/lordec-correct --sam2psl /path/to/sam2psl.py --fusioncatcher /path/to/fusioncatcher 
Example 18 (all):
Run all pipeline steps (long-read example) 
run_rnacocktail.py all --outdir out --workdir work --threads 10 --U seq_short.fa --long seq_long.fa --sample all_C --ref_gtf genes.GRCh37.gtf --ref_genome genome.GRCh37.fa --align_idx hisat2-idx --quantifier_idx salmon_fmd_idx --unzip --file_format fasta --CleanSam --knownsites dbsnp_138.b37.vcf --strand_pos test/GRCh37_strand_pos.bed --genes_pos test/GRCh37_genes_pos.bed --data_dir /path/to/fusioncatcher/ensembl/data/directory/ --giremi_dir /path/to/giremi/directory/ --gatk /path/to/gatk.jar --htslib_dir /path/to/htslib/directory/ --star_genome_dir /path/to/STAR/genome_idx/ --genome_bowtie2_idx genome.bt2_idx --transcriptome_bowtie2_idx genes.bt2_idx --uniqueness_bedgraph uniqueness.bedGraph --gmap_idx gmap_idx --ref_all_gpd hg19.all.refSeq_gencode_ensemble_EST_known.gpd --ref_gpd genes.GRCh37.refFlat.txt --read_length 100 --picard /path/to/picard.jar --hisat2_opts \"-f\" --idp /path/to/idp/runIDP.py --idpfusion /path/to/idpfusion/runIDP.py --samtools /path/to/samtools --hisat2 /path/to/hisat2 --hisat2_sps /path/to/hisat2_extract_splice_sites.py --stringtie /path/to/stringtie --salmon /path/to/salmon --featureCounts /path/to/featureCounts --oases /path/to/oases --velveth /path/to/velveth --velvetg /path/to/velvetg --lordec /path/to/lordec-correct --sam2psl /path/to/sam2psl.py --fusioncatcher /path/to/fusioncatcher 
Command line options
General options
Option
Definition
--sample STRING
Sample name
--threads INT
Number of threads to use (default: 1)
--start INT
It re-starts executing the workflow/pipeline from the given step number. This can be used when the pipeline has crashed/stopped and one wants to re-run it from from the step where it stopped without re-running from the beginning the entire pipeline. 0 is for restarting automatically and 1 is the first step. (default is '0').
--timeout INT
Maximum run time for commands (in seconds) (default 10000000)
Short-read alignment options
run_rnacocktail.py align 
Option
Definition
--sr_aligner STRING
Short-read alignment tool (default: HISAT2)
--align_idx STRING
The basename of the index generated by the alignment tool for the reference genome
--1 STRING
Comma-separated list of files containing mate 1s (filename usually includes _1), e.g. --1 A_1.fq,B_1.fq.
--2 STRING
Comma-separated list of files containing mate 2s (filename usually includes _2), e.g. --2 A_2.fq,B_2.fq.
--U STRING
Comma-separated list of files containing unpaired reads to be aligned, e.g. --U A.fq,B.fq.
--sra STRING
Comma-separated list of SRA accession numbers, e.g. --sra SRR353653,SRR353654. Information about read types is available at here, where sra is SRA accession number.
--ref_gtf STRING
The reference transcriptome annotation file (in GTF or GFF3 format) to guide the analysis. ( --known-splicesite-infile option for HISAT will be created based on this file)
--hisat2 STRING
Path to HISAT2 executable (Optional. Can be on PATH or defined on defaults.py)
--hisat2_sps STRING
Path to hisat2_extract_splice_sites.py script. Can be found in HISAT2 package. (Optional. Can be on PATH or defined on defaults.py)
--samtools STRING
Path to samtools executable. (Optional. Can be on PATH or defined on defaults.py)
--hisat2_opts STRING
Other options used for HISAT2 aligner. (should be put between " ") (For HISAT2 check here).
Short-read transcriptome reconstruction options
run_rnacocktail.py reconstruct 
Option
Definition
--reconstructor STRING
The transcriptome reconstruction tool to use (default: StringTie)
--alignment_bam STRING
A BAM file with RNA-Seq read mappings which must be sorted by their genomic location (e.g. The output BAM file generated in align mode).
--ref_gtf STRING
The reference transcriptome annotation file (in GTF or GFF3 format) to guide the analysis.
--stringtie STRING
Path to StringTie executable (Optional. Can be on PATH or defined on defaults.py)
--samtools STRING
Path to samtools executable. (Optional. Can be on PATH or defined on defaults.py)
--stringtie_opts STRING
Other options used for StringTie transcriptome reconstruction. (should be put between " ") (For StringTie check here).
Alignment-free transcript quantification options
run_rnacocktail.py quantify 
Option
Definition
--quantifier STRING
The quantification tool to use (default: Salmon-SMEM)
--quantifier_idx STRING
The index generated for the reference transcriptome. (FMD-based index for Salmon-SMEM)
--1 STRING
Comma-separated list of files containing mate 1s (filename usually includes _1), e.g. --1 A_1.fq,B_1.fq.
--2 STRING
Comma-separated list of files containing mate 2s (filename usually includes _2), e.g. --2 A_2.fq,B_2.fq.
--U STRING
Comma-separated list of files containing unpaired reads to be aligned, e.g. --U A.fq,B.fq.
--salmon_k INT
SMEM's smaller than this size will not be considered by Salmon. (default 19).
--libtype STRING
Format string describing the library type. (For Salmon check here).
--unzip
The sequence files are zipped. So unzip them first
--salmon STRING
Path to Salmon executable (Optional. Can be on PATH or defined on defaults.py)
--salmon_smem_opts STRING
Other options used for Salmon-SMEM quantifications. (should be put between " ") (For Salmon check here).
Differential Analysis options
run_rnacocktail.py diff 
Option
Definition
--difftool STRING
The differential analysis tool to use. (default: DESeq2)
--quant_files STRING
Quantification files for each sample (e.g. Salmon's quant.sf outputs). Replicates in same sample should be listed comma separated. e.g --quant_files A1/quant.sf,A2/quant.sf B1/quant.sf,B2/quant.sf
--transcripts_gtfs STRING
Reconstructed transcript GTF files (for instance StringTie's transcripts.gtf output). Replicates in same sample should be listed comma separated. e.g --transcripts_gtfs A1/transcripts.gtf,A2/transcripts.gtf B1/transcripts.gtf,B2/transcripts.gtf
--alignments STRING
Alignment BAM files for each sample (for instance HISAT2's output). Replicates in same sample should be listed comma separated. e.g --alignments A1/alignments.bam,A2/alignments.bam B1/alignments.bam,B2/alignments.bam
--ref_gtf STRING
The reference transcriptome annotation file (in GTF or GFF3 format) to guide the analysis.
--sample STRING
Sample names. Number of samples and replicates should match the input quantification (--quant_files) or alignemnt (--alignments). Replicates in same sample should be listed comma separated. e.g --sample A1,A2 B1,B2
--mincount INT
Minimum read counts per transcripts. Differential analysis pre-filtering step removes transcripts that have less than this number of reads. (default 2)
--alpha FLOAT
Adjusted p-value significance level for differential analysis. (default 0.05)
--R STRING
Path to R executable (DESeq2, readr, tximport should have been installed in R) (Optional. Can be on PATH or defined on defaults.py)
--featureCounts STRING
Path to featureCounts executable. (Optional. Can be on PATH or defined on defaults.py)
--stringtie STRING
Path to StringTie executable. (Optional. Can be on PATH or defined on defaults.py)
--stringtie_merge_opts STRING
Other options used for StringTie merge. Can be set when the reconstructed transcript GTFs are used.(should be put between " ") (For StringTie check here).
--featureCounts_opts STRING
Other options used for featureCounts. (should be put between " ") (For options check here).
De novo assembly options
run_rnacocktail.py denovo 
Option
Definition
--assembler STRING
The de novo assembler to use. (default Oases)
--assmebly_hash INT
Odd integer, or a comma separated list of odd integers that specify the assembly has length (for Oases/Velvet).
--file_format STRING
Input file format for de novo assembly Options: fasta, fastq, raw, fasta.gz, fastq.gz, raw.gz, sam, bam, fmtAuto. (default fasta)
--read_type STRING
Input sequence read type for de novo assembly Options: short, shortPaired, short2, shortPaired2, long, longPaired, reference. (Check here for description) (default short)
--1 STRING
Comma-separated list of files containing mate 1s (filename usually includes _1), e.g. --1 A_1.fq,B_1.fq.
--2 STRING
Comma-separated list of files containing mate 2s (filename usually includes _2), e.g. --2 A_2.fq,B_2.fq.
--U STRING
Comma-separated list of files containing unpaired reads to be aligned, e.g. --U A.fq,B.fq.
--I STRING
Comma-separated list of files containing interleaved paired-end reads to be assembled, e.g. --I A.fq,B.fq.
--oases STRING
Path to oases executable. (Optional. Can be on PATH or defined on defaults.py)
--velvetg STRING
Path to velvetg executable. (Optional. Can be on PATH or defined on defaults.py)
--velveth STRING
Path to velveth executable. (Optional. Can be on PATH or defined on defaults.py)
--velveth_opts STRING
Other options used for assembly by velveth. (For velvet options check herehttps://github.com/dzerbino/velvet/blob/master/Manual.pdf).
--velvetg_opts STRING
Other options used for assembly by velvetg. (should be put between " ") (For velvet options check herehttps://github.com/dzerbino/velvet/blob/master/Manual.pdf).
--velveth_opts STRING
Other options used for assembly by velveth. (For velvet options check herehttps://github.com/dzerbino/velvet/blob/master/Manual.pdf).
--oases_opts STRING
Other options used for assembly by Oases. (should be put between " ") (For Oases options check herehttps://github.com/dzerbino/oases).
Long read error correction options
run_rnacocktail.py long_correct 
Option
Definition
--long_corrector STRING
The long-read error correction tool to use. (default LoRDEC).
--kmer INT
LoRDEC k-mer length
--solid INT
LoRDEC solidity abundance threshold for k-mers
--long STRING
The FASTA file containing long reads
--short STRING
The FASTA or FASTQ file containing short reads. (can be compressed .gz file)
--lordec STRING
Path to LoRDEC executable (Optional. Can be on PATH or defined on defaults.py)
--lordec_opts STRING
Other options used for LoRDEC. (should be put between " ") (For LoRDEC check here).
Long read alignment options
run_rnacocktail.py long_align 
Option
Definition
--long_aligner STRING
The long-read alignment tool to use. (default STARlong).
--long STRING
The FASTA file containing long reads
--genome_dir STRING
Specifies path to the genome directory where STAR genome indices where generated
--ref_gtf STRING
The reference transcriptome annotation file (in GTF or GFF3 format) to guide the analysis.
--starlong STRING
Path to STARlong executable (version 2.5.0a or later) (Optional. Can be on PATH or defined on defaults.py)
--sam2psl STRING
Path to the sam2psl.py script. Can be found in FusionCatcher package. (Optional. Can be on PATH or defined on defaults.py)
--samtools STRING
Path to samtools executable. (Optional. Can be on PATH or defined on defaults.py)
--starlong_opts STRING
Other options used for LoRDEC. (should be put between " ") (For LoRDEC check here). As the default we use the following options as advised in here: --outSAMattributes NH HI NM MD --readNameSeparator space --outFilterMultimapScoreRange 1 --outFilterMismatchNmax 2000 --scoreGapNoncan -20 --scoreGapGCAG -4 --scoreGapATAC -8 --scoreDelOpen -1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 --alignEndsType Local --seedSearchStartLmax 50 --seedPerReadNmax 100000 --seedPerWindowNmax 1000 --alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000 .
Long read transcriptome reconstruction options
run_rnacocktail.py long_reconstruct 
Option
Definition
--long_reconstructor STRING
The long-read transcriptome reconstruction tool to use. (default IDP).
--alignment STRING
A BAM/SAM file with short RNA-Seq read mappings (e.g. The output BAM file generated in align mode). If BAM file is given, it will be converted to SAM.
--short_junction STRING
A BED file with short RNA-Seq read junctions (e.g. The bed junction file generated in align mode.) For bed file format check here.
--long_alignment STRING
A PSL file with long RNA-Seq read mappings (e.g. The output PSL file generated in long_align mode)
--mode_number INT
You can run IDP in two steps. If for a reason IDP finished isoform candidate construction step but was terminated in candidate selection step, you can restart the candidate selection step without re-running the isoform candidate construction. mode 0 (default): end-to-end IDP run. mode 1: generates isoform candidate pool (file: isoform_construction.NisoXX.gpd). mode 2: runs candidate selection step. Note: make sure isoform candidate pool (file: isoform_construction.NisoXX.gpd) file is already generated in temp folder.
--ref_genome STRING
The reference genome FASTA file
--ref_all_gpd STRING
GPD format annotation file for the whole genome splicing data from multiple sources including ESTs and reference genome databases. For hg19 you may use the full genome example in here.
--ref_gpd STRING
The reference transcriptome annotation file (in GPD format) to guide the analysis.
--ref_genome STRING
The reference genome FASTA file
--read_length INT
The short-read length. (default: 100).
--samtools STRING
Path to samtools executable. (Optional. Can be on PATH or defined on defaults.py)
--idp STRING
Path to runIDP.py script. (Optional. Can be on PATH or defined on defaults.py)
--idp_cfg STRING
the .cfg file that include other options used for IDP long read transcriptome reconstruction. (For IDP check here) These options will be used to generate the .cfg file.
Long read fusion prediction options
run_rnacocktail.py long_reconstruct 
Option
Definition
--long_fusion STRING
The long-read fusion detection tool to use. (default IDP-fusion).
--alignment STRING
A BAM/SAM file with short RNA-Seq read mappings (e.g. The output BAM file generated in align mode). If BAM file is given, it will be converted to SAM.
--short_junction STRING
A BED file with short RNA-Seq read junctions (e.g. The bed junction file generated in align mode.) For bed file format check here.
--long_alignment STRING
A PSL file with long RNA-Seq read mappings (The output PSL by GMAP). This is optional, if you don't provide it the code will automatically call gmap to align.
--long_fasta STRING
The FASTA file containing long reads.
--short_fasta STRING
The FASTA file containing short reads.
--mode_number INT
You can run IDP in two steps. If for a reason IDP finished isoform candidate construction step but was terminated in candidate selection step, you can restart the candidate selection step without re-running the isoform candidate construction. mode 0 (default): end-to-end IDP run. mode 1: generates isoform candidate pool (file: isoform_construction.NisoXX.gpd). mode 2: runs candidate selection step. Note: make sure isoform candidate pool (file: isoform_construction.NisoXX.gpd) file is already generated in temp folder.
--ref_genome STRING
The reference genome FASTA file
--ref_all_gpd STRING
GPD format annotation file for the whole genome splicing data from multiple sources including ESTs and reference genome databases. For hg19 you may use the full genome example in here.
--ref_gpd STRING
The reference transcriptome annotation file (in GPD format) to guide the analysis.
--ref_genome STRING
The reference genome FASTA file
--uniqueness_bedgraph STRING
File with the uniqueness scores in bedgraph format. Used to annotate the uniqueness of regions flanking fusions sites. Duke Uniqueness track from UCSC genome browser in bedGraph format can be used.
--genome_bowtie2_idx STRING
The reference genome bowtie2 index file.
--transcriptome_bowtie2_idx STRING
The reference transcriptome bowtie2 index file.
--read_length INT
The short-read length. (default: 100).
--star_dir STRING
Path to the directory with STAR executable.
--bowtie2_dir STRING
Path to the directory with bowtie2 executable.
--gmap_idx STRING
Path to the directory with GMAP index for the reference genome.
--samtools STRING
Path to samtools executable. (Optional. Can be on PATH or defined on defaults.py)
--idpfusion STRING
Path to runIDP.py script. (Optional. Can be on PATH or defined on defaults.py)
--idpfusion_cfg STRING
the .cfg file that include other options used for IDP-fusion long read fusion detection. (For IDP-fusion check here) These options will be used to generate the .cfg file.
--gmap STRING
Path to GMAP executable.
RNA-Seq variant calling options
run_rnacocktail.py variant 
Option
Definition
--variant_caller STRING
The variant caller to use. For GATK's general approach used for calling variants in RNAseq check here (default GATK).
--alignment STRING
A BAM/SAM file with RNA-Seq read mappings. (e.g. The output BAM file generated in align mode).
--CleanSam
Use Picard's CleanSam command to clean the input alignment.
--no_BaseRecalibrator
Don't run BaseRecalibrator step.
--ref_genome
The reference genome FASTA file
--knownsites
A database of known polymorphic sites (e.g. dbSNP). Used in GATK BaseRecalibrator. NOTE: to run BaseRecalibrator step knownsites should be provided.
--picard
Path to picard executable. (Optional. Can be on PATH or defined on defaults.py)
--gatk
Path to GATK executable. (Optional. Can be on PATH or defined on defaults.py)
--java
Path to JAVA executable. (Optional. Can be on PATH or defined on defaults.py)
--java_opts
Java options used for picard and GATK commands. (should be put between " ") (default: -Xms1g -Xmx5g) 
--AddOrReplaceReadGroups_opts
Other options used for picard AddOrReplaceReadGroups command. (should be put between " ") (For Picard check here) (default: "SO=coordinate RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample") 
--MarkDuplicates_opts
Other options used for picard MarkDuplicates command. (should be put between " ") (For Picard check here) (default: "CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT") 
--SplitNCigarReads_opts
Other options used for GATK SplitNCigarReads command. (should be put between " ") (For GATK SplitNCigarReads check here) (default: "-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS") 
--BaseRecalibrator_opts
Other options used for GATK BaseRecalibrator command. (should be put between " ") (For GATK BaseRecalibrator check here).
--ApplyBQSR_opts
Other options used for GATK ApplyBQSR command. (should be put between " ") (For GATK ApplyBQSR check here).
--HaplotypeCaller_opts
Other options used for GATK HaplotypeCaller command. (should be put between " ") (For GATK HaplotypeCaller check here). (default: "-stand-call-conf --dont-use-soft-clipped-bases" 
--VariantFiltration_opts
Other options used for GATK VariantFiltration command. (should be put between " ") (For GATK VariantFiltration check here). (default: "-window 35 -cluster 3 --filter-name FS -filter 'FS > 30.0' --filter-name QD -filter 'QD < 2.0'" 
RNA editing detection options
run_rnacocktail.py editing 
Option
Definition
--editing_caller STRING
The RNA Editing caller to use. (default GIREMI).
--alignment STRING
A BAM/SAM file with RNA-Seq read mappings. (e.g. The output BAM file generated in align mode).
--variant STRING
A VCF file with variants. (e.g. The output VCF file generated in variant calling mode).
--strand_pos STRING
A BED file which specifies the strand of the genes/transcripts. Each row should have 5 columns: chromosome,start,end,name,score(can be .),strand (+ or -). Examples for Human on GRCh37 can be found in test directory. You can generate this file using reference transcript annotations.
--genes_pos STRING
A BED file which specifies the positions in the genome that genes reside. Each row should have 3 columns: chromosome,start,end,name. Examples for Human on GRCh37 can be found in test directory. You can generate this file using reference transcript annotations.
--ref_genome STRING
The reference genome FASTA file
--knownsites
A database of known polymorphic sites (e.g. dbSNP) in VCF format.
--giremi_dir
Path to giremi directory that include giremi executable and giremi.r R script. (required)
--htslib_dir
Path to HTSlib library directory. (Optional. Can be on LD_LIBRARY_PATH or defined on defaults.py)
--samtools STRING
Path to samtools executable. (Optional. Can be on PATH or defined on defaults.py)
--gatk
Path to GATK executable. (Optional. Can be on PATH or defined on defaults.py)
--java
Path to JAVA executable. (Optional. Can be on PATH or defined on defaults.py)
--java_opts
Java options used for GATK commands. (should be put between " ") (default: -Xms1g -Xmx5g) 
--giremi_opts
Other options used for GIREMI. (should be put between " ").(For GIREMI check here) 
--VariantAnnotator_opts
Other options used for GATK VariantAnnotator command. (should be put between " ") (For GATK VariantAnnotator check here). 
RNA fusion detection options
run_rnacocktail.py fusion 
Option
Definition
--fusion_caller STRING
The RNA fusion caller to use.(default FusionCatcher).
--data_dir STRING
The data directory where all the annotations files from Ensembl database are placed. This directory should be built using 'fusioncatcher-build'.
--input STRING
The input file(s) or directory. The files should be in FASTQ or SRA format and may be or not compressed using gzip or zip. A list of files can be specified by given the filenames separated by comma. If a directory is given then it will analyze all the files found with the following extensions: .sra, .fastq, .fastq.zip, .fastq.gz, .fastq.bz2, fastq.xz, .fq, .fq.zip, .fq.gz, .fq.bz2, fz.xz, .txt, .txt.zip, .txt.gz, .txt.bz2 .
--fusioncatcher
Path to FusionCatcher executable. (Optional. Can be on PATH or defined on defaults.py)
--fusioncatcher_opts
Other options used for FusionCatcher. (should be put between " ") (For FusionCatcher check here).
Run all pipeline steps options
run_rnacocktail.py all 
Option
Definition
--sr_aligner STRING
Short-read alignment tool (default: HISAT2)
--reconstructor STRING
The transcriptome reconstruction tool to use (default: StringTie)
--quantifier STRING
The quantification tool to use (default: Salmon-SMEM)
--difftool STRING
The differential analysis tool to use. (default: DESeq2)
--assembler STRING
The de novo assembler to use. (default Oases)
--long_corrector STRING
The long-read error correction tool to use. (default LoRDEC).
--long_reconstructor STRING
The long-read transcriptome reconstruction tool to use. (default IDP).
--editing_caller STRING
The RNA Editing caller to use. (default GIREMI).
--variant_caller STRING
The variant caller to use. For GATK's general approach used for calling variants in RNAseq check here (default GATK).
--fusion_caller STRING
The RNA fusion caller to use.(default FusionCatcher).
--long_fusion STRING
The long-read fusion detection tool to use. (default IDP-fusion).
--long_aligner STRING
The long-read alignment tool to use. (default STARlong).
--1 STRING
List of files containing mate 1s (filename usually includes _1). Replicates in same sample should be listed comma-separated : e.g. --1 A1_1.fq,A2_1.fq B1_1.fq,B2_1.fq.
--2 STRING
List of files containing mate 2s (filename usually includes _2). Replicates in same sample should be listed comma-separated : e.g. --1 A1_2.fq,A2_2.fq B1_2.fq,B2_2.fq.
--U STRING
List of files containing unpaired reads to be aligned. Replicates in same sample should be listed comma-separated : e.g. --U A1.fq,A2.fq B1.fq,B2.fq.
--long STRING
List of FASTA files containing long reads. Replicates in same sample should be listed comma-separated : e.g. --long A1.fasta,A2.fasta B1.fasta,B2.fasta.
--exclude STRING
List of exlcluded steps.
--sample STRING
Sample names. Number of samples and replicates should match the input sequences. Replicates in same sample should be listed comma separated. e.g --sample A1,A2 B1,B2
--align_idx STRING
The basename of the index generated by the alignment tool for the reference genome
--quantifier_idx STRING
The index generated for the reference transcriptome. (FMD-based index for Salmon-SMEM)
--star_genome_dir STRING
Specifies path to the genome directory where STAR genome indices where generated
--genome_bowtie2_idx STRING
The reference genome bowtie2 index file.
--transcriptome_bowtie2_idx STRING
The reference transcriptome bowtie2 index file.
--gmap_idx STRING
Path to the directory with GMAP index for the reference genome.
--read_length INT
The short-read length. (default: 100).
--salmon_k INT
SMEM's smaller than this size will not be considered by Salmon. (default: 19).
--libtype STRING
Format string describing the library type. (For Salmon check here). (default: IU)
--unzip
The sequence files are zipped. So unzip them first
--mincount INT
Minimum read counts per transcripts. Differential analysis pre-filtering step removes transcripts that have less than this number of reads. (default 2)
--alpha FLOAT
Adjusted p-value significance level for differential analysis. (default 0.05)
--assmebly_hash INT
Odd integer, or a comma separated list of odd integers that specify the assembly has length (for Oases/Velvet).
--file_format STRING
Input file format for de novo assembly Options: fasta, fastq, raw, fasta.gz, fastq.gz, raw.gz, sam, bam, fmtAuto. (default fasta)
--read_type STRING
Input sequence read type for de novo assembly Options: short, shortPaired, short2, shortPaired2, long, longPaired, reference. (Check here for description) (default short)
--kmer INT
LoRDEC k-mer length
--solid INT
LoRDEC solidity abundance threshold for k-mers
--mode_number INT
You can run IDP in two steps. If for a reason IDP finished isoform candidate construction step but was terminated in candidate selection step, you can restart the candidate selection step without re-running the isoform candidate construction. mode 0 (default): end-to-end IDP run. mode 1: generates isoform candidate pool (file: isoform_construction.NisoXX.gpd). mode 2: runs candidate selection step. Note: make sure isoform candidate pool (file: isoform_construction.NisoXX.gpd) file is already generated in temp folder.
--CleanSam
Use Picard's CleanSam command to clean the input alignment.
--no_BaseRecalibrator
Don't run BaseRecalibrator step.
--data_dir STRING
The data directory where all the annotations files from Ensembl database are placed. This directory should be built using 'fusioncatcher-build'.
--input STRING
The input file(s) or directory. The files should be in FASTQ or SRA format and may be or not compressed using gzip or zip. A list of files can be specified by given the filenames separated by comma. If a directory is given then it will analyze all the files found with the following extensions: .sra, .fastq, .fastq.zip, .fastq.gz, .fastq.bz2, fastq.xz, .fq, .fq.zip, .fq.gz, .fq.bz2, fz.xz, .txt, .txt.zip, .txt.gz, .txt.bz2 .
--ref_gtf STRING
The reference transcriptome annotation file (in GTF or GFF3 format) to guide the analysis. ( --known-splicesite-infile option for HISAT will be created based on this file)
--ref_genome STRING
The reference genome FASTA file
--ref_all_gpd STRING
GPD format annotation file for the whole genome splicing data from multiple sources including ESTs and reference genome databases. For hg19 you may use the full genome example in here.
--ref_gpd STRING
The reference transcriptome annotation file (in GPD format) to guide the analysis.
--uniqueness_bedgraph STRING
File with the uniqueness scores in bedgraph format. Used to annotate the uniqueness of regions flanking fusions sites. Duke Uniqueness track from UCSC genome browser in bedGraph format can be used.
--knownsites
A database of known polymorphic sites (e.g. dbSNP). Used in GATK BaseRecalibrator. NOTE: to run BaseRecalibrator step knownsites should be provided.
--strand_pos STRING
A BED file which specifies the strand of the genes/transcripts. Each row should have 5 columns: chromosome,start,end,name,score(can be .),strand (+ or -). Examples for Human on GRCh37 can be found in test directory. You can generate this file using reference transcript annotations.
--genes_pos STRING
A BED file which specifies the positions in the genome that genes reside. Each row should have 3 columns: chromosome,start,end,name. Examples for Human on GRCh37 can be found in test directory. You can generate this file using reference transcript annotations.
--hisat2 STRING
Path to HISAT2 executable (Optional. Can be on PATH or defined on defaults.py)
--hisat2_sps STRING
Path to hisat2_extract_splice_sites.py script. Can be found in HISAT2 package. (Optional. Can be on PATH or defined on defaults.py)
--samtools STRING
Path to samtools executable. (Optional. Can be on PATH or defined on defaults.py)
--hisat2_opts STRING
Other options used for HISAT2 aligner. (should be put between " ") (For HISAT2 check here).
--stringtie STRING
Path to StringTie executable (Optional. Can be on PATH or defined on defaults.py)
--samtools STRING
Path to samtools executable. (Optional. Can be on PATH or defined on defaults.py)
--stringtie_opts STRING
Other options used for StringTie transcriptome reconstruction. (should be put between " ") (For StringTie check here).
--salmon STRING
Path to Salmon executable (Optional. Can be on PATH or defined on defaults.py)
--salmon_smem_opts STRING
Other options used for Salmon-SMEM quantifications. (should be put between " ") (For Salmon check here).
--R STRING
Path to R executable (DESeq2, readr, tximport should have been installed in R) (Optional. Can be on PATH or defined on defaults.py)
--featureCounts STRING
Path to featureCounts executable. (Optional. Can be on PATH or defined on defaults.py)
--stringtie_merge_opts STRING
Other options used for StringTie merge. Can be set when the reconstructed transcript GTFs are used.(should be put between " ") (For StringTie check here).
--featureCounts_opts STRING
Other options used for featureCounts. (should be put between " ") (For options check here).
--oases STRING
Path to oases executable. (Optional. Can be on PATH or defined on defaults.py)
--velvetg STRING
Path to velvetg executable. (Optional. Can be on PATH or defined on defaults.py)
--velveth STRING
Path to velveth executable. (Optional. Can be on PATH or defined on defaults.py)
--velveth_opts STRING
Other options used for assembly by velveth. (For velvet options check herehttps://github.com/dzerbino/velvet/blob/master/Manual.pdf).
--velvetg_opts STRING
Other options used for assembly by velvetg. (should be put between " ") (For velvet options check herehttps://github.com/dzerbino/velvet/blob/master/Manual.pdf).
--velveth_opts STRING
Other options used for assembly by velveth. (For velvet options check herehttps://github.com/dzerbino/velvet/blob/master/Manual.pdf).
--oases_opts STRING
Other options used for assembly by Oases. (should be put between " ") (For Oases options check herehttps://github.com/dzerbino/oases).
--lordec STRING
Path to LoRDEC executable (Optional. Can be on PATH or defined on defaults.py)
--lordec_opts STRING
Other options used for LoRDEC. (should be put between " ") (For LoRDEC check here).
--starlong STRING
Path to STARlong executable (version 2.5.0a or later) (Optional. Can be on PATH or defined on defaults.py)
--sam2psl STRING
Path to the sam2psl.py script. Can be found in FusionCatcher package. (Optional. Can be on PATH or defined on defaults.py)
--starlong_opts STRING
Other options used for LoRDEC. (should be put between " ") (For LoRDEC check here). As the default we use the following options as advised in here: --outSAMattributes NH HI NM MD --readNameSeparator space --outFilterMultimapScoreRange 1 --outFilterMismatchNmax 2000 --scoreGapNoncan -20 --scoreGapGCAG -4 --scoreGapATAC -8 --scoreDelOpen -1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 --alignEndsType Local --seedSearchStartLmax 50 --seedPerReadNmax 100000 --seedPerWindowNmax 1000 --alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000 .
--idp STRING
Path to runIDP.py script. (Optional. Can be on PATH or defined on defaults.py)
--idp_cfg STRING
the .cfg file that include other options used for IDP long read transcriptome reconstruction. (For IDP check here) These options will be used to generate the .cfg file.
--star_dir STRING
Path to the directory with STAR executable.
--idpfusion STRING
Path to runIDP.py script. (Optional. Can be on PATH or defined on defaults.py)
--idpfusion_cfg STRING
the .cfg file that include other options used for IDP-fusion long read fusion detection. (For IDP-fusion check here) These options will be used to generate the .cfg file.
--bowtie2_dir STRING
Path to the directory with bowtie2 executable.
--gmap STRING
Path to GMAP executable.
--picard
Path to picard executable. (Optional. Can be on PATH or defined on defaults.py)
--gatk
Path to GATK executable. (Optional. Can be on PATH or defined on defaults.py)
--java
Path to JAVA executable. (Optional. Can be on PATH or defined on defaults.py)
--java_opts
Java options used for picard and GATK commands. (should be put between " ") (default: -Xms1g -Xmx5g) 
--AddOrReplaceReadGroups_opts
Other options used for picard AddOrReplaceReadGroups command. (should be put between " ") (For Picard check here) (default: "SO=coordinate RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample") 
--MarkDuplicates_opts
Other options used for picard MarkDuplicates command. (should be put between " ") (For Picard check here) (default: "CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT") 
--SplitNCigarReads_opts
Other options used for GATK SplitNCigarReads command. (should be put between " ") (For GATK SplitNCigarReads check here) (default: "-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS") 
--BaseRecalibrator_opts
Other options used for GATK BaseRecalibrator command. (should be put between " ") (For GATK BaseRecalibrator check here).
--ApplyBQSR_opts
Other options used for GATK ApplyBQSR command. (should be put between " ") (For GATK ApplyBQSR check here).
--HaplotypeCaller_opts
Other options used for GATK HaplotypeCaller command. (should be put between " ") (For GATK HaplotypeCaller check here). (default: "-stand-call-conf --dont-use-soft-clipped-bases" 
--VariantFiltration_opts
Other options used for GATK VariantFiltration command. (should be put between " ") (For GATK VariantFiltration check here). (default: "-window 35 -cluster 3 --filter-name FS -filter 'FS > 30.0' --filter-name QD -filter 'QD < 2.0'" 
--giremi_dir
Path to giremi directory that include giremi executable and giremi.r R script. (required)
--htslib_dir
Path to HTSlib library directory. (Optional. Can be on LD_LIBRARY_PATH or defined on defaults.py)
--gatk
Path to GATK executable. (Optional. Can be on PATH or defined on defaults.py)
--giremi_opts
Other options used for GIREMI. (should be put between " ").(For GIREMI check here) 
--VariantAnnotator_opts
Other options used for GATK VariantAnnotator command. (should be put between " ") (For GATK VariantAnnotator check here). 
--fusioncatcher
Path to FusionCatcher executable. (Optional. Can be on PATH or defined on defaults.py)
--fusioncatcher_opts
Other options used for FusionCatcher. (should be put between " ") (For FusionCatcher check here).
Preparing genome index files
RNACocktail requires the user to separately build the indexes for the genomic and/or transcriptomic references. We show below how this can be done based on the tools to be invoked in RNACocktail.
HISAT2
hisat2-build [options] reference.fa hisat2_index_basename 
Check here for more information.
Salmon-SMEM
salmon index -t reference.fa -i salmon_index_basename --type fmd 
Check here for more information.
STAR
STAR --runMode genomeGenerate --genomeDir STAR_index_basename --genomeFastaFiles reference.fa --runThreadN 4 
Check here for more information.
GMAP
gmap_build -d gmap_index_basename reference.fa 
Check here for more information.
bowtie2
bowtie2-build reference.fa bowtie2_index_basename 
Check here for more information.
