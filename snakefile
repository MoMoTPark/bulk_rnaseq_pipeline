## Author: Moe
## Bulk RNA-seq analysis pipeline

# Note that relative paths are used and snakefile should be located in a directory (e.g., code) within root analysis directory
# Note that threads are adjusted for a 32 core EC2 machine; modify based on your environment.

# Sample wildcards (modify to reflect your data)
IDS, = glob_wildcards("../data/{id}_R1.fastq.gz")

rule all:
    input:
        expand("../salmon_quant_raw/{id}/quant.sf", id=IDS),
        expand("../qorts/{id}/{id}.QC.summary.txt", id=IDS),
        expand("../salmon_quant/{id}/quant.sf", id=IDS),

ruleorder: samtools_import > star_align_1 > splice_junction_filter > star_align_2

# Quantify with Salmon Raw Fastq (for fast turn-around jobs)
rule salmon_quant_raw:
    input:
        read_1 = "../data/{id}_R1.fastq.gz", # Modify to refect your data
        read_2 = "../data/{id}_R2.fastq.gz", # Modify to refect your data
        index = '../ref/salmon_index/'
    output:
        quant = '../salmon_quant_raw/{id}/quant.sf'
    params:
        outfolder = '../salmon_quant_raw/{id}'
    log: "../logs/{id}_salmon_quant_raw.log"
    benchmark: "../benchmarks/{id}_salmon_quant_raw.benchmark"
    priority: 100
    resources:
        load = 33
    shell:
        '''salmon quant -i {input.index} -l A -1 {input.read_1} -2 {input.read_2} --validateMappings --seqBias --gcBias --posBias --threads 8  -o {params.outfolder} 2> {log}'''


# Generate QC html from raw reads for agent trim
# -A: Disable adaptor trimming
# -G: Disable trim poly-G
# -Q: Disable quality read filtering
# -L: Disable length filtering
# -w: number of threads
# -h: html output
rule fastp_qc:
    input:
        read_1 = "../data/{id}_R1.fastq.gz", # Modify to refect your data
        read_2 = "../data/{id}_R2.fastq.gz"  # Modify to refect your data
    output:
        html = "../fastp_qc/{id}_fastp_qc.html",
        json = "../fastp_qc/{id}_fastp_qc.json"
    log: "../logs/{id}_fastp_qc.log"
    benchmark: "../benchmarks/{id}_fastp_qc.benchmark"
    shell:
        '''fastp -i {input.read_1} -I {input.read_2} -A -G -Q -L -w 4 -h {output.html} -j {output.json} 2> {log}'''

# -v2 specifies HS2 library prep protocol
rule agent_trim:
    input:
        read_1 = "../data/{id}_R1.fastq.gz", # Modify to refect your data
        read_2 = "../data/{id}_R2.fastq.gz", # Modify to refect your data
        html = "../fastp_qc/{id}_fastp_qc.html"
    output:
        agent_out = "../agent_trim/{id}_MBC.txt.gz",
        read_1 = "../agent_trim/{id}_R1.fastq.gz",
        read_2 = "../agent_trim/{id}_R2.fastq.gz"
    params:
        out_prefix = "../agent_trim/{id}"
    log: "../logs/{id}_agent_trim.log"
    benchmark: "../benchmarks/{id}_agent_trim.benchmark"
    shell:
        '''agent trim -v2 -fq1 {input.read_1} -fq2 {input.read_2} -polyG 3 -qualityTrimming 15 -out {params.out_prefix} 2> {log}'''

# Trim reads with fastp with the following parameters
# -q, the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
# -u, how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
# -y, enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).
# -Y, the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. (int [=30])
rule fastp_trim:
    input:
        read_1 = "../agent_trim/{id}_R1.fastq.gz",
        read_2 = "../agent_trim/{id}_R2.fastq.gz"
    output:
        html = "../fastp/{id}_fastp.html",
        json = "../fastp/{id}_fastp.json",
        trimmed_read_1 = "../trimmed_reads/{id}_R1_qctrim.fastq.gz",
        trimmed_read_2 = "../trimmed_reads/{id}_R2_qctrim.fastq.gz"
    log: "../logs/{id}_fastp_trim.log"
    benchmark: "../benchmarks/{id}_fastp_trim.benchmark"
	shell:
		'''fastp -i {input.read_1} -I {input.read_2} -q 15 -u 40 -y -Y 30 -h {output.html} -j {output.json} -o {output.trimmed_read_1} -O {output.trimmed_read_2} 2> {log}'''

# Preserve fastq headers with barcode information included and output as bam with BC tags included 
rule samtools_import:
    input:
        read_1 = "../trimmed_reads/{id}_R1_qctrim.fastq.gz",
        read_2 = "../trimmed_reads/{id}_R2_qctrim.fastq.gz"
    output:
        bam = "../bam_files/{id}.bam"
    log: "../logs/{id}_samtools_import.log"
    benchmark: "../benchmarks/{id}_samtools_import.benchmark"
    shell:
        '''samtools import -i -1 {input.read_1} -2 {input.read_2} -o {output.bam} -O BAM -@ 4 -T '*' 2> {log}'''

# First alignment with STAR
rule star_align_1:
    input:
        bam_reads = "../bam_files/{id}.bam",
        genome_index = '../ref/star_index'
    output:
        star_out = '../star_align1/{id}.SJ.out.tab'
    params:
        genome_index = '../ref/star_index',
        outfile_prefix = '../star_align1/{id}.',
        outfile_tmp_dir = '../tmp/{id}'
    log: "../logs/{id}_star_align1.log"
    benchmark: "../benchmarks/{id}_star_align1.benchmark"
    resources:
        load = 25,
        tmpdir = '../tmp/'
    shell:
        '''	STAR \
        --runThreadN 4 \
        --genomeDir {params.genome_index} \
        --readFilesIn {input.bam_reads} \
        --readFilesType SAM PE \
        --readFilesCommand samtools view \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFileNamePrefix {params.outfile_prefix} \
        --outSAMtype None	\
        2> {log}'''	    

        # ENCODE alignment parameters
        # --outFilterMultimapNmax 20 \
        # --alignSJoverhangMin 8 \
        # --alignSJDBoverhangMin 1 \
        # --outFilterMismatchNmax 999 \
        # --outFilterMismatchNoverReadLmax 0.04 \
        # --alignIntronMin 20 \
        # --alignIntronMax 1000000 \
        # --alignMatesGapMax 1000000 \

# Prepare for realigning reads to assembled contigs
rule splice_junction_filter:
    input:
        splice_junctions = '../star_align1/{id}.SJ.out.tab'
    output:
        splice_junctions_filtered = '../star_align1/{id}.SJ.out.tab.filter'
    log: "../logs/{id}_splice_junction_filter.log"
    benchmark: "../benchmarks/{id}_splice_junction_filter.benchmark"
    shell:
        '''awk '$1~/chr[1-2XY]/ && $6==0 && $5>0 && $7>0' {input.splice_junctions} > {output.splice_junctions_filtered} 2> {log}'''

# Align reads to splice junctions and reference  
rule star_align_2:
    input:
        bam_reads = "../bam_files/{id}.bam",
        genome_index = '../ref/star_index',
        splice_junctions_filtered = '../star_align1/{id}.SJ.out.tab.filter'
    output:
        star_aligned2 = '../star_align2/{id}.Aligned.out.bam'
    params:
        genome_index = '../ref/star_index',
        outfile_prefix = '../star_align2/{id}.'
    log: "../logs/{id}_star_align2.log"
    benchmark: "../benchmarks/{id}_star_align2.benchmark"
    resources:
        load = 33
    priority: 50
    shell:
        '''STAR \
        --runThreadN 8 \
        --genomeDir {params.genome_index} \
	    --readFilesIn {input.bam_reads} \
        --readFilesType SAM PE \
        --readFilesCommand samtools view \
	    --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
	    --sjdbFileChrStartEnd {input.splice_junctions_filtered} \
	    --quantMode TranscriptomeSAM GeneCounts \
        --limitSjdbInsertNsj 20000000 \
	    --outBAMcompression 9 \
	    --outFileNamePrefix {params.outfile_prefix} \
	    --outSAMtype BAM Unsorted \
        2> {log}'''

        # ENCODE alignment parameters
        # --outFilterMultimapNmax 20 \
        # --alignSJoverhangMin 8 \
        # --alignSJDBoverhangMin 1 \
        # --outFilterMismatchNmax 999 \
        # --outFilterMismatchNoverReadLmax 0.04 \
        # --alignIntronMin 20 \
        # --alignIntronMax 1000000 \
        # --alignMatesGapMax 1000000 \

# Sort BAM files 
rule sort_bam:
	input: '../star_align2/{id}.Aligned.out.bam'
	output: '../star_align2/{id}.sorted.bam'
	log: "../logs/{id}_pos_sort_bam.log"
	benchmark: "../benchmarks/{id}_pos_sort_bam.benchmark"
	shell:
		'''samtools sort -@ 4 -o {output} {input} 2> {log}'''	

# Query sort aligned bam file
rule qsort_bam:
    input: "../star_align2/{id}.sorted.bam"
    output: "../star_align2/{id}.qsorted.bam"
    log: "../logs/{id}_qsort_bam.log"
    benchmark: "../benchmarks/{id}_qsort_bam.benchmark"
    shell:
        '''samtools sort -@ 4 -n -o {output} -O bam {input}'''

# Sort MBC output from agent trim and generate fastq.gz output
rule fgbio_sort_fastq:
    input:
        jar_file = "/home/dvargas/bin/fgbio/share/fgbio/fgbio.jar",
        mbc = "../agent_trim/{id}_MBC.txt.gz"
    output:
        mbc_fastq = "../dedup/{id}_MBC.fastq.gz"
    log: "../logs/{id}_fgbio_sort_fastq.log"
    benchmark: "../benchmarks/{id}_fgbio_sort_fastq.benchmark"
    shell:
        '''java -jar {input.jar_file} SortFastq -i {input.mbc} -o {output.mbc_fastq} 2> {log}'''

# De-duplicate reads and generate deduplicated bam files
rule agent_locatit:
    input:
        queryname_sorted_bam = "../star_align2/{id}.qsorted.bam",
        mbc_fastq = "../dedup/{id}_MBC.fastq.gz"
    output: "../dedup/{id}_dedup.bam"
    params:
        tmp = "../tmp/{id}"
    log: "../logs/{id}_agent_locatit.log"
    benchmark: "../benchmarks/{id}_agent_locatit.benchmark"
    shell:
        '''agent -Xmx20G locatit -i -R -C -q 0 -c 2500 -U -L -m 1 -d 1 \
        -X {params.tmp} -PM:xm,Q:xq,q:nQ,r:nR \
        -IB -OB -o {output} {input.queryname_sorted_bam} {input.mbc_fastq} 2> {log}'''

# Sort bam files
rule sort_dedup_bam:
	input: "../dedup/{id}_dedup.bam"
	output: "../dedup/{id}_dedup_sorted.bam"
	log: "../logs/{id}_dedup_pos_sort_bam.log"
	benchmark: "../benchmarks/{id}_dedup_pos_sort_bam.benchmark"
	shell:
		'''samtools sort -@ 4 -o {output} {input} 2> {log}'''

# Query sort bam files
rule qsort_dedup_bam:
    input: "../dedup/{id}_dedup.bam"
    output: "../dedup/{id}_dedup.qsorted.bam"
    log: "../logs/{id}_qsort_dedup_bam.log"
    benchmark: "../benchmarks/{id}_qsort_dedup_bam.benchmark"
    shell:
        '''samtools sort -m 8G -@ 4 -n -o {output} {input} 2> {log}'''

# Convert BAM to FASTQ
rule bam_to_fastq:
    input:
        bam = "../dedup/{id}_dedup.qsorted.bam"
    output:
        read_1 = "../fastq_files/{id}_dedup_R1.fastq",
        read_2 = "../fastq_files/{id}_dedup_R2.fastq"
    log: "../logs/{id}_bam2fastq.log"
	benchmark: "../benchmarks/{id}_bam2fastq.benchmark"
	shell:
		'''bedtools bamtofastq -i {input.bam} -fq {output.read_1} -fq2 {output.read_2} 2> {log}'''

# Compress fastq files
rule compress_fq:
    input:
        fq1 = "../fastq_files/{id}_dedup_R1.fastq",
        fq2 = "../fastq_files/{id}_dedup_R2.fastq"
    output:
        gz1 = "../fastq_files/{id}_dedup_R1.fastq.gz",
        gz2 = "../fastq_files/{id}_dedup_R2.fastq.gz"
    log: "../logs/{id}_compress_fq.log"
    benchmark: "../benchmarks/{id}_compress_fq.benchmark"
    priority: 50
    shell:
        '''pigz -8 {input.fq1} {input.fq2}'''


# Quantify with Salmon
rule salmon_quant:
    input:
        read_1 = "../fastq_files/{id}_dedup_R1.fastq.gz",
        read_2 = "../fastq_files/{id}_dedup_R2.fastq.gz",
        index = '../ref/salmon_index/'
    output:
        quant = '../salmon_quant/{id}/quant.sf'
    params:
        outfolder = '../salmon_quant/{id}'
    log:
        "../logs/{id}_salmon_quant.log"
    benchmark:
        "../benchmarks/{id}_salmon_quant.benchmark"
    resources:
        load = 33
    shell:
        '''salmon quant -i {input.index} -l A -1 {input.read_1} -2 {input.read_2} --validateMappings --seqBias --gcBias --posBias --threads 8 -o {params.outfolder} 2> {log}'''

# Index BAM files
rule index_bam:
	input: "../dedup/{id}_dedup_sorted.bam"
	output: "../dedup/{id}_dedup_sorted.bam.bai"
	log: "../logs/{id}_index_bam.log"
	benchmark: "../benchmarks/{id}_index_bam.benchmark"
	shell:
		'''samtools index -@ 8 {input} 2> {log}'''	

# Create statistics report 
rule samtools_idx_chrD:
	input:
		bam = "../dedup/{id}_dedup_sorted.bam",
		index = "../dedup/{id}_dedup_sorted.bam.bai"
	output:
		chrD = '../qorts/{id}.chrD.txt'
	log: "../logs/{id}_chrD.log"
	benchmark: "../benchmarks/{id}_chrD.benchmark"
	shell:
		'''
        samtools idxstats {input.bam} | cut -f1 | grep  -v '^chr[0-9]' | grep -v '^chrX\|chrY\|chrM\|ERCC\|\*' | perl -pe 's/\n/,/gi'| perl -pe 's/,$//gi' > {output.chrD} 2> {log}
        '''
	
# Run QoRTs
rule run_qorts:
    input:
        bam = "../dedup/{id}_dedup_sorted.bam",
        chrD = '../qorts/{id}.chrD.txt',
        gtf = "../ref/gencode.v44.primary_assembly.annotation.sorted.biotype.gtf"
    output:
        qorts = '../qorts/{id}/{id}.QC.summary.txt',
    params:
        outfile_prefix = '{id}',
        outfile_path = '../qorts/{id}/',
    log: "../logs/{id}_qorts.log",
    benchmark: "../benchmarks/{id}_qorts.benchmark",
    priority: 10
    shell:
        '''java -Xmx10G -jar /home/dvargas/bin/QoRTs.jar QC --stranded --outfilePrefix {params.outfile_prefix}. --dropChrom $(cat {input.chrD}) \
        {input.bam} {input.gtf} {params.outfile_path} 2> {log}'''