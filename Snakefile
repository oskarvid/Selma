# Import modules for tsv file handling and globbing of certain output files
import pandas as pd
import glob
import sys
import os
import shutil

# Change the working directory to 'workspace' to separate code from input and output files
workdir: 'workspace'

# Read functions.smk file to define the functions that are used to create the read group and select input files
include: 'scripts/functions.smk'

# Define path for the config file
configfile: 'config.yaml'

# This variable is define when the workflow is started by using "--config version=hg38|b37"
refversion = config['version']

# This variable is defined when the workflow is started by using '--config interval="${INTERVAL}"'
interval = config['interval']

# The following variables are used to select interval files and create output file names
SCATTERCOUNT = config['scattercount']
DIRECTORIES=[]
for i in range(SCATTERCOUNT):
	DIRECTORIES.append(str(i).zfill(4))

# Create variables to select sample names, lane numbers and flowcell names from the sample.tsv file
samples = pd.read_csv(config["samples"], sep='\t', dtype=str).set_index(["flowcell", "sample", "lane"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels]) # enforce str in index

# The following lines are used to dynamically set the CONTIGS variable
# The CONTIGS variable is only a range of values that are used to run BaseRecalibrator, ApplyBQSR and 
# GenotypeGVCFs in scatter/gather mode
ref_dict = config[refversion]['dict']

# Make a list of all contigs, extract the lengths, find the longest one
with open(ref_dict, "r") as ref_dict_file:
    sequence_tuple_list = []
    longest_sequence = 0
    for line in ref_dict_file:
        if line.startswith("@SQ"):
            line_split = line.split("\t")
            sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
    longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]

# Initialize the tsv string
string = " "

# Initialize variable for determination of total length of combined contig lengths
temp_size = sequence_tuple_list[0][1]

# For loop and conditional that goes through each contig and checks the combined length of the 
# contig lengths to create groups that are roughly the same length
for sequence_tuple in sequence_tuple_list[1:]:
    if temp_size + sequence_tuple[1] <= longest_sequence:
        temp_size += sequence_tuple[1]
    else:
        string += "\n"
        temp_size = sequence_tuple[1]

# add a final last line to "add" the unmapped contig job as well
string += "\n"

# the CONTIGS variable is only for genotypegvcfs because it doesn't have the last bed file which contains the unmapped regions
CONTIGS = range(0, len(string.splitlines()))

# the CONTIGSWUNMAPPED variable contains the unmapped regions and is used for basercalibrator and applybqsr
CONTIGSWUNMAPPED = range(0, len(string.splitlines()) + 1)


# Determine if the mode is CPU or GPU for NVScoreVariants
def mode():
	if config['mode'] == "gpu":
		mode = "gpu"
		return mode
	else:
		mode = "cpu"
		return mode

rule all:
	input:
		expand("Outputs/Stats/VcfPlots/{sample}_VcfPlots",
			sample=samples['sample']),
		expand("Outputs/Stats/BamPlots/{sample}/quals.gp",
			sample=samples['sample']),
		expand("Outputs/ApplyVqsrSnp/{sample}_SnpApplyVQSR.g.vcf.gz", 
			sample=samples['sample']),
		expand("Outputs/ApplyVqsrIndel/{sample}_IndelApplyVQSR.g.vcf.gz", 
			sample=samples['sample']),
#                expand("Outputs/DoppelMark/{sample}_doppelMark.bam",
#                        sample=samples['sample']),

# Split the interval list for HaplotypeCaller into sub intervals for scatter gather execution
rule MakeIntervalLists:
	input:
		interval = interval,
		fasta = config[refversion]['fasta'],
	output:
		"Outputs/MakeIntervalLists/{split}-scattered.interval_list",
	priority:
		30
	conda: "conda/gatk4.yaml"
	shell:
		"touch {output} && gatk \
		SplitIntervals \
		-L {input.interval} \
		-R {input.fasta} \
		--scatter-count {SCATTERCOUNT} \
		-O Outputs/MakeIntervalLists/"

# Map fastq files to reference genome
rule BwaMem:
	input:
		fastq1 = get_fastq1,
		fastq2 = get_fastq2,
		fasta = config[refversion]['fasta'],
	params:
		rgs = get_BwaRG,
	output:
#		temp("Outputs/BwaMem/{sample}_{lane}_{flowcell}.mapped.bam"),
		"Outputs/BwaMem/{sample}_{lane}_{flowcell}.mapped.bam",
	benchmark:
		"Outputs/benchmarks/{sample}_{lane}_{flowcell}.bwa.tsv",
	threads:
		12
	priority:
		0
	conda: "conda/bwa.yaml"
	shell:
		r"bwa mem -t {threads} \
		-R '{params.rgs}' \
		-M {input.fasta} \
		{input.fastq1} \
		{input.fastq2} \
		| samtools view -Sb - > {output}"

# Create unmapped bam files from the fastq files
rule FastqtoSam:
	input:
		fastq1 = get_fastq1,
		fastq2 = get_fastq2,
		fasta = config[refversion]['fasta'],
	output:
#		bam = temp("Outputs/FastqToSam/{sample}_{lane}_{flowcell}.unmapped.bam"),
#		tmp = directory(temp("Outputs/FastqToSam/{sample}_{lane}_{flowcell}.tmp")),
		bam = "Outputs/FastqToSam/{sample}_{lane}_{flowcell}.unmapped.bam",
		tmp = directory("Outputs/FastqToSam/{sample}_{lane}_{flowcell}.tmp"),
	benchmark:
		"Outputs/benchmarks/{sample}_{lane}_{flowcell}.FastqToSam.tsv",
	params:
		lane = get_FQLN,
		sample = get_FQSM,
		flowcell = get_FQFC,
		library = get_FQLIB,
	priority: 2
	conda: "conda/gatk4.yaml"
	shell:
		r"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		FastqToSam \
		-O {output.bam} \
		--PLATFORM ILLUMINA \
		--FASTQ {input.fastq1} \
		--FASTQ2 {input.fastq2} \
		--SAMPLE_NAME '{params.sample}' \
		--LIBRARY_NAME '{params.library}' \
		--READ_GROUP_NAME '{params.flowcell}.{params.lane}' \
		--PLATFORM_UNIT '{params.flowcell}.{params.lane}.{params.sample}' \
		--TMP_DIR {output.tmp}"

# Merge output files from bwa and FastqToSam
rule MergeBamAlignment:
	input:
		fasta = config[refversion]['fasta'],
		mapped = "Outputs/BwaMem/{sample}_{lane}_{flowcell}.mapped.bam",
		unmapped = "Outputs/FastqToSam/{sample}_{lane}_{flowcell}.unmapped.bam",
	output:
		bam = "Outputs/MergeBamAlignment/{sample}_{lane}_{flowcell}.merged.bam",
#		tmp = directory(temp("Outputs/MergeBamAlignment/{sample}_{lane}_{flowcell}.tmp")),
		tmp = directory("Outputs/MergeBamAlignment/{sample}_{lane}_{flowcell}.tmp"),
	benchmark:
		"Outputs/benchmarks/{sample}_{lane}_{flowcell}.MergeBamAlignments.tsv",
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options '-Xmx14G -Djava.io.tempdir=$(pwd)/tmp' \
		MergeBamAlignment \
		-O {output.bam} \
		--ADD_MATE_CIGAR true \
		--CLIP_ADAPTERS false \
		--SORT_ORDER coordinate \
		--ATTRIBUTES_TO_RETAIN X0 \
		--ALIGNED_READS_ONLY false \
		--EXPECTED_ORIENTATIONS FR \
		--MAX_RECORDS_IN_RAM 2000000 \
		--PROGRAM_RECORD_ID 'bwamem' \
		--ALIGNED_BAM {input.mapped} \
		--PROGRAM_GROUP_NAME 'bwamem' \
		--IS_BISULFITE_SEQUENCE false \
		--VALIDATION_STRINGENCY SILENT \
		--UNMAPPED_BAM {input.unmapped} \
		--MAX_INSERTIONS_OR_DELETIONS -1 \
		--REFERENCE_SEQUENCE {input.fasta} \
		--PROGRAM_GROUP_VERSION '0.7.12-r1039' \
		--PRIMARY_ALIGNMENT_STRATEGY MostDistant \
		--PROGRAM_GROUP_COMMAND_LINE 'bwa mem -t 15 -R -M Input1 Input2 > output.sam' \
		--TMP_DIR {output.tmp}"

# Checkpoint so that MarkDuplicates can find the output files from MergeBamAlignment
checkpoint MarkDupCheckpoint:
	input:
		expand("Outputs/MergeBamAlignment/{sample}_{lane}_{flowcell}.merged.bam", zip,
			sample=samples['sample'],
			lane=samples['lane'],
			flowcell=samples['flowcell']),
	output:
		touch("Outputs/MergeBamAlignment/placeholder"),
	shell:
		"echo 'Running placeholder checkpoint rule to create correct dependency for MarkDuplicates to start after MergeBamAlignment and be able to find the output files correctly'"

rule DoppelMark:
	input:
		flag = "Outputs/MergeBamAlignment/placeholder",
		files = lambda wcs: glob.glob('Outputs/MergeBamAlignment/%s*.bam' % wcs.sample),
		index = lambda wcs: glob.glob('Outputs/MergeBamAlignment/%s*.bai' % wcs.sample),
	singularity: "docker://oskarv/doppelmark"
	output:
		bam = "Outputs/DoppelMark/{sample}_doppelMark.bam",
		bai = "Outputs/DoppelMark/{sample}_doppelMark.bai",
	threads: 1
	shell:
		"doppelmark \
		--bam {input.files} \
		--index {input.index} \
		--clip-padding 250 \
		--output {output.bam} \
		--parallelism {threads}"

# Mark duplicates in the output files from MergeBamAlignment
rule MarkDup:
	input:
		flag = "Outputs/MergeBamAlignment/placeholder",
		files = lambda wcs: glob.glob('Outputs/MergeBamAlignment/%s*.bam' % wcs.sample),
	output:
#		tmp = directory(temp("Outputs/MarkDuplicates/{sample}_tmp")),
#		bam = temp("Outputs/MarkDuplicates/{sample}_markedDuplicates.bam"),
#		bai = temp("Outputs/MarkDuplicates/{sample}_markedDuplicates.bai"),
		tmp = directory("Outputs/MarkDuplicates/{sample}_tmp"),
		bam = "Outputs/MarkDuplicates/{sample}_markedDuplicates.bam",
		bai = "Outputs/MarkDuplicates/{sample}_markedDuplicates.bai",
		metrics = "Outputs/MarkDuplicates/{sample}_markedDuplicates.metrics",
	benchmark:
		"Outputs/benchmarks/{sample}.MarkDup.tsv",
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options '-Xmx14G -Djava.io.tempdir=$(pwd)/tmp' \
		MarkDuplicates \
		-O {output.bam} \
		--CREATE_INDEX true \
		--VALIDATION_STRINGENCY LENIENT \
		--METRICS_FILE {output.metrics} \
		--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 200000 \
		$(echo ' {input.files}' | sed 's/ / --INPUT /g') \
		--TMP_DIR {output.tmp}"

# This rule creates bed files with the contigs grouped into roughly equal lengths for each file.
rule MakeSequenceGroupings:
	input:
		config[refversion]['dict'],
	output:
		touch("Outputs/MakeContigBeds/flag"),
	priority:
		30
	shell:
		"python scripts/split-bedfile.py {input} Outputs/MakeContigBeds/"

# Run checkpoint so that glob() in BaseRecalibrator finds the input files
checkpoint BaseRecalibratorCheckpoint:
	input:
		"Outputs/MakeContigBeds/flag",
	output:
		touch("Outputs/MakeContigBeds/placeholder"),
	shell:
		"echo 'Running checkpoint rule to create correct dependency for BaseRecalibrator to start after MakeContigBeds and be able to find the bed files correctly'"

# Do base quality score recalibration
rule BaseRecalibrator:
	input:
		fasta = config[refversion]['fasta'],
		dbsnp = config[refversion]['dbsnp'],
		mills = config[refversion]['mills'],
		v1000g = config[refversion]['v1000g'],
		flag = "Outputs/MakeContigBeds/placeholder",
#                bam = "Outputs/DoppelMark/{sample}_doppelMark.bam",
		bam = "Outputs/MarkDuplicates/{sample}_markedDuplicates.bam",
		bai = "Outputs/MarkDuplicates/{sample}_markedDuplicates.bai",
		contigs = lambda wcs: glob.glob('Outputs/MakeContigBeds/contigs_%s.bed' % wcs.contigs),
	threads:
		1
	output:
		grp = "Outputs/BaseRecalibrator/{sample}_BQSR_{contigs}.grp",
	benchmark:
		"Outputs/benchmarks/{sample}_{contigs}.BaseRecalibrator.tsv",
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options '-Xmx1330M -Djava.io.tempdir=$(pwd)/tmp' \
		BaseRecalibrator \
		-O {output.grp} \
		--input {input.bam} \
		--reference {input.fasta} \
		--known-sites {input.mills} \
		--known-sites {input.dbsnp} \
		--known-sites {input.v1000g} \
		--tmp-dir Outputs/BaseRecalibrator \
		$(cat {input.contigs})"

# Merge bqsr files from BaseRecalibrator into one
rule GatherBQSRReports:
	input:
		expand("Outputs/BaseRecalibrator/{{sample}}_BQSR_{directory}.grp",
			directory=CONTIGSWUNMAPPED),
	output:
		"Outputs/GatherBQSR/{sample}_GatheredBQSR.grp"
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		GatherBQSRReports \
		-O {output} \
		$(echo ' {input}' | sed 's/ / --input /g')"

# Apply base quality score recalibration 
rule ApplyBQSR:
	input:
		fasta = config[refversion]['fasta'],
		flag = "Outputs/MakeContigBeds/placeholder",
		grp = "Outputs/GatherBQSR/{sample}_GatheredBQSR.grp",
		bam = "Outputs/MarkDuplicates/{sample}_markedDuplicates.bam",
		bai = "Outputs/MarkDuplicates/{sample}_markedDuplicates.bai",
		contigs = lambda wcs: glob.glob('Outputs/MakeContigBeds/contigs_%s.bed' % wcs.con),
	output:
		bam = "Outputs/ApplyBQSR/{sample}_{con}_recalibrated.bam",
		bai = "Outputs/ApplyBQSR/{sample}_{con}_recalibrated.bai",
	benchmark:
		"Outputs/benchmarks/{sample}_{con}.ApplyBQSR.tsv",
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options '-Xmx1330M -Djava.io.tempdir=$(pwd)/tmp' \
		ApplyBQSR \
		-O {output.bam} \
		-bqsr {input.grp} \
		--input {input.bam} \
		--reference {input.fasta} \
		$(cat {input.contigs}) \
		--create-output-bam-index true \
		--tmp-dir Outputs/ApplyBQSR"

# Merge bam files from ApplyBQSR into one
rule GatherApplyBQSRbams:
	input:
		bai = expand("Outputs/ApplyBQSR/{{sample}}_{contigs}_recalibrated.bai",
			contigs=CONTIGSWUNMAPPED),
		bam = expand("Outputs/ApplyBQSR/{{sample}}_{contigs}_recalibrated.bam",
			contigs=CONTIGSWUNMAPPED),
	output:
		bam = "Outputs/GatherBamFiles/{sample}_GatheredABQSRFiles.bam",
		bai = "Outputs/GatherBamFiles/{sample}_GatheredABQSRFiles.bai",
	benchmark:
		"Outputs/benchmarks/{sample}.GatheredBamFiles.tsv",
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		GatherBamFiles \
		-O {output.bam} \
		$(echo ' {input.bam}' | sed 's/ / --INPUT /g') \
		--CREATE_INDEX true"


rule BamStats:
	input:
		bam = "Outputs/GatherBamFiles/{sample}_GatheredABQSRFiles.bam",
		bai = "Outputs/GatherBamFiles/{sample}_GatheredABQSRFiles.bai",
	output:
		"Outputs/Stats/BamPlots/{sample}_BamStats",
	benchmark:
		"Outputs/benchmarks/{sample}.BamStats.tsv",
	threads:
		1
	priority: 1
	conda: "conda/samtools.yaml"
	shell:
		"samtools stats \
		--threads {threads} \
		{input.bam} > {output}"

rule BamPlot:
	input:
		"Outputs/Stats/BamPlots/{sample}_BamStats",
	output:
		gp = "Outputs/Stats/BamPlots/{sample}/quals.gp",
	params:
		prefix = "Outputs/Stats/BamPlots/{sample}/",
	benchmark:
		"Outputs/benchmarks/{sample}.BamPlot.tsv",
	priority: 1
	conda: "conda/samtools.yaml"
	shell:
		"plot-bamstats \
		-p {params.prefix} \
		{input}"

# Call germline SNPs and indels
# The default is to use an interval list for wgs (that excludes regions such as centromeres) and for manual parallelization.
# {threads} is set to 2 because Snakemake interprets the core count as cores*2=threads, but that's not right. Therefore it's
# necessary to make it think HaplotypeCaller consumes 2 threads per process so that only 16 jobs are started if Snakemake is
# given 16 cores, otherwise it would try to start 32 jobs which would consume way too much RAM.
rule HaplotypeCaller:
	input:
		fasta = config[refversion]['fasta'],
		bam = "Outputs/GatherBamFiles/{sample}_GatheredABQSRFiles.bam",
		bai = "Outputs/GatherBamFiles/{sample}_GatheredABQSRFiles.bai",
		intervals = "Outputs/MakeIntervalLists/{split}-scattered.interval_list",
	output:
		vcf = "Outputs/HaplotypeCaller/{sample}_{split}_rawVariants.g.vcf.gz",
		tbi = "Outputs/HaplotypeCaller/{sample}_{split}_rawVariants.g.vcf.gz.tbi",
	benchmark:
		"Outputs/benchmarks/{sample}_{split}.HaplotypeCaller.tsv",
	threads:
		1
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options '-Xmx1100M -Djava.io.tempdir=$(pwd)/tmp' \
		HaplotypeCaller \
		-ERC GVCF \
		-I {input.bam} \
		-O {output.vcf} \
		-R {input.fasta} \
		-L {input.intervals} \
		--tmp-dir Outputs/HaplotypeCaller"

# The output files from HaplotypeCaller need to be gathered into one file
# This tool is currently unable to index the output file, the indexing tool below handles that
rule GatherHTCVCFs:
	input:
		vcf = expand("Outputs/HaplotypeCaller/{{sample}}_{split}_rawVariants.g.vcf.gz", 
			split=DIRECTORIES),
		tbi = expand("Outputs/HaplotypeCaller/{{sample}}_{split}_rawVariants.g.vcf.gz.tbi",
			split=DIRECTORIES),
	output:
		vcf = "Outputs/GatherHTCVCFs/{sample}_GatheredHTCVCFs.g.vcf.gz",
	benchmark:
		"Outputs/benchmarks/{sample}.GatherHTCVCFs.tsv",
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		GatherVcfs \
		-O {output.vcf} \
		$(echo ' {input.vcf}' | sed 's/ / --INPUT /g')"

# Index the output file from GatherHTCVCFs
rule IndexGatheredHTCVCFs:
	input:
		"Outputs/GatherHTCVCFs/{sample}_GatheredHTCVCFs.g.vcf.gz",
	output:
		"Outputs/GatherHTCVCFs/{sample}_GatheredHTCVCFs.g.vcf.gz.tbi",
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options '-Xmx14G -Djava.io.tempdir=$(pwd)/tmp' \
		IndexFeatureFile \
		-I {input} \
		-O {output}"

# Perform joint genotyping
rule GenotypeGVCFs:
	input:
		fasta = config[refversion]['fasta'],
		flag = "Outputs/MakeContigBeds/placeholder",
		vcf = "Outputs/GatherHTCVCFs/{sample}_GatheredHTCVCFs.g.vcf.gz",
		tbi = "Outputs/GatherHTCVCFs/{sample}_GatheredHTCVCFs.g.vcf.gz.tbi",
		contigs = lambda wcs: glob.glob('Outputs/MakeContigBeds/contigs_%s.bed' % wcs.con),
	output:
		vcf = "Outputs/GenotypeGVCFs/{sample}_{con}_genotypes.g.vcf.gz",
		tbi = "Outputs/GenotypeGVCFs/{sample}_{con}_genotypes.g.vcf.gz.tbi",
	benchmark:
		"Outputs/benchmarks/{sample}_{con}.GenotypeGVCFs.tsv",
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options '-Xmx1330M -Djava.io.tempdir=$(pwd)/tmp' \
		GenotypeGVCFs \
		-V {input.vcf} \
		-O {output.vcf} \
		-R {input.fasta} \
		$(cat {input.contigs}) \
		--tmp-dir Outputs/GenotypeGVCFs"

# The output files from GenotypeGVCFs need to be gathered into one file
# This tool is currently unable to index the output file, the indexing tool below handles that
rule GatherGenotypeGVCFs:
	input:
		tbi = expand("Outputs/GenotypeGVCFs/{{sample}}_{contigs}_genotypes.g.vcf.gz.tbi",
			contigs=CONTIGS),
		vcf = expand("Outputs/GenotypeGVCFs/{{sample}}_{contigs}_genotypes.g.vcf.gz", 
			contigs=CONTIGS),
	output:
		vcf = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz",
	benchmark:
		"Outputs/benchmarks/{sample}.GatherGenotypeGVCFs.tsv",
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		GatherVcfs \
		-O {output.vcf} \
		$(echo ' {input.vcf}' | sed 's/ / --INPUT /g')"

# Index the output file from GatherGenotypeGVCFs
rule IndexGatheredGVCFs:
	input:
		"Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz",
	output:
		"Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz.tbi",
	benchmark:
		"Outputs/benchmarks/{sample}.IndexGatheredGenotypeGVCFs.tsv",
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options '-Xmx14G -Djava.io.tempdir=$(pwd)/tmp' \
		IndexFeatureFile \
		-I {input} \
		-O {output}"

rule VcfStats:
	input:
		vcf = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz",
		tbi = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz.tbi"
	output:
		"Outputs/Stats/VcfPlots/{sample}/BcftoolsStats",
	benchmark:
		"Outputs/benchmarks/{sample}.BcftoolsStats.tsv",
	threads:
		8
	priority: 1
	conda: "conda/bcftools.yaml"
	shell:
		"bcftools stats \
		--threads {threads} \
		{input.vcf} > {output}"

rule VcfPlot:
	input:
		"Outputs/Stats/VcfPlots/{sample}/BcftoolsStats",
	output:
		directory("Outputs/Stats/VcfPlots/{sample}_VcfPlots")
	benchmark:
		"Outputs/benchmarks/{sample}.VcfPlot.tsv",
	priority: 1
	conda: "conda/bcftools.yaml"
	shell:
		"plot-vcfstats \
		-p {output} \
		{input} || true"

rule NVScoreVariants:
	input:
		vcf = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz",
		bam = "Outputs/GatherBamFiles/{sample}_GatheredABQSRFiles.bam",
		bai = "Outputs/GatherBamFiles/{sample}_GatheredABQSRFiles.bai",
		fasta = config[refversion]['fasta'],
		model = "NVScoreVariants_models/",
	output:
		vcf = "Outputs/NVScoreVariants/{sample}_NVScoreVariants.vcf"
	params:
		mode = mode()
#		mode = config['mode'],
	conda: "conda/nvscorevariants_environment.yml"
	shell:
		"""
		if [[ {params.mode} == "gpu" ]]; then
			python3 scripts/nvscorevariants.py \
			--output-file {output.vcf} \
			--vcf-file {input.vcf} \
			--ref-file {input.fasta} \
			--tensor-type read_tensor \
			--batch-size 16 \
			--seed 724 \
			--tmp-file /tmp/nvscoretemp.txt \
			--model-directory {input.model} \
			--input-file {input.bam} \
			--gpus 0 
		else
			java -jar ../gatk-package-NVSCOREVARIANTS-PREVIEW-SNAPSHOT-local.jar NVScoreVariants \
			-V {input.vcf} \
			-R {input.fasta} \
			--tensor-type read_tensor \
	                -I {input.bam} \
		        --batch-size 16 \
			-O {output}
		fi
		"""

#		else:

rule VcfStats2:
	input:
		vcf = "Outputs/NVScoreVariants/{sample}_NVScoreVariants.vcf",
	output:
		"Outputs/Stats/VcfPlots/{sample}/BcftoolsStats-NVScore",
#	benchmark:
#		"Outputs/benchmarks/{sample}.BcftoolsStats.tsv",
	threads:
		8
	priority: 1
	conda: "conda/bcftools.yaml"
	shell:
		"bcftools stats \
		--threads {threads} \
		{input.vcf} > {output}"

rule VcfPlot2:
	input:
		"Outputs/Stats/VcfPlots/{sample}/BcftoolsStats-NVScore",
	output:
		directory("Outputs/Stats/VcfPlots/{sample}_VcfPlots-NVScore")
#	benchmark:
#		"Outputs/benchmarks/{sample}.VcfPlot.tsv",
	priority: 1
	conda: "conda/bcftools.yaml"
	shell:
		"plot-vcfstats \
		-p {output} \
		{input} || true"

# Build a recalibration model to score variant quality for filtering purposes
rule VariantRecalibratorSNP:
	input:
		omni = config[refversion]['omni'],
		fasta = config[refversion]['fasta'],
		dbsnp = config[refversion]['dbsnp'],
		v1000g = config[refversion]['v1000g'],
		hapmap = config[refversion]['hapmap'],
		vcf = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz",
		tbi = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz.tbi",
	output:
#		recal = temp("Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.recal"),
#		idx = temp("Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.recal.idx"),
#		tranches = temp("Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.tranches"),
		recal = "Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.recal",
		idx = "Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.recal.idx",
		tranches = "Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.tranches",
	benchmark:
		"Outputs/benchmarks/{sample}.VariantRecalibratorSNP.tsv",
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options '-Xmx7G -Djava.io.tempdir=$(pwd)/tmp' \
		VariantRecalibrator \
		--mode SNP \
		-V {input.vcf} \
		-R {input.fasta} \
		--max-gaussians 6 \
		--output {output.recal} \
		--tranches-file {output.tranches} \
		-an QD -an MQ -an DP -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
		-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 \
		-tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 \
		-tranche 97.0 -tranche 90.0 \
		--resource:omni,known=false,training=true,truth=true,prior=12.0 {input.omni} \
		--resource:dbsnp,known=true,training=false,truth=false,prior=7.0 {input.dbsnp} \
		--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} \
		--resource:v1000G,known=false,training=true,truth=false,prior=10.0 {input.v1000g} \
		--tmp-dir Outputs/VariantRecalibratorSNP"

# Build a recalibration model to score variant quality for filtering purposes
rule VariantRecalibratorINDEL:
	input:
		fasta = config[refversion]['fasta'],
		dbsnp = config[refversion]['dbsnp'],
		mills = config[refversion]['mills'],
		vcf = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz",
		tbi = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz.tbi",
	output:
#		recal = temp("Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.recal"),
#		idx = temp("Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.recal.idx"),
#		tranches = temp("Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.tranches"),
		recal = "Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.recal",
		idx = "Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.recal.idx",
		tranches = "Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.tranches",
	benchmark:
		"Outputs/benchmarks/{sample}.VariantRecalibratorINDEL.tsv",
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options '-Xmx7G -Djava.io.tempdir=$(pwd)/tmp' \
		VariantRecalibrator \
		--mode INDEL \
		-V {input.vcf} \
		-R {input.fasta} \
		--max-gaussians 4 \
		--output {output.recal} \
		--tranches-file {output.tranches} \
		-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
		--resource:mills,known=false,training=true,truth=true,prior=12.0 {input.mills} \
		--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} \
		-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 \
		-tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
		--tmp-dir Outputs/VariantRecalibratorINDEL"

# Apply a score cutoff to filter variants based on a recalibration table
rule ApplyVqsrSnp:
	input:
		fasta = config[refversion]['fasta'],
		recal = "Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.recal",
		idx = "Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.recal.idx",
		vcf = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz",
		tranches = "Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.tranches",
		tbi = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz.tbi",
	output:
		vcf = "Outputs/ApplyVqsrSnp/{sample}_SnpApplyVQSR.g.vcf.gz",
	benchmark:
		"Outputs/benchmarks/{sample}.ApplyVqsrSnp.tsv",
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options '-Xmx7G -Djava.io.tempdir=$(pwd)/tmp' \
		ApplyVQSR \
		--mode SNP \
		-V {input.vcf} \
		-O {output.vcf} \
		-R {input.fasta} \
		-ts-filter-level 99.7 \
		-recal-file {input.recal} \
		-tranches-file {input.tranches} \
		--tmp-dir Outputs/ApplyVqsrSnp"

# Apply a score cutoff to filter variants based on a recalibration table
rule ApplyVqsrIndel:
	input:
		fasta = config[refversion]['fasta'],
		vcf = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz",
		recal = "Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.recal",
		idx = "Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.recal.idx",
		tbi = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz.tbi",
		tranches = "Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.tranches",
	output:
		vcf = "Outputs/ApplyVqsrIndel/{sample}_IndelApplyVQSR.g.vcf.gz",
	benchmark:
		"Outputs/benchmarks/{sample}.ApplyVqsrIndel.tsv",
	conda: "conda/gatk4.yaml"
	shell:
		"gatk --java-options '-Xmx7G -Djava.io.tempdir=$(pwd)/tmp' \
		ApplyVQSR \
		--mode INDEL \
		-V {input.vcf} \
		-O {output.vcf} \
		-R {input.fasta} \
		-ts-filter-level 99.7 \
		-recal-file {input.recal} \
		-tranches-file {input.tranches} \
		--tmp-dir Outputs/ApplyVqsrIndel"
