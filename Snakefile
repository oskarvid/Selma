# Import modules for tsv file handling and globbing of certain output files
import pandas as pd
import glob

# Change the working directory to '/data/workspace' to separate code from input and output files
workdir: '/data/workspace'

# Read functions.smk file to define the functions that are used to create the read group and select input files
include: '/data/functions.smk'

# Define path for the config file
configfile: '/data/config.yaml'

# This variable is defined when the workflow is started by using "--config version=hg38|b37"
refversion = config['version']

# The following variables are used to select interval files and create output file names
CONTIGCOUNT = config[refversion]['contigfiles']
CONTIGS = range(1, CONTIGCOUNT)
GATHERCONTIGS = range(1, (( CONTIGCOUNT + 1 )))

SCATTERCOUNT = config['scattercount']
DIRECTORIES = range(1, (( SCATTERCOUNT + 1 )))

# Create variables to select sample names, lane numbers and flowcell names from the sample.tsv file
samples = pd.read_csv(config["samples"], sep='\t', dtype=str).set_index(["flowcell", "sample", "lane"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index

# Define variable names and collect output files from the final rule
rule all:
	input:
		expand("Outputs/BwaMem/{sample}_{lane}_{flowcell}.mapped.bam", zip,
			sample=samples['sample'], 
			lane=samples['lane'],
			flowcell=samples['flowcell']),
		expand("Outputs/FastqToSam/{sample}_{lane}_{flowcell}.unmapped.bam", zip,
			sample=samples['sample'], 
			lane=samples['lane'],
			flowcell=samples['flowcell']),
		expand("Outputs/ApplyVqsrSnp/{sample}_SnpApplyVQSR.g.vcf.gz", 
			sample=samples['sample']),
		expand("Outputs/ApplyVqsrIndel/{sample}_IndelApplyVQSR.g.vcf.gz", 
			sample=samples['sample']),

# This rule creates a tsv file with the contigs grouped into roughly equal combined lengths.
# The script should be improved to create a separate file per group so as to avoid the need \
# for the MakeTSVs.sh file which performs this function.
rule MakeSequenceGroupings:
	input:
		config[refversion]['dict'],
	output:
		"Outputs/MakeContigBeds/sequence_grouping_with_unmapped.tsv",
	priority:
		30
	shell:
		"python2 /data/scripts/split-bedfile.py {input} Outputs/MakeContigBeds/"

# Split the sequence_grouping_with_unmapped.tsv file into one file per line.
# This script should be made obsolete by performing its function in the split-bedfile.py instead.
rule MakeContigBeds:
	input:
		"Outputs/MakeContigBeds/sequence_grouping_with_unmapped.tsv"
	output:
		flag = touch("Outputs/MakeContigBeds/flag"),
	priority:
		30
	shell:
		"bash /data/scripts/MakeTSVs.sh {input} Outputs/MakeContigBeds/"

# Split the interval list for HaplotypeCaller into sub intervals for scatter gather execution
rule MakeIntervalLists:
	input:
		config[refversion]['intervals'],
	output:
		expand("Outputs/MakeIntervalLists/{directory}_of_{total}/scattered.interval_list", 
			directory=DIRECTORIES, 
			total=SCATTERCOUNT),
	priority:
		30
	shell:
		"python2 /data/scripts/create_scatter_intervals.py {input} {SCATTERCOUNT} 4 'Outputs/MakeIntervalLists' 'This is a placeholder that should or could contain information about something useful regarding the interval lists'"

# Map fastq files to reference genome
rule BwaMem:
	input:
		fasta = config[refversion]['fasta'],
		fastq1 = get_fastq1,
		fastq2 = get_fastq2,
	params:
		rgs = get_BwaRG,
	output:
		temp("Outputs/BwaMem/{sample}_{lane}_{flowcell}.mapped.bam"),
	threads:
		12
	priority:
		0
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
		fasta = config[refversion]['fasta'],
		fastq1 = get_fastq1,
		fastq2 = get_fastq2,
	output:
		bam = temp("Outputs/FastqToSam/{sample}_{lane}_{flowcell}.unmapped.bam"),
		tmp = directory(temp("Outputs/FastqToSam/{sample}_{lane}_{flowcell}.tmp")),
	params:
		sample = get_FQSM,
		flowcell = get_FQFC,
		lane = get_FQLN,
		library = get_FQLIB,
	priority: 1
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

# Due to the hack used to define the input files for MarkDuplicates, the bam \
# file from MergeBamAlignment cannot be flagged as temp()
rule MergeBamAlignment:
	input:
		fasta = config[refversion]['fasta'],
		mapped = "Outputs/BwaMem/{sample}_{lane}_{flowcell}.mapped.bam",
		unmapped = "Outputs/FastqToSam/{sample}_{lane}_{flowcell}.unmapped.bam",
	output:
		bam = "Outputs/MergeBamAlignment/{sample}_{lane}_{flowcell}.merged.bam",
		tmp = directory(temp("Outputs/MergeBamAlignment/{sample}_{lane}_{flowcell}.tmp")),
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		MergeBamAlignment \
		-O {output.bam} \
		--ADD_MATE_CIGAR true \
		--CLIP_ADAPTERS false \
		--SORT_ORDER coordinate \
		--ATTRIBUTES_TO_RETAIN X0 \
		--ALIGNED_READS_ONLY false \
		--EXPECTED_ORIENTATIONS FR \
		--MAX_RECORDS_IN_RAM 200000 \
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
		--PROGRAM_GROUP_COMMAND_LINE 'bwa mem -t 18 -R -M Input1 Input2 > output.sam' \
		--TMP_DIR {output.tmp}"

# This checkpoint placeholder was the only way I could find that \
# let me run MarkDuplicates correctly, it's a hack, but it works, \
# feel free to solve it cleanly, I don't know how to.
# A clean solution would allow the current glob() solution but without \
# the need for the checkpoint. Unfortunately snakemake needs the output files \
# to exist for the glob() function to work, clearly this is impossible when \
# the workflow is started, and this is therefore a limitation of snakemake. 
# A clean solution would probably create a string with the input file path as \
# input for MarkDuplicates, but for now the checkpoint solution does the trick.
# It would have been convenient to make MergeBamAlignment the checkpoint \
# but it wouldn't recognize it for some reason, hence the checkpoint rule hack below.

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

# Mark duplicates in the output files from MergeBamAlignment
rule MarkDup:
	input:
		flag = "Outputs/MergeBamAlignment/placeholder",
		files = lambda wcs: glob.glob('Outputs/MergeBamAlignment/%s*.bam' % wcs.sample),
	output:
		tmp = directory(temp("Outputs/MarkDuplicates/{sample}_tmp")),
		bam = temp("Outputs/MarkDuplicates/{sample}_markedDuplicates.bam"),
		metrics = "Outputs/MarkDuplicates/{sample}_markedDuplicates.metrics",
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		MarkDuplicates \
		-O {output.bam} \
		--CREATE_INDEX true \
		--VALIDATION_STRINGENCY LENIENT \
		--METRICS_FILE {output.metrics} \
		--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 200000 \
		$(echo ' {input.files}' | sed 's/ / --INPUT /g') \
		--TMP_DIR {output.tmp}"

# Ideally the intervals in the contig bed file would be stored in a list or something like that, \
# that way it would probably be possible to use a cleaner solution than "$(cat {input.contigs})" to \
# supply the "--intervals chr1:1+" etc. 
rule BaseRecalibrator:
	input:
		fasta = config[refversion]['fasta'],
		dbsnp = config[refversion]['dbsnp'],
		mills = config[refversion]['mills'],
		v1000g = config[refversion]['v1000g'],
		flag = "Outputs/MakeContigBeds/flag",
		bam = "Outputs/MarkDuplicates/{sample}_markedDuplicates.bam",
		contigs = lambda wcs: glob.glob('Outputs/MakeContigBeds/contigs_%s.bed' % wcs.contigs),
	threads:
		1
	output:
		grp = temp("Outputs/BaseRecalibrator/{sample}_BQSR_{contigs}.grp"),
	shell:
		"gatk --java-options '-Xmx3375M -Djava.io.tempdir=$(pwd)/tmp' \
		BaseRecalibrator \
		-O {output.grp} \
		--input {input.bam} \
		--reference {input.fasta} \
		--known-sites {input.mills} \
		--known-sites {input.dbsnp} \
		--known-sites {input.v1000g} \
		--tmp-dir Outputs/BaseRecalibrator \
		$(cat {input.contigs})"

# The bqsr files from BaseRecalibrator need to be merged into one
# The input substitution with the sed pipe is a bit hacky, the previous solution \
# relied on a run directive with python code to create the {input} definition, \
# but due to the use of singularity this isn't possible, therefore the current bash \
# command substitution is used.
rule GatherBQSRReports:
	input:
		expand("Outputs/BaseRecalibrator/{{sample}}_BQSR_{directory}.grp",
			directory=CONTIGS),
	output:
		temp("Outputs/GatherBQSR/{sample}_GatheredBQSR.grp")
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		GatherBQSRReports \
		-O {output} \
		$(echo ' {input}' | sed 's/ / --input /g')"

# This rule applies the recalibration that BaseRecalibrator calculated, due to the usage of glob() \
# it's necessary to include the flag from MakeContigBeds, otherwise glob() won't find the input files.
# The cat input substitution can probably be done smoother somehow, it's necessary due to the use of glob()
rule ApplyBQSR:
	input:
		fasta = config[refversion]['fasta'],
		flag = "Outputs/MakeContigBeds/flag",
		grp = "Outputs/GatherBQSR/{sample}_GatheredBQSR.grp",
		bam = "Outputs/MarkDuplicates/{sample}_markedDuplicates.bam",
		contigs = lambda wcs: glob.glob('Outputs/MakeContigBeds/contigs_%s.bed' % wcs.con),
	output:
		bam = temp("Outputs/ApplyBQSR/{sample}_{con}_recalibrated.bam"),
	shell:
		"gatk --java-options '-Xmx3375M -Djava.io.tempdir=$(pwd)/tmp' \
		ApplyBQSR \
		-O {output.bam} \
		-bqsr {input.grp} \
		--input {input.bam} \
		--reference {input.fasta} \
		$(cat {input.contigs}) \
		--create-output-bam-index true \
		--tmp-dir Outputs/ApplyBQSR"

# This rule merges the bam files from ApplyBQSR into one.
# The input substitution with the sed pipe is a bit hacky, the previous solution \
# relied on a run directive with python code to create the {input} definition, \
# but due to the use of singularity this isn't possible, therefore the current bash \
# command substitution is used.
rule GatherApplyBQSRbams:
	input:
		expand("Outputs/ApplyBQSR/{{sample}}_{directory}_recalibrated.bam",
			directory=GATHERCONTIGS),
	output:
		bam = temp("Outputs/GatherBamFiles/{sample}_GatheredBamFiles.bam"),
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		GatherBamFiles \
		-O {output.bam} \
		$(echo ' {input}' | sed 's/ / --INPUT /g') \
		--CREATE_INDEX true"

# Call germline SNPs and indels, use interval list for wgs that excludes regions such as centromeres and \
# also for manual parallelization
rule HaplotypeCaller:
	input:
		fasta = config[refversion]['fasta'],
		bam = "Outputs/GatherBamFiles/{sample}_GatheredBamFiles.bam",
		intervals = "Outputs/MakeIntervalLists/{directory}_of_30/scattered.interval_list",
	output:
		vcf = temp("Outputs/HaplotypeCaller/{sample}_{directory}_rawVariants.g.vcf.gz"),
	threads: 1
	shell:
		"gatk --java-options '-Xmx3375M -Djava.io.tempdir=$(pwd)/tmp' \
		HaplotypeCaller \
		-ERC GVCF \
		-I {input.bam} \
		-O {output.vcf} \
		-R {input.fasta} \
		-L {input.intervals} \
		--use-new-qual-calculator TRUE \
		--native-pair-hmm-threads {threads} \
		--tmp-dir Outputs/HaplotypeCaller"

# The output files from HaplotypeCaller need to be merged
# The input substitution with the sed pipe is a bit hacky, the previous solution \
# relied on a run directive with python code to create the {input} definition, \
# but due to the use of singularity this isn't possible, therefore the current bash \
# command substitution is used.
rule GatherHTCVCFs:
	input:
		expand("Outputs/HaplotypeCaller/{{sample}}_{directory}_rawVariants.g.vcf.gz", 
			directory=DIRECTORIES)
	output:
		vcf = temp("Outputs/GatherVCFs/{sample}_GatheredVCFs.g.vcf.gz"),
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		MergeVcfs \
		-O {output.vcf} \
		$(echo ' {input}' | sed 's/ / --INPUT /g') \
		--CREATE_INDEX true"

# Perform joint genotyping
# The cat input substitution can probably be done smoother somehow, it's necessary due to the use of glob()
rule GenotypeGVCFs:
	input:
		fasta = config[refversion]['fasta'],
		flag = "Outputs/MakeContigBeds/flag",
		vcf = "Outputs/GatherVCFs/{sample}_GatheredVCFs.g.vcf.gz",
		contigs = lambda wcs: glob.glob('Outputs/MakeContigBeds/contigs_%s.bed' % wcs.contigs),
	output:
		vcf = temp("Outputs/GenotypeGVCFs/{sample}_{contigs}_genotypes.g.vcf.gz"),
	shell:
		"gatk --java-options '-Xmx3375M -Djava.io.tempdir=$(pwd)/tmp' \
		GenotypeGVCFs \
		-V {input.vcf} \
		-O {output.vcf} \
		-R {input.fasta} \
		$(cat {input.contigs}) \
		--tmp-dir Outputs/GenotypeGVCFs"

# Gather output files from GenotypeGVCFs
# The input substitution with the sed pipe is a bit hacky, the previous solution \
# relied on a run directive with python code to create the {input} definition, \
# but due to the use of singularity this isn't possible, therefore the current bash \
# command substitution is used.
rule GatherGenotypeGVCFs:
	input:
		expand("Outputs/GenotypeGVCFs/{{sample}}_{contigs}_genotypes.g.vcf.gz", 
			contigs=CONTIGS)
	output:
		vcf = temp("Outputs/GatherVCFs2/{sample}_GatheredVCFs2.g.vcf.gz"),
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		MergeVcfs \
		-O {output.vcf} \
		$(echo ' {input}' | sed 's/ / --INPUT /g') \
		--CREATE_INDEX true"

# Build a recalibration model to score variant quality for filtering purposes
rule VariantRecalibratorSNP:
	input:
		omni = config[refversion]['omni'],
		fasta = config[refversion]['fasta'],
		dbsnp = config[refversion]['dbsnp'],
		v1000g = config[refversion]['v1000g'],
		hapmap = config[refversion]['hapmap'],
		vcf = "Outputs/GatherVCFs2/{sample}_GatheredVCFs2.g.vcf.gz",
	output:
		recal = temp("Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.recal"),
		tranches = temp("Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.tranches"),
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		VariantRecalibrator \
		--mode SNP \
		-V {input.vcf} \
		-R {input.fasta} \
		--max-gaussians 4 \
		--output {output.recal} \
		-tranche 97.0 -tranche 90.0 \
		--tranches-file {output.tranches} \
		-an QD -an MQ -an DP -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
		-tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 \
		-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 \
		--resource:omni,known=false,training=true,truth=true,prior=12.0 {input.omni} \
		--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} \
		--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} \
		--resource:v1000G,known=false,training=true,truth=false,prior=10.0 {input.v1000g} \
		--tmp-dir Outputs/VariantRecalibratorSNP"

# Build a recalibration model to score variant quality for filtering purposes
rule VariantRecalibratorINDEL:
	input:
		fasta = config[refversion]['fasta'],
		dbsnp = config[refversion]['dbsnp'],
		mills = config[refversion]['mills'],
		vcf = "Outputs/GatherVCFs2/{sample}_GatheredVCFs2.g.vcf.gz",
	output:
		recal = "Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.recal",
		tranches = "Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.tranches",
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
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
		vcf = "Outputs/GatherVCFs2/{sample}_GatheredVCFs2.g.vcf.gz",
		recal = "Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.recal",
		tranches = "Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.tranches"
	output:
		vcf = "Outputs/ApplyVqsrSnp/{sample}_SnpApplyVQSR.g.vcf.gz",
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		ApplyVQSR \
		--mode SNP \
		-V {input.vcf} \
		-O {output.vcf} \
		-R {input.fasta} \
		-ts-filter-level 99.6 \
		-recal-file {input.recal} \
		-tranches-file {input.tranches} \
		--tmp-dir Outputs/ApplyVqsrSnp"

# Apply a score cutoff to filter variants based on a recalibration table
rule ApplyVqsrIndel:
	input:
		fasta = config[refversion]['fasta'],
		vcf = "Outputs/GatherVCFs2/{sample}_GatheredVCFs2.g.vcf.gz",
		recal = "Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.recal",
		tranches = "Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.tranches"
	output:
		vcf = "Outputs/ApplyVqsrIndel/{sample}_IndelApplyVQSR.g.vcf.gz",
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		ApplyVQSR \
		--mode INDEL \
		-V {input.vcf} \
		-O {output.vcf} \
		-R {input.fasta} \
		-ts-filter-level 95.0 \
		-recal-file {input.recal} \
		-tranches-file {input.tranches} \
		--tmp-dir Outputs/ApplyVqsrIndel"
