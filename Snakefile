# Import modules for tsv file handling and globbing of certain output files
import pandas as pd
import glob

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
CONTIGCOUNT = config[refversion]['contigfiles']
CONTIGS = range(1, CONTIGCOUNT)
GATHERCONTIGS = range(1, (( CONTIGCOUNT + 1 )))

SCATTERCOUNT = config['scattercount']
DIRECTORIES=[]
for i in range(SCATTERCOUNT):
	DIRECTORIES.append(str(i).zfill(4))

# Create variables to select sample names, lane numbers and flowcell names from the sample.tsv file
samples = pd.read_csv(config["samples"], sep='\t', dtype=str).set_index(["flowcell", "sample", "lane"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels]) # enforce str in index

rule all:
	input:
		expand("Outputs/ApplyVqsrSnp/{sample}_SnpApplyVQSR.g.vcf.gz", 
			sample=samples['sample']),
		expand("Outputs/ApplyVqsrIndel/{sample}_IndelApplyVQSR.g.vcf.gz", 
			sample=samples['sample']),

# This rule creates a tsv file with the contigs grouped into roughly equal combined lengths.
# The script should be improved to create a separate file per group so as to avoid the need
# for the MakeTSVs.sh script which performs this function
rule MakeSequenceGroupings:
	input:
		config[refversion]['dict'],
	output:
		temp("Outputs/MakeContigBeds/sequence_grouping_with_unmapped.tsv"),
	priority:
		30
	shell:
		"python2 scripts/split-bedfile.py {input} Outputs/MakeContigBeds/"

# Split the sequenc_grouping_with_unmapped.tsv file into one file per line
# This script should be made osbolete by performin its function in the split-bedfile.py instead
rule MakeContigBeds:
	input:
		"Outputs/MakeContigBeds/sequence_grouping_with_unmapped.tsv"
	output:
		flag = touch("Outputs/MakeContigBeds/flag"),
	priority:
		30
	shell:
		"bash scripts/MakeTSVs.sh {input} Outputs/MakeContigBeds/"

# Split the interval list for HaplotypeCaller into sub intervals for scatter gather execution
rule MakeIntervalLists:
	input:
		interval = interval,
		fasta = config[refversion]['fasta'],
	output:
		temp("Outputs/MakeIntervalLists/{split}-scattered.interval_list"),
	priority:
		30
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
		temp("Outputs/BwaMem/{sample}_{lane}_{flowcell}.mapped.bam"),
	benchmark:
		"Outputs/benchmarks/{sample}_{lane}_{flowcell}.bwa.tsv",
	threads:
		15
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
		fastq1 = get_fastq1,
		fastq2 = get_fastq2,
		fasta = config[refversion]['fasta'],
	output:
		bam = temp("Outputs/FastqToSam/{sample}_{lane}_{flowcell}.unmapped.bam"),
		tmp = directory(temp("Outputs/FastqToSam/{sample}_{lane}_{flowcell}.tmp")),
	benchmark:
		"Outputs/benchmarks/{sample}_{lane}_{flowcell}.FastqToSam.tsv",
	params:
		lane = get_FQLN,
		sample = get_FQSM,
		flowcell = get_FQFC,
		library = get_FQLIB,
	priority: 2
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

# Due to the usage of the MarkDupCheckpoint below it is not possible to flag the bam
# file from MergeBamAlignment as temp(), this seems like a bug
rule MergeBamAlignment:
	input:
		fasta = config[refversion]['fasta'],
		mapped = "Outputs/BwaMem/{sample}_{lane}_{flowcell}.mapped.bam",
		unmapped = "Outputs/FastqToSam/{sample}_{lane}_{flowcell}.unmapped.bam",
	output:
		bam = "Outputs/MergeBamAlignment/{sample}_{lane}_{flowcell}.merged.bam",
		tmp = directory(temp("Outputs/MergeBamAlignment/{sample}_{lane}_{flowcell}.tmp")),
	benchmark:
		"Outputs/benchmarks/{sample}_{lane}_{flowcell}.MergeBamAlignments.tsv",
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
		--PROGRAM_GROUP_COMMAND_LINE 'bwa mem -t 15 -R -M Input1 Input2 > output.sam' \
		--TMP_DIR {output.tmp}"

# It would have been convenient to make MergeBamAlignment the checkpoint
# but it wouldn't recognize it for some reason, hence the checkpoint rule below.
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
		bai = temp("Outputs/MarkDuplicates/{sample}_markedDuplicates.bai"),
		metrics = "Outputs/MarkDuplicates/{sample}_markedDuplicates.metrics",
	benchmark:
		"Outputs/benchmarks/{sample}.MarkDup.tsv",
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
		bam = "Outputs/MarkDuplicates/{sample}_markedDuplicates.bam",
		bai = "Outputs/MarkDuplicates/{sample}_markedDuplicates.bai",
		contigs = lambda wcs: glob.glob('Outputs/MakeContigBeds/contigs_%s.bed' % wcs.contigs),
	threads:
		1
	output:
		grp = temp("Outputs/BaseRecalibrator/{sample}_BQSR_{contigs}.grp"),
	benchmark:
		"Outputs/benchmarks/{sample}_{contigs}.BaseRecalibrator.tsv",
	shell:
		"gatk --java-options '-Xmx3500M -Djava.io.tempdir=$(pwd)/tmp' \
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
			directory=CONTIGS),
	output:
		temp("Outputs/GatherBQSR/{sample}_GatheredBQSR.grp")
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
		bam = temp("Outputs/ApplyBQSR/{sample}_{con}_recalibrated.bam"),
		bai = temp("Outputs/ApplyBQSR/{sample}_{con}_recalibrated.bai"),
	benchmark:
		"Outputs/benchmarks/{sample}_{con}.ApplyBQSR.tsv",
	shell:
		"gatk --java-options '-Xmx3500M -Djava.io.tempdir=$(pwd)/tmp' \
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
		bai = expand("Outputs/ApplyBQSR/{{sample}}_{con}_recalibrated.bai",
			con=GATHERCONTIGS),
		bam = expand("Outputs/ApplyBQSR/{{sample}}_{con}_recalibrated.bam",
			con=GATHERCONTIGS),
	output:
		bam = "Outputs/GatherBamFiles/{sample}_GatheredABQSRFiles.bam",
		bai = "Outputs/GatherBamFiles/{sample}_GatheredABQSRFiles.bai",
	benchmark:
		"Outputs/benchmarks/{sample}.GatheredBamFiles.tsv",
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		GatherBamFiles \
		-O {output.bam} \
		$(echo ' {input.bam}' | sed 's/ / --INPUT /g') \
		--CREATE_INDEX true"

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
		vcf = temp("Outputs/HaplotypeCaller/{sample}_{split}_rawVariants.g.vcf.gz"),
		tbi = temp("Outputs/HaplotypeCaller/{sample}_{split}_rawVariants.g.vcf.gz.tbi"),
	benchmark:
		"Outputs/benchmarks/{sample}_{split}.HaplotypeCaller.tsv",
	threads:
		2
	shell:
		"gatk --java-options '-Xmx3500M -Djava.io.tempdir=$(pwd)/tmp' \
		HaplotypeCaller \
		-ERC GVCF \
		-I {input.bam} \
		-O {output.vcf} \
		-R {input.fasta} \
		-L {input.intervals} \
		--tmp-dir Outputs/HaplotypeCaller"

# The output files from HaplotypeCaller need to be gathered into one file
# This tool is currently unable to index the output file, the indexing tool handles that
# The tbi input files are only used as input so that they are deleted at the right time using the temp directive
rule GatherHTCVCFs:
	input:
		vcf = expand("Outputs/HaplotypeCaller/{{sample}}_{split}_rawVariants.g.vcf.gz", 
			split=DIRECTORIES),
		tbi = expand("Outputs/HaplotypeCaller/{{sample}}_{split}_rawVariants.g.vcf.gz.tbi",
			split=DIRECTORIES),
	output:
		vcf = temp("Outputs/GatherHTCVCFs/{sample}_GatheredHTCVCFs.g.vcf.gz"),
	benchmark:
		"Outputs/benchmarks/{sample}.GatherHTCVCFs.tsv",
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
		temp("Outputs/GatherHTCVCFs/{sample}_GatheredHTCVCFs.g.vcf.gz.tbi"),
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		IndexFeatureFile \
		-F {input} \
		-O {output}"
		
# Perform joint genotyping
rule GenotypeGVCFs:
	input:
		fasta = config[refversion]['fasta'],
		flag = "Outputs/MakeContigBeds/placeholder",
		vcf = "Outputs/GatherHTCVCFs/{sample}_GatheredHTCVCFs.g.vcf.gz",
		tbi = "Outputs/GatherHTCVCFs/{sample}_GatheredHTCVCFs.g.vcf.gz.tbi",
		contigs = lambda wcs: glob.glob('Outputs/MakeContigBeds/contigs_%s.bed' % wcs.contigs),
	output:
		vcf = temp("Outputs/GenotypeGVCFs/{sample}_{contigs}_genotypes.g.vcf.gz"),
		tbi = temp("Outputs/GenotypeGVCFs/{sample}_{contigs}_genotypes.g.vcf.gz.tbi"),
	benchmark:
		"Outputs/benchmarks/{sample}_{contigs}.GenotypeGVCFs.tsv",
	shell:
		"gatk --java-options '-Xmx3500M -Djava.io.tempdir=$(pwd)/tmp' \
		GenotypeGVCFs \
		-V {input.vcf} \
		-O {output.vcf} \
		-R {input.fasta} \
		$(cat {input.contigs}) \
		--tmp-dir Outputs/GenotypeGVCFs"

# The output files from GenotypeGVCFs need to be gathered into one file
# This tool is currently unable to index the output file, the indexing tool handles that
# The tbi input files are only used as input so that they are deleted at the right time using the temp directive
rule GatherGenotypeGVCFs:
	input:
		tbi = expand("Outputs/GenotypeGVCFs/{{sample}}_{contigs}_genotypes.g.vcf.gz.tbi",
			contigs=CONTIGS),
		vcf = expand("Outputs/GenotypeGVCFs/{{sample}}_{contigs}_genotypes.g.vcf.gz", 
			contigs=CONTIGS),
	output:
		vcf = temp("Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz"),
	benchmark:
		"Outputs/benchmarks/{sample}.GatherGenotypeGVCFs.tsv",
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
		temp("Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz.tbi"),
	benchmark:
		"Outputs/benchmarks/{sample}.IndexGatheredGenotypeGVCFs.tsv",
	shell:
		"gatk --java-options -Djava.io.tempdir=$(pwd)/tmp \
		IndexFeatureFile \
		-F {input} \
		-O {output}"

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
		recal = temp("Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.recal"),
		idx = temp("Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.recal.idx"),
		tranches = temp("Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.tranches"),
	benchmark:
		"Outputs/benchmarks/{sample}.VariantRecalibratorSNP.tsv",
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
		vcf = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz",
		tbi = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz.tbi",
	output:
		recal = temp("Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.recal"),
		idx = temp("Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.recal.idx"),
		tranches = temp("Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.tranches"),
	benchmark:
		"Outputs/benchmarks/{sample}.VariantRecalibratorINDEL.tsv",
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
		recal = "Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.recal",
		idx = "Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.recal.idx",
		vcf = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz",
		tranches = "Outputs/VariantRecalibratorSNP/{sample}_SnpVQSR.tranches",
		tbi = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz.tbi",
	output:
		vcf = "Outputs/ApplyVqsrSnp/{sample}_SnpApplyVQSR.g.vcf.gz",
	benchmark:
		"Outputs/benchmarks/{sample}.ApplyVqsrSnp.tsv",
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
		vcf = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz",
		recal = "Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.recal",
		idx = "Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.recal.idx",
		tbi = "Outputs/GatherGenotypeGVCFs/{sample}_GatheredGVCFs.g.vcf.gz.tbi",
		tranches = "Outputs/VariantRecalibratorINDEL/{sample}_IndelVQSR.tranches",
	output:
		vcf = "Outputs/ApplyVqsrIndel/{sample}_IndelApplyVQSR.g.vcf.gz",
	benchmark:
		"Outputs/benchmarks/{sample}.ApplyVqsrIndel.tsv",
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
