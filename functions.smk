# These functions are used by the snakefile to select the input files and create \
# the readgroup
# It's probably possible to simplify the functions somehow

# get_fastq1 and 2 select the correct pair of fastq files from the samples.tsv file
def get_fastq1(wildcards):
	return samples.loc[(wildcards.flowcell, wildcards.sample, wildcards.lane), ["R1"]].dropna()

def get_fastq2(wildcards):
	return samples.loc[(wildcards.flowcell, wildcards.sample, wildcards.lane), ["R2"]].dropna()

# This function creates the complete readgroup for bwa mem
def get_BwaRG(wc):
	library = samples.loc[(wc.flowcell, wc.sample, wc.lane), ["library"]].unique()[0]
	rgs = r'@RG:\tID:%(flowcell)s.%(lane)s\tSM:%(sample)s\tLB:%(library)s\tPU:%(flowcell)s.%(lane)s.%(sample)s\tPL:Illumina' \
		% {'flowcell': wc.flowcell,
			'lane': wc.lane,
			'sample': wc.sample,
			'library': library}
	return rgs

# get_FQ{SM,FC,LN,LIB} select each respective column and row to construct the \
# readgroup for FastqToSam
def get_FQSM(wc):
	library = samples.loc[(wc.flowcell, wc.sample, wc.lane), ["library"]].unique()[0]
	rgs = r'%(sample)s' \
		% {'sample': wc.sample}
	return rgs

def get_FQFC(wc):
	library = samples.loc[(wc.flowcell, wc.sample, wc.lane), ["library"]].unique()[0]
	rgs = r'%(flowcell)s' \
		% {'flowcell': wc.flowcell}
	return rgs

def get_FQLN(wc):
	library = samples.loc[(wc.flowcell, wc.sample, wc.lane), ["library"]].unique()[0]
	rgs = r'%(lane)s' \
		% {'lane': wc.lane}
	return rgs

def get_FQLIB(wc):
	library = samples.loc[(wc.flowcell, wc.sample, wc.lane), ["library"]].unique()[0]
	rgs = r'%(library)s' \
		% {'library': library}
	return rgs
