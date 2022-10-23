# Germline Variant Calling Pipeline built in Snakemake  

<h1 align="center">
  <br>
  <a href="https://github.com/oskarvid/Selma"><img src="https://raw.githubusercontent.com/oskarvid/Selma/master/.selma.svg?sanitize=true" alt="Selma" width="300"></a>
</h1>

[![Travis Build Status](https://api.travis-ci.com/oskarvid/.selma.svg?branch=master)](https://travis-ci.com/oskarvid/Selma?branch=master)  

#### Graphical visualization of the workflow steps
![Graphical visualization of the pipeline steps](https://raw.githubusercontent.com/oskarvid/Selma/master/.simplifieddag.png)

## Instructions
**Raw metal (or Guix shell environment) execution**  
**N.B** No matter which method you use you always need to edit the `workspace/config.yaml` and `workspace/samples.tsv` files with correct paths, samples etc.
**N.B**

The raw metal method is if you have taken care of the dependency installation yourself.  
```
snakemake -j --config version=hg38 interval=/path/to/hg38/interval_list
```
The following instructions guides you how to set up the `Selma` environment using `docker`, `conda` or `guix`.

**Docker execution**  
Either download the image manually with `docker pull oskarv/selma` or run 
```
docker run --rm -ti -v $PWD:/data -w /data selma snakemake -j --config version=hg38 interval=/path/to/hg38/interval_list
```
and it'll get downloaded automatically. Alternatively build it manually with the method described farther down.

**Conda execution**  
Selma can also run in a virtualized environment using conda as such:  
```
conda env create -n selma --file conda/env.yaml
snakemake -j --config version=hg38/or/b37 interval=/path/to/interval_list
```
Alternatively you can use the `--use-conda` flag:  
```
snakemake -j --config version=hg38 interval=/path/to/hg38/interval_list --use-conda
```

**Guix execution**  
[Guix](https://guix.gnu.org/) can be used to install all dependencies apart from `gatk4` as such:  
```
guix shell -m manifest.scm
```
`gatk4` can be installed in a location of your choice with this:  
```
wget --no-check-certificate https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip -O $PWD/gatk4.zip && \
unzip -q $PWD/gatk4.zip -d $PWD/ && \
mv $PWD/gatk*/gatk* $PWD/ && \
rm -r $PWD/gatk*/ $PWD/gatk4.zip && \
export PATH="$PATH:$PWD/" && \
export GATK_LOCAL_JAR=$PWD/gatk-package-4.3.0.0-local.jar
```
Running `gatk` should by now work as expected.

**Guix docker build**  
The docker image was built as follows:  
```
guix pack -f docker -S python=python3 -S /usr/bin/env=bin/env -S /bin=bin samtools gnuplot bwa snakemake bcftools python-pandas openjdk nss-certs bash wget unzip coreutils python2-minimal python-minimal sed python-matplotlib tectonic texlive-base
```
This produces a `tar.gz` named `/gnu/store/a-long-hash-samtools-gnuplot-bwa-snakemake-bcftools-docker-pack.tar.gz` and this file is then loaded as a docker image like this:  
```
docker load < /gnu/store/a-long-hash-samtools-gnuplot-bwa-snakemake-bcftools-docker-pack.tar.gz
```
Then it's on to building the `Selma` docker image to add and configure `gatk4`:  
```
docker build -t oskarv/selma .
```
And now you can run `Selma` with docker!
```
docker run --rm -ti -v $PWD:/data -w /data selma snakemake -j --config version=hg38 interval=/path/to/hg38/interval_list
```

**Reference files**
The default reference files are the hg38 reference files from the Broad Institute, they host them at their public ftp server here:  
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle  
There is no password. You can automatically download the hg38 folder with this command:  
`wget -m ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38`  

If you haven't indexed the fasta file with bwa you must do that before you run the pipeline.  

## Hardware requirements and optimizations  
At the current state the pipeline is highly optimized for use on a single server with 16 threads, 64GB RAM and at least 500GB storage assuming that there are 8 
fastq.gz files totalling 51GB with ~30x coverage. But when using the test files 
in the fastq folder it should run on any laptop using 2 threads and 8GB RAM, but 
preferrably 4 threads and 16GB RAM, the storage requirements apart from the 
reference files is negligible.  
The run time on my current test machine that has 16 threads and 64 GB RAM has 
been between 16 hours and 14 minutes to 16 hours and 25 minutes with 8 fastq.gz 
file pairs totalling ~51GB/30x coverage.  
The execution time on a server with 16 threads and 16 GB RAM is roughly 18 hours 
and 30 minutes if each scatter gather tool is given 2GB RAM each and using the 
same input files as above.  

=======

## About Selma
Selma is a whole genome germline variant calling workflow initially developed at the University of Bergen heavily inspired by the GATK best practices workflow. The guiding philosophy behind it is that it should be easy to setup, easy to use and that it utilizes system resources efficiently. The workflow is based on [Snakemake](https://snakemake.readthedocs.io/en/stable/) and supports [Conda](https://anaconda.org/), [Guix](https://guix.gnu.org/), [Docker](https://www.docker.com/) and (soon to be tested) [Singularity](https://singularity.lbl.gov/) execution modes.  
Selma is named after the mythical Norwegian sea serpent that supposedly lives in [Lake Seljord](https://en.wikipedia.org/wiki/Selma_(lake_monster))


###### This is a simplified graph portraying the key steps that the workflow goes through, [this](https://raw.githubusercontent.com/elixir-no-nels/Selma/master/.completedag.png) is a complete overview including every single step. The steps that have been left out only perform "administrative" functions and don't add to the data analysis per se.

### Tools
[bwa](http://bio-bwa.sourceforge.net/bwa.shtml) version 0.7.17 - Maps fastq file to reference genome  
[samtools](http://www.htslib.org/doc/samtools.html) version 1.14 - bwa pipes its output to samtools to make a bam output file  
The following tools are all [gatk](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.3.0.0/) version 4.3.0.0  
[SplitIntervals](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.3.0.0/org_broadinstitute_hellbender_tools_walkers_SplitIntervals.php) - Splits interval list for scatter gather parallelization  
[FastqToSam](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.3.0.0/picard_sam_FastqToSam.php) - Converts fastq files to unmapped bam files  
[MergeBamAlignment](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.3.0.0/picard_sam_MergeBamAlignment.php) - Merge aligned BAM file from bwa with the unmapped BAM file from FastqToSam  
[MarkDuplicates](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.3.0.0/picard_sam_markduplicates_MarkDuplicates.php) - Identifies duplicate reads  
[BaseRecalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.3.0.0/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php) - Generates recalibration table for Base Quality Score Recalibration  
[GatherBQSRReports](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.3.0.0/org_broadinstitute_hellbender_tools_walkers_bqsr_GatherBQSRReports.php) - Gather base recalibration files from BaseRecalibrator  
[ApplyBQSR](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.3.0.0/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php) - Apply base recalibration from BaseRecalibrator  
[GatherBamFiles](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.3.0.0/picard_sam_GatherBamFiles.php) - Concatenate efficiently BAM files from ApplyBQSR  
[HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.3.0.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php) - Call germline SNPs and indels via local re-assembly of haplotypes  
[GenotypeGVCFs](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.3.0.0/org_broadinstitute_hellbender_tools_walkers_GenotypeGVCFs.php) - Perform genotyping on one pre-called sample from HaplotypeCaller  
[VariantRecalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.3.0.0/org_broadinstitute_hellbender_tools_walkers_vqsr_VariantRecalibrator.php) - Build a recalibration model to score variant quality for filtering purposes  
[ApplyVQSR](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.3.0.0/org_broadinstitute_hellbender_tools_walkers_vqsr_ApplyVQSR.php) -  Apply a score cutoff to filter variants based on a recalibration table
