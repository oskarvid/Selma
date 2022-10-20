# Germline Variant Calling Pipeline built in Snakemake  

## Introduction
**Docker**  
The Snakefile is adapted to run inside a docker container that I prepared for the pipeline. Either download it manually with `docker pull oskarv/snakemake-germline-tools`
or run the start script and it'll get downloaded automatically if it isn't already downloaded. Alternatively build it manually with the Dockerfile.  

**Conda**  
Selma can also run in a virtualized environment using conda as such:  
```
cd Selma
conda env create -n selma --file conda/env.yaml
snakemake -j --config version=hg38/or/b37 interval=/path/to/interval_list
```
After it's done building the environment you can run it straight on the commandline without the need for docker or the helper scripts. You still need to edit the `workspace/config.yaml` and `workspace/samples.tsv` files with correct paths, samples etc. 

![Graphical visualization of the pipeline steps](https://github.com/oskarvid/snakemake_germline/blob/master/dag.png)

# Instructions  
Edit `scripts/start-pipeline.sh` and change the file path for `REFERENCES` to the filepath where you keep your reference files. The default reference files 
are the hg38 reference files from the Broad Institute, they host them at their public ftp server here:  
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle  
There is no password. You can automatically download the hg38 folder with this command:  
`wget -m ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38`  

If you haven't indexed the fasta file with bwa you must do that before you run the pipeline.  

Run the pipeline with `sh scripts/start-pipeline.sh` to run it in the oskarv/snakemake-germline-tools docker container with snakemake, bwa, samtools and 
gatk installed.  
You can also run it locally with `snakemake -j`, just edit the relevant paths in the script and make sure all tools are installed locally.  
Singularity is not supported due to the use of "run:", the Singularity directive is only allowed with shell, script or wrapper directives.

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
<h1 align="center">
  <br>
  <a href="https://github.com/elixir-no-nels/Selma"><img src="https://raw.githubusercontent.com/elixir-no-nels/Selma/master/.selma.svg?sanitize=true" alt="Selma" width="300"></a>
</h1>

[![Travis Build Status](https://api.travis-ci.com/elixir-no-nels/Selma.svg?branch=milestone2)](https://travis-ci.com/elixir-no-nels/Selma?branch=milestone2)  

## About Selma
Selma is a whole genome (germline) variant calling workflow developed at the University of Bergen based on the GATK suite of tools. The guiding philosophy behind it is that it should be easy to setup, easy to use and that it utilizes system resources efficiently. This is achieved by adopting a user centric frame of mind that aims to simplify complex tasks without sacrificing functionality. The workflow itself is based on [Snakemake](https://snakemake.readthedocs.io/en/stable/) and all dependencies are handled by using [Docker](https://www.docker.com/) and [Singularity](https://singularity.lbl.gov/) container technology. The current intended platform is [TSD](https://www.uio.no/tjenester/it/forskning/sensitiv/) but support for [HUNT-cloud](https://www.ntnu.edu/mh/huntcloud) as well as local execution is planned for future releases.  
Selma is named after the mythical Norwegian sea serpent that supposedly lives in [Lake Seljord](https://en.wikipedia.org/wiki/Selma_(lake_monster))

The workflow development is currently supported by [Elixir2](https://elixir-europe.org/), [NorSeq](https://www.norseq.org/) and [Tryggve2](https://neic.no/tryggve/), and in the past also by [BioBank Norway](https://bbmri.no/). 

#### Graphical visualization of the workflow steps
![Graphical visualization of the workflow steps](https://raw.githubusercontent.com/elixir-no-nels/Selma/master/.simplifieddag.png)
###### This is a simplified graph portraying the key steps that the workflow goes through, [this](https://raw.githubusercontent.com/elixir-no-nels/Selma/master/.completedag.png) is a complete overview including every single step. The steps that have been left out only perform "administrative" functions and don't add to the data analysis per se.

### Documentation
* [TSD-instructions](https://github.com/elixir-no-nels/Selma/blob/master/docs/TSD-instructions.md)  
* [Instructions for local use](https://github.com/elixir-no-nels/Selma/blob/master/docs/instructions-for-local-use.md)  
* [Developer-instructions](https://github.com/elixir-no-nels/Selma/blob/master/docs/developer-instructions.md)  

### Tools
[bwa](http://bio-bwa.sourceforge.net/bwa.shtml) version 0.7.15-2+deb9u1 - Maps fastq file to reference genome  
[samtools](http://www.htslib.org/doc/samtools.html) version 1.3.1-3 - bwa pipes its output to samtools to make a bam output file  
The following tools are all [gatk](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/) version 4.1.2.0  
[SplitIntervals](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_SplitIntervals.php) - Splits interval list for scatter gather parallelization  
[FastqToSam](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/picard_sam_FastqToSam.php) - Converts fastq files to unmapped bam files  
[MergeBamAlignment](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/picard_sam_MergeBamAlignment.php) - Merge aligned BAM file from bwa with the unmapped BAM file from FastqToSam  
[MarkDuplicates](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/picard_sam_markduplicates_MarkDuplicates.php) - Identifies duplicate reads  
[BaseRecalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php) - Generates recalibration table for Base Quality Score Recalibration  
[GatherBQSRReports](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_bqsr_GatherBQSRReports.php) - Gather base recalibration files from BaseRecalibrator  
[ApplyBQSR](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php) - Apply base recalibration from BaseRecalibrator  
[GatherBamFiles](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/picard_sam_GatherBamFiles.php) - Concatenate efficiently BAM files from ApplyBQSR  
[HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php) - Call germline SNPs and indels via local re-assembly of haplotypes  
[GenotypeGVCFs](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_GenotypeGVCFs.php) - Perform genotyping on one pre-called sample from HaplotypeCaller  
[VariantRecalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_vqsr_VariantRecalibrator.php) - Build a recalibration model to score variant quality for filtering purposes  
[ApplyVQSR](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_vqsr_ApplyVQSR.php) -  Apply a score cutoff to filter variants based on a recalibration table


### Credits  
**Supervisor**  
[Kjell Petersen](mailto:kjell.petersen@uib.no)

**Main developer**  
[Oskar Vidarsson](mailto:oskar.vidarsson@uib.no)
