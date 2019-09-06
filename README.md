<h1 align="center">
  <br>
  <a href="https://github.com/elixir-no-nels/Selma"><img src="https://github.com/elixir-no-nels/Selma/blob/master/.Selma.svg?sanitize=true" alt="Selma" width="300"></a>
</h1>

## About Selma
Selma is a germline variant calling workflow developed at the University of Bergen. The guiding philosophy behind it is that it should be easy to setup, easy to use and that it utilizes system resources efficiently. This is achieved by adopting a user centric frame of mind that aims to simplify complex tasks without sacrificing functionality. The workflow itself is based on [Snakemake](https://snakemake.readthedocs.io/en/stable/) and all dependencies are handled by using [Docker](https://www.docker.com/) and [Singularity](https://singularity.lbl.gov/) container technology. The current intended platform is [TSD](https://www.uio.no/tjenester/it/forskning/sensitiv/) but support for [HUNT-cloud](https://www.ntnu.edu/mh/huntcloud) as well as local execution is planned for future releases.  
Selma is named after the mythical Norwegian sea serpent that supposedly lives in [Lake Seljord](https://en.wikipedia.org/wiki/Selma_(lake_monster))

The workflow development is funded by [Elixir2](https://elixir-europe.org/), [NorSeq](https://www.norseq.org/) and [Tryggve2](https://neic.no/tryggve/). 

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
