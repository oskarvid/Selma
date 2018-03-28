# snakemake_germline
Germline Variant Calling Pipeline built in Snakemake

# This repository is still in progress and can not be regarded as a finished pipeline, although it should run 100% with a proper set of input files.

The Snakefile is adapted to run inside a docker container that I prepared for the pipeline. Either download it manually with `docker pull oskarv/snakemake-germline-tools`
or run the start script and it'll get downloaded automatically if it isn't already downloaded. 

# Instructions  
Edit scripts/start-pipeline.sh and change the file path for REFERENCES to the filepath where you keep your reference files. The default reference files 
are the hg38 reference files from the Broad Institute, they host them at their public ftp server here: 
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle
There is no password. You can automatically download the hg38 folder with this command:
wget -m ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38

If you haven't indexed the fasta file with bwa you must do that before you run the pipeline.  

Run the pipeline with "sh scripts/start-pipeline.sh" to run it in the docker container with snakemake, bwa, samtools and gatk installed.
Or run it with singularity with "snakemake -j --use-singularity", it will use another docker container that I normally use for my wdl germline pipeline.
Or run it locally and make sure all tools are installed already.
