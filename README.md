## The current version of this readme is outdated, incomplete and incorrect, it's only a placeholder until it gets updated for real.
# Germline Variant Calling Pipeline built in Snakemake  

The Snakefile is adapted to run inside a docker container that I prepared for the pipeline. Either download it manually with `docker pull oskarv/snakemake-germline-tools`
or run the start script and it'll get downloaded automatically if it isn't already downloaded. Alternatively build it manually with the Dockerfile.  

![Graphical visualization of the pipeline steps](https://github.com/elixir-no-nels/snakemake_germline/blob/master/dag.png)

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
At the current state the pipeline is highly optimized for use on a single server with 16 thread, 64GB RAM and at least 500GB storage assuming that there are 8 
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
