# Germline Variant Calling Pipeline built in Snakemake  
test
The Snakefile is adapted to run inside a docker container that I prepared for the pipeline. Either download it manually with `docker pull oskarv/snakemake-germline-tools`
or run the start script and it'll get downloaded automatically if it isn't already downloaded. Alternatively build it manually with the Dockerfile.  

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

`<rant>` Compared with my WDL 
based pipeline for germline variant calling, this is 5-6 hours faster. The reason for this speed increase is due to parallelization options that aren't
available in WDL. In WDL you are not able to manually limit a scatter/gather process to loop over each input file for one tool, this causes an inefficiency
for bwa since all input files must run at the same time, as well as all FastqToSam processes, meaning that you must either choose between not overloading 
the system and not parallelize bwa, which would mean that you run e.g 8 input pairs for bwa and FastqToSam until FastqToSam is finished, which takes ~25 
minutes, and thus temporarily creating a system load of 16, but then only use 8 threads once FastqToSam is finished, or temporarily overload the system 
and parallelize bwa with at least 3 threads, since I expect using 2 threads won't actually parallelize anything since one thread is usually used to control
the rest of the threads, meaning you would use one thread to control one thread if only two threads are used to parallelize bwa.  

Thus using 3 threads would
create a system load of 3x8 + 8 until FastqToSam is finished, and then a consistent system load of 3x8 until bwa is finished. On a 16 thread machine this 
is suboptimal, a better solution, which Snakemake enables, is to loop over the input files with bwa, allowing you to use 16 threads per pair which means 
that each pair takes roughly 50 minutes to finish, and then run 8 parallel processes for FastqToSam, which takes roughly 25 minutes. That way you don't 
overload the system and gain in speed.  

After correspondence with The Broad Institute, the organization that develops WDL, their stance is that WDL should
 rather be used on e.g Google cloud, and that the shards should be distributed to a compute node each. This is not always possible, hence this feature is
 sorely needed in WDL since the lack of it causes unecessary inefficiens. To be fair there are optimizations in this pipeline that could be implemented
in my current WDL germline pipeline that should decrease the execution time by at least ~70 minutes. `</rant>`

## Planned features and testing  
I am still learning Snakemake, and so far I am planning to enable the use of a config file to define input paths and variables.  

