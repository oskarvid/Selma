## The following instructions are not up to date! - Basic setup instructions for local execution  
The Snakefile is adapted to run inside a singularity container built with singularity version 2.5.1, build the image with the [BuildSingularityImage.sh](https://github.com/elixir-no-nels/Selma/blob/master/singularity/BuildSingularityImage.sh) script. It uses the `oskarv/snakemake-germline-tools` docker image to build the singularity image. It is not possible to download a prebuilt singularity image.  

Edit `./start-workflow.sh` and change the file path for `REFERENCES` to the filepath where the reference files are stored. In order to support more than one reference file version the file path needs to be of this format: `/storage/workflows/references/`, and in the `references` folder there should be a directory called `hg38` and a directory called `b37`. The [config file](https://github.com/elixir-no-nels/Selma/blob/master/config.yaml) has file paths for. Use these links to download [b37](https://console.cloud.google.com/storage/browser/gatk-legacy-bundles/b37) and [hg38](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/) reference files from Googles cloud, a Google account is necessary to access the links.  
There is also a public ftp server here: [ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle)  
There is no password.  
Use this command to automatically download the hg38 folder:  
`wget -m ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38`  

The wgs interval list for b37 isn't available on the ftp so use the google link above to download that as well. The hg38 bundle on the ftp is complete.

Index the fasta file with bwa before running the pipeline like so: `bwa index -a bwtsw hg38.fasta`.  
