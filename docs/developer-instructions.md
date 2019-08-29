## Setting up the Selma workflow on TSD from scratch
These instructions describe how to install a clean copy of Selma in the `/tsd/shared/bioinformatics/workflows` directory. The assumption is that everything needs to be done completely from scratch.  
You will need to go through these steps:  
**[Clone this repository and build the Singularity image locally](file:///home/oskar/01-workspace/04-pipelines/Selma/docs/developer-instructions.html#cloning-this-repository-and-building-the-singularity-image)**  
**[Download and index reference files](file:///home/oskar/01-workspace/04-pipelines/Selma/docs/developer-instructions.html#downloading-and-indexing-reference-files)**  

### Cloning this repository and building the Singularity image
**Install Singularity**  
The current version of the workflow has been tested using Singularity version 2.6.1, it may or may not function as expected using any other version. You will find Singularity installation instructions here: [https://www.sylabs.io/docs/](https://www.sylabs.io/docs/)  

**Clone the repository and build the Singularity image**  
Now run this:
```
git clone https://github.com/elixir-no-nels/Selma
cd Selma
./singularity/BuildSingularityImage.sh
```

### Copy the Selma directory to TSD
When the Singularity image has been built, make a `.tar.gz` archive of the Selma directory like so:
```bash
cd ../
tar -zcvf Selma-$(date "+%F").tar.gz Selma/
```
Then use your [preferred method](https://www.uio.no/english/services/it/research/sensitive-data/use-tsd/import-export/) to upload the `Selma-2019-01-01.tar.gz` archive to TSD.  
Once it has been copied to TSD, run the following steps:  
```bash
mv /tsd/shared/bioinformatics/workflows/Selma /tsd/shared/bioinformatics/workflows/Selma_backup-$(date "+%F")
tar -zxvf Selma-2019-01-01.tar.gz -C /tsd/shared/bioinformatics/workflows/
```
This will only work if you have write permission to the `/tsd/shared/bioinformatics/` directory.
### Downloading and indexing reference files  
There is a public ftp server with the reference files here: [ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle)  
To download the supported reference files you need to run these commands.  

```bash
mkdir hg38 && cd hg38 && wget -i /directory/with/Selma/scripts/hg38-reference-dl-list
mkdir b37 && cd b37 && wget -i /directory/with/Selma/scripts/b37-reference-dl-list
```
N.B: One or more files may fail to download if more than 25 users are connected to the ftp at the same time.

Unfortunately the b37 wgs interval file isn't available on the ftp site so you need to download it from the Google cloud bucket with this:
```bash
cd /directory/with/downloaded/b37/files/
wget https://storage.googleapis.com/gatk-legacy-bundles/b37/wgs_calling_regions.v1.interval_list
```

There is a second method to download the reference files, I will document it here, but the files aren't the same as on the ftp server so this is only here just in case. Using these files would require fixing the paths in the workflow too so don't use them unless you _really_ need to.  
<details><summary>Click here for the instructions</summary>
<p>
The files are hosted in Google buckets here: https://console.cloud.google.com/storage/browser/gatk-legacy-bundles/b37/  
And here: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/  
A Google account is necessary to browse the file buckets, but you can still download them with gsutils without having a Google account.  
The current Selma version uses the ftp hosted reference files, and a difference between the two is that e.g `dbsnp_146.vcf` is **not** available in the hg38 version hosted on Google, maybe this will change, or maybe one day Selma will use another dbsnp.vcf version, but until then this why you should download the ftp hosted reference files.  
So, just in case this becomes relevant in the future; To download all reference files in one go you should use `gsutils`. Either install it following the appropriate instructions here: https://cloud.google.com/storage/docs/gsutil_install or use docker like this:  
For b37 do:
```bash
docker run --rm -ti -v $(pwd):/data -w /data google/cloud-sdk:latest gsutil -m cp -r gs://gatk-legacy-bundles/b37 .
```
And for hg38 do:
```bash
docker run --rm -ti -v $(pwd):/data -w /data google/cloud-sdk:latest gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38/v0/ .
mv v0 hg38
```

If you use this docker method you need to run `sudo chown -R username:username b37orhg38` to get ownership of the files.

</p>
</details>

#### Indexing the reference files
There is a script that will unzip and index the fasta file and run tabix to index all vcf files.
Once you have downloaded the files; copy the script from `Selma/utilities/reference-preparation.sh` to the b37 and/or hg38 directories and run it with `./reference-preparation.sh`  
It will run the following steps:  
  * Tabix index on bgzipped vcf files  
  * Samtools faidx on the fasta files  
  * bwa index on the fasta files  

### Copy the reference files to TSD
Make a `tar.gz` archive of the reference file directories before you upload them to TSD.
```bash
cd /path/to/directory/with/hg38orb37
tar -zcvf hg38.tar.gz hg38
tar -zcvf b37.tar.gz b37
```
Then you can upload the reference file archives to TSD with your [preferred method](https://www.uio.no/english/services/it/research/sensitive-data/use-tsd/import-export/).  
Once the reference file archives have been uploaded to TSD, untar them like so:
```bash
cd /tsd/pXX/data/durable/file-import/pXX-member-group
tar -zxvf hg38.tar.gz -C /tsd/shared/bioinformatics/reference-data/
tar -zxvf b37.tar.gz -C /tsd/shared/bioinformatics/reference-data/
rm hg38.tar.gz b37.tar.gz
```

## Hardware requirements and optimizations  
The workflow should be able to run on a single server with 16 thread, 64GB RAM and at least 500GB storage assuming that there are 8 fastq.gz files totalling 51GB with ~30x coverage. When using the test files in the fastq folder it should run on any laptop using 2 threads and 8GB RAM, but preferrably 4 threads and 16GB RAM, the storage requirements apart from the reference files is negligible.  
The run time on my current test machine that has 16 threads and 32 GB RAM has been between 16 hours and 14 minutes to 16 hours and 25 minutes with 8 fastq.gz file pairs totalling ~51GB/30x coverage.  
The execution time on a server with 16 threads and 16 GB RAM is roughly 18 hours and 30 minutes if each scatter gather tool is given 2GB RAM each and using the same input files as above.  
