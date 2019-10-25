# Running the workflow on TSD
These instructions have been verified to work with singularity version 2.6.1 on the Colossus 3.0 slurm cluster. Any other slurm setups and/or singularity versions may or may not work.

## Short instructions  
Assuming that Selma is already installed, begin by running `cd /path/to/Selma/directory/`

Then create a _tab separated file_ using the the header below and add your sample information in a new row below it:  
```bash
flowcell	sample	library	lane	R1	R2
``` 
Or use [this](https:/raw.githubusercontent.com/elixir-no-nels/Selma/master/samples.tsv) as a template.  
Populate the columns with appropriate information, then save the file and name it `my-samples.tsv` or something suitable. Remember to tab separate the columns.  
Assuming you already have the input files ready, and that the output directory exists, you can now start the workflow as such:  
```bash 
./scripts/start-workflow.sh -i /tsd/p172ncspm/data/durable/input-data/ -t /tsd/p172ncspm/data/durable/input-data/my-samples.tsv -o /tsd/p172ncspm/data/durable/Selma-outputs -r hg38
```
This will use hg38 reference files, you can also use b37 reference files.

<details><summary>Optimizing the execution of Selma</summary>
<p>

## The default mode
Out of the box Selma will run all samples on one execute node, this is acceptable if time is no issue or if you have only one or a few samples, but it is very slow for multiple samples. If you have multiple samples to run the default mode makes you wait for _all samples_ to finish running before getting your hands on the output files.

## One sample per node
To optimize the execution you can use the `-s` flag to make Selma run one sample per execute node. This means that if you have 200 samples to run, each sample is put in the slurm queue and will start as soon as there is a slot free in the slurm queue. The output files are transferred back to the final storage as soon as the job on the node is finished. This is the fastest mode Selma offers.

## Select subset of samples from "master" sample.tsv file
Let's assume you have a `big-collection-of-samples.tsv` file that has 200 rows with 13 different samples in it, you start Selma with this tsv file as input and unfortunately there was an issue with 5 of the samples. You now want to rerun these 5 samples and the question is what is the best way of doing this? Since you are smart and don't want to spend lots of time copy/pasting a bunch of rows from one file to another, you simply create a new file that _only contains the sample names_ in one row each like this:
```markdown
MELANO1
ENDOMET37
HEPATO8
LYMPHO15
SARCO3
```
Now name the new tsv file something clever like `5samples.tsv` and use this _in addition to_ your master `big-collection-of-samples.tsv` file like this: `./scripts/start-workflow.sh -t big-collection-of-samples.tsv -o /path/to/output/disk ... -e 5samples.tsv`.  
Hint: Also use the `-s` flag to run each sample on a separate node for optimal execution times.

</p>
</details>

## Long instructions  

### Making a clean copy of Selma and the reference files with the setup script
Before you can run Selma for the first time you need to make a clean copy first. The suggested method is to do it with the setup script that is located in `/tsd/shared/bioinformatics/workflows/Selma/utilities/Selma-setup.sh`  
It needs the directory path to where you want to put your own installation of Selma and another path for where to store the reference directories. Let's assume you want to put Selma in `/cluster/projects/p172/UiO-Cancer/` and the reference files in `/cluster/projects/p172/Selma-references/`, simply run the following command:
```bash
/tsd/shared/bioinformatics/workflows/Selma/utilities/Selma-setup.sh -s /cluster/projects/p172/UiO-Cancer/ -b /cluster/projects/p172/Selma-references/ -g /cluster/projects/p172/Selma-references/
```
This will copy the b37 (-b) and hg38 (-g) to the `/cluster/projects/p172/Selma-references/` directory.

### Configuring directory paths
#### Reference files
You need to edit the `settings/settings.conf` file with the correct directory path to the reference file directories as well as the staging directory. You need to use a directory that is readable by Colossus, a suggestion is `/cluster/projects/p172` because any directory on the `/cluster/projects/pXX` disk is readable by Colossus. The hg38 and b37 reference file directories need to be located in the same directory like this:
```bash
/cluster/projects/p172/Selma-references/
├── b37
│   └── files
└── hg38
    └── files
```
In this case you would put `REFERENCES=/cluster/projects/p172/Selma-references/` in the `settings/settings.conf` file.

#### File staging directory
Next up is setting the file staging directory in the `settings/settings.conf` file. This is where Selma will do all the preparation steps before starting the actual workflow on Colossus, and this is also where the output files from the finished Colossus data analysis will end up temporarily before being sent to the final storage directory that you define with the `-o` option when you start the workflow. The directory needs to be on a disk that is writeable by Colossus, so using something like `FILESTAGING=/cluster/projects/p172/Selma-staging` is a suggestion, run `mkdir /cluster/projects/p172/Selma-staging` to create it.  

#### Editing the sbatch file
Now you need to edit the `scripts/RunOnNode.sbatch` file and change the `#SBATCH --account=p172` line and put your slurm account name there.

Let's continue by using a thought experiment to understand how to supply the workflow with correct options.

### Locating your input files
Your input data in this thought experiment is located in `/tsd/p172ncspm/data/durable/input-data/`, this directory has two files and one directory that also contains two files like this:  
```
/tsd/p172ncspm/data/durable/input-data/
├── breast_cancer
│   ├── ductal_carcinoma_R1_L001.fastq.gz
│   └── ductal_carcinoma_R2_L001.fastq.gz
├── human_adenoma_R1_L001.fastq.gz
└── human_adenoma_R2_L001.fastq.gz
    
```

The first flag that we can set based on this information is the `-i` flag, the `-i` flag takes the input file directory as argument, in this case it will look like this `./scripts/start-workflow.sh -i /tsd/p172ncspm/data/durable/input-data/ [...]`.  
The next step is to prepare the tsv file.

### Preparing the tsv file
The tsv file is used to tell the workflow which files belong to the same sample. If you want to analyse many samples at once you simply put all information per sample and file in the tsv file. The sample information is also used to create proper read groups, this will help you identify the files since this information is put in the header of the bam files by the workflow. This information is also used by the workflow to name the output files.  
Either use the sample file [here](https:/github.com/elixir-no-nels/Selma/blob/master/samples.tsv) or copy/paste the header below (remember to separate the columns with tabs!) into a new file and add your sample information below:  
```bash
flowcell	sample	library	lane	R1	R2
```  
If you don't have information about flowcell id, sample id, library id or lane number, you still need to put something there, just pretend that you have the correct information if you don't have it. What matters most to get things to run is that you write the file paths correctly and that you use the same sample name for each sample input file pair.  

Let's continue the thought experiment and use our pretend input files from above in our new tsv file like this:  

```bash
flowcell	sample	library	lane	R1	R2
FlowcellX	HA001	libHA	L001	human_adenoma_R1_L001.fastq.gz	human_adenoma_R2_L001.fastq.gz
FlowcellX	BCD001	libBCD	L001	breast_cancer/ductal_carcinoma_R1_L001.fastq.gz	breast_cancer/ductal_carcinoma_R2_L001.fastq.gz
```
Now save the file as `my-samples.tsv` and add `-t /tsd/p172ncspm/data/durable/input-data/my-samples.tsv` in the start command like so: 
```bash
./scripts/start-workflow.sh -i /tsd/p172ncspm/data/durable/input-data/ -t /tsd/p172ncspm/data/durable/input-data/my-samples.tsv [...]
```

Now let's move on to the output directory.  

### Set the output directory
Begin by creating a new directory like so: `mkdir /tsd/p172ncspm/data/durable/Selma-outputs` and simply add `-o /tsd/p172ncspm/data/durable/Selma-outputs` to the command line like so: 
```bash
./scripts/start-workflow.sh -i /tsd/p172ncspm/data/durable/input-data/ -t /tsd/p172ncspm/data/durable/input-data/my-samples.tsv -o /tsd/p172ncspm/data/durable/Selma-outputs [...]
```
The workflow will write the outputs to an automatically generated new uniqely named directory using the naming pattern `Selma-yyyy-mm-dd-HH-MM-SS` inside the `/tsd/p172ncspm/data/durable/Selma-outputs` directory. 

The final required step is to select the reference file version. 

### Selecting reference file version
You have a choice of two reference file versions, either the `b37` decoy version, or the `hg38` version.
[This article](https:/arxiv.org/pdf/1404.0929.pdf) has some interesting comparisons of the b37/GRCh37 and hg38/GRCh38 reference files.

If you don't know which one to choose you should probably use hg38, it's generally more complete compared to b37 according to the article above.

The flag for reference version selection is `-r`, so the resulting command line so far looks like this:  
```bash
./scripts/start-workflow.sh -i /tsd/p172ncspm/data/durable/input-data/ -t /tsd/p172ncspm/data/durable/input-data/my-samples.tsv -o /tsd/p172ncspm/data/durable/Selma-outputs -r hg38
```

And that's it! You should be able to run the workflow now by running the following:  
```bash
cd /cluster/projects/p172/UiO-Cancer/
./scripts/start-workflow.sh -i /tsd/p172ncspm/data/durable/input-data/ -t /tsd/p172ncspm/data/durable/input-data/my-samples.tsv -o /tsd/p172ncspm/data/durable/Selma-outputs -r hg38
```
This will run Selma on Colossus using the Singularity image that was built with [this](https:/github.com/elixir-no-nels/Selma/blob/master/singularity/BuildSingularityImage.sh) script.

### Optional custom interval file
This feature has not been quality tested yet but should in principle work as it should.  
If you want to do use a custom interval list you need to be certain that your interval file is compatible with the default reference files.  
The hg38 reference files have chromosome names of the format: `chr1	chr2	chr3	etc`  
For a complete list of contigs in the hg38 reference fasta file you can check out this .dict file: [ftp:/gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict](ftp:/gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict)  
In case it asks for password there is none, just hit enter.

The b37 reference files have chromosome names like this: `1	2	3	etc`.  
For a complete list of contigs in the b37 reference fasta file you can check out this .dict file: [ftp:/gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.dict.gz](ftp:/gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.dict.gz)  
In case it asks for password there is none, just hit enter.

There's a default hg38 and a b37 interval list that is always used for wgs calling, the intervals are made to exclude e.g centromeric regions since these regions don't add any useful information for variant calling. Here's some more information: [https:/software.broadinstitute.org/gatk/documentation/article?id=11009](https:/software.broadinstitute.org/gatk/documentation/article?id=11009)  
Here's the most relevant part from that link:  
> **Whole genomes (WGS)**  
	For whole genome sequence, the intervals lists don’t depend on the prep (since in principle you captured the “whole genome”) so instead it depends on what regions of the genome you want to blacklist (e.g. centromeric regions that waste your time for nothing) and how the reference genome build enables you to cut up regions (separated by Ns) for scatter-gather parallelizing.  
	We make our WGS interval lists available, and the good news is that, as long as you're using the same genome reference build as us [The genome reference build on TSD is the same as build as they are using], you can use them with your own data even if it comes from somewhere else -- assuming you agree with our decisions about which regions to blacklist! Which you can examine by looking at the intervals themselves. However, we don't currently have documentation on their provenance, sorry -- baby steps.

If you want to use a custom interval list the flag is `-l` and you need to provide a complete file path to it like so: `-l /put/the/actual/absolute/path/here/hg38_custom_interval.list`  
And the command line:

```bash
./scripts/start-workflow.sh -i /tsd/p172ncspm/data/durable/input-data/ -t /tsd/p172ncspm/data/durable/input-data/my-samples.tsv -o /tsd/p172ncspm/data/durable/Selma-outputs -r hg38 -l /put/the/actual/absolute/path/here/hg38_custom_interval.list
```
Supported interval list formats include `.bed`, `.intervals`, `.list` and `.interval_list`
