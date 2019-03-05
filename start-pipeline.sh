#!/bin/bash
######################
# This script is both a place holder for the singularity start command \
# as well as a literal start script.
######################

# Change the path to where your reference files are located
REFERENCES=/home/oskar/01-workspace/01-data/refdata

# Checks if an argument has been supplied and if that argument is either \
# hg38 or b37, if no argument is given, or if it's not hg38 or b37 the script will \
# exit.
if [[ ! -z $1 && $1 == hg38 || $1 == b37 ]]; then
	REF=$1
else
	printf "Reference version has not been set or has an incorrect value, either use hg38 or b37.\n"
	printf "Usage: ./start-pipeline.sh reference-version (hg38 or b37)\n"
	exit 1
fi

# The workflow assumes that the top level directories for the reference files \
# is named /references and that the workflow top level directory is /data.
# Hardcoding these filepaths simplify any further assumptions that are needed \
# for the workflow.
# The working directory is set to /data because that's where the snakefile is.
# -j enables running rules in parallel.
# --config version="${REF}" is used to select the reference file version from the config file.
singularity exec \
-B $REFERENCES:/references \
-B $(pwd):/data \
-W /data \
singularity/snakemake-germline-tools.simg snakemake -j --config version="${REF}"
