#!/bin/bash
REFERENCES=/home/oskar/01-workspace/01-data/refdata

if [[ ! -z $1 && $1 == hg38 || $1 == b37 ]]; then
	REF=$1
else
	printf "Reference version has not been set or has an incorrect value, either use hg38 or b37.\n"
	printf "Usage: ./start-pipeline.sh reference-version (hg38 or b37)\n"
	exit 1
fi

singularity exec \
-B $REFERENCES:/references \
-B $(pwd):/data \
-W $(pwd) \
singularity/snakemake-germline-tools.simg snakemake -j --config version="${REF}"

