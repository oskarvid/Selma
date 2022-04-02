#!/bin/bash

# Uncomment for debugging
#set -o xtrace

SDEST=false
BDEST=false
GDEST=false

usage () {
        echo "Usage: $0 [OPTION] DESTINATION"
        echo "-s, Destination for Selma workflow installation"
        echo "-b, Destination for b37 reference files"
        echo "-g, Destination for hg38 reference files"
        echo "-h, this help message"
}

rprog () {
	rsync -ahv --progress $@
}

if [[ $# -eq 0 ]]; then
        usage
fi

while getopts 's:b:g:h' flag; do
        case "${flag}" in
        s)
                s=${OPTARG}
                if [[ ! -d $s ]]; then
                        echo "ERROR: $s directory can not be found"
                        echo "Create the directory or select another destination"
                        exit
                fi

                if [[ -d $s/Selma ]]; then
                        echo "ERROR: Selma directory already exists"
                        echo "Remove Selma directory in the submitted destination directory or select a different installation directory"
                else
                        SDEST="${s%/}"
                fi
                ;;
        b)
                b=${OPTARG}
                if [[ ! -d $b ]]; then
                        echo "ERROR: $b directory can not be found"
                        echo "Create the directory or select another destination"
                        exit
                fi

                if [[ -d $b/b37 ]]; then
                        echo "ERROR: b37 directory already exists"
                        echo "Remove b37 directory in the submitted destination directory or select a different destination directory"
                else
                        BDEST="${b%/}"
                fi
                ;;
        g)
                g=${OPTARG}
                if [[ ! -d $g ]]; then
                        echo "ERROR: $g directory can not be found"
                        echo "Create the directory or select another destination"
                        exit
                fi

                if [[ -d $g/hg38 ]]; then
                        echo "ERROR: hg38 directory already exists"
                        echo "Remove hg38 directory in the submitted destination directory or select a different destination directory"
                else
                        GDEST="${g%/}"
                fi
                ;;
        h)
                h=${OPTARG}
                usage
		;;
        ?)
                usage
                ;;
        esac
done

if [[ $GDEST == false ]]; then
        :
elif [[ $GDEST == "${g%/}" ]]; then
        rprog /tsd/shared/bioinformatics/reference-data/selma-references/hg38/ $GDEST/hg38
fi

if [[ $BDEST == false ]]; then
        :
elif [[ $BDEST == "${b%/}" ]]; then
        rprog /tsd/shared/bioinformatics/reference-data/selma-references/b37/ $BDEST/b37
fi

if [[ $SDEST == false ]]; then
        :
elif [[ $SDEST == "${s%/}" ]]; then
        rprog --exclude=.git /tsd/shared/bioinformatics/workflows/Selma $SDEST/
	rprog /tsd/shared/bioinformatics/containers/singularity/snakemake-germline-tools.simg $SDEST/Selma/singularity
fi
