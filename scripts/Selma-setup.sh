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
                        echo "ERROR: $s directory does not exist"
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
                        echo "ERROR: $b directory does not exist"
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
                        echo "ERROR: $g directory does not exist"
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
        rprog /tsd/p172ncspmdata/data/no-backup/oskar/01-workspace/01-data/references/hg38/ $GDEST/hg38
#       rsync -ah --progress /tsd/shared/bioinformatics/reference-data/b37/intervals/Nimblegen_SeqCap_EZ_Exome_v3/ $GDEST/hg38
fi

if [[ $BDEST == false ]]; then
        :
elif [[ $BDEST == "${b%/}" ]]; then
        rprog /tsd/p172ncspmdata/data/no-backup/oskar/01-workspace/01-data/references/b37/ $BDEST/b37
#       rsync -ah --progress /tsd/shared/bioinformatics/reference-data/b37/intervals/Nimblegen_SeqCap_EZ_Exome_v3/ $BDEST/b37
fi

if [[ $SDEST == false ]]; then
        :
elif [[ $SDEST == "${s%/}" ]]; then
        rprog /tsd/p172ncspmdata/data/no-backup/oskar/01-workspace/04-pipelines/Selma/ $SDEST/Selma
#       rsync -ah --progress /tsd/shared/bioinformatics/reference-data/b37/intervals/Nimblegen_SeqCap_EZ_Exome_v3/ $SDEST/Selma
fi




## Fun fact: If the reference files are stored as tar.gz archives it's possible to save 9.5GB storage on the source disk and at least up to 2.5 minutes in transfer time compared with using rsync to transfer the untarred archives
## To be more specific, in the hg38 benchmark I did, I used the line below to untar the pigzhg38.tar.gz, which is 9GB, to the destination directory in 13 minutes and 57 seconds
## tar -xvzf /tsd/p172ncspmdata/data/no-backup/oskar/01-workspace/01-data/references/pigzhg38.tar.gz -C $GDEST/hg38
## The rsync line below transfered the 14GB large directory that is created by untarring the pigzhg38 archive, in 14 minutes and 42 seconds
## rprog /tsd/p172ncspmdata/data/no-backup/oskar/01-workspace/01-data/references/hg38/ $GDEST/hg38
## This saves 40ish seconds, or about 5%, in a 14-15 minute long transfer and that isn't much, but on the other hand, saving 5GB storage, or about 35%, in a 14GB file, is much.

## I did the same test with the b37 reference files and it took 10 minutes and 32 seconds with tar and 12 minutes and 38 seconds with rsync.
## So the main benefit is disk savings rather than transfer durations because it's not super critical since the transfer is only meant to be done once very rarely, if ever twice in the lifetime of a project.
## But the source files, archived or not, need to be permanently available, so storing them in a tar.gz format isn't unreasonable if it's worth the price in inconvenience if you at some point need access to a single file from the tarred directory
