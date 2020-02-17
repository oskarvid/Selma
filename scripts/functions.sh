#!/bin/bash

#set -o xtrace

# Function that prints help message
usage () {
echo "${GREEN}""Usage:${RESET} $0 [ -h ] [ -i input directory ] [ -t tsv file ] [ -o output directory ] [ -r genome version] [ -l interval file ] [ -c singularity or docker ] [ -m execution mode ] [ -s one sample per node ] [ -e sample subset file ] \
-h : Print this help message \

-r : Human genome reference: b37 or hg38, required \

-i : Path of the input file directory, path must be absolute, required \

-t : Tab separated file with fastq info, path must be absolute, required \

-o : Path of the output folder, path must be absolute, required \

-l : Optional interval file, path must be absolute \

-c : Container version, valid options are singularity or docker, only for local execution \

-m : Mode, valid options are local or tsd \

-s : Split samples on one node each, improves execution time considerably for two or more samples \

-e : Select sample subset from master sample file \
" 1>&2; exit 0; }

# Make your terminal text output beautiful <3
# The following variables are used to define the color in the functions below
BLACK=$(tput setaf 0)
RED=$(tput setaf 1)
GREEN=$(tput setaf 2)
YELLOW=$(tput setaf 3)
BLUE=$(tput setaf 4)
MAGENTA=$(tput setaf 5)
CYAN=$(tput setaf 6)
WHITE=$(tput setaf 7)
RESET=$(tput sgr 0)

# Show date and time in "YYYY-MM-DD HH:MM:SS" format
showdate () {
	local DATE
	DATE=$(date "+%F %T")
	printf "$DATE"
}

# Print red "ERROR" message
err () {
	local DATE
	DATE=$(showdate)
	printf "[$DATE][${RED}ERROR${RESET}]: $1\n" 1>&2
}

# Print yellow "WARNING" message
warn () {
	local DATE
	DATE=$(showdate)
	printf "[$DATE][${YELLOW}WARNING${RESET}]: $1\n"
}

# Print neutral "INFO" message
inf () {
	local DATE
	DATE=$(showdate)
	printf "[$DATE][INFO]: $1\n"
}

# Print message with "SUCCESS" in green
succ () {
	local DATE
	DATE=$(showdate)
	printf "[$DATE][${GREEN}SUCCESS${RESET}]: $1\n"
}

# Nice rsync function that shows file transfer progress
rprog () {
	rsync -ah --progress $@
}

# Nice Selma ascii art
selma-ascii () {
inf "Starting Selma"
cat << "EOF"
                 _..--+~/@-~--.
             _-=~      (  .   "}
          _-~     _.--=.\ \""""
        _~      _-       \ \_\
       =      _=          '--'
      '      =                             .
     :      :       ____                   '=_. ___
___  |      ;                            ____ '~--.~.
     ;      ;                               _____  } |
  ___=       \ ___ __     __..-...__           ___/__/__
     :        =_     _.-~~          ~~--.__
_____ \         ~-+-~                   ___~=_______
     ~@#~~ == ...______ __ ___ _--~~--_
EOF
}

# remove files with no stdout output
clean () {
	rm $@ >/dev/null
}

# Count columns in tsv file, this is only for the select_samples function
count_cols () {
	awk '{ print NF }' $1 | sort -nrk1,1 | head -1
}

# Select individual samples from a mother sample file using a custom tsv file
select_samples () {
	SSAM=$1
	SAMPLEFILE=$2
	DATE="${3}"
	COLS=$(count_cols $SSAM)

	mkdir -p .tmpSamples-$DATE
	RMTMP=".tmpSamples-$DATE"

	if [[ "${COLS}" -eq 1 ]]; then
		mapfile -t SAMPLES < "${SSAM}"
		printf "flowcell\tsample\tlibrary\tlane\tR1\tR2\n" > .tmpSamples-$DATE/cust-list-$DATE.tsv
		for sample in ${SAMPLES[@]}; do
			if [[ ! $(grep $sample $SAMPLEFILE) ]]; then
				err "Sample '$sample' from $SSAM does not exist in $SAMPLEFILE"
				err "Make sure that the custom samples file is correct"
				die
			else
				grep $sample $SAMPLEFILE >> .tmpSamples-$DATE/cust-list-$DATE.tsv
			fi
		done
	else
		err "Custom input file has $COLS columns"
		err "The custom sample selection file cannot have more than one column"
	fi

	CUSTTSV=".tmpSamples-$DATE/cust-list-$DATE.tsv"
	export CUSTTSV
	export RMTMP
}

# Run one sample per node
one_sample_per_node () {
# Put the TSV file and the selected samples file in named variables
	SAMPLEFILE=$1
	DATE="${2}"

	# Clear the TSVFILES variable
	TSVFILES=
	mapfile -t SAMPLES < <(awk '{ print $2 }' $SAMPLEFILE | uniq | tail -n +2)
	mkdir -p .tmpSamples-$DATE
	RMTMP=".tmpSamples-$DATE"

	for sample in ${SAMPLES[@]}; do
		printf "flowcell\tsample\tlibrary\tlane\tR1\tR2\n" > .tmpSamples-$DATE/$sample.tsv
		grep -w $sample $SAMPLEFILE >> .tmpSamples-$DATE/$sample.tsv
		TSVFILES+=("$(pwd)/.tmpSamples-$DATE/$sample.tsv")
	done

	export TSVFILES
	export RMTMP
}

grdelstag () {

	LOGFILE=$1
	STAGINGDIR=$2
	WFDIRNAME=$3

	read -a RMSTAGING < <(strings $LOGFILE | grep -ow "$STAGINGDIR/$WFDIRNAME-[0-9]*-[0-9]*-[0-9]*\.[0-9]*\.[0-9]*.[0-9]*" | sort -u | tr '\n' ' ')
        RMSTAGING=${RMSTAGING% }
        if [[ -a ${RMSTAGING[0]} ]]; then
                echo "[$(showdate)][${YELLOW}WARNING${RESET}]: Will now attempt to delete ${RMSTAGING[@]}"
                rm -fr ${RMSTAGING[@]} && \
                echo "[$(showdate)][${YELLOW}WARNING${RESET}]: Deleted ${RMSTAGING[@]}" || echo "[$(showdate)][${RED}ERROR${RESET}]: Deleting ${RMSTAGING[@]} failed. Try deleting manually"
        fi
}

# Delete the file staging directory from the checker script
delstag () {

	STAGINGDIR=$1
	WFDIRNAME=$2
	DATE=$3

        if [[ -a $STAGINGDIR/$WFDIRNAME-$DATE ]]; then
                warn "Will now attempt to delete $STAGINGDIR/$WFDIRNAME-$DATE"
                rm -r $STAGINGDIR/$WFDIRNAME-$DATE && \
                succ "$STAGINGDIR/$WFDIRNAME-$DATE has been deleted." || err "Deleting $STAGINGDIR/$WFDIRNAME-$DATE failed. Try deleting it manually."
        fi
}

grslcan () {

	LOGFILE=$1

        read -a SLURMID < <(strings $LOGFILE | grep -ow "ID [0-9]*" | sed -e 's/[A-Z]*//g' -e 's/^ //' -e 's/ $//' | sort -u | tr '\n' ' ')
        if [[ ${SLURMID[0]} == [0-9]* ]]; then
                echo "[$(showdate)][${YELLOW}WARNING${RESET}]: Stopping slurm job(s) ${SLURMID[@]}"
                scancel ${SLURMID[@]} && \
                echo "[$(showdate)][${YELLOW}WARNING${RESET}]: Slurm job(s) ${SLURMID[@]} has/have successfully been stopped"
        fi
}

# Cancel slurm job if it's still running
slcan () {

	SLURMID=$1
	STATUS=$(squeue -h -t all -j $SLURMID -o %t 2>&1)
	if [[ $SLURMID == [0-9]* && $STATUS != '' ]]; then
		warn "Stopping slurm job $SLURMID"
		scancel $SLURMID && \
		warn "Slurm job $SLURMID has successfully been stopped."
	fi
}

# Move slurm-123456.out file to the final output directory
mvslot () {

	ALL=("$@")
	WFDIR=${ALL[0]}
	OUTPUTDIR=${ALL[1]}
	WFDIRNAME=${ALL[2]}
	SLURMID=${ALL[@]:3}

        if [[ -e $WFDIR/slurm-${SLURMID[0]}.out && -d ${OUTPUTDIR}/$WFDIRNAME-$DATE ]]; then
                for ID in ${SLURMID[@]}; do
                        inf "Copying slurm-$ID.out to ${OUTPUTDIR}/$WFDIRNAME-$DATE"
                        rprog $WFDIR/slurm-$ID.out ${OUTPUTDIR}/$WFDIRNAME-$DATE/slurm-$ID.out
                        inf "Here is the slurm job file: ${OUTPUTDIR}/$WFDIRNAME-$DATE/slurm-$ID.out"
		done
        fi
}

# Remove temporary directory for split up tsv files
rmtmp () {

	RMTMP=$1

	if [[ -d $RMTMP ]]; then
		warn "Removing $RMTMP directory"
		rm -r $RMTMP
		warn "The $RMTMP directory has been removed" || warn "Deleting $RMTMP failed. Try deleting manually"
	fi
}

# Print benchmarking data
benchmark () {

	START=$1
	LOGFILE=$2

	FINISH=$(date +%s)
	EXECTIME=$(( $FINISH-$START ))
	inf "Exited on $(date)"
	printf "[$(showdate)][INFO]: Duration: %dd:%dh:%dm:%ds\n" $((EXECTIME/86400)) $((EXECTIME%86400/3600)) $((EXECTIME%3600/60)) $((EXECTIME%60))
}
