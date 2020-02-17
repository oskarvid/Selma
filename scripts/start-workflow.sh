#!/bin/bash

set -Eeo pipefail
# Uncomment to enable debugging
#set -vo xtrace

source scripts/functions.sh

trap 'die && exec 2>&4 1>&3' SIGTERM SIGINT SIGKILL ERR

# Set default variables
START=$(date +%s)
DATE=$(date +%F.%H.%M.%S)
SPLIT=FALSE

# Function that gives the user troubleshooting information and attempts to clean up stray files so as not to litter the intermediary workflow directory
die () {
	printf "#########################################\n"
	err "$0 failed at line $BASH_LINENO"
	sleep 1
	grdelstag $LOGFILE $STAGINGDIR $WFDIRNAME
	grslcan $LOGFILE
	mvslot $WFDIR $OUTPUTDIR $WFDIRNAME ${SLURMID[@]}
	rmtmp $RMTMP
	benchmark $START $LOGFILE
	if [[ -f $LOGFILE ]]; then
		err "Execution stopped unexpectedly, here's the logfile: $LOGFILE"
	fi
	printf "#########################################\n"
	kill -s 15 -$(ps opgid= $PID | tr -d ' ')
	exec 2>&4 1>&3
	exit 1
}

# Print the Selma ascii serpent
selma-ascii

# Print usage message and exit if no arguments are given
if [[ ! $@ ]]; then
	usage
fi

while getopts 'm:c:i:t:o:r:l:e:hs' flag; do
	case "${flag}" in
	m)
		m=${OPTARG}
		if [[ $m == local ]]; then
			MODE="${m}"
			inf "Using $MODE mode"
			source settings/settings.conf
			continue
		elif [[ $m == tsd ]]; then
			if [[ $(which sbatch) ]]; then
				MODE="${m}"
				inf "Using $MODE mode"
				source settings/settings.conf
				continue
			else
				err "tsd mode selected but sbatch is not found, are you sure you are using the correct mode?"
				die
			fi
		elif [[ $m != local || ! $m == tsd ]]; then
			err "Selected mode \"$m\" is not supported, valid options are tsd and local"
			die
		fi
		;;

	c)
		c=${OPTARG}
		if [[ $c == docker || $c == singularity ]]; then
			CONT="${c}"
		else
			err "\"${c}\" container is not supported, valid options are docker and singularity"
		fi
		;;
	s)
		s=${OPTARG}
		SPLIT=TRUE
		inf "Running one sample per node"
		;;
	e)
		e=${OPTARG}
		if [[ -a $e ]]; then
			SSAM=$(realpath $e)
			inf "Using custom sample selection from file: ${e}"
		elif [[ ! -a $e ]]; then
			err "The specified custom sample file does not seem to exist"
			usage
		fi
		;;
	i)
		i=${OPTARG}
		if [[ ! -d $(realpath "${i}") ]]; then
			err "Folder path \"$(realpath "${i}")\" cannot be found"
			usage
		else
			INPUTS=$(realpath "${i}")
		fi
		;;
	t)
		t=${OPTARG}
		if [[ ! -f $(realpath "${t}") ]]; then
			err "Samples file \"$(realpath "${t}")\" cannot be found"
			usage
		else
			TSV=$(realpath "${t}")
		fi
		;;
	o)
		o=${OPTARG}
		if [[ ! -d $(realpath "${o}") ]]; then
			err "Output folder path \"$(realpath "${o}")\" cannot be found"
			usage
		else
			OUTPUTDIR=$(realpath "${o}")
		fi
		;;
	r)
		r=${OPTARG}
		if [[ "${r}" == "hg38" ]]; then
			HGVER="${r}"
		elif [[ "${r}" == "b37" ]]; then
			HGVER="${r}"
		else
			err "Reference version must be set to hg38 or b37, \"${r}\" is not a valid option"
			usage
		fi
		;;
	l)
		# Optional interval list
		l=${OPTARG}
		if [[ ! -f $(realpath "${l}") ]]; then
			err "Interval file \"$(realpath "${l}")\" cannot be found"
			usage
		else
			INTERVAL=$(realpath "${l}")
		fi
		;;
	h)
		h=${OPTARG}
		usage
		;;
	\?)
		usage
		;;
	esac
done

# Check that container mode has been set correctly
if [[ -z "${m}" ]]; then
	err "-m is required and must be set to local or tsd"
	usage
elif [[ "${m}" == tsd && "${c}" == docker ]]; then
	err "Cannot select docker mode when running on TSD, only singularity is possible"
	die
elif [[ $MODE == local ]]; then
	if [[ -z "${c}" ]]; then
		err "-c is required, valid options are docker or singularity"
		usage
	elif [[ $CONT == singularity || $CONT == docker ]]; then
		inf "Using $CONT profile"
	else
		err "-c must be set to docker or singularity, \"$c\" is not a valid option"
		usage
	fi
fi

# make a log file in $OUTPUTDIR/logfiles named logfile-$DATE.log
mkdir -p $OUTPUTDIR/logfiles
LOGFILE=$OUTPUTDIR/logfiles/logfile-$DATE.log
exec 3>&1 4>&2
exec &> >(tee $LOGFILE) 2>&1

# Check that input file path was set
if [[ -z "${i}" ]]; then
	err "-i is required and path must be absolute"
	usage
fi

# Check that input tsv file is set
if [[ -z "${t}" ]]; then
	err "-t is required and path must be absolute"
	usage
fi

# Check that output directory is set
if [[ -z "${o}" ]]; then
	err "-o is required and path must be absolute"
	usage
fi

# Check that reference file version is set
if [[ -z "${r}" ]]; then
	err "-r is required, valid options are hg38 or b37"
	usage
fi

# Check that the supplied input file folder contains supported input files before proceeding
if [[ ! $(find $INPUTS -type f -name *.fq) && ! $(find $INPUTS -type f -name *.fq.gz) && ! $(find $INPUTS -type f -name *.fastq) && ! $(find $INPUTS -type f -name *.fastq.gz) ]]; then
	err "The provided input directory $INPUTS does not contain supported input files"
	err "Supported file formats are *.fq, *.fq.gz, *.fastq or *.fastq.gz"
	err "Please check if the correct folder path is used"
	warn "Here is the content of the supplied input file folder: \n$(ls -lh $INPUTS)"
	usage
fi

# Check if a custom interval file has been set
# Then check if the file exists
# If no custom interval has been set it will use the default interval file for the selected reference file version
# If a custom interval file has been set and it exists it will be used
# If a custom interval file has been set but it does not exist, or if it has the wrong file suffix the scstart script will stop with an error message
if [[ -n "${INTERVAL:-}" ]]; then
	if [[ -a "${INTERVAL:-}" ]]; then
		if [[ "${INTERVAL:-}" == *.bed || "${INTERVAL:-}" == *.interval_list || "${INTERVAL:-}" == *.list || "${INTERVAL:-}" == *.intervals ]]; then
			inf "Using custom interval file: ${INTERVAL:-}"
				if [[ $MODE == local ]]; then
					mkdir -p $(pwd)/workspace/staging/
					rprog $INTERVAL $(pwd)/workspace/staging/
					INTERVAL=staging/$(basename $INTERVAL)
				elif [[ $MODE == tsd ]]; then
					rprog ${INTERVAL:-} $STAGINGDIR/$WFDIRNAME-$DATE/custom-interval.${INTERVAL##*.}
					INTERVAL=$STAGINGDIR/$WFDIRNAME-$DATE/custom-interval.${INTERVAL##*.}
				fi
		else
			err "The custom interval file must be a .bed, .interval_list, .list or .intervals file. The detected file format is \'.${INTERVAL##*.}\'"
			usage
		fi
	else
		err "A custom interval file was set but it does not seem to exist, check the file path and try again."
		usage
	fi
else
	if [[ "${HGVER}" == hg38 ]]; then
		INTERVAL=/references/"${HGVER}"/wgs_calling_regions.hg38.interval_list
	elif [[ "${HGVER}" == b37 ]]; then
		INTERVAL=/references/"${HGVER}"/b37_wgs_calling_regions.v1.interval_list
	fi
fi

## The $TSV variable needs to be kept outside of the array, otherwise it can't be set by the for loop in slurm mode
BASIC=($INPUTS $HGVER $INTERVAL $OUTPUTDIR)
SLURM=($JOBNAME $ACCOUNT $TIME $MEMPERCPU $CPUSPERTASK)

# This logic checks what the different variables are set to and sets the correct mode based on the values
# Either it runs both custom sample selection and one sample per node or just one of them
if [[ -n "${SSAM:-}" && "${SPLIT:-}" == TRUE ]]; then
	select_samples "${SSAM:-}" "${TSV}" "${DATE}"
	one_sample_per_node "${CUSTTSV:-}" "${DATE}"
elif [[ -n "${SSAM:-}" ]]; then
	select_samples "${SSAM:-}" "${TSV}" "${DATE}"
elif [[ "${SPLIT:-}" == TRUE ]]; then
	one_sample_per_node "${TSV}" "${DATE}"
fi

# This loop starts the slurm job(s).
SAMPLESFILES="${TSVFILES[@]:-"${CUSTTSV:-$TSV}"}"
for TSV in ${SAMPLESFILES[@]:-$TSV}; do
	sleep 2
	./scripts/checker.sh $TSV ${BASIC[@]} ${SLURM[@]} $RMTMP $LOGFILE &> >(tee $LOGFILE) 2>&1 & #This ampersand forks each iteration of the for loop to the background
	PID=$$
done
wait
