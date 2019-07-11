#!/bin/bash

set -Ee
# Uncomment to enable debugging
#set -vo xtrace

#source $(dirname $0)/functions.sh
#source $(dirname $0)/config.sh
source scripts/functions.sh
source settings/settings.conf

# Run the die function to clean up files and give some troubleshooting help in case the workflow is killed or fails for some reason
trap die ERR SIGTERM SIGINT

# Print the Selma ascii serpent
selma-ascii

# Function that gives the user troubleshooting information and attempts to clean up stray files so as not to litter the intermediary workflow directory
die () {
	printf "#########################################\n"
	err "$0 failed at line $BASH_LINENO"
	warn "Will now attempt to delete $STAGINGDIR/$WFDIRNAME-$DATE"
	rm -r $STAGINGDIR/$WFDIRNAME-$DATE && \
	warn "$STAGINGDIR/$WFDIRNAME-$DATE has been deleted." || err "Deleting $INTERMEDSTOR/$WFDIRNAME-$DATE failed. Try deleting it manually."

	if [[ $SLURMID == [0-9]* ]]; then
		warn "Stopping slurm job $SLURMID"
		scancel $SLURMID && \
		warn "Slurm job $SLURMID has successfully been stopped."
	fi

	if [[ -e $WFDIR/slurm-$SLURMID.out ]]; then
		inf "Copying slurm-$SLURMID.out to ${OUTPUTDIR}/$WFDIRNAME-$DATE"
		rprog $WFDIR/slurm-$SLURMID.out ${OUTPUTDIR}/$WFDIRNAME-$DATE/slurm-$SLURMID.out
		inf "Here is the slurm job file: ${OUTPUTDIR}/$WFDIRNAME-$DATE/slurm-$SLURMID.out"
	fi

	err "Something seems to have gone wrong, here's the logfile: $LOGFILE"
	printf "#########################################\n"
	exit 1
}

# Print usage message and exit if no arguments are given
if [[ ! $@ ]]; then
	usage
fi

# Set default variables
START=$(date +%s)
DATE=$(date +%F.%H.%M.%S)

while getopts 'i:t:o:r:l:h' flag; do
	case "${flag}" in
	i)
		i=${OPTARG}
		if [[ ! -d /"${i}" ]]; then
			err "Folder path for the -i flag must be absolute"
			usage
		else
			INPUTS="${i%/}"
		fi
		;;
	t)
		t=${OPTARG}
		if [[ ! -f /"${t}" ]]; then
			err "File path for the -t flag must be absolute"
			usage
		else
			TSV="${t}"
		fi
		;;
	o)
		o=${OPTARG}
		if [[ ! -d /"${o}" ]]; then
			err "Output folder for the -o flag must exist and the path must be absolute"
			usage
		else
			OUTPUTDIR="${o%/}"
		fi
		;;
	r)
		r=${OPTARG}
		if [[ "${r}" == "hg38" ]]; then
			HGVER="${r}"
			inf "Reference file set to "${HGVER}""
		elif [[ "${r}" == "b37" ]]; then
			HGVER="${r}"
			inf "Reference file set to "${HGVER}""
		else
			err "Reference version must be set to hg38 or b37, "${r}" is not a valid option"
			usage
		fi
		;;
	l)
		l=${OPTARG}

		# Optional interval list
		if [[ ! -f /"${l}" ]]; then
			err "Interval file "${l}" does not exist or file path is not absolute"
			usage
		else
			INTERVAL="${l}"
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

#### Verify that all mandatory flags were set
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

# make a log file at $OUTPUTDIR/$WFDIRNAME-$DATE/ named logfile.log
mkdir $OUTPUTDIR/$WFDIRNAME-$DATE
LOGFILE=$OUTPUTDIR/$WFDIRNAME-$DATE/logfile.log
exec &> >(tee $LOGFILE)

# Check that the supplied input file folder contains supported input files before proceeding
if [[ ! $(compgen -G $INPUTS/*.fq) && ! $(compgen -G $INPUTS/*.fq.gz) && ! $(compgen -G $INPUTS/*.fastq) && ! $(compgen -G $INPUTS/*.fastq.gz) ]]; then
	err "The provided input directory $INPUTS does not contain supported input files"
	err "Supported file formats are *.fq, *.fq.gz, *.fastq or *.fastq.gz"
	err "Please check if the correct folder path is used"
	warn "Here is the content of the supplied input file folder: \n$(ls -lh $INPUTS)"
	usage
fi

# Make uniquely named temporary directory on /cluster/projects/p172/ so there can be no mistakes if the pipeline is started again when a 
# file transfer is already ongoing so there is no confusion about which input and output files belong to which execution
mkdir $STAGINGDIR/$WFDIRNAME-$DATE

# Check if a custom interval file has been set
# Then check if the file exists
# If no custom interval has been set it will use the default interval file for the selected reference file version
# If a custom interval file has been set and it exists it will be used
# If a custom interval file has been set but it does not exist, or if it has the wrong file suffix the scstart script will stop with an error message

if [[ -n "${INTERVAL:-}" ]]; then
	if [[ -a "${INTERVAL:-}" ]]; then
		if [[ "${INTERVAL:-}" == *.bed || "${INTERVAL:-}" == *.interval_list || "${INTERVAL:-}" == *.list || "${INTERVAL:-}" == *.intervals ]]; then
			inf "Using custom interval file: ${INTERVAL:-}"
			rprog ${INTERVAL:-} $STAGINGDIR/$WFDIRNAME-$DATE/custom-interval.${INTERVAL##*.}
			INTERVAL=$STAGINGDIR/$WFDIRNAME-$DATE/custom-interval.${INTERVAL##*.}
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
		INTERVAL="${REFERENCES}"/"${HGVER}"/wgs_calling_regions.hg38.interval_list
	elif [[ "${HGVER}" == b37 ]]; then
		INTERVAL="${REFERENCES}"/"${HGVER}"/b37_wgs_calling_regions.v1.interval_list
	fi
fi

# Make fresh copy of the workflow directory in case multiple executions of the workflow are started simultaneously
inf "Copying $WFDIR to staging directory: $STAGINGDIR/$WFDIRNAME-$DATE\n"
rprog $WFDIR/ $STAGINGDIR/$WFDIRNAME-$DATE

# Copy TSV file to $STAGINGDIR/$WFDIR-$DATE/workspace/samples.tsv
inf "Copying $TSV to $STAGINGDIR/$WFDIRNAME-$DATE/workspace/samples.tsv\n"
rprog $TSV $STAGINGDIR/$WFDIRNAME-$DATE/workspace/samples.tsv

# Say that the input files/folder is being copied to $STAGINGDIR/$WFDIRNAME-$DATE/workspace/ on the /cluster disk
inf "Copying input files to $STAGINGDIR/$WFDIRNAME-$DATE/workspace/ this may take a while if your input files are large.\n"
rprog $INPUTS/* "$STAGINGDIR/$WFDIRNAME-$DATE"/workspace/

# Run sbatch with all variables set and store the slurm id in $SLURMID
SLURMID=$(sbatch \
--export=\
HGVER=$HGVER,\
REFERENCES=$REFERENCES,\
INTERVAL=$INTERVAL,\
WFDIR=$STAGINGDIR/$WFDIRNAME-$DATE \
scripts/RunOnNode.sbatch | awk '{ print $4 }')

if [[ $SLURMID == '' ]]; then
	err "The job did not get submitted to slurm, check the sbatch command."
	err "Maybe the logfile has some answers too: $LOGFILE"
	die
fi

succ "Workflow with slurm ID $SLURMID has been submitted to the Colossus slurm queue\n"

# Set to 0 so that the $STATUS variable is checked the first time the $STATUS variable is created
PREVSTAT=0

# These *NOCHECK variables will be set to 1 when their respective trigger flag has been detected once
# This is to avoid printing the information message over and over
IPININOCHECK=0
IPSUCCNOCHECK=0
SGSTARTNOCHECK=0
SGNOCHECK=0

# Create variables with file transfer flags instead of using the complete paths in the checks below
JOBSUCCESS=$STAGINGDIR/$WFDIRNAME-$DATE/job-success.txt
JOBFAILED=$STAGINGDIR/$WFDIRNAME-$DATE/job-failed.txt
INPTRINI=$STAGINGDIR/$WFDIRNAME-$DATE/input-file-transfer-ini.txt
INPTRSUC=$STAGINGDIR/$WFDIRNAME-$DATE/input-file-transfer-succ.txt
SINGSTART=$STAGINGDIR/$WFDIRNAME-$DATE/sing-start.txt
SINGDONE=$STAGINGDIR/$WFDIRNAME-$DATE/sing-done.txt

#######################
# The following if else statements check the slurm job status and informs the user what the current slurm job status is
while true; do
	if [[ -e $JOBFAILED && ! -e $JOBSUCCESS ]]; then
		err "The job failure trigger has been detected, will now remove the intermediate files and stop slurm job $SLURMID"
		die

	elif [[ $IPININOCHECK == 0 && -e $INPTRINI ]]; then
		inf "Input files transfer to Colossus initiated, this may take a while if the input files are large"
		IPININOCHECK=1

	elif [[ $IPSUCCNOCHECK == 0 && -e $INPTRSUC ]]; then
		inf "Input files have been successfully transfered to Colossus, starting data analysis"
		IPSUCCNOCHECK=1

	elif [[ $SGSTARTNOCHECK == 0 && -e $SINGSTART ]]; then
		inf "Data analysis has started"
		SGSTARTNOCHECK=1

	elif [[ $SGNOCHECK == 0 && -e $SINGDONE ]]; then
		succ "Data analysis has finished, initiating output file transfer to staging directory"
		SGNOCHECK=1

	elif [[ -e $JOBSUCCESS ]]; then
		succ "File transfer trigger from slurm job with ID $SLURMID has been detected, the file transfer has completed."
		sleep 1
		inf "Will now attempt to copy output files from $STAGINGDIR/$WFDIRNAME-$DATE/workspace/Outputs/ to ${OUTPUTDIR}/$WFDIRNAME-$DATE\n"
		sleep 5
		break
	fi

	# Run squeue, store the current job state code and use it to give information to the user as well as act according to the status
	STATUS=$(squeue -h -t all -j $SLURMID -o %t 2>&1)

	# Begin checking the $STATUS variable by checking if the job has been canceled
	if [[ $STATUS == CA ]]; then
		err "Slurm job with ID $SLURMID has been canceled, will now remove intermediary files"
		die

	# If the new status does not change from the previous status the rest of the checks do not need to be performed
	# OR if the $STATUS variable is empty there is no reason to run the rest of the checks again
	# The $STATUS variable becomes empty when the squeue command fails, it fails quite often due to various connection errors that are due to TSD problems
	# So an empty $STATUS variable does not mean that the status has changed or that the job should get terminated, it should be ignored until squeue establishes contact with slurm again
	# Or in case the job is done $STATUS will remain empty forever, so the trigger files in the if statements above will determine what to do next.
	elif [[ $STATUS == $PREVSTAT ]] || [[ $STATUS == '' ]]; then
		sleep 10
		continue

	# PD means that the slurm job is pending
	elif [[ $STATUS == PD ]]; then
		inf "Slurm job with ID $SLURMID is pending"

	# CF means that the execute node is being configured
	elif [[ $STATUS == CF ]]; then
		inf "Execute node for slurm job with ID $SLURMID is being configured"

	# R means that the sbatch script has started running, it does not mean that the singularity job has started
	# $SINGSTART means that the singularity job has started
	elif [[ $STATUS == R ]]; then
		inf "Slurm job with ID $SLURMID is running"

	# CG means that the slurm job is completing
	elif [[ $STATUS == CG ]]; then
		inf "Slurm job with ID $SLURMID is completing"

	# CD means that the job has completed
	elif [[ $STATUS == CD ]]; then
		inf "Slurm job with ID $SLURMID has completed"

	# F means that the slurm job failed
	elif [[ $STATUS == F ]]; then
		err "Slurm job with ID $SLURMID has failed, will now remove intermediary files"
		die
	fi

	PREVSTAT=$STATUS
	sleep 10
done

# Only copy the final bam file that goes into HaplotypeCaller as well as the VCF files from ApplyVqsr
############ And temporarily also the benchmark files
rprog $STAGINGDIR/$WFDIRNAME-$DATE/workspace/Outputs/benchmarks $INTERMEDSTOR/$WFDIRNAME-$DATE/workspace/Outputs/GatherBamFiles $INTERMEDSTOR/$WFDIRNAME-$DATE/workspace/Outputs/ApplyVqsr* ${OUTPUTDIR}/$WFDIRNAME-$DATE/ && \
succ "Files have been copied from $STAGINGDIR/$WFDIRNAME-$DATE/workspace/Outputs/ to ${OUTPUTDIR}/$WFDIRNAME-$DATE"

# Copy the slurm file to the final output directory
#################### Will this work if the wf is started from any directory? Does the code need to be made specific to something like $(pwd) or something?
inf "Copying slurm-$SLURMID.out to ${OUTPUTDIR}/$WFDIRNAME-$DATE"
rprog $WFDIR/slurm-$SLURMID.out ${OUTPUTDIR}/$WFDIRNAME-$DATE/slurm-$SLURMID.out
#####################

# Clean up any intermediary files
inf "Attempting to delete staging directory: $STAGINGDIR/$WFDIRNAME-$DATE"
rm -r $STAGINGDIR/$WFDIRNAME-$DATE && \
succ "Staging directory $STAGINGDIR/$WFDIRNAME-$DATE has been deleted\n"

# All done!
if [[ ! -z "$(ls "${OUTPUTDIR}/$WFDIRNAME-$DATE")" ]]; then
        succ "Your output files from slurm job $SLURMID are in ${OUTPUTDIR}/$WFDIRNAME-$DATE"
	succ "Workflow execution is all done!"
else
        err "The output directory "${OUTPUTDIR}/$WFDIRNAME-$DATE" from slurm job $SLURMID seems to be empty, something is broken, sorry. Please contact oskar.vidarsson@uib.no for support. \n"
fi

# Print benchmarking data
FINISH=$(date +%s)
EXECTIME=$(( $FINISH-$START ))
printf "The entire workflow, including input file transfer, job execution on Colossus, and output file transfer back to the p172ncspmdata disk, took %dd:%dh:%dm:%ds\n" $((EXECTIME/86400)) $((EXECTIME%86400/3600)) $((EXECTIME%3600/60)) $((EXECTIME%60))
