#!/bin/bash

set -Eeo pipefail
# Uncomment to enable debugging
#set -vo xtrace

source scripts/functions.sh
source settings/settings.conf

trap 'die && exec 2>&4 1>&3' ERR SIGTERM SIGINT SIGKILL

START=$(date +%s)
DATE=$(date +%F.%H.%M.%S)

# Set variables
TSV=$1
INPUTS=$2
HGVER=$3
INTERVAL=$4
OUTPUTDIR=$5
JOBNAME=$6
ACCOUNT=$7
TIME=$8
MEMPERCPU=$9
CPUSPERTASK=${10}
RMTMP=${11}
LOGFILE=${12}

# Function that gives the user troubleshooting information and attempts to clean up stray files so as not to litter the intermediary workflow directory
die () {
	printf "#########################################\n"
	err "$0 failed at line $BASH_LINENO"
	sleep 1

	delstag $STAGINGDIR $WFDIRNAME $DATE
	slcan $SLURMID
	mvslot $WFDIR $OUTPUTDIR $WFDIRNAME $SLURMID
	rmtmp $RMTMP
	benchmark $START $LOGFILE

	if [[ -f $LOGFILE ]]; then
		err "Execution stopped unexpectedly, here's the logfile: $LOGFILE"
	fi

	printf "#########################################\n"
	exec 2>&4 1>&3
	exit 1
}

# Function that gives the user useful information when the workflow has finished successfully
# This function also attempts to clean up stray files so as not to litter the intermediary workflow directory
cleanexit () {
	printf "#########################################\n"

	delstag $STAGINGDIR $WFDIRNAME $DATE
	mvslot $WFDIR $OUTPUTDIR $WFDIRNAME $SLURMID
	rmtmp $RMTMP
	benchmark $START $LOGFILE

	# All done!
	if [[ "$(ls -A $OUTPUTDIR/$WFDIRNAME-$DATE)" ]]; then
		succ "Your output files for slurm job $SLURMID are in $OUTPUTDIR/$WFDIRNAME-$DATE"
	else
		err "The output directory $OUTPUTDIR/$WFDIRNAME-$DATE for slurm job $SLURMID seems to be empty, something is broken, sorry. Please contact oskar.vidarsson@uib.no for support."
	fi

	printf "#########################################\n"
	exec 2>&4 1>&3
	exit 1
}

# Make fresh copy of the workflow directory in case multiple executions of the workflow are started simultaneously
inf "Copying $WFDIR to staging directory: $STAGINGDIR/$WFDIRNAME-$DATE\n"
rprog --exclude=.git $WFDIR/ $STAGINGDIR/$WFDIRNAME-$DATE

# Copy TSV file to $STAGINGDIR/$WFDIR-$DATE/workspace/samples.tsv
inf "Copying $TSV to $STAGINGDIR/$WFDIRNAME-$DATE/workspace/samples.tsv\n"
rprog $TSV $STAGINGDIR/$WFDIRNAME-$DATE/workspace/samples.tsv

# Say that the input files/folder is being copied to $STAGINGDIR/$WFDIRNAME-$DATE/workspace/ on the /cluster disk
inf "Copying input files to $STAGINGDIR/$WFDIRNAME-$DATE/workspace/ this may take a while if your input files are large.\n"
mapfile -t FILES < <(awk '{ print $5, $6}' $TSV | tail -n +2)
rprog $(echo " ${FILES[@]}" | sed -e 's, , '"$INPUTS"'/,g') "$STAGINGDIR/$WFDIRNAME-$DATE"/workspace/

# Add time stamp to make event trigger files unique
STAMP=$(date +%s)

# Submit job to slurm
SLURMID=$(sbatch \
-A $ACCOUNT \
-c $CPUSPERTASK \
-J $JOBNAME \
-t $TIME \
--mem-per-cpu=$MEMPERCPU \
--export=\
WFOUTPUTS=$WFOUTPUTS,\
HGVER=$HGVER,\
REFERENCES=$REFERENCES,\
INTERVAL=$INTERVAL,\
WFDIR=$STAGINGDIR/$WFDIRNAME-$DATE,\
STAMP=$STAMP \
scripts/RunOnNode.sbatch | awk '{ print $4 }')

# Check if job was submitted to slurm successfully
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
JOBSUCCESS=$STAGINGDIR/$WFDIRNAME-$DATE/job-success-$STAMP.txt
JOBFAILED=$STAGINGDIR/$WFDIRNAME-$DATE/job-failed-$STAMP.txt
INPTRINI=$STAGINGDIR/$WFDIRNAME-$DATE/input-file-transfer-ini-$STAMP.txt
INPTRSUC=$STAGINGDIR/$WFDIRNAME-$DATE/input-file-transfer-succ-$STAMP.txt
SINGSTART=$STAGINGDIR/$WFDIRNAME-$DATE/sing-start-$STAMP.txt
SINGDONE=$STAGINGDIR/$WFDIRNAME-$DATE/sing-done-$STAMP.txt

#######################
# The following if else statements check the slurm job status and informs the user what the current slurm job status is
while true; do
	if [[ -e $JOBFAILED && ! -e $JOBSUCCESS ]]; then
		err "The job failure trigger for slurm job $SLURMID has been detected, will now remove the intermediate files and stop slurm job $SLURMID"
		die

	elif [[ $IPININOCHECK == 0 && -e $INPTRINI ]]; then
		inf "Input files transfer to Colossus for slurm job $SLURMID initiated, this may take a while if the input files are large"
		IPININOCHECK=1

	elif [[ $IPSUCCNOCHECK == 0 && -e $INPTRSUC ]]; then
		inf "Input files for slurm job $SLURMID have been successfully transfered to Colossus, starting data analysis"
		IPSUCCNOCHECK=1

	elif [[ $SGSTARTNOCHECK == 0 && -e $SINGSTART ]]; then
		inf "Data analysis for slurm job $SLURMID has started"
		SGSTARTNOCHECK=1

	elif [[ $SGNOCHECK == 0 && -e $SINGDONE ]]; then
		succ "Data analysis for slurm job $SLURMID has finished, initiating output file transfer to staging directory"
		SGNOCHECK=1

	elif [[ -e $JOBSUCCESS ]]; then
		succ "File transfer trigger from slurm job $SLURMID has been detected, the file transfer has completed."
		sleep 1
		inf "Will now attempt to copy output files from $STAGINGDIR/$WFDIRNAME-$DATE/workspace/Outputs/ to ${OUTPUTDIR}/$WFDIRNAME-$DATE\n for slurm job $SLURMID"
		sleep 5
		break
	fi

	# Run squeue, store the current job state code and use it to give information to the user as well as act according to the status
	# Use "|| true" because sometimes the squeue command fails and that kills the workflow prematurely
	STATUS=$(squeue -h -t all -j $SLURMID -o %t) || true #2>&1)

	# Begin checking the $STATUS variable by checking if the job has been canceled
	if [[ $STATUS == CA ]]; then
		err "Slurm job $SLURMID has been canceled, will now remove intermediary files"
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
		inf "Slurm job $SLURMID is pending"

	# CF means that the execute node is being configured
	elif [[ $STATUS == CF ]]; then
		inf "Execute node for slurm job $SLURMID is being configured"

	# R means that the sbatch script has started running, it does not mean that the singularity job has started
	# $SINGSTART means that the singularity job has started
	elif [[ $STATUS == R ]]; then
		inf "Slurm job $SLURMID is running"

	# CG means that the slurm job is completing
	elif [[ $STATUS == CG ]]; then
		inf "Slurm job $SLURMID is completing"

	# CD means that the job has completed
	elif [[ $STATUS == CD ]]; then
		inf "Slurm job $SLURMID has completed"

	# F means that the slurm job failed
	elif [[ $STATUS == F ]]; then
		err "Slurm job $SLURMID has failed, will now remove intermediary files"
		die
	fi

	PREVSTAT=$STATUS
	sleep 10
done

# Make output directory
mkdir -p ${OUTPUTDIR}/$WFDIRNAME-$DATE

# Copy the slurm file to the final output directory
inf "Copying slurm-$SLURMID.out to ${OUTPUTDIR}/$WFDIRNAME-$DATE"
rprog $WFDIR/slurm-$SLURMID.out ${OUTPUTDIR}/$WFDIRNAME-$DATE/slurm-$SLURMID.out

# Only copy the final bam file that goes into HaplotypeCaller as well as the VCF files from ApplyVqsr
rprog $STAGINGDIR/$WFDIRNAME-$DATE/workspace/Outputs/GatherBamFiles $STAGINGDIR/$WFDIRNAME-$DATE/workspace/Outputs/ApplyVqsr* ${OUTPUTDIR}/$WFDIRNAME-$DATE/ && \
succ "Files for slurm job $SLURMID have been copied from $STAGINGDIR/$WFDIRNAME-$DATE/workspace/Outputs/ to ${OUTPUTDIR}/$WFDIRNAME-$DATE"

# Run the clean exit function when everything worked as it should!
cleanexit
