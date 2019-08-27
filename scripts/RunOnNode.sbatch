#!/bin/bash
#
#SBATCH --job-name=Selma
#SBATCH --account=p172
#SBATCH --time=00-25:00:00
#SBATCH --mem-per-cpu=3750
#SBATCH --cpus-per-task=16

set -Eeou #xtrace

source /cluster/bin/jobsetup
source $WFDIR/settings/slurm-settings.conf

module purge
module load singularity/2.6.1

trap die ERR SIGINT SIGTERM

die () {
	printf "\n\n#########################################\n"
	printf "$0 failed at line $BASH_LINENO\n"
	printf "Will now attempt to stop the start-workflow.sh script\n"
	touch $WFDIR/job-failed.txt
	printf "Created job-failed.txt flag in $WFDIR\n"
	printf "#########################################\n"
	exit 0
}

WFDIRBASENAME=$(basename $WFDIR)
cd $WFDIR

# Say that the workflow is being copied over to the scratch disk
# And then create trigger that signals that the file copy was successful
printf "Copying $WFDIR directory to $SCRATCH\n"
touch $WFDIR/input-file-transfer-ini.txt
rsync -ah --progress $WFDIR $SCRATCH && \
printf "File transfer from $WFDIR to $SCRATCH was successful, creating file transfer trigger in $WFDIR to signal the file transfer success\n"
touch $WFDIR/input-file-transfer-succ.txt

## cd to the execute node to run the job from /scratch
cd $SCRATCH/$WFDIRBASENAME

# The $DECOY variable is only there to swallow the exit signal, which, as a side note, is 1 even when the workflow starts with no issues.
# If DECOY= isn't used the exit signal will be interpreted as a command, so its only purpose is to not make this happen
DECOY=$(singularity exec \
-B /cluster \
-B /tsd \
-B /net \
--home $(pwd):/srv \
-B $REFERENCES:/references \
-W $(pwd) \
singularity/snakemake-germline-tools.simg snakemake -j --config version=$HGVER interval=$INTERVAL) &

# Signal that snakemake started
touch $WFDIR/sing-start.txt
printf "Snakemake workflow started\n"

# Wait until the singularity process has exited
wait

# Create trigger that signals that the singularity job has exited
touch $WFDIR/sing-done.txt
printf "Snakemake workflow is finished\n"

# Copy files back to the output directory
printf "Copying all output files back to $WFDIR\n"
rsync -ah --progress $SCRATCH/$WFDIRBASENAME/$WFOUTPUTS/ $WFDIR/$WFOUTPUTS/ && \
printf "File transfer to $SCRATCH/$WFDIRBASENAME/$WFOUTPUTS was successful, creating file transfer trigger in $WFDIR to initiate file transfer to final storage disk\n"

# Create file transfer trigger so that the start-workflow.sh script can move the
# files from $INTERMEDSTOR/$WFDIRNAME-$DATE/$WFOUTPUTS to $OUTPUTDIR/$WFDIRNAME-$DATE
# These variables are found in the start-workflow.sh script
touch $WFDIR/job-success.txt

# All done!
printf "Slurm execution all done!\n"
