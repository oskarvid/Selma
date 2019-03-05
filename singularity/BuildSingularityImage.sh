# Run this script from the workflow base directory, not the directory containing this script.
# It builds the singularity image from the docker image that is stored on docker hub.

singularity build \
singularity/snakemake-germline-tools.simg \
docker://oskarv/snakemake-germline-tools
