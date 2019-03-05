#!/bin/bash
# This script can be used to quickly restore the workflow to a clean state, \
# it supresses any output to the terminal by default
rm -r workspace/Outputs/* workspace/.snakemake &>/dev/null