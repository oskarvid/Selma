#!/bin/bash
# count the number of lines in the input file to create a variable that is \
# used to split the input file into that many individual files.
LINES=$(wc -l "${1}" | awk '{ print $1 }')

# split the input file into individual bed files
for line in $(seq 1 "${LINES}"); do
	awk -v TI="$line" "NR==TI" "${1}" > "${2}"/contigs_"${line}".bed
done
