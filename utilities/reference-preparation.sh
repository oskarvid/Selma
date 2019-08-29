#!/bin/bash

if [[ ! $(which tabix) ]]; then
	echo "tabix is missing"
	exit
elif [[ ! $(which samtools) ]]; then
	echo "samtools is missing"
	exit
elif [[ ! $(which bwa) ]]; then
	echo "bwa is missing"
	exit
fi
# Bgzip all uncompressed vcf files
for file in *.vcf.gz; do
#	echo "Running bgzip on $file"
#	bgzip -c $file > $file.gz
#	echo "Running tabix on $file.gz"
	echo "Running tabix on $file"
	tabix -p vcf $file
#	rm $file
done

# Unpack all compressed vcf files
echo "Unpacking compressed FASTA file(s)"
for file in *.fasta.gz; do
	echo "Unpacking $file"
	zcat $file > $(basename $file .gz)
	rm $file
done

# Index the fasta files, it's a loop because b37 has a decoy and non decoy version
FASTA=(*.fasta)
for fasta in ${FASTA[@]}; do
	echo "Running bwa index on $fasta"
	bwa index -a bwtsw $fasta

	echo "Running samtools idx on $fasta"
	samtools faidx $fasta
done


echo "All done!"
