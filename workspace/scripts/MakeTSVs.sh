LINES=$(wc -l "${1}" | awk '{ print $1 }')

for line in $(seq 1 "${LINES}"); do
	echo $line
	awk -v TI="$line" "NR==TI" "${1}" > "${2}"/contigs_"${line}".bed
done
