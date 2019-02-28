# usr/bin/python
import os
import sys
import shutil

# CLI arguments
ref_dict = sys.argv[1]
dir_name = sys.argv[2]

with open(ref_dict, "r") as ref_dict_file:
    sequence_tuple_list = []
    longest_sequence = 0
    for line in ref_dict_file:
        if line.startswith("@SQ"):
            line_split = line.split("\t")
            # (Sequence_Name, Sequence_Length)
            sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
    longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
# We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in 
# some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
hg38_protection_tag = ":1+"
# initialize the tsv string with the first sequence
tsv_string = "--intervals " + sequence_tuple_list[0][0] + hg38_protection_tag
temp_size = sequence_tuple_list[0][1]
for sequence_tuple in sequence_tuple_list[1:]:
    if temp_size + sequence_tuple[1] <= longest_sequence:
        temp_size += sequence_tuple[1]
        tsv_string += "\t" + "--intervals " + sequence_tuple[0] + hg38_protection_tag
    else:
        tsv_string += "\n" + "--intervals " + sequence_tuple[0] + hg38_protection_tag
        temp_size = sequence_tuple[1]
# add the unmapped sequences as a separate line to ensure that they are recalibrated as well
path = dir_name
if os.path.exists(path):
    shutil.rmtree(path)
os.mkdir(path)

tsv_string += '\n' + "--intervals unmapped\n"

with open(path+"sequence_grouping_with_unmapped.tsv","w") as tsv_file_with_unmapped:
  tsv_file_with_unmapped.write(tsv_string)
  tsv_file_with_unmapped.close()
