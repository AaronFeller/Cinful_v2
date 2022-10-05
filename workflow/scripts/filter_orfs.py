"""
This script will filter the output from find_orfs.py
to pull out sub-orfs contained inside full-size orfs.
"""

import re
from Bio import SeqIO

# Note: Alternative start codons can occur within the larger open reading frames.

# Since we are filtering by size, we might remove all ORFs where the larger frame is too big.
# The smaller alternative sequences might be under the size limit and we don't want to miss them.
# This script deals with this issue.

smk = snakemake #type: ignore
input_path = smk.input[0]
output_path = smk.output[0]

orf_list = list(SeqIO.parse(input_path, "fasta"))

with open(output_path, 'w', encoding="UTF-8") as output_file:
    for rec in orf_list:

        # extracting ORF information from the label/name
        ORF_SEQUENCE = str(rec.seq)
        ORF_DESCRIPTION = str(rec.description)
        split_orf = re.split(' ', ORF_DESCRIPTION)
        location = split_orf[1].strip('(+)')
        split_location = location.strip('[]').split('-')
        start_loc = int(split_location[0])
        end_loc = int(split_location[1].strip(']('))

        if len(ORF_SEQUENCE) <= 150 and len(ORF_SEQUENCE) > 30: #if ORF in size limit, include it
            output_file.write(">" + ORF_DESCRIPTION + "\n")
            output_file.write(ORF_SEQUENCE + "\n")

        if len(ORF_SEQUENCE) > 30: # excluding ORFs that are too small
            SUB_ORF_COUNT = 0 # how many sub-ORFs I find
            for position, aa in enumerate(ORF_SEQUENCE):
                if position != 0 and aa == 'M': # if "M" found within an ORF (not first residue)
                    SUB_ORF_COUNT += 1
                    new_start = start_loc + position + 1 # recalculate the new start (same end)
                    NEW_LOCATION = "[" + str(new_start) + "-" + str(end_loc) + "]"

                    # next five lines create a new fasta label/name for the new sub-ORF:
                    split_orf = re.split(r'[\[\]]', ORF_DESCRIPTION)
                    split_orf[0] = split_orf[0].replace(' ', '-' + str(SUB_ORF_COUNT))
                    split_orf[1] = NEW_LOCATION
                    NEW_LABEL = ' '.join(split_orf)
                    NEW_LABEL = NEW_LABEL.replace('] (', '](')

                    # getting the sub-ORF sequence:
                    new_orf = ORF_SEQUENCE[position:]

                    # if sub-ORF is within size limit, then add it to the new fasta:
                    if len(new_orf) <= 150 and len(new_orf) >= 30:
                        output_file.write(">" + NEW_LABEL + "\n")
                        output_file.write(new_orf + "\n")

print("Finished filtering ORFs.")
