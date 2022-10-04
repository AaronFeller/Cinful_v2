from hashlib import new
from itertools import count
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import re
import string

# filter ORFs by size and find smaller frames

#Note: Alternative start codons can occur within the larger open reading frames. 
# Since we are filtering by size, we might remove all ORFs where the larger frame is too big. 
# The smaller alternative sequences might be under the size limit and we don't want to miss them.
# This script deals with this issue.


# input_path = "colicinV_Ecoli/GCF_000026285.1_ORFer.fasta"
# output_path = "colicinV_Ecoli/GCF_000026285.1_ORFer_Filtered.fasta"
input_path = snakemake.input[0]
output_path = snakemake.output[0]


orf_list = list(SeqIO.parse(input_path, "fasta"))

with open(output_path, 'w') as output_file:
    for rec in orf_list:
        
        # extracting ORF information from the label/name
        orf_sequence = str(rec.seq)
        orf_description = str(rec.description)
        split_orf = re.split(' ', orf_description)
        location = split_orf[1].strip('(+)')
        split_location = location.strip('[]').split('-')
        start_loc = int(split_location[0])
        end_loc = int(split_location[1].strip(']('))

        if len(orf_sequence) <= 150 and len(orf_sequence) > 30: #if the ORF is already within size limit, we include it.
            output_file.write(">" + orf_description + "\n")
            output_file.write(orf_sequence + "\n")

        if len(orf_sequence) > 30: # excluding ORFs that are too small
            sub_ORF_count = 0 # how many sub-ORFs I find
            for position, aa in enumerate(orf_sequence):
                if position != 0 and aa == 'M': # if "M" is found within an ORF (not the first residue)
                    sub_ORF_count += 1
                    new_start = start_loc + position + 1 # recalculate the new start (end will be the same)
                    new_location = "[" + str(new_start) + "-" + str(end_loc) + "]"

                    # next five lines create a new fasta label/name for the new sub-ORF:
                    split_orf = re.split(r'[\[\]]', orf_description)
                    split_orf[0] = split_orf[0].replace(' ', '-' + str(sub_ORF_count))
                    split_orf[1] = new_location
                    new_label = ' '.join(split_orf)
                    new_label = new_label.replace('] (', '](')

                    # getting the sub-ORF sequence:
                    new_orf = orf_sequence[position:]
                    
                    # if sub-ORF is within size limit, then add it to the new fasta:
                    if len(new_orf) <= 150 and len(new_orf) >= 30:
                        output_file.write(">" + new_label + "\n")
                        output_file.write(new_orf + "\n")

print("Finished filtering ORFs! :-)")
                

           
