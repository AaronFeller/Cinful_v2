"""
Script finds all ORFs inside of a genome
"""

from Bio import SeqIO

# Find all potential open reading frames in both strands.
smk = snakemake #type: ignore
genome_path = smk.input[0]
output_path = smk.output[0]

genome = SeqIO.parse(genome_path, "fasta")

noncanonical = ['R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N']

gencode = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'STOP', 'TAG':'STOP',
    'TGC':'C', 'TGT':'C', 'TGA':'STOP', 'TGG':'W' }

ORF_COUNT = 0
ORF = ""
strands = ["(+)", "(-)"]

with open(output_path, 'w', encoding="UTF-8") as output_file:
    for rec in genome:
        seq = rec.seq
        name = rec.description.split(" ")
        ORF_type = name[5][0:-1]

        for strand in strands:
            if strand == "(-)": # looking in the (-) strand, sequence is reverse complement
                seq = seq.reverse_complement()
            for i in range(0,3):   # 3 possible frames
                FOUND_START = False  # Initally, we have not found a start
                START_LOC = 0
                END_LOC = 0
                for j in range(i, len(seq), 3):  # moving along the genome sequence every codon
                    codon = seq[j:j+3]

                    # if 'N' in codon:
                    #     continue

                    if len(codon) == 3:

                        if FOUND_START is False and any(x in noncanonical for x in codon):
                            continue

                        if FOUND_START is True and any(x in noncanonical for x in codon):
                            FOUND_START = False
                            END_LOC = j + 3
                            LOCATION = "[" + str(START_LOC) + "-" + str(END_LOC) + "]"
                            LENGTH = len(ORF)
                            START = 'ATG'
                            STOP = 'Noncanonical base in sequencing'
                            frame = i + 1
                            continue

                        aa = gencode[codon]

                        if aa == 'M' and FOUND_START is False:  # If we found a start
                            FOUND_START = True
                            START_LOC = j + 1
                            ORF += 'M'
                            ORF_COUNT += 1

                        elif FOUND_START is True and aa != 'STOP': # Translating/elongating our ORF
                            ORF += aa

                        elif FOUND_START is True and aa == 'STOP':  #if we run into stop after start
                            FOUND_START = False
                            END_LOC = j + 3
                            LOCATION = "[" + str(START_LOC) + "-" + str(END_LOC) + "]"
                            LENGTH = len(ORF)
                            START = 'ATG'
                            STOP = codon
                            frame = i + 1

                            LABEL = str(">" + name[0] + "_ORF." + str(ORF_COUNT) + " " + LOCATION +
                                strand + " type:" + ORF_type + " length:" + str(LENGTH) +
                                " frame:" + str(frame) + " start:" + START + " stop:" + STOP)
                            output_file.write(LABEL + "\n")
                            output_file.write(ORF + "\n")
                            ORF = ""

                    else:
                        break

print("Finished extracting ORFs.")
