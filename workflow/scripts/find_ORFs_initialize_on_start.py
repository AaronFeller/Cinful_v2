from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid

# Find all potential open reading frames in both strands.
genome_path = snakemake.input[0]
# This is used for BLAST search to full bacteria database
output_protein = snakemake.output[0]

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

orf_count = 1
orf = ""
dna = ""
strands = ["(+)", "(-)"]

with open (output_protein, 'w') as output_file_protein:

    for rec in genome:
        seq = rec.seq
        name = rec.description.split(" ")
        sample = name[0]

        for strand in strands:  
            if strand == "(-)": # if we are looking in the (-) strand, our sequence is the reverse complement
                seq = seq.reverse_complement()
            for i in range(0,3):   # 3 possible frames
                found_start = True  # Initally, we have not found a start
                first_codon = True # initially we're on the 'first codon'
                start = '*'
                start_loc = 0
                end_loc = 0
                for j in range(i, len(seq), 3):  # moving along the genome sequence every codon
                    codon = seq[j:j+3]
                    if first_codon == True:
                        first_codon = False
                        if codon in ['TTA', 'TAG', 'TGA']:
                            found_start = False
                        continue
                    # ensure full length codon
                    if len(codon) == 3:

                        # if we haven't found a start yet, and we find a start codon, we can start translating
                        if found_start == False and any(x in noncanonical for x in codon):
                            continue

                        if found_start == True and any(x in noncanonical for x in codon):
                            found_start = False
                            end_loc = j + 3
                            location = "[" + str(start_loc) + "-" + str(end_loc) + "]"
                            length = len(orf)
                            start = 'ATG'
                            stop = 'Noncanonical base in sequencing'
                            frame = i + 1
                            continue

                        aa = gencode[codon]
                        
                        #Include codons as a check, rather than the AA.. then add M (as all start codons will use fMet)
                        if codon in ['ATG', 'GTG', 'TTG'] and found_start == False:  # If we found a start for the first time after last ORF
                            found_start = True
                            start_loc = j + 1
                            orf += 'M'
                            dna += codon
                            start = codon
                            orf_count += 1

                        elif found_start == True and aa != 'STOP': # Translating/elongating our ORF
                            orf += aa
                            dna += codon

                        elif found_start == True and aa == 'STOP':  # if we run into stop after start
                            found_start = False
                            end_loc = j + 3
                            location = "[" + str(start_loc) + "-" + str(end_loc) + "]"
                            length = len(dna)

                            if length < 63:
                                orf = ""
                                dna = ""
                                continue

                            dna += codon
                            orf += '*'  # add stop codon to AA sequence
                            stop = codon
                            frame = i + 1

                            pephash = seguid(orf)
                            dnahash = seguid(dna)
                            
                            label = (
                                        f">{name[0]}_ORF.{orf_count}|{location}{strand}|{sample}|{length}|"
                                        f"{int(int(length)/3)}|{start}|{stop}|{pephash}|{dnahash}|{str(dna)}|{str(orf)}"
                                    )

                            output_file_protein.write(label + "\n")
                            output_file_protein.write(str(orf) + "\n")
                            
                            #reset variables for next ORF
                            orf = ""
                            dna = ""
                            
                    else:
                        break

print("Finished extracting ORFs.")
