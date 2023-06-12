from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import seqhash
import pandas as pd
import csv

#sets wildcards
#SAMPLES, = glob_wildcards("{sample}.fna")

# rule merge_fasta:
#     input:
#         expand(config["outdir"] + "/results/ORFs/{sample}_ORFs.faa", sample=SAMPLES)
#     output:
#         fasta = protected(config["outdir"] + "/results/ORFs/all_ORFs.faa")
#     run:
#         for file in input:
#             # Open the input text file in read mode
#             with open(file, 'r') as input_file:
#                 # Read the lines from the input file
#                 lines = input_file.readlines()
#                 # Open the output file in append mode
#                 with open(output.fasta, 'a') as output_file:
#                     # Write each line from the input file to the output file
#                     for line in lines:
#                         output_file.write(line)

# rule fasta_to_csv:
#     input:
#         #rules.ORF_filter.output[0]
#         fasta=config["outdir"] + "/results/ORFs/{sample}_ORFs_filtered.faa"
#     output:
#         csv=config["outdir"] + "/results/ORFs/{sample}_ORFs_filtered.csv"
#     priority: 50
#     threads: 1
#     run:
#         for fasta_file, csv_file in zip(input.fasta, output.csv):
#             with open(fasta_file) as handle:
#                 sample = str(fasta_file).split(config["outdir"] + "/results/ORFs/")[1].strip(".faa")
#                 # Create an empty dictionary
#                 idDict = {}

#                 for seq_record in SeqIO.parse(handle, "fasta"):
#                     # pull out the sequence for dna and protein
#                     sequence = str(seq_record.seq)
#                     protein = str(seq_record.seq.transla2023-05-26T164013.227395.snakemake.logte())    
#                     pephash = seqhash.seqhash(protein, dna_type='PROTEIN')
#                     # separate out the description
#                     descriptionParts = seq_record.description.split(" ")
#                     start = descriptionParts[1].split("]")[0].split("-")[0].strip("[")
#                     stop = descriptionParts[1].split("]")[0].split("-")[1]
#                     strand = descriptionParts[1].split("]")[1].strip("()")
#                     contig = descriptionParts[0].strip(">")

#                     seqID = f"{sample}|{contig}|{start}:{stop}:{strand}"
#                     idDict[seqID] = [pephash, sample, contig, start, stop, strand, sequence, protein]
#                 idDF = pd.DataFrame.from_dict(idDict, orient="index", columns=["pephash","sample","contig","start","stop","strand","dna","seq"]).reset_index()
#                 idDF.rename(columns={'index': 'id'}, inplace=True)
#                 idDF.to_csv(csv_file, index=None)

# rule merge_csv:
#     input:
#         expand(config["outdir"] + "/results/ORFs/{sample}_ORFs.csv", sample=SAMPLES)
#     output:
#         csv = protected(config["outdir"] + "/results/ORFs/all_ORFs.csv")
#     run:
#         # Open the output file in append mode
#         with open(output.csv, 'a') as output_file:
#             output_file.write("id,pephash,sample,contig,start,stop,strand,dna,seq\n")

#             # Open the input csv file in read mode
#             for file in input:
#                 with open(file, 'r') as input_file:
#                     # Create a csv reader object
#                     reader = csv.reader(input_file)
#                     # Skip the header row
#                     next(reader)
#                     # Create a csv writer object
#                     writer = csv.writer(output_file)
#                     # Write each row from the input file to the output file
#                     for row in reader:
#                         writer.writerow(row)

# rule dna_to_protein:
#     input:
#         #rules.ORF_filter.output[0]
#         config["outdir"] + "/results/ORFs/{sample}_ORFs_filtered.faa"
#     output:
#         config["outdir"] + "/results/ORFs/{sample}_proteinSeq_filtered.faa"
#     priority: 51
#     threads: 1
#     run:
#         def dna_to_protein(input_file, output_file):
#             # Read the DNA sequences from the input FASTA file
#             dna_records = list(SeqIO.parse(input_file, "fasta"))

#             # Translate DNA sequences into protein sequences
#             protein_records = []
#             for dna_record in dna_records:
#                 # Translate the DNA sequence into protein sequence
#                 protein_seq = dna_record.seq.translate()
#                 # Create a new record with the same ID and description as the original one
#                 protein_record = SeqIO.SeqRecord(protein_seq, id=dna_record.id, description=dna_record.description)
#                 protein_records.append(protein_record)

#             # Write the protein sequences to the output FASTA file
#             SeqIO.write(protein_records, output_file, "fasta")

#         # Provide input and output file paths
#         input_file = input[0]  # Path to the input DNA FASTA file
#         output_file = output[0]  # Path to save the output protein FASTA file

#         # Convert DNA to protein and save the result
#         dna_to_protein(input_file, output_file)

rule merge_protein_fasta:
    input:
        expand(config["outdir"] + "/results/ORFs/{sample}_aminoAcid.faa", sample=SAMPLES)
    output:
        fasta = config["outdir"] + "/results/ORFs/all_AAs.faa"
    run:
        for file in input:
            # Open the input text file in read mode
            with open(file, 'r') as input_file:
                # Read the lines from the input file
                lines = input_file.readlines()
                # Open the output file in append mode
                with open(output.fasta, 'a') as output_file:
                    # Write each line from the input file to the output file
                    for line in lines:
                        output_file.write(line)

# rule protein_ORF:
#     input:
#         dna = config["outdir"] + "/results/ORFs/ORFs_filtered_all.faa"
#     output:
#         fasta = protected(config["outdir"] + "/results/ORFs/proteins_filtered_all.faa")
#     threads:
#         workflow.cores * 0.9
#     script:
#         "../scripts/dna_to_protein.py"

# rule protein_ORF:
# 	input:
# 		dna = config["outdir"] + "/results/ORFs/ORFs_filtered_all.faa"
# 	output:
# 		fasta = config["outdir"] + "/results/ORFs/proteins_filtered_all.faa"
# 	run:
# 		# Open the output file for writing
# 		with open(input.dna, 'r') as input_file:
# 		    with open(output.fasta, "w") as out_handle:
#         	    # Loop over the records in the input file
#                 for record in SeqIO.parse(input_file, "fasta"):
#                     # Get the DNA sequence as a Seq object
#                     dna_seq = record.seq
#                     # Translate the DNA sequence into protein sequence using the standard genetic code
#                     protein_seq = dna_seq.translate()
#                     # Create a new record with the same id and description as the original one
#                     protein_record = SeqIO.SeqRecord(protein_seq, id=record.id, description=record.description)
#                     # Write the protein record to the output file in fasta format
#                     SeqIO.write(protein_record, out_handle, "fasta")