from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import seqhash
import pandas as pd
from csv import writer


# def hasAllStandardAA(seq, alphabet="ACDEFGHIKLMNPQRSTVWY",ignore="*"):
# 	return (set(seq) - set(alphabet+ignore)) == set()


# SAMPLES, = glob_wildcards("../resources/genomes/{sample}.fna")
# if SAMPLES == []:
# 	SAMPLES, = glob_wildcards(config["outdir"] + "/results/prodigal/{sample}.faa")


rule merge_ORF_Finder:
    input:
        expand(config["outdir"] + "/results/ORFs/{sample}_ORFs_filtered.faa", sample=SAMPLES)
    output:
        fasta = config["outdir"] + "/results/ORFs/ORFs_filtered_all.faa",
        csv = config["outdir"] + "/results/ORFs/ORFs_filtered_all.csv"
    run:
        #idDict = {}
        hashDict = {}
        csvDF = pd.DataFrame.from_dict(idDict, orient="index", columns=["pephash","sample","contig","start","stop","strand","dna","seq"]).reset_index()
        csvDF.rename(columns = {'index': 'cinful_id'}, inplace = True)
        csvDF.to_csv(output.csv, index = None)
        fastaDF = pd.DataFrame()
        fastaDF.to_csv(output.fasta, index=False, header=0)

        with open(output.csv, "a") as csv_file, open(output.fasta, "w") as fasta_file:
            for file in input:
                sample = file.split(config["outdir"] + "/results/ORFs/")[1].strip(".faa")

                with open(file) as handle:

                    for seq_record in SeqIO.parse(handle, "fasta"):
                        sequence = str(seq_record.seq)
                        protein = str(seq_record.seq.translate())    
                        pephash = seqhash.seqhash(sequence,dna_type='PROTEIN')
                        hashDict[pephash] = sequence
                        descriptionParts = seq_record.description.split(" ")
                        start = descriptionParts[1].split("]")[0].split("-")[0].strip("[")
                        stop = descriptionParts[1].split("]")[0].split("-")[1]
                        strand = descriptionParts[1].split("]")[1].strip("()")
                        contig = descriptionParts[0].strip(">")
                        #allStandardAA = hasAllStandardAA(sequence)
                        seqID = f"{sample}|{contig}|{start}:{stop}:{strand}"
                        csv_input = [pephash, sample, contig, start, stop, strand, sequence, protein]
                        csv_writer = writer(csv_file)
                        csv_writer.writerow(csv_input)

            for pephash in hashDict:
                outRecord = SeqRecord(
                    Seq(hashDict[pephash]),
                    id=pephash,
                    description=""
                    )
                SeqIO.write(outRecord, fasta_file, "fasta")

        # idDF = pd.DataFrame.from_dict(idDict, orient="index", columns=["pephash","sample","contig","start","stop","strand","dna","seq"]).reset_index()
        # idDF.rename(columns={'index': 'id'}, inplace = True)
        # idDF.to_csv(output.csv, index = None)


        # for file in input:
        #     sample = file.split(config["outdir"] + "/results/ORFs/")[1].strip(".faa")
        #     with open(file) as handle:
        #         for seq_record in SeqIO.parse(handle, "fasta"):
        #             sequence = str(seq_record.seq)
        #             protein = str(seq_record.seq.translate())    
        #             pephash = seqhash.seqhash(sequence,dna_type='PROTEIN')
        #             hashDict[pephash] = sequence
        #             descriptionParts = seq_record.description.split(" ")
        #             start = descriptionParts[1].split("]")[0].split("-")[0].strip("[")
        #             stop = descriptionParts[1].split("]")[0].split("-")[1]
        #             strand = descriptionParts[1].split("]")[1].strip("()")
        #             contig = descriptionParts[0].strip(">")
        #             #allStandardAA = hasAllStandardAA(sequence)
        #             seqID = f"{sample}|{contig}|{start}:{stop}:{strand}"
        #             idDict[seqID] = [pephash, sample, contig, start, stop, strand, sequence, protein]
        # with open(output.fasta,"w") as fasta_file:
        #     for pephash in hashDict:
        #         outRecord = SeqRecord(
        #             Seq(hashDict[pephash]),
        #             id=pephash,
        #             description="")
        #         SeqIO.write(outRecord, fasta_file, "fasta")
        # idDF = pd.DataFrame.from_dict(idDict, orient="index", columns=["pephash","sample","contig","start","stop","strand","dna","seq"]).reset_index()
        # idDF.rename(columns={'index': 'id'}, inplace = True)
        # idDF.to_csv(output.csv, index = None)


rule protein_ORF:
	input:
		dna = config["outdir"] + "/results/ORFs/ORFs_filtered_all.faa"
	output:
		fasta = config["outdir"] + "/results/ORFs/proteins_filtered_all.faa"
	run:
		# Open the output file for writing
		with open(input.dna, 'r') as input_file:
		    with open(output.fasta, "w") as out_handle:
        	    # Loop over the records in the input file
                for record in SeqIO.parse(input_file, "fasta"):
                    # Get the DNA sequence as a Seq object
                    dna_seq = record.seq
                    # Translate the DNA sequence into protein sequence using the standard genetic code
                    protein_seq = dna_seq.translate()
                    # Create a new record with the same id and description as the original one
                    protein_record = SeqIO.SeqRecord(protein_seq, id=record.id, description=record.description)
                    # Write the protein record to the output file in fasta format
                    SeqIO.write(protein_record, out_handle, "fasta")