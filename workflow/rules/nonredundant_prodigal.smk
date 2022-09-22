from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import seqhash
import pandas as pd


def hasAllStandardAA(seq, alphabet="ACDEFGHIKLMNPQRSTVWY",ignore="*"):
	return (set(seq) - set(alphabet+ignore)) == set()


SAMPLES, = glob_wildcards("../../resources/genomes/{sample}.fna")
if SAMPLES == []:
	SAMPLES, = glob_wildcards(config["outdir"] + "/results/prodigal/{sample}.faa")


rule nonredundant_prodigal:
    input:
        expand(config["outdir"] + "/results/prodigal/{sample}.faa", sample=SAMPLES)
    output:
        fasta = config["outdir"] + "/results/prodigal/prodigal_out_all_nr.faa",
        csv = config["outdir"] + "/results/prodigal/prodigal_out_all_nr_expanded.csv"
    run:
        hashDict = {}
        idDict = {}
        print("INPUT:",input)
        for file in input:
            sample = file.split(config["outdir"] + "/results/prodigal/")[1].strip(".faa")
            with open(file) as handle:
                for seq_record in SeqIO.parse(handle, "fasta"):
                    sequence = str(seq_record.seq)
                    pephash = seqhash.seqhash(sequence.strip("*"),dna_type='PROTEIN')
                    hashDict[pephash] = sequence
                    descriptionParts = seq_record.description.split("#")
                    start = descriptionParts[1].strip()
                    stop = descriptionParts[2].strip()
                    strand = descriptionParts[3].strip()
                    contig = '_'.join(seq_record.id.split("_")[:-1])
                    allStandardAA = hasAllStandardAA(sequence)
                    seqID = f"{sample}|{contig}|{start}:{stop}:{strand}"
                    idDict[seqID] = [pephash, sample, contig, start, stop, strand, allStandardAA, sequence]
        with open(output.fasta,"w") as fasta_file:
            for pephash in hashDict:
                outRecord = SeqRecord(
                    Seq(hashDict[pephash]),
                    id=pephash,
                    description="")
                SeqIO.write(outRecord, fasta_file, "fasta")

        idDF = pd.DataFrame.from_dict(idDict, orient="index", columns=["pephash","sample","contig","start","stop","strand","allStandardAA","seq"]).reset_index()
        idDF.rename(columns={'index': 'cinful_id'}, inplace = True)
        idDF.to_csv(output.csv, index = None)

