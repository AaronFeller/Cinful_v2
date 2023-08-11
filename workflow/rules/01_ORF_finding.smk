from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import seqhash
import pandas as pd
import csv

#sets wildcards
SAMPLES, = glob_wildcards("{sample}.fna")

rule ORF_finder:
	input:
		genomes = "{sample}.fna"
	output:
		temp(config["outdir"] + "/results/temp/ORFs/{sample}_aminoAcid.faa")#,
#		config["outdir"] + "/results/ORFs/{sample}_dna.faa",
#		config["outdir"] + "/results/ORFs/{sample}_ORFs.csv"
	script:	
		"../scripts/find_ORFs_initialize_on_start.py"


# rule ORF_filter:
# 	input:
# 	    config["outdir"] + "/results/ORFs/{sample}_ORFs.faa"
# 	output:
# 		temp(config["outdir"] + "/results/ORFs/{sample}_ORFs_filtered.faa")
# 	script:
# 		"../scripts/filter_ORFs.py"


rule edit_contig_names:
    input:
        genomes = config["outdir"] + "/results/temp/ORFs/{sample}_aminoAcid.faa"
    output:
        checkfiles = temp(config["outdir"] + "/results/temp/ORFs/{sample}_successfully_edited.txt")
    shell:
        """
        filename=$(basename {input.genomes}) &&
        sed -i 's/ /_/g' {input.genomes} &&
        sed -i "s/>/>{{$filename}}_/g" {input.genomes} &&
        touch {output.checkfiles}
        """

rule merge_protein_fasta:
    input:
        files = expand(config["outdir"] + "/results/temp/ORFs/{sample}_aminoAcid.faa", sample=SAMPLES),
        checkfiles = expand(config["outdir"] + "/results/temp/ORFs/{sample}_successfully_edited.txt", sample=SAMPLES)
    output:
        fasta = config["outdir"] + "/results/temp/ORFs/all_AAs.faa"
    run:
        for file in input.files:
            # Open the input text file in read mode
            with open(file, 'r') as input_file:
                # Read the lines from the input file
                lines = input_file.readlines()
                # Open the output file in append mode
                with open(output.fasta, 'a') as output_file:
                    # Write each line from the input file to the output file
                    for line in lines:
                        output_file.write(line)
