from collections import defaultdict
from Bio import SearchIO
import pandas as pd
import re
import os


rule extract_hits:
    input:
        tblout = config["outdir"] + "/results/hmmsearch/microcin_hmmsearch_tblout.txt"
    output:
        extracted_hits = config["outdir"] + "/results/hmmsearch/extracted_hits.csv"
    run:
        attribs = ['id', 'bitscore', 'evalue']
        hits = defaultdict(list)

        with open(input.tblout) as handle:
            for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
                for hit in queryresult.hits:
                    for attrib in attribs:
                        hits[attrib].append(getattr(hit, attrib))
        
        extracted_df = pd.DataFrame.from_dict(hits)
        extracted_df.to_csv(output.extracted_hits, index=False)

rule extract_hits_domain:
    input:
        domtblout = config["outdir"] + "/results/hmmsearch/microcin_hmmsearch_domtblout.txt"
    output:
        extracted_hits = config["outdir"] + "/results/hmmsearch/extracted_hits_domain.csv"
    run:
        attributes = ['hit_id', 'evalue', 'hit_start', 'query_start']
        hits = defaultdict(list)

        with open(input.domtblout) as handle:
            for queryresult in SearchIO.parse(handle, 'hmmsearch3-domtab'):
                for hit in queryresult.hits:
                    for hsp in hit:
                        for attrib in attributes:
                            hits[attrib].append(getattr(hsp, attrib))
        extracted_df = pd.DataFrame.from_dict(hits)
        extracted_df.to_csv(output.extracted_hits, index=False)

rule merge_hmmscan:
    input:
        extracted_hits = config["outdir"] + "/results/hmmsearch/extracted_hits.csv"
    output:
        hmm_hits = config["outdir"] + "/results/hmmsearch/hmm_hits.csv"
    run:
        # Load hits
        extracted_df = pd.read_csv(input.extracted_hits).rename(columns={'id': 'contig'})

        # Add blank columns
        added_columns = ['id', 'pephash', 'dnahash', 'sample', 'start', 'stop', 'strand', 'dna', 'seq']
        for column in added_columns:
            extracted_df[column] = ''

        # Iterate through each row of extracted_df
        for i in range(len(extracted_df)):
            # Find the genome that the contig belongs to
            filename = extracted_df.loc[i, 'contig'].split('_aminoAcid.faa')[0].strip('{') + "_ORFs.csv"
            genome = extracted_df.loc[i, 'contig'].split('_')[0].strip('{') + '_' + extracted_df.loc[i, 'contig'].split('_')[1]

            # Generate the string for ORF file for the genome found
            orf_file = ""

            # Check the first file
            file_path = config["outdir"] + "/results/ORFs/" + filename
            if os.path.isfile(file_path):
                orf_file = file_path

            # If the first file doesn't exist, check the next one
            if not orf_file:
                file_path = config["outdir"] + "/results/ORFs/genbank/bacteria/" + genome + "/" + filename
                if os.path.isfile(file_path):
                    orf_file = file_path

            # If the second file doesn't exist, check the last one
            if not orf_file:
                file_path = config["outdir"] + "/results/ORFs/refseq/bacteria/" + genome + "/" + filename
                if os.path.isfile(file_path):
                    orf_file = file_path

            # Read ORF file and assign it to orfs_df
            orfs_df = pd.read_csv(orf_file)

            # Find the row in orfs_df that matches the 'contig' (should be only 1)
            matching_row = orfs_df[orfs_df['contig'] == extracted_df.loc[i, 'contig'].split('}_')[1].split('_[')[0]]

            # Add values from columns to extracted_df
            for column in added_columns:
                extracted_df.loc[i, column] = matching_row.iloc[0][column]

        # Write the final completed DataFrame to the output file
        extracted_df.to_csv(output.hmm_hits, index=False)

