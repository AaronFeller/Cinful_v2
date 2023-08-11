from collections import defaultdict
from Bio import SearchIO
import pandas as pd
import re
import os


rule extract_hits:
    input:
        tblout = config["outdir"] + "/results/temp/hmmsearch/microcin_hmmsearch_tblout.txt"
    output:
        extracted_hits = config["outdir"] + "/results/temp/hmmsearch/hmm_hits.csv"
    run:
        attribs = ['id', 'bitscore', 'evalue']
        hits = defaultdict(list)

        with open(input.tblout) as handle:
            for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
                for hit in queryresult.hits:
                    for attrib in attribs:
                        hits[attrib].append(getattr(hit, attrib))
        
        extracted_df = pd.DataFrame.from_dict(hits)

        # split id for each row into seqID, pephash, dnahash, sample, contig, start, stop, strand, str(dna), orf
        extracted_df[['seqID', 'location', 'sample', 'dna_length', 'protein_length', 'start', 'stop', 'pephash', 'dnahash', 'dna', 'seq']] = extracted_df['id'].str.split('|', expand=True)
        # seqID trimming after first '|'
        extracted_df['seqID'] = extracted_df['seqID'].str.split('|').str[0]

        extracted_df.to_csv(output.extracted_hits, index=False)


rule extract_hits_domain:
    input:
        domtblout = config["outdir"] + "/results/temp/hmmsearch/microcin_hmmsearch_domtblout.txt"
    output:
        extracted_hits = config["outdir"] + "/results/temp/hmmsearch/extracted_hits_domain.csv"
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
