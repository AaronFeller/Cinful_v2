from collections import defaultdict
from Bio import SearchIO
import pandas as pd

# Extract hits from the hmmsearch output
rule extract_hits:
    input:
        tblout = config["outdir"] + "/results/temp/hmmsearch/{sample}_microcin_hmmsearch_tblout.txt"
    output:
        extracted_hits = temp(config["outdir"] + "/results/temp/hmmsearch/{sample}_hmm_hits.csv")
    run:
        attribs = ['id', 'bitscore', 'evalue']
        hits = defaultdict(list)

        with open(input.tblout) as handle:
            verify = False
            lines = handle.readlines()
            # Check if there are any lines representing hits
            for line in lines:
                if not line.startswith('#'):
                    verify = True
                    break
            if verify:
                # reset the file pointer to the beginning of the file
                handle.seek(0)
                for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
                    if queryresult.hits:
                        for hit in queryresult.hits:
                            for attrib in attribs:
                                hits[attrib].append(getattr(hit, attrib))
        
                extracted_df = pd.DataFrame.from_dict(hits)
                # OLD: split id for each row into seqID, pephash, dnahash, sample, contig, start, stop, strand, str(dna), orf
                extracted_df[['seqID', 'location', 'sample', 'dna_length', 'protein_length', 'start', 'stop', 'pephash', 'dnahash', 'dna', 'seq']] = extracted_df['id'].str.split('|', expand=True)
                extracted_df['seqID'] = extracted_df['seqID'].str.split('|').str[0]
                extracted_df.to_csv(output.extracted_hits, index=False)
            else:
                extracted_df = pd.DataFrame(columns=['seqID', 'location', 'sample', 'dna_length', 'protein_length', 'start', 'stop', 'pephash', 'dnahash', 'dna', 'seq'])
                extracted_df.to_csv(output.extracted_hits, index=False)

# Merge all the hits into a single file
rule merge_hits:
    input:
        extracted_hits = expand(config["outdir"] + "/results/temp/hmmsearch/{sample}_hmm_hits.csv", sample=SAMPLES)
    output:
        merged_hits = config["outdir"] + "/results/temp/hmmsearch/hmm_hits.csv"
    run:
        hits = pd.concat([pd.read_csv(hit) for hit in input.extracted_hits if os.path.getsize(hit) > 20])
        hits.to_csv(output.merged_hits, index=False)
        # if the file is empty, write an empty file
        if os.path.getsize(output.merged_hits) <= 20:
            with open(output.merged_hits, 'w') as f:
                f.write('seqID,location,sample,dna_length,protein_length,start,stop,pephash,dnahash,dna,seq\n')


# Extract hits from the hmmsearch output
rule extract_hits_domain:
    input:
        domtblout = config["outdir"] + "/results/temp/hmmsearch/{sample}_microcin_hmmsearch_domtblout.txt"
    output:
        extracted_hits = temp(config["outdir"] + "/results/temp/hmmsearch/{sample}_extracted_hits_domain.csv")
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

# Merge all the hits into a single file
rule merge_hits_domain:
    input:
        extracted_hits = expand(config["outdir"] + "/results/temp/hmmsearch/{sample}_extracted_hits_domain.csv", sample=SAMPLES)
    output:
        merged_hits = config["outdir"] + "/results/temp/hmmsearch/extracted_hits_domain.csv"
    run:
        # concatenate all the files that are not empty
        hits = pd.concat([pd.read_csv(hit) for hit in input.extracted_hits if os.path.getsize(hit) > 20])
        hits.to_csv(output.merged_hits, index=False)
        # if the file is empty, write an empty file
        if os.path.getsize(output.merged_hits) <= 20:
            with open(output.merged_hits, 'w') as f:
                f.write('hit_id,evalue,hit_start,query_start\n')
