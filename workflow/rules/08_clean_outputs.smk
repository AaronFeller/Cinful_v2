import pandas as pd
import numpy as np

rule merge_domains_with_hits:
    input:
        hits = config["outdir"] + "/results/process_files/hmm_hits_with_pident.csv",
        ss_domains = config["outdir"] + "/results/hmmsearch/ss_domains.csv",
        gp_domains = config["outdir"] + "/results/hmmsearch/ss_gram_positive_domains.csv",
        full_domains = config["outdir"] + "/results/hmmsearch/extracted_hits_domain.csv"
    output:
        final = config["outdir"] + "/results/process_files/hmm_hits_with_domain_hits.csv"
    run:
        # Load data
        hits_df = pd.read_csv(input.hits)
        ss_df = pd.read_csv(input.ss_domains)
        gp_df = pd.read_csv(input.gp_domains)
        microcin_domain_hmm_df = pd.read_csv(input.full_domains)

        #filter hits_df where column 'reg' is greater than 1
        ss_df = ss_df[ss_df['evalue'] < 1]
        gp_df = gp_df[gp_df['evalue'] < 1]
        microcin_domain_hmm_df = microcin_domain_hmm_df[microcin_domain_hmm_df['evalue'] < 1]

        #rename evalue column to evalue_signalsequence in ss_df
        ss_df = ss_df.rename(columns={'evalue': 'evalue_signal_sequence'})
        gp_df = gp_df.rename(columns={'evalue': 'evalue_gram_positive'})
        microcin_domain_hmm_df = microcin_domain_hmm_df.rename(columns={'evalue': 'evalue_microcin_domain'})

        #rename column 'hit_id' to 'pephash' in ss_df and gp_df
        ss_df = ss_df.rename(columns={'hit_id': 'pephash'})
        gp_df = gp_df.rename(columns={'hit_id': 'pephash'})
        microcin_domain_hmm_df = microcin_domain_hmm_df.rename(columns={'hit_id': 'contig'})

        #merge left hits_df and right hmm_dfs on pephash and target name
        merged_df = pd.merge(hits_df, ss_df, how='left', on='pephash')
        merged_df = pd.merge(merged_df, gp_df, how='left', on='pephash')
        merged_df = pd.merge(merged_df, microcin_domain_hmm_df, how='left', on='contig')

        #rename columns hit_start_x to hit_start_signalsequence and hit_start_y to hit_start_gram_positive
        merged_df = merged_df.rename(columns={'hit_start_x': 'hit_start_signal_sequence', 'hit_start_y': 'hit_start_gram_positive', 'hit_start': 'hit_start_microcin_domain'})

        #rename columns query_start_x to query_start_signalsequence and query_start_y to query_start_gram_positive
        merged_df = merged_df.rename(columns={'query_start_x': 'query_start_signal_sequence', 'query_start_y': 'query_start_gram_positive', 'query_start': 'query_start_microcin_domain'})

        #write csv merged_df
        merged_df.to_csv(output.final, index=False)

rule final_filter:
    input:
        config["outdir"] + "/results/process_files/hmm_hits_with_domain_hits.csv"
    output:
        config["outdir"] + "/results/final/output.csv"
    run:
        #load df from csv of file hmm_hits_with_ss_hit.csv
        df = pd.read_csv(input[0])

        #fill empty cells with nan
        df = df.fillna(np.nan)

        # Remove all rows that have dont have 'hit_start' and have an 'evalue' < 1
        sub_df = df[(df['evalue'] < 1) | (df['hit_start_signal_sequence'].notna()) | (df['hit_start_gram_positive'].notna())]

        # Shorten seq if signal sequence hmm matched
        for index, row in sub_df.iterrows():
            if row.hit_start_signal_sequence > 0:
                position = int(row.hit_start_signal_sequence)+int(row.query_start_signal_sequence)
                while all(elem not in ['M'] for elem in row.seq[position]):
                    position -= 1
                sub_df.at[index, 'seq'] = row.seq[position:]
                
        # Shorten seq if gram positive hmm matched
        for index, row in sub_df.iterrows():
            if row.hit_start_gram_positive > 0:
                position = int(row.hit_start_gram_positive)+int(row.query_start_gram_positive)
                while all(elem not in ['M'] for elem in row.seq[position]):
                    position -= 1
                sub_df.at[index, 'seq'] = row.seq[position:]

        # remove all rows with query_start_microcin_domain > 140
        sub_df = sub_df[(sub_df['query_start_microcin_domain'] <= 5) | sub_df['query_start_microcin_domain'].isna()]

        # Shorten seq with original full sequence HMM if neither gram+ or ss matched
        for index, row in sub_df.iterrows():
            if row.hit_start_microcin_domain > 0:
                if not pd.notna(row.hit_start_signal_sequence) and not pd.notna(row.hit_start_gram_positive):
                    position = int(row.hit_start_microcin_domain)+int(row.query_start_microcin_domain)
                    while all(elem not in ['M'] for elem in row.seq[position]):
                        position -= 1
                    # Replace seq in sub_df with row.seq[position:]
                    sub_df.at[index, 'seq'] = row.seq[position:]

        # remove rows with length over 140
        sub_df = sub_df[sub_df['seq'].str.len() <= 140]

        # Write CSV
        sub_df.to_csv(output[0], index=False)
