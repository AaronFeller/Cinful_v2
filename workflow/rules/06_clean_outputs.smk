import pandas as pd
import numpy as np

rule merge_domains_with_hits:
    input:
        hits = config["outdir"] + "/results/temp/process_files/hmm_hits_with_pident.csv",
        ss_domains = config["outdir"] + "/results/temp/hmmsearch/ss_domains.csv",
        gp_domains = config["outdir"] + "/results/temp/hmmsearch/ss_gram_positive_domains.csv",
        full_domains = config["outdir"] + "/results/temp/hmmsearch/extracted_hits_domain.csv"
    output:
        final = config["outdir"] + "/results/temp/process_files/hmm_hits_with_domain_hits.csv"
    run:
#Load hits
        hits_df = pd.read_csv(input.hits)

# Load microcin signal sequence domains from HMM
        try:
            ss_df = pd.read_csv(input.ss_domains)
        except pd.errors.EmptyDataError:
            ss_df = pd.DataFrame()
            #add 'hit_id' 'evalue' 'hit_start' & 'query_start'
            ss_df['hit_id'] = np.nan
            ss_df['evalue'] = np.nan
            ss_df['hit_start'] = np.nan
            ss_df['query_start'] = np.nan

# Load gram positive domains from HMM
        try:
            gp_df = pd.read_csv(input.gp_domains)
        except pd.errors.EmptyDataError:
            gp_df = pd.DataFrame()
            #add 'hit_id' 'evalue' 'hit_start' & 'query_start'
            gp_df['hit_id'] = np.nan
            gp_df['evalue'] = np.nan
            gp_df['hit_start'] = np.nan
            gp_df['query_start'] = np.nan

# Load microcin domains from HMM
        try:
            microcin_domain_hmm_df = pd.read_csv(input.full_domains)
        except pd.errors.EmptyDataError:
            microcin_domain_hmm_df = pd.DataFrame()
            #add 'hit_id' 'evalue' 'hit_start' & 'query_start'
            microcin_domain_hmm_df['hit_id'] = np.nan
            microcin_domain_hmm_df['evalue'] = np.nan
            microcin_domain_hmm_df['hit_start'] = np.nan
            microcin_domain_hmm_df['query_start'] = np.nan

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
        microcin_domain_hmm_df = microcin_domain_hmm_df.rename(columns={'hit_id': 'id'})

        #merge left hits_df and right hmm_dfs on pephash and target name
        merged_df = pd.merge(hits_df, ss_df, how='left', on='pephash')
        merged_df = pd.merge(merged_df, gp_df, how='left', on='pephash')
        merged_df = pd.merge(merged_df, microcin_domain_hmm_df, how='left', on='id')

        #rename columns hit_start_x to hit_start_signalsequence and hit_start_y to hit_start_gram_positive
        merged_df = merged_df.rename(columns={'hit_start_x': 'hit_start_signal_sequence', 'hit_start_y': 'hit_start_gram_positive', 'hit_start': 'hit_start_microcin_domain'})

        #rename columns query_start_x to query_start_signalsequence and query_start_y to query_start_gram_positive
        merged_df = merged_df.rename(columns={'query_start_x': 'query_start_signal_sequence', 'query_start_y': 'query_start_gram_positive', 'query_start': 'query_start_microcin_domain'})

        #write csv merged_df
        merged_df.to_csv(output.final, index=False)


rule final_filter:
    input:
        config["outdir"] + "/results/temp/process_files/hmm_hits_with_domain_hits.csv"
    output:
        config["outdir"] + "/results/final/unprocessed_output.csv"
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

rule clean_outputs:
    input:
        config["outdir"] + "/results/final/unprocessed_output.csv"
    output:
        config["outdir"] + "/results/final/final_output.csv",
        config["outdir"] + "/results/temp/temp_file_check.txt"
    run:
        # load df from csv of file preprocessed_output.csv
        df = pd.read_csv(input[0])

        # rename columns
        df = df.rename(columns={'id': 'input_file', 'seqID': "contig_ORF"})

        # remove before } and then first 1 character in 'contig_ORF'
        df['contig_ORF'] = df['contig_ORF'].str.split('}').str[1].str[1:]

        # select text between { and } in column 'contig'
        df['input_file'] = df['input_file'].str.extract(r'{(.*)}')

        # strip * from seq
        df['seq'] = df['seq'].str.strip('*')

        # change column 'protein_length' to 'length'
        df = df.rename(columns={'protein_length': 'protein'})
        # make 'length' len(seq)
        df['length'] = df['seq'].str.len()

        # remove columns
        df = df.drop(columns=['pephash', 'dnahash', 'dna', 'start', 'stop', 'length', 'mismatch', 'sample', 'dna_length',
                              'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue_signal_sequence', 
                              'hit_start_signal_sequence', 'query_start_signal_sequence', 'evalue_gram_positive', 
                              'hit_start_gram_positive', 'query_start_gram_positive', 'evalue_microcin_domain', 
                              'hit_start_microcin_domain', 'query_start_microcin_domain'])

        # write csv
        df.to_csv(output[0], index=False)

        # Touch temp file
        shell("touch {output[1]}")

del_temp = config["del_temp"]

rule delete_temp_files:
    input:
        config["outdir"] + "/results/temp/temp_file_check.txt"
    output:
        temp(config["outdir"] + "/results/remove_temp_files.txt")
    run:
        # Touch temp file
        shell("touch {output[0]}")

        # check if config[del_temp] == True and delete if so
        if del_temp:
            # Generate text in python for bash command to delete temp files
            delete_command = f"rm -r {config['outdir']}/results/temp"
            shell(delete_command)
