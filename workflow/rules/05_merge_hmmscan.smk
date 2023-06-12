from collections import defaultdict
from Bio import SearchIO
import pandas as pd

# rule merge_hmmscan:
#     input:
#         tblout = config["outdir"] + "/results/hmmsearch/ss_hmmsearch_tblout.txt",
#         ORFs = config["outdir"] + "/results/ORFs/all_ORFs.csv"
#     output:
#         hmm_hits = protected(config["outdir"] + "/results/hmmsearch/hmm_hits.csv")
#     run:

#         attribs = ['id', 'bitscore', 'evalue']
#         hits = defaultdict(list)

#         with open(input.tblout) as handle:
#             for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
#                 for hit in queryresult.hits:
#                     for attrib in attribs:
#                         hits[attrib].append(getattr(hit, attrib))
#         pd.DataFrame.from_dict(hits).merge(pd.read_csv(input.ORFs), on='id', how='left').to_csv(output.hmm_hits, index=False)

rule extract_hits:
    input:
        tblout = config["outdir"] + "/results/hmmsearch/ss_hmmsearch_tblout.txt"
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

rule build_contig_to_genome_table:
    input: 
        expand(config["outdir"] + "/results/ORFs/{sample}_ORFs.csv", sample=SAMPLES)
    output:
        contig_to_genome_table = config["outdir"] + "/results/ORFs/contig_to_genome_table.csv"
    run:
        # Create empty DataFrame
        contig_to_genome_table = pd.DataFrame()
        subset_contig_to_genome = pd.DataFrame(columns=['genome', 'contig'])

        # Iterate through each sample
        count = -1
        for i in range(len(SAMPLES)):    
            # Assign sample to the current sample
            sample = SAMPLES[i]
            
            # Read ORF file and assign it to orfs_df
            orfs_df = pd.read_csv(str(config["outdir"] + "/results/ORFs/" + sample + "_ORFs.csv"))

            # Create a new column called 'sample' and assign it to the sample name
            for j in range(len(orfs_df['sample'].unique())):
                count += 1
                
                # Make DataFrame with two columns
                line = pd.DataFrame([[sample, orfs_df['sample'].unique()[j]]], columns=['genome', 'contig'])
                
                # Concatenate line with contig_to_genome_table
                contig_to_genome_table = pd.concat([contig_to_genome_table, line], ignore_index=True)

        # Write the final completed DataFrame to the output file
#        print(contig_to_genome_table)
        contig_to_genome_table.to_csv(output.contig_to_genome_table, index=False)

rule merge_hmmscan:
    input:
        extracted_hits = config["outdir"] + "/results/hmmsearch/extracted_hits.csv",
        contig_to_genome_table = config["outdir"] + "/results/ORFs/contig_to_genome_table.csv"
    output:
        hmm_hits = protected(config["outdir"] + "/results/hmmsearch/hmm_hits.csv")
    run:
        # Load hits
        extracted_df = pd.read_csv(input.extracted_hits).rename(columns={'id': 'contig'})

        # Add blank columns
        added_columns = ['id', 'pephash', 'sample', 'start', 'stop', 'strand', 'dna', 'seq']
        for column in added_columns:
            extracted_df[column] = ''

        # Iterate through each row of extracted_df
        for i in range(len(extracted_df)):
            # Find the genome that the contig belongs to
            contig_to_genome_table = pd.read_csv(input.contig_to_genome_table)
            genome_row = contig_to_genome_table[contig_to_genome_table['contig'] == extracted_df.loc[i, 'contig'].split('_ORF')[0]]
            if not genome_row.empty:
                genome = genome_row.iloc[0]['genome']

                # Generate the string for ORF file for the genome found
                orf_file = config["outdir"] + "/results/ORFs/" + genome + "_ORFs.csv"

                # Read ORF file and assign it to orfs_df
                orfs_df = pd.read_csv(orf_file)

                # Find the row in orfs_df that matches the 'contig' (should be only 1)
                matching_row = orfs_df[orfs_df['contig'] == extracted_df.loc[i, 'contig']]

                # Add values from columns to extracted_df
                for column in added_columns:
                    extracted_df.loc[i, column] = matching_row.iloc[0][column]

        # Write the final completed DataFrame to the output file
        extracted_df.to_csv(output.hmm_hits, index=False)


# rule merge_hmmscan:
#     input:
#         extracted_hits = config["outdir"] + "/results/hmmsearch/extracted_hits.csv",
#         orfs_csv = config["outdir"] + "/results/ORFs/all_ORFs.csv"
#     output:
#         hmm_hits = protected(config["outdir"] + "/results/hmmsearch/hmm_hits.csv")
#     run:
#         extracted_df = pd.read_csv(input.extracted_hits).rename(columns={'id': 'contig'})
#         orfs_df = pd.read_csv(input.orfs_csv)
#         merged_df = extracted_df.merge(orfs_df, on='contig', how='left')
#         merged_df.to_csv(output.hmm_hits, index=False)
