from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Generate blast db from input
rule makeblastdb:
    input:
        "../../resources/input/pblast_database/pblast_database.fa"
    output:
        "../../resources/input/pblast_database/pblast_database.fa.phr"
    shell:
        "makeblastdb -in {input} -dbtype prot -parse_seqids"

# Generate FASTA file from hmmsearch results
rule hmmsearch2fasta:
    input:
        config["outdir"] + "/results/temp/hmmsearch/hmm_hits.csv"
    output:
        config["outdir"] + "/results/temp/hmmsearch/hmm_hits.protein.fasta",
        config["outdir"] + "/results/temp/hmmsearch/hmm_hits.dna.fasta"
    run:

        def pandas_df_to_fasta(dataframe, header_column, sequence_column, output_file):
            records = []

            for _, row in dataframe.iterrows():
                seq_id = row[header_column]
                sequence = row[sequence_column]
                sequence = Seq(sequence)
                record = SeqRecord(sequence, id=seq_id, description="")
                records.append(record)

            with open(output_file, "w") as output_handle:
                SeqIO.write(records, output_handle, "fasta")

        df = pd.read_csv(input[0])
        protein_df = df[["pephash", "seq"]]
        #Trim * from seq
        protein_df["seq"] = protein_df["seq"].str.replace("*", "", regex=False)

        #select unique rows based on pephash
        protein_df = protein_df.drop_duplicates(subset=["pephash"])

        dna_df = df[["dnahash", "dna"]]
        #select unique rows based on dnahash
        dna_df = dna_df.drop_duplicates(subset=["dnahash"])

        pandas_df_to_fasta(protein_df, "pephash", "seq", output[0])
        pandas_df_to_fasta(dna_df, "dnahash", "dna", output[1])

#Run blastp
rule blastp:
    input:
        query = config["outdir"] + "/results/temp/hmmsearch/hmm_hits.protein.fasta",
        database = "../../resources/input/pblast_database/pblast_database.fa",
        blastdb = "../../resources/input/pblast_database/pblast_database.fa.phr"
    output:
        config["outdir"] + "/results/temp/blast/blastp_results.txt"
    threads:
        workflow.cores * 0.9
    shell:
        "blastp -num_threads {threads} -db {input.database} -query {input.query} -out {output} -outfmt 6 -evalue 1 -max_target_seqs 1"

# Merge pid and hit with hmm_hits.csv
rule merge_HMM_and_blastp:
    input:
        blast = config["outdir"] + "/results/temp/blast/blastp_results.txt",
        hmm = config["outdir"] + "/results/temp/hmmsearch/hmm_hits.csv"
    output:
        config["outdir"] + "/results/temp/process_files/hmm_hits_with_pident.csv"
    run:
        import pandas as pd
        blast = pd.read_csv(input.blast, sep='\t', header=None)
        blast.columns = ["pephash", "subject", "pid", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue_blastp", "bitscore_blastp"]
        hmm = pd.read_csv(input.hmm)
        hmm = hmm.merge(blast, on="pephash", how='left')
        hmm.to_csv(output[0], index=False)
