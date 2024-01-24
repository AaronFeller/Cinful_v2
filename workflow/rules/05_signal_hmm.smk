import Bio.SearchIO.HmmerIO as HmmerIO
from Bio import SearchIO

"""
This section was removed for timesaving, as the alignments and hmm are now precomputed and included in resources

However, this code can be used if the user would like to create new alignments and hmm models

### THIS PORTION IS FOR THE 41 VALIDATED MICROCIN SIGNAL SEQUENCES ###

# Create a Mafft alignment of signal sequences
rule mafft_ss:
    input: 
        "../../resources/input/ss_verified.fa"
    output:
        config["outdir"] + "/results/temp/mafft/ss_mafft.aln"
    threads:
        workflow.cores * 0.9
    shell:
        "mafft --auto --globalpair --maxiterate 1000 --reorder --amino {input} > {output}"

# Build the pHMM using microcin sequences
rule buildhmm_signal_sequence:
    input:
        config["outdir"] + "/results/temp/mafft/ss_mafft.aln"
    output:
        config["outdir"] + "/results/temp/hmmsearch/ss.hmm"
    threads:
        workflow.cores * 0.9
    shell:
        "hmmbuild --cpu {threads} --amino {output} {input}"# && hmmcalibrate {output}"
"""

#Run hmmsearch to find signal sequences in the ORFs
rule signal_sequence_hmmsearch:
    input:
        hmmfile = "../../resources/input/hmm_models/ss.hmm",
        seqdb = config["outdir"] + "/results/temp/hmmsearch/hmm_hits.protein.fasta"
    output:
        tblout = config["outdir"] + "/results/temp/hmmsearch/ss_hmmsearch_tblout.txt",
        domtblout = config["outdir"] + "/results/temp/hmmsearch/ss_hmmsearch_domtblout.txt",
        fa = config["outdir"] + "/results/temp/hmmsearch/ss_hmmsearch.fa"
    threads:
        workflow.cores * 0.9
    shell:
        "threads=$(({threads}-1)) && hmmsearch -E 10 --cpu $threads --tblout {output.tblout} --domtblout {output.domtblout} {input.hmmfile} {input.seqdb} > {output.fa}"

rule extract_ss_hits:
    input:
        domtblout = config["outdir"] + "/results/temp/hmmsearch/ss_hmmsearch_domtblout.txt"
    output:
        extracted_hits = config["outdir"] + "/results/temp/hmmsearch/ss_domains.csv"
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

"""
This section was removed for timesaving, as the alignments and hmm are now precomputed and included in resources

However, this code can be used if the user would like to create new alignments and hmm models

### THIS PORTION IS FOR GRAM POSITIVE SIGNALS ###

# Create a Mafft alignment of signal sequences
rule mafft_ss_gram_positive:
    input: 
        "../../resources/input/ss_gram_positive.fa"
    output:
        config["outdir"] + "/results/temp/mafft/ss_gram_positive_mafft.aln"
    threads:
        workflow.cores * 0.9
    shell:
        "mafft --auto --globalpair --maxiterate 1000 --reorder --amino {input} > {output}"

# Build the pHMM using microcin sequences
rule buildhmm_ss_gram_positive:
    input:
        config["outdir"] + "/results/temp/mafft/ss_gram_positive_mafft.aln"
    output:
        config["outdir"] + "/results/temp/hmmsearch/ss_gram_positive.hmm"
    threads:
        workflow.cores * 0.9
    shell:
        "hmmbuild --cpu {threads} --amino {output} {input}"# && hmmcalibrate {output}"

"""


#Run hmmsearch to find signal sequences in the ORFs
rule ss_gram_positive_hmmsearch:
    input:
        hmmfile = "../../resources/input/hmm_models/ss_gram_positive.hmm",
        seqdb = config["outdir"] + "/results/temp/hmmsearch/hmm_hits.protein.fasta"
    output:
        tblout = config["outdir"] + "/results/temp/hmmsearch/ss_gram_positive_hmmsearch_tblout.txt",
        domtblout = config["outdir"] + "/results/temp/hmmsearch/ss_gram_positive_hmmsearch_domtblout.txt",
        fa = config["outdir"] + "/results/temp/hmmsearch/ss_gram_positive_hmmsearch.fa"
    threads:
        workflow.cores * 0.9
    shell:
        "threads=$(({threads}-1)) && hmmsearch -E 10 --cpu $threads --tblout {output.tblout} --domtblout {output.domtblout} {input.hmmfile} {input.seqdb} > {output.fa}"

rule extract_ss_gram_positive_hits:
    input:
        domtblout = config["outdir"] + "/results/temp/hmmsearch/ss_gram_positive_hmmsearch_domtblout.txt"
    output:
        extracted_hits = config["outdir"] + "/results/temp/hmmsearch/ss_gram_positive_domains.csv"
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
