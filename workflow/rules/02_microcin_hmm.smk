from Bio import SeqIO

# Inputs from snakemake config
num_queries = config["num_queries"]
biased_composition_filter = config["biased_composition_filter"]
evalue = config["evalue"]

# Generate string for running hmmsearch
start = 'threads=$(({threads}-1)) && '
command = f'hmmsearch -E {evalue} '
if num_queries == 0:
    Z = ''
else:
    Z = f'-Z {num_queries} '
if biased_composition_filter:
    bias = ''
else:
    bias = '--nobias '
T = f'--cpu $threads '
end = '--tblout {output.tblout} --domtblout {output.domtblout} {input.hmmfile} {input.seqdb} > {output.fa}'

# Produce HMMER command
hmmsearch = start + command + Z + bias + T + end

# Create a Mafft alignment of microcins
rule mafft_microcin:
    input:
        "../../resources/input/microcins_verified.fa"
    output:
        config["outdir"] + "/results/temp/mafft/microcin_mafft.aln"
    threads:
        workflow.cores * 0.9
    shell:
        "mafft --auto --globalpair --maxiterate 1000 --reorder --amino {input} > {output}"

# Build the pHMM using microcin sequences
rule buildhmm_microcin:
    input:
        config["outdir"] + "/results/temp/mafft/microcin_mafft.aln"
    output:
        # config["outdir"] + "/results/temp/hmmsearch/microcin.hmm"
        "../../resources/input/hmm_models/microcin.hmm"
    threads:
        workflow.cores * 0.9
    shell:
        "hmmbuild --cpu {threads} --amino {output} {input}"# && hmmcalibrate {output}"

# Remove sequences longer than 10,000 amino acids, avoiding memory issues
rule remove_long_seqs:
    input:
        fasta=config["outdir"] + "/results/temp/ORFs/{sample}_aminoAcid.faa",
        checkfiles=config["outdir"] + "/results/temp/ORFs/{sample}_successfully_edited.txt"
    output:
        temp(config["outdir"] + "/results/temp/ORFs/{sample}_aminoAcid_clean.faa")
    run:
        with open(input.fasta, "r") as f:
            records = list(SeqIO.parse(f, "fasta"))
        with open(output[0], "w") as f:
            for record in records:
                if len(record.seq) <= 10000:
                    SeqIO.write(record, f, "fasta")
        # delete checkfiles
        os.remove(input.checkfiles)

# Run hmmsearch to find signal sequences in the ORFs
rule hmmsearch:
    input:
        hmmfile = "../../resources/input/hmm_models/microcin.hmm",
        seqdb = config["outdir"] + "/results/temp/ORFs/{sample}_aminoAcid_clean.faa"
    output:
        tblout = config["outdir"] + "/results/temp/hmmsearch/{sample}_microcin_hmmsearch_tblout.txt",
        domtblout = config["outdir"] + "/results/temp/hmmsearch/{sample}_microcin_hmmsearch_domtblout.txt",
        fa = config["outdir"] + "/results/temp/hmmsearch/{sample}_microcin_hmmsearch.fa"
    threads:
        3
    shell:
        hmmsearch
