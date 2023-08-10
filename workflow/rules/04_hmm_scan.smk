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

hmmsearch = start + command + Z + bias + T + end

# OLD: "threads=(({threads}-1)) && hmmsearch --cpu $threads --tblout {output.tblout} {input.hmmfile} {input.seqdb} > {output.fa}"


#Create a Mafft alignment of signal sequences
rule mafft_microcin:
    input: 
        "../../resources/input/microcins_verified_40.fa"
    output:
        config["outdir"] + "/results/mafft/microcin_mafft.aln"
    threads:
        workflow.cores * 0.9
    shell:
        "mafft --auto --globalpair --maxiterate 1000 --reorder --amino {input} > {output}"

#Build the pHMM using microcin sequences
rule buildhmm_microcin:
    input:
        config["outdir"] + "/results/mafft/microcin_mafft.aln"
    output:
        config["outdir"] + "/results/hmmsearch/microcin.hmm"
    threads:
        workflow.cores * 0.9
    shell:
        "hmmbuild --cpu {threads} --amino {output} {input}"# && hmmcalibrate {output}"

#Run hmmsearch to find signal sequences in the ORFs
rule hmmsearch:
    input:
        hmmfile = config["outdir"] + "/results/hmmsearch/microcin.hmm",
        seqdb = config["outdir"] + "/results/ORFs/all_AAs.faa"
    output:
        tblout = config["outdir"] + "/results/hmmsearch/microcin_hmmsearch_tblout.txt",
        domtblout = config["outdir"] + "/results/hmmsearch/microcin_hmmsearch_domtblout.txt",
        fa = config["outdir"] + "/results/hmmsearch/microcin_hmmsearch.fa"
    threads:
        workflow.cores * 0.9
    shell:
        hmmsearch
