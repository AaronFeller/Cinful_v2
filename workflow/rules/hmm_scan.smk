#Create a Mafft alignment of signal sequences
rule mafft_ss:
    input: 
        "../resources/input/ss_verified.fa"
    output:
        config["outdir"] + "/results/mafft/ss_mafft.aln"
    threads:
        3
    shell:
        "mafft --auto --globalpair --maxiterate 1000 --reorder --amino {input} > {output}"

#Build the pHMM using microcin double glycine signal sequence
rule buildhmm_microcin_ss:
    input:
        config["outdir"] + "/results/mafft/ss_mafft.aln"
    output:
        config["outdir"] + "/results/signalSeq/ss.hmm"
    threads:
        3
    shell:
        "hmmbuild --amino {output} {input}"


rule hmmsearch:
    input:
        hmmfile = config["outdir"] + "/results/signalSeq/ss.hmm",
        seqdb = config["outdir"] + "/results/ORFs/ORFs_filtered_all.faa"
    output:
        tblout = config["outdir"] + "/results/hmmsearch/ss_hmmsearch_tblout.txt",
        fa = config["outdir"] + "/results/hmmsearch/ss_hmmsearch.fa"
    threads:
        3
    shell:
        "hmmsearch --cpu {threads} -A {output.fa} --tblout {output.tblout} {input.hmmfile} {input.seqdb}"

# This rule will run phmmer for each input sequence individually
# rule phmmer_ss:
#     input:
#         seqfile = config["outdir"] + "/results/signalSeq/ss_mafft.aln",
#         seqdb = config["outdir"] + "/results/microcins/filtered_nr.fa"
#     output:
#         tblout = config["outdir"] + "/results/phmmer/ss_hmm_out.txt",
#     threads:375
#     shell:
#         "phmmer --cpu {threads} --tblout {output.tblout} {input.seqfile} {input.seqdb}"