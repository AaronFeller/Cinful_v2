rule mafft_ss:
    input: 
        "resources/input/ss_verified.fa"
    output:
        "results/mafft/ss_mafft.aln"
    shell:
        "mafft --auto --globalpair --maxiterate 10 --reorder --nuc {input} > {output}"

# rule phmmer_ss:
#     input:
#         seqfile = "results/signalSeq/ss_mafft.aln",
#         seqdb = "results/microcins/filtered_nr.fa"
#     output:
#         tblout = "results/phmmer/ss_hmm_out.txt",
#     threads:workflow.cores * 0.75
#     shell:
#         "phmmer --cpu {threads} --tblout {output.tblout} {input.seqfile} {input.seqdb}"

rule hmmsearch:
    input:
        hmmfile = "results/signalSeq/SS.verified.hmm",
        seqdb = "results/microcins/filtered_nr.fa"
    output:
        "results/hmmsearch/ss_hmmsearch_tblout.txt"
    threads:workflow.cores * 0.75
    shell:
        "hmmsearch --cpu {threads} --tblout {output} {input.hmmfile} {input.seqdb}"