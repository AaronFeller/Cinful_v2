threads_max = 4

rule phmmer_ss:
    input:
        seqfile = "results/signalSeq/SS.verified.aln",
        seqdb = "results/prodigal/prodigal_out.all.nr.faa"
    output:
        hmm_txt = "results/phmmer/ss_hmm_hits.txt",
        hmm_A = "results/phmmer/ss_hmm_hit2.txt",
        tblout = "results/phmmer/ss_hmm_tblout.txt",
        domtblout = "results/phmmer/ss_hmm_domtblout.txt",
        pfamtblout = "results/phmmer/ss_hmm_pfamtblout.txt"

    threads:threads_max
    shell:
        "phmmer --cpu {threads} -o {output.hmm_txt} -A {output.hmm_A} --tblout {output.tblout} --domtblout {output.domtblout} --pfamtblout {output.pfamtblout} {input.seqfile} {input.seqdb}"


rule hmmsearch:
    input:
        hmmfile = "results/signalSeq/SS.verified.hmm",
        seqdb = "results/prodigal/prodigal_out.all.nr.faa"
    output:
        "results/hmmsearch/hmmsearch_tblout.txt"
    threads:threads_max
    shell:
        "hmmsearch --cpu {threads} --tblout {output} {input.hmmfile} {input.seqdb} "