#Create a Mafft alignment of signal sequences
rule mafft_ss:
    input: 
        "../../resources/input/microcins_verified_40.fa"
    output:
        config["outdir"] + "/results/mafft/ss_mafft.aln"
    threads:
        workflow.cores * 0.9
    shell:
        "mafft --auto --globalpair --maxiterate 1000 --reorder --amino {input} > {output}"

#Build the pHMM using microcin double glycine signal sequence
rule buildhmm_microcin_ss:
    input:
        config["outdir"] + "/results/mafft/ss_mafft.aln"
    output:
        config["outdir"] + "/results/signalSeq/ss.hmm"
    threads:
        workflow.cores * 0.9
    shell:
        "hmmbuild --cpu {threads} --amino {output} {input}"# && hmmcalibrate {output}"

#Run hmmsearch to find signal sequences in the ORFs
rule hmmsearch:
    input:
        hmmfile = config["outdir"] + "/results/signalSeq/ss.hmm",
        seqdb = config["outdir"] + "/results/ORFs/all_AAs.faa"
    output:
        tblout = config["outdir"] + "/results/hmmsearch/ss_hmmsearch_tblout.txt",
        fa = config["outdir"] + "/results/hmmsearch/ss_hmmsearch.fa"
    threads:
        workflow.cores * 0.9
    shell:
        "threads=$(({threads}-1)) && hmmsearch --cpu $threads --tblout {output.tblout} {input.hmmfile} {input.seqdb} > {output.fa}"
