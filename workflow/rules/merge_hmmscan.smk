rule merge_hmmscan:
    input:
        tblout = config["outdir"] + "/results/hmmsearch/ss_hmmsearch_tblout.txt",
        prodigal_all = config["outdir"] + "/results/prodigal/prodigal_out_all_nr_expanded.csv"
    output:
        out = config["outdir"] + "/results/hmmsearch/hmm_hits.csv"
    script:
        "../scripts/merge_hmmsearch_w_peptides.R"