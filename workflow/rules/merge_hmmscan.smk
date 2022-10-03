rule merge_hmmscan:
    input:
        tblout = config["outdir"] + "/results/hmmsearch/ss_hmmsearch_tblout.txt",
        prodigal_all = config["outdir"] + "/results/ORFs/ORFs_filtered_all.csv"
    output:
        hmm_hits = config["outdir"] + "/results/hmmsearch/hmm_hits.csv"
    run:
        from collections import defaultdict
        from Bio import SearchIO
        import pandas as pd

        attribs = ['accession', 'bias', 'bitscore', 'description', 'cluster_num', 'domain_exp_num',  'domain_included_num', 'domain_obs_num', 'domain_reported_num', 'env_num', 'evalue', 'id', 'overlap_num', 'region_num']
        hits = defaultdict(list)

        with open(input.tblout) as handle:
            for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
                for hit in queryresult.hits:
                    for attrib in attribs:
                        hits[attrib].append(getattr(hit, attrib))
        pd.DataFrame.from_dict(hits).rename(columns={'id':'pephash'}).merge(pd.read_csv(input.prodigal_all), on='pephash', how='left').to_csv(output.hmm_hits, index=False)
