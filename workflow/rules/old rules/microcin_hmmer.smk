import os
import sys
import pandas as pd
import pyhmmer
from functools import reduce

def hasAllStandardAA(seq, alphabet="ACDEFGHIKLMNPQRSTVWY",ignore="*"):
	return (set(seq) - set(alphabet+ignore)) == set()

def load_hmm(hmmFile):
    with pyhmmer.plan7.HMMFile(hmmFile) as h:
        hmm = next(h)
    return hmm

def run_hmmsearch(queryFile, hmmFile):
    hmm = load_hmm(hmmFile)
    with pyhmmer.easel.SequenceFile(queryFile) as seq_file:
        sequences = [ seq.digitize(hmm.alphabet) for seq in seq_file ]
    pipeline = pyhmmer.plan7.Pipeline(hmm.alphabet)
    hits = pipeline.search_hmm(hmm, sequences) # Has lots of goodies!

    return hits , hmm.name.decode() # [hit.name.decode() for hit in hmmerOut]

def hmmsearch(queryFile, hmm):

	with pyhmmer.easel.SequenceFile(queryFile) as seq_file:
		sequences = [ seq.digitize(hmm.alphabet) for seq in seq_file ]
	pipeline = pyhmmer.plan7.Pipeline(hmm.alphabet)
	hits = pipeline.search_hmm(hmm, sequences)
	return hits

def build_hmm(alnFile):
	abc = pyhmmer.easel.Alphabet.amino()
	builder = pyhmmer.plan7.Builder(alphabet=abc)

	with pyhmmer.easel.MSAFile(alnFile) as msa_file:
		msa_file.set_digital(abc)
		msa = next(msa_file)
  # MSA must have a name, otherwise building will fail
	if msa.name is None:
		msa.name = b"alignment"
	builder = pyhmmer.plan7.Builder(abc)
	background = pyhmmer.plan7.Background(abc)
	hmm, _, _ = builder.build_msa(msa, background)

	return hmm


rule signalSeqHMM:
    input:
        input_seqs = "results/microcins/filtered_nr.fa",
        signalSeqAln = "results/database/SP.verified.aln"
    output:
        ss_hit = "results/microcins/signalSeq.hit.csv"
    run:
        signalSeqHMM = build_hmm(input.signalSeqAln)

        signalSeqHits = hmmsearch(input.input_seqs, signalSeqHMM)
        signalSeqHitStr = [hit.name.decode('utf-8') for hit in signalSeqHits]
        matchDF = pd.DataFrame.from_dict({"signalMatch":signalSeqHitStr})
        matchDF.to_csv(output.ss_hit)
