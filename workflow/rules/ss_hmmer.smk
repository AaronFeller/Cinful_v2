threads_max = 4

rule msa_microcin_ss:
	input:
		"results/database/SP.verified.aln"
	output:
		"results/signalSeq/SS.verified.aln"
	threads:threads_max
	shell:
		"mafft --thread {threads} --auto {input} > {output}"


rule buildhmm_microcin_ss:
	input:
		"results/signalSeq/SS.verified.aln"
	output:
		"results/signalSeq/SS.verified.hmm"
	shell:
		"hmmbuild {output} {input}"