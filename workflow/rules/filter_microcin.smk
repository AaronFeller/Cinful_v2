rule filter_microcin:
	input:
		"results/prodigal/prodigal_out.all.nr.faa"
	output:
		"results/microcins/filtered_nr.fa"
	shell:
		"seqkit seq -m 30 -M 150 {input} | seqkit rmdup -s > {output}"