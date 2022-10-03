rule filter_microcin:
	input:
		config["outdir"] + "/results/prodigal/prodigal_out_all_nr.faa"
	output:
		config["outdir"] + "/results/microcins/filtered_nr.fa"
	shell:
		"seqkit seq -m 30 -M 150 {input} | seqkit rmdup -s > {output}"