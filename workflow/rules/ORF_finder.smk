#sets wildcards
SAMPLES, = glob_wildcards("{sample}.fna")

rule ORF_finder:
	input:
	    genomes = "{sample}.fna"
	output:
		config["outdir"] + "/results/ORFs/{sample}_ORFs.faa"
	script:	
		"../scripts/find_ORFs.py"

rule ORF_filter:
	input:
	    config["outdir"] + "/results/ORFs/{sample}_ORFs.faa"
	output:
		config["outdir"] + "/results/ORFs/{sample}_ORFs_filtered.faa"
	script:
		"../scripts/filter_ORFs.py"
