#sets wildcards
SAMPLES, = glob_wildcards("{sample}.fna")

rule ORF_finder:
	input:
		genomes = "{sample}.fna"
	output:
#		temp(config["outdir"] + "/results/ORFs/{sample}_ORFs.faa"),
		config["outdir"] + "/results/ORFs/{sample}_aminoAcid.faa",
		config["outdir"] + "/results/ORFs/{sample}_ORFs.csv"
	script:	
		"../scripts/find_ORFs_initialize_on_start.py"

# rule ORF_filter:
# 	input:
# 	    config["outdir"] + "/results/ORFs/{sample}_ORFs.faa"
# 	output:
# 		temp(config["outdir"] + "/results/ORFs/{sample}_ORFs_filtered.faa")
# 	script:
# 		"../scripts/filter_ORFs.py"