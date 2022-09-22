threads_max = 4 

#sets wildcards
SAMPLES, = glob_wildcards("{sample}.fna")

# runs prodigal on normal mode with inputs of all .fna files and outputs .gff3, .cds, and .faa

rule prodigal:
	input:
	    fna = "../../resources/genomes/{sample}.fna"
	output:
		gff3 = config["outdir"] + "/results/prodigal/{sample}.gff3",
		cds = config["outdir"] + "/results/prodigal/{sample}.cds",
		aa = config["outdir"] + "/results/prodigal/{sample}.faa"
	shell:
		"prodigal -i {input.fna} -o {output.gff3} -a {output.aa} -d {output.cds}"
