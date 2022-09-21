
#sets wildcards
SAMPLES, = glob_wildcards("{sample}.fna")

# runs prodigal on normal mode with inputs of all .fna files and outputs .gff3, .cds, and .faa

rule prodigal:
	input:
	    fna = "resources/genomes/{sample}.fna"
	output:
		gff3 = "results/prodigal/{sample}.gff3",
		cds = "results/prodigal/{sample}.cds",
		aa = "results/prodigal/{sample}.faa"
	shell:
		"prodigal -i {input.fna} -o {output.gff3} -a {output.aa} -d {output.cds}"
