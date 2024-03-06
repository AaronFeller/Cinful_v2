# Run orf finder script on all genomes
rule ORF_finder:
    input:
        genomes = "{sample}.fna"
    output:
        temp(config["outdir"] + "/results/temp/ORFs/{sample}_aminoAcid.faa")
    script:
        "../scripts/find_ORFs_initialize_on_start.py"

# Edit the contig names in the protein fasta files
rule edit_contig_names:
    input:
        genomes = config["outdir"] + "/results/temp/ORFs/{sample}_aminoAcid.faa"
    output:
        checkfiles = temp(config["outdir"] + "/results/temp/ORFs/{sample}_successfully_edited.txt")
    shell:
        """
        filename=$(basename {input.genomes}) &&
        sed -i 's/ /_/g' {input.genomes} &&
        sed -i "s/>/>{{$filename}}_/g" {input.genomes} &&
        touch {output.checkfiles}
        """
