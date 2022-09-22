# import dependencies

from io import StringIO
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
import seqhash

# define a function to produce a hash dataframe from the fasta input file, this will create a peptide hash of each input line 

def fa2hashDF(fasta_file):
    outDict = {"header":[],"pephash":[],"sequence":[]}
    with open(fasta_file) as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            pephash = seqhash.seqhash(record.seq,dna_type='PROTEIN')
            outDict["header"].append(record.id)
            outDict["pephash"].append(pephash)
            outDict["sequence"].append(record.seq)
        return pd.DataFrame.from_dict(outDict)

# input initial fasta with verified proteins and generate pephash

rule write_CvaB:
	input:
		verified_CvaB = "resources/input/CvaB.verified.pep"
	output:
		CvaB_fasta = config["outdir"] + "/results/database/CvaB.verified.pep",
		CvaB_pepHash = config["outdir"] + "/results/database/CvaB.verified.pephash.csv"
	run:
		verified_CvaB = input.verified_CvaB
		verified_CvaB_SeqIO = SeqIO.parse(StringIO(open(verified_CvaB).read()), "fasta")
		verified_CvaB_pepHashDF = fa2hashDF(verified_CvaB)
		verified_CvaB_pepHashDF.to_csv(output.CvaB_pepHash)
		with open(output.CvaB_fasta,"w") as seq_out:
			SeqIO.write(verified_CvaB_SeqIO, seq_out,"fasta")

rule write_MFP:
	input:
		verified_MFP = "resources/input/MFP.verified.pep"
	output:
		MFP_fasta = config["outdir"] + "/results/database/MFP.verified.pep",
		MFP_pepHash = config["outdir"] + "/results/database/MFP.verified.pephash.csv"
	run:
		verified_MFP = input.verified_MFP
		verified_MFP_SeqIO = SeqIO.parse(StringIO(open(verified_MFP).read()), "fasta")
		verified_MFP_pepHashDF = fa2hashDF(verified_MFP)
		verified_MFP_pepHashDF.to_csv(output.MFP_pepHash)
		with open(output.MFP_fasta,"w") as seq_out:
			SeqIO.write(verified_MFP_SeqIO, seq_out,"fasta")

rule write_microcins:
	input:
		verified_microcins = "resources/input/microcins.verified.pep"
	output:
		microcin_fasta = config["outdir"] + "/results/database/Microcins.verified.pep",
		microcin_pepHash = config["outdir"] + "/results/database/Microcins.verified.pephash.csv"
	run:
		verified_microcins = input.verified_microcins
		verified_microcins_SeqIO = SeqIO.parse(StringIO(open(verified_microcins).read()), "fasta")
		verified_microcins_pepHashDF = fa2hashDF(verified_microcins)
		verified_microcins_pepHashDF.to_csv(output.microcin_pepHash)
		with open(output.microcin_fasta,"w") as seq_out:
			SeqIO.write(verified_microcins_SeqIO, seq_out,"fasta")

rule write_immunity_proteins:
	input:
		verified_immunity_proteins = "resources/input/immunity_proteins.verified.pep"
	output:
		immunity_protein_fasta = config["outdir"] + "/results/database/immunity_proteins.verified.pep",
		immunity_protein_pepHash = config["outdir"] + "/results/database/immunity_proteins.verified.pephash.csv"
	run:
		verified_immunity_proteins = input.verified_immunity_proteins
		verified_immunity_proteins_SeqIO = SeqIO.parse(StringIO(open(verified_immunity_proteins).read()), "fasta")
		verified_immunity_proteins_pepHashDF = fa2hashDF(verified_immunity_proteins)
		verified_immunity_proteins_pepHashDF.to_csv(output.immunity_protein_pepHash)
		with open(output.immunity_protein_fasta,"w") as seq_out:
  			SeqIO.write(verified_immunity_proteins_SeqIO, seq_out,"fasta")

rule microcin_signal:
	input:
		verified_SP = "resources/input/SP.verified.pep"
	output:
		SP_alignment = config["outdir"] + "/results/database/SP.verified.aln"
	run:
		verified_SP = open(input.verified_SP)
		verified_SP_msa = AlignIO.read(StringIO(verified_SP.read()), "fasta")
		with open(output.SP_alignment,"w") as alignment_out:
  			AlignIO.write(verified_SP_msa,alignment_out,"fasta")
