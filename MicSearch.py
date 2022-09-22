import os
import sys
import subprocess as sp
import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import errno

#Confirms the directory for input and output are true directories
class readable_dir(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        prospective_dir=values
        if not os.path.isdir(prospective_dir):
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace,self.dest,prospective_dir)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

#Confirms the input files are true files
class valid_file(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        prospective_file=values
        if not os.path.exists(prospective_file):
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_file))


def subDF2Fasta(df, subName, baseDir,seqCol="seq", idCol = "cinful_id"):
	recordDict = {}
	for row in df.to_dict(orient="records"):
		cinful_id = row["cinful_id"]
		sample = cinful_id.split("|")[0]
		header = cinful_id.split(sample)[1].strip("|")
		seq = row["seq"]

		record = SeqRecord(
    		Seq(seq),
    		id=header,
    		description=""
		)

		if sample not in recordDict:
			recordDict[sample] = []
		recordDict[sample].append(record)

	for sample in recordDict:
		outFile = Path(baseDir)/ "03_best_hits"/"fastas"/ subName / f"{sample}.{subName}.fa"
		if not os.path.exists(os.path.dirname(outFile)):
			try:
				os.makedirs(os.path.dirname(outFile))
			except OSError as exc: # Guard against race condition
				if exc.errno != errno.EEXIST:
					raise
		with open(outFile, "w") as fastaOut:
			SeqIO.write(recordDict[sample], fastaOut, "fasta")

def main():
    # set up command line arguments
    parser = argparse.ArgumentParser(description='cinful')
    parser.add_argument('-d', '--directory', action=readable_dir,required=True,
        help='Must be a directory containing uncompressed FASTA formatted genome assemblies with .fna extension. Files within nested directories are fine')
    parser.add_argument('-o', '--outDir', type=str,default="default_output",
        help='This directory will contain all output files. It will be nested under the input directory.')
    parser.add_argument('-t', '--threads', type=int, default=1,
        help = 'This specifies how many threads to allow snakemake to have access to for parallelization')
    # parser.add_argument('--snakemake_params', type=list, nargs='+')
    args = parser.parse_args()

    threads = args.threads
    workdir = args.directory
    outdir = args.outDir

    # absolute path of this file, needed to find snakefile
    currentAbsPath = os.path.dirname(os.path.abspath(__name__))
    snakefile = f"{currentAbsPath}/workflow/Snakefile"

    # snakemake command to run
    cmd = [ "python3",
            "-m",
            "snakemake",
            "--snakefile",
            snakefile,
            "-j",
            str(threads),
            "--directory",
            workdir,
            "--config",
            f"outdir={outdir}", 
            ]
    print("Running the following command:")
    print(" ".join(cmd))

    try:
        sp.check_output(cmd)
        print("MicSearch finished succesfuly!\n")
        # all_hits = Path(workdir)/ outdir / "02_homology_results/all_merged.csv"
        # print(f"Now checking that any hits were identified in {all_hits}")
        # all_hitsDF = pd.read_csv(all_hits).sort_values(['component','bitscore'],ascending=False)
        # microcinDF = all_hitsDF[all_hitsDF["component"] == "microcins.verified" ]
        # CvaBDF = all_hitsDF[all_hitsDF["component"] == "CvaB.verified" ]
        # MFPDF = all_hitsDF[all_hitsDF["component"] == "MFP.verified" ]
        # print("\nFinal fasta results will be written to:",Path(workdir)/outdir/ "03_best_hits"/"fastas")
        # subDF2Fasta(microcinDF,"microcin",Path(workdir)/outdir)
        # subDF2Fasta(CvaBDF,"CvaB",Path(workdir)/outdir)
        # subDF2Fasta(MFPDF,"MFP",Path(workdir)/outdir)
    except sp.CalledProcessError as e:
        print(e.output)

if __name__ == "__main__":
    main()