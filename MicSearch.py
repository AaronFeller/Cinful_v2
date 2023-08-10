import os
import subprocess as sp
import argparse
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
    parser.add_argument('-d', '--directory', action=readable_dir, required=True,
        help= 'Must be a directory containing uncompressed FASTA formatted genome assemblies with .fna extension. Files within nested directories from refseq/genbank are acceptable.')
    parser.add_argument('-o', '--outDir', type=str, required=True,
        help= 'This directory will contain all output files. It will be nested under the input directory.')
    parser.add_argument('-t', '--threads', type=int, default=4,
        help = 'This specifies how many threads to allow snakemake to have access to for parallelization. Default is 1.')
    parser.add_argument('-b', '--biased_composition_filter', type=bool, default=True,
        help = 'This specifies boolian for whether to use the biased composition filter. Default is True.')
    parser.add_argument('-e', '--evalue', type=int, default="100",
        help = 'This specifies the evalue threshold for initial HMMER search. Best_hits filters to an evalue of 1 Default is 100.')
    parser.add_argument('-n', '--num_queries', type=int, default="0",
        help = 'This specifies the number of queries to use for HMMER calculation of evalue. Default is 0 (or true database size).')
    
    # parser.add_argument('--snakemake_params', type=list, nargs='+')
    args = parser.parse_args()

    threads = args.threads
    workdir = args.directory
    outdir = args.outDir
    biased_composition_filter = args.biased_composition_filter
    evalue = args.evalue
    num_queries = args.num_queries

    # absolute path of this file, needed to find snakefile
    currentAbsPath = os.path.dirname(os.path.abspath(__name__))
    snakefile = f"{currentAbsPath}/workflow/Snakefile"

    # snakemake command to run
    cmd = ["python3",
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
        f"biased_composition_filter={biased_composition_filter}",
        f"evalue={evalue}",
        f"num_queries={num_queries}",
        "--rerun-triggers",
        "mtime",
        ]

    print("Running the following command:")
    print(" ".join(cmd))

    try:
        sp.check_output(cmd)
        print("MicSearch finished successfully!\n")
    except sp.CalledProcessError as error:
        print('Encountered error while running snakemake.')
        #print(error.output)

if __name__ == "__main__":
    main()
