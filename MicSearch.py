import os
import subprocess as sp
import argparse

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

def main():
    # set up command line arguments
    parser = argparse.ArgumentParser(description='cinful')
    parser.add_argument('-d', '--directory', action=readable_dir, required=True,
        help= 'Must be a directory containing uncompressed FASTA formatted genome assemblies with .fna extension. Files within nested directories from refseq/genbank are acceptable.')
    parser.add_argument('-o', '--outDir', type=str, required=True,
        help= 'This directory will contain all output files. It will be nested under the input directory.')
    parser.add_argument('-t', '--threads', type=int, default=4,
        help = 'This specifies how many threads to allow snakemake to have access to for parallelization. Default is 1.')
    parser.add_argument('-f', '--biased_composition_filter', type=bool, default=True,
        help = 'This specifies boolian for whether to use the biased composition filter. Default is True.')
    parser.add_argument('-e', '--evalue', type=int, default="100",
        help = 'This specifies the evalue threshold for initial HMMER search. Best_hits filters to an evalue of 1 Default is 100.')
    parser.add_argument('-n', '--num_queries', type=int, default="0",
        help = 'This specifies the number of queries to use for HMMER calculation of evalue. Default is 0 (or true database size).')
    parser.add_argument('-x', '--del_temp', action='store_true',
        help = 'Including this will delete temporary files.')
    parser.add_argument('-b', '--use_blast', type=bool, default=False,
        help = 'This specifies boolian for whether to use blast instead of HMMER. Default is False.')
    
    # parser.add_argument('--snakemake_params', type=list, nargs='+')
    args = parser.parse_args()

    threads = args.threads
    workdir = args.directory
    outdir = args.outDir
    biased_composition_filter = args.biased_composition_filter
    evalue = args.evalue
    num_queries = args.num_queries
    del_temp = args.del_temp
    use_blast = args.use_blast

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
        f"del_temp={del_temp}",
        f"use_blast={use_blast}",
        "--rerun-triggers",
        "mtime"
        ]

    print("Running the following command:")
    print(" ".join(cmd))

    try:
        sp.check_output(cmd)

        if args.del_temp:
            print("Deleting temporary files...")
        else: 
            print("Temporary files not deleted.")        

        print("MicSearch finished successfully!\n")
    except sp.CalledProcessError as error:
        print('Encountered error while running snakemake.')
        #print(error.output)

if __name__ == "__main__":
    main()
