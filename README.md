# MicSearch
A fully automated workflow for identification of microcins using HMMER.

MicSearch is developed by the [Wilke lab](https://wilkelab.org/) at the [Department of Integrative Biology](https://integrativebio.utexas.edu/) in collaboration with the [Davies lab](https://bwdaviesutaustin.org/) at the [Department of Molecular Biosciences](https://molecularbiosci.utexas.edu/), both at [The University of Texas at Austin](https://www.utexas.edu/).

## Installation

### Setting up conda environment:
```bash
$ conda create --name <environment_name> python=3.8 pip
$ conda activate <environment_name>
```

### Installing dependencies:
```bash
$ conda install mamba -c conda-forge
$ mamba install biopython=1.79 blast=2.14 hmmer=3.3 mafft=7.508 numpy=1.24 pandas=1.5 snakemake=7.18 -c conda-forge -c bioconda 
$ pip install seqhash==1.0 blake3==0.2
```

### Downloading repository:
```bash
$ git clone https://github.com/AaronFeller/MicSearch.git
```
## Running MicSearch:

### Testing your install:
If installed properly, running `python MicSearch.py -h` should give you a list of instructions.

Basic run:
```bash
conda activate <environment_name>
python MicSearch.py -d resources/test_genome -o test_run -t <number_of_threads>
```

Output from this run is contained in the test_genome folder.

### Deleting temp files:
You are free to manually delete temp files all contained in the temp folder inside of your output directory.

Additionally, including the argument `-x` will remove the temp folder. 

After analysis has completed without the `-x` command and you would like to remove the temp fiels, you can rerun the same command with `-x` and they will be removed.
```bash
python MicSearch.py -d resources/test_genome -o test_run -t <number_of_threads> -x
```

# Contributing

MicSearch currently exists as a wrapper to a series of snakemake subroutines, so adding functionality to it is as simple as adding additional subroutines. If there are any subroutines that you see are needed feel free to raise an issue.
