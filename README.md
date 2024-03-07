# MicSearch
A fully automated workflow for identification of microcins using HMMER.

MicSearch (Cinful v2) is developed by the [Wilke lab](https://wilkelab.org/) at the [Department of Integrative Biology](https://integrativebio.utexas.edu/) in collaboration with the [Davies lab](https://bwdaviesutaustin.org/) at the [Department of Molecular Biosciences](https://molecularbiosci.utexas.edu/), both at [The University of Texas at Austin](https://www.utexas.edu/).

## Installation

### Setting up conda environment:
```bash
conda create --name <environment_name> python=3.8
conda activate <environment_name>
```

### Installing dependencies:
```bash
conda install mamba -c conda-forge
mamba install biopython=1.79 blast=2.14 hmmer=3.3 mafft=7.508 numpy=1.24 pandas=1.5 snakemake=7.18 -c conda-forge -c bioconda
```

### Downloading repository:
```bash
git clone https://github.com/AaronFeller/MicSearch.git
```
## Running MicSearch:

### Testing your install:
If downloaded properly, running `python MicSearch.py -h` should give you a list of instructions.

Basic run:
```bash
cd MicSearch
conda activate <environment_name>
python MicSearch.py -d resources/test_genome -o test_run -t <number_of_threads>
```

The precomputed results from this test are contained in the test_genome folder for your comparison.

### File naming required
This pipeline will only identify files that end in the extension '.fna'.
If you are working with another extension, the simplest fix is to rename these files and add the .fna extension.

### Deleting temp files:
You are able to manually delete temp files contained in the temp folder inside of your output directory.

Additionally, including the argument `-x` will remove the temp folder at runtime. 

After analysis has completed without the `-x` command and you would like to remove the temp files, you can rerun the same command with `-x` and they will be removed.
```bash
python MicSearch.py -d resources/test_genome -o test_run -t <number_of_threads> -x
```

# Contributing

MicSearch currently exists as a wrapper to a series of snakemake rules, so adding functionality to it is as simple as adding additional rules. If there are any additions that you feel are needed, feel free to raise an issue.
