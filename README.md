# MicSearch
A fully automated workflow for identification of microcins using HMMER.

is developed by the [Wilke lab](https://wilkelab.org/) at the [Department of Integrative Biology](https://integrativebio.utexas.edu/) in collaboration with the [Davies lab](https://bwdaviesutaustin.org/) at the [Department of Molecular Biosciences](https://molecularbiosci.utexas.edu/), both at [The University of Texas at Austin](https://www.utexas.edu/).

## Installation

### Setup conda environment
```bash
$ conda create --name <environment_name> python=3.8 pip
$ conda activate <environment_name>
```

### Install dependencies:
```bash
$ conda install mamba -c conda-forge
$ mamba install biopython=1.79 blast=2.14 hmmer=3.3 mafft=7.508 numpy=1.24 pandas=1.5 snakemake=7.18 -c conda-forge -c bioconda 
$ pip install seqhash==1.0
```


#seqkit=2.3