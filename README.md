# GRMT

GRMT (Generative Reconstruction of Mutation Tree from scratch) is a tool for tumor mutation tree inference from single-cell DNA sequencing data. GRMT is able to efficiently recover the chronological order of mutations and scale well to large datasets. It exploits [BayesOpt] (https://github.com/rmcantin/bayesopt) to estimate error rates in mutation data based on the Bayesian optimization algorithm.

## Requirements

* Linux systems.
* CMake3.0+.
* g++.

## Installation

To build binary, do as follows:

```
tar -zxvf GRMT.tar.gz
cd GRMT
cmake .
make
```

After the installation, the main programs of GRMT are generated in “bin” directory. Type following command if you want to add GRMT to system path:
```
make install
```

## Usage

GRMT uses the genotype matrix derived from single-cell DNA sequencing data to infer mutation tree.

Example:

```
grmt -i testdata/example.txt -o testdata/example -a 0.01 -b 0.2
```

## Input Files

### 1. Genotype Matrix

The mutational profile of single cells are denoted as a genotype matrix. Each row defines the mutation states of a single cell, and each column represents one mutation. Columns are separated by tabs. The genotype matrix can be binary/ternary.

#### a) Binary: Only presence/absence of a mutation is distinguished.
The entry at position [i,j] should be

* 0 if mutation i is not observed in cell j,
* 1 if mutation i is observed in cell j, or
* 3 if the genotype information is missing

#### b) Ternary: Heterozygous and homozygous mutations are distinguished.
The entry at position [i,j] should be

* 0 if mutation i is not observed in cell j,
* 1 if heterozygous mutation i is observed in cell j
* 2 if homozygous mutation i is observed in cell j
* 3 if the genotype information is missing

### 2. Cell names (optional)

A file listing the names of the single cells. Each row specifies the name of a single cell.
If no such file is specified, the cells are numbered from 1 to n.

### 3. Mutation names (optional)

A file listing the names of the mutations. Each row specifies the name of a mutation.
If no such file is specified, the mutations are numbered from 1 to n.

## Output File

### Mutation tree

The mutation tree is written to a file in GraphViz format. The base name of the output file is provided by users.

## Arguments

* `-i, --input <filename>` Replace \<filename\> with the file containing the genotype matrix.

* `-o, --output <string>` Replace \<string\> with the base name of the output file.

## Optional arguments

* `-c, --clabel <filename>` Replace \<filename\> with the path to the file containing the names of the cells.

* `-m, --mlabel <filename>` Replace \<filename\> with the path to the file containing the names of the mutations.

* `-l, --lambda <Double>` Set \<Double\> to a value between 0.5 and 1. This specifies the value of the hyper-parameter \lambda. The default value is 0.7.

* `-K, --kappa <Double>` Set \<Double\> to a positive value. This specifies the value of the hyper-parameter \kappa. The default value is 1.0.

* `-k, --maxl <INT>`  Set \<INT\> to the desired maximum number of times that a mutation can be lost. The default value is 0.

* `-t, --threads <INT>`  Set \<INT\> to the number of threads to use. The default value is 1.

* `-a, --alpha <Double>` Set \<Double\> to estimated false positive rate of the single-cell sequencing experiment.

* `-b, --beta <Double>` Set \<Double\> to estimated false negative rate of the single-cell sequencing experiment.

* `-n, --n_init <INT>`  Set \<INT\> to the desired number of initial data points to sample in the BO algorithm. The default value is 100.

* `-N, --n_iter <INT>`  Set \<INT\> to the desired number of data points to evaluate in the BO algorithm. The default value is 30.

## Contact

If you have any questions, please contact zhyu@nxu.edu.cn.