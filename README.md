# EMERALD manual
[![Anaconda-Server Badge](https://anaconda.org/bioconda/emerald/badges/version.svg)](https://anaconda.org/bioconda/emerald)
[![License](https://img.shields.io/badge/licence-GPLv3-blue)](https://www.gnu.org/licenses/gpl-3.0.html)

1. [Introduction](#sec1)
2. [Installation](#sec2)\
2.1. [Installing from conda](#sec2.1)\
2.2. [Compile from source](#sec2.2)
3. [Running EMERALD](#sec3)\
3.1. [EMERALD input](#sec3.1)\
3.2. [Command line options](#sec3.2)\
3.3. [EMERALD output](#sec3.3)\
3.4. [Example](#sec3.4)

<a name="sec1"></a>
# Introduction
EMERALD is a command line protein sequence aligner that explores the suboptimal space and calculates **$\alpha$-safety windows**: partial alignments that are contained in an $\alpha$ proportion of all suboptimal alignments.
EMERALD takes FASTA cluster files and aligns one selected representative sequence to all the other sequences.\
EMERALD's features include
- using custom substitution matrices (by default: BLOSUM62) and affine-linear gap score
- approximating large integers with floats for gain of computational speed but loss of accurate output
- multi threading
- selecting a custom representative sequence

<a name="overview"></a>
![EMERALD overview](./figs/emerald_overview.png)
Schematic representation of EMERALD’s safety window calculation

<a name="sec2"></a>
# Installation
EMERALD is already compiled for Linux and Mac OS silicon.
You can download the EMERALD binary in the ![Releases page](https://github.com/algbio/emerald/releases) and run it on the command line.

<a name="sec2.1"></a>
### Conda installation
EMERALD can be installed via conda:
```
 conda install -c bioconda emerald 
 ```

<a name="sec2.2"></a>
### Compile from source
EMERALD is written in C++ and uses the [gmp library](https://gmplib.org/) for the representation of big integers.
Additionally, [cmake](https://cmake.org/install/) is needed for the compilation.
After installing gmp and downloading the source, navigate to its main directory and run
```
cmake .
```
followed by
```
make
```
to compile.

<a name="sec3"></a>
# Running EMERALD
Use ``` --help ``` for a first overview of the commands.

<a name="sec3.1"></a>
### EMERALD input
EMERALD expects `.fasta` cluster files of protein sequences.\
EMERALD defines two kinds of sequences: the singular `representative sequence` and `cluster members` for all the other sequences. The representative sequence is aligned with all the cluster members, resulting in $n-1$ alignments for a cluster of size $n$.

<a name="sec3.2"></a>
### Command line options
##### The basic options are the following
- ```-f, --file {FILE}``` Path to FASTA file, mandatory argument.
- ```-a, --alpha {value}``` $\alpha$ value for safety, $0.5 < \alpha \leq 1$, by default: 0.75. The safety windows will be partial alignments contained in an $\alpha$ proportion of all alignments. If $\alpha$ is chosen outside this range, a warning will be displayed. EMERALD will keep running but it can crash.
- ```-d, --delta {value}``` $\Delta$ value for the size of the suboptimal space, any positive integer, by default: 0. If $\Delta$ is larger, more alignments will be considered suboptimal, which will decrease the number and lengths of the safety windows.
- ```-i, --threads {value}``` How many threads to use. By default 1 thread is used.
- ```-r, --reference {sequence}``` Select a specific sequence as representative sequence by some unique identitifer in the sequence description. By default the first sequence in the cluster will be the representative.
##### More advanced options
- ```-c, --costmat {file}``` This file is a lower triangular matrix C which for which C[a][b] is the aligning score of the amino acids a and b. The amino acids are given in the following order:
```Ala  Arg  Asn  Asp  Cys  Gln  Glu  Gly  His  Ile  Leu  Lys  Met  Phe  Pro  Ser  Thr  Trp  Tyr  Val```. Examples are given in the [utils](./utils) directory.
- ```-s, --special {value}``` is an integer assigned to the score of aligned amino acids in which one of the two is not included in the list above.
- ```-g, --gapcost {value}``` and ```-e, --startgap {value}``` Defines the affine-linear gap score function, by default -1 and -11, respectively.

<a name="sec3.3"></a>
### EMERALD output
EMERALD's output is in stdout. The first part of the output is the following.
```
Representative sequence description
Representative sequence
Number of aligned sequence pairs
```
Following for every aligned sequence pair:
```
Cluster sequence description
Cluster sequence
Number of safety windows
```
Finally, every safety window will be printed in a separate line: $L_0\,R_0\,L_1\,R_1$, first for the representative sequence $[L_0, R_0)$ and then for the cluster sequence $[L_1, R_1)$.\
Safety windows are half open intervals, the left index is inclusive and the right index is exclusive, and indexing starts at 0.

<a name="sec3.4"></a>
### Example
`examples/ex1.fasta` (same as in the [Overview](#overview)):
```
>Representative sequence
MSFDLKSKFLG
>Cluster member 1
MSKLKDFLFKS
>Cluster member 2
MSLGSFKDKFL
>Cluster member 3
MSLKDKKFLKS
>Cluster member 4
MSFLKKKFDSL
```
Output:
```
$ ./emerald -f examples/ex1.fasta -a 0.75 -d 8
>Representative sequence
MSFDLKSKFLG
5
>Cluster member 1
MSKLKDFLFKS
3
0 2 0 2
4 6 3 5
8 11 8 11
>Cluster member 2
MSLGSFKDKFL
2
0 3 0 3
4 9 5 10
>Cluster member 3
MSLKDKKFLKS
2
0 2 0 2
7 10 6 9
>Cluster member 4
MSFLKKKFDSL
2
0 3 0 3
5 9 4 8
```
