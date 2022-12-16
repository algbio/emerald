# Introduction
Emerald is a commandline sequence aligner that explores the suboptimal space and calculates $\alpha$-safety windows: parts of a suboptimal alignment that are contained in an $\alpha$ proportion of all suboptimal alignments.
Emerald takes FASTA cluster files and aligns one selected representative sequence to all the other sequences.
It is possible to
- control the length of the safety windows ( $\alpha$ parameter) and control the size of the suboptimal space ( $\Delta$ parameter)
- use a custom substitution matrix (by default: BLOSUM62) and affine-linear gap score
- approximate large integers with floats for gain of computational speed but loss of accurate output
- use multi threading (only use as many threads as 
- select a custom representative sequence

# First steps
### Downloading the binary
You can download Emerald in the Releases page and run it on the command line.

### Alternatively: Compile from source
Emerald is written in C++ and uses the [gmp library](https://gmplib.org/) for the representation of big integers.
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
Alternatively, it is possible to link the source files and the gmo library directly yourself.

### How to run
Use ``` --help ``` for a first overview of the commands.
#### Some quick examples
Run with $\alpha$ = 0.75 and $\Delta$ = 8.
```
./emerald -f example.fasta -a 0.75 -d 8
```
Use edit distance with only optimal alignments and 100%-safety and multi threading:
```
./emerald -f example.fasta -a 1 -d 0 --costmat edit_distance.txt --gapcost 1 --startgap 0
```
