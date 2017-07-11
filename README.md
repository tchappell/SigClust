# SigClust
k-means clustering for biological sequences using topological signatures

## Usage
`./SigClust (options) [fasta input]`

## Options
* -sw [signature width (default = 256)]
* -k [kmer length (default = 5)]
* -d [signature density (default = 0.0476..)]
* -c [cluster count (default = 1000)]
* -i [k-means iterations (default = 4)]
* --fasta-output

## Requirements

A version of gcc (or compatible compiler) with support for OpenMP and `__builtin_popcountll()`. 

## Installation

Type `make` to install the software. If make or gcc is not available, this software only consists of a single source file which can be compiled manually. If OpenMP is not available the `#pragma` directives can be ignored to build a single-threaded version of the software. If `__builtin_popcountll()` is not available, it can be replaced with whatever builtin is needed to emit a 64-bit `POPCNT` instruction with your compiler and architecture.

## Operation

Run SigClust, passing it the fasta file containing the sequences to be clustered and any options needed. It will then produce a series of clusters as output to standard output, which can then be redirected as necessary. Note that both the sequences and the signatures generated for each sequence are stored in memory during clustering.

SigClust is multithreaded with OpenMP. The number of threads used can hence be controlled with OpenMP environment variables such as `OMP_NUM_THREADS`.

## Detailed description of options

### -sw [signature width]

This is the number of bits to use to store the signatures that are generated to represent each sequence. This value must be a multiple of 64.

The signature width affects the representational capacity of the signatures, and when sequences are long and/or long k-mer lengths are used to represent the signatures, a larger signature width can be beneficial. This comes with a cost to signature generation and clustering speed.

### -k [kmer length]

This is the number of characters making up each k-mer. k-mers of the given length are created for every position within each sequence to represent that sequence for the purposes of signature construction. For example, with a k-mer length of 5, the sequence `ACAAGATGCCATTG` results in the following k-mers: `ACAAG CAAGA AAGAT AGATG GATGC ATGCC TGCCA GCCAT CCATT CATTG`.

The k-mer length plays an important role in both the representational capacity of the generated signatures and the fragility of those sequences, and should only be tweaked with testing to ensure that the new value improves performance.

### -d [signature density]

When signatures are created from sequences, each k-mer is hashed into a k-mer signature of the density provided to this parameter. The k-mer signatures are then combined to create a signature for the sequence. The signature density determines the portion of bits set in each k-mer signature. This value may need to be tweaked based on the length of the signatures and/or the k-mer length.

### -c [cluster count]

The number of clusters to be created by SigClust. This value should be lower than the number of unique sequences provided to SigClust or the clustering will not work. Larger cluster counts will require more clustering time.

### -i [k-means iterations]

The number of iterations of k-means to run. Each iteration will consume additional time and k-means rarely improves much after four iterations; however, based on your available computational resources you may want to increase or decrease this number.

### --fasta-output

By default SigClust will produce a two-column CSV consisting of the sequence ID and cluster ID of each sequence. An alternative output is available by passing in this parameter; instead, SigClust will produce a fasta-format file containing the same sequences passed in, but with the name of each sequence replaced with the cluster number that sequence is a part of.

## Examples of operation

Given the fasta file `test.fasta`, containing the following data:
```
>a
AAAAAAAAAACAAAAAACAAAAACAACAACAA
>b
ACAAAAAACAAAAAAACAAAAAAAACAAACAAAAA
>c
GGGGTGGGGGAGGGGGATGGGGGGGGTGGGGGGG
>d
GGAGGGGGGGCGGGGGGGTGGGGGGAGGGGCGGGGTGG
>e
AGAGAAAGAGCGAGGGAGAGGAGGAAGGAGCAGGATGG
```

Running SigClust to split these five sequences into two clusters with the following command:
```
./SigClust -c 2 test.fasta > test.clusters
```
...should produce the following output in the `test.clusters` file:
```
0,1
1,1
2,0
3,0
4,0
```
The format of the default output is a two-column CSV, the first column containing the ID of the sequence (starting from 0 and based on the order the sequences were in in the original file) and the second column containing the ID of the cluster that sequence was placed int (again, starting from 0). Here, the first two sequences have been placed into one cluster (#1) and the latter three sequences into another cluster (#0). The numbering of the individual clusters is unimportant; all that is important is the grouping.

To instead create three clusters, the following command is used instead:
```
./SigClust -c 3 test.fasta > test.clusters
```
...producing this output:
```
0,1
1,1
2,2
3,2
4,0
```
If fasta output is desired:
```
./SigClust -c 3 --fasta-output test.fasta > test.clusters
```
...producing this output:
```
>1
AAAAAAAAAACAAAAAACAAAAACAACAACAA
>1
ACAAAAAACAAAAAAACAAAAAAAACAAACAAAAA
>2
GGGGTGGGGGAGGGGGATGGGGGGGGTGGGGGGG
>2
GGAGGGGGGGCGGGGGGGTGGGGGGAGGGGCGGGGTGG
>0
AGAGAAAGAGCGAGGGAGAGGAGGAAGGAGCAGGATGG
```
