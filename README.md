# VCF-RLZ

## What is this about?

### TL;DR

Implementation of a VCF to RLZ conversion module, which is intended for searching patterns inside genomic collections, reporting all the occurrences, characterized by (Sample ID, Chromosome, Allele, Position within allele).

### Summary

In order to add searches on genomic collections in VCF (Variant Calling Format), we present a conversion module which takes the information contained in VCF and, without decompressing it, builds RLZ [1] (Relative Lempel-Ziv), a compressed index for repetitive collections.

The implemented solution consists of three VCF processing stages: edits characterization, characterization clustering and interpretation for RLZ build process. In parallel, from the interpretation, we generate an extra compact data structure that will allow us to support the transformation of PLZ occurrence's position to VCF. This module was validated on different datasets generated from the VCFs published by the 1000 Genomes project, with the objective of evaluate the conversion time and volume of the products. Processing up to 72 human genomes in 6 hours, and leaving the information ready to be consumed by RLZ.

However, it was found that RLZ still requires the original sequence to be decompressed, which disqualifies it for use with massive collections. Therefore, we also propose an alternative RLZ construction format, which does not require access to the original text.

Full detailed explanation is [here](https://repositorio.uchile.cl/handle/2250/191843) (Spanish).

## How does it work?

As is described in the previous section, this process consist of three main steps, those are PARSE, SORT and BUILD. In order to follow these steps, please open the [Construction script](https://github.com/SanquirinoB/VCF_RLZ/blob/main/experiments/ConstructionTime.cpp).

### [PARSE](https://github.com/SanquirinoB/VCF_parsing/tree/master)

Takes as input the reference (FASTA format) and the VCF collection. From the reference, recovers all the bases for each chromosome described, then the main reference is the concatenation of each chromosome (Named R), and save the relative position within the main reference of each one. 
After that, for each VCF, we discard the meta-information lines, jumping directly into the registries. 
For each line, we identify the edit and for each individual which has this edit, we create a phrase with the form: (C, A, ID, P, L, PE, LE), where: C as the chromosome number, A as the number of allele, ID the number of sample, P the position where the edit will be applied, L the length of the section to be replaced, PE the position within the reference where the edit appears and LE the length of the edit.

### [SORT](https://github.com/SanquirinoB/VCF_RLZ/blob/main/src/VCFParsingSorter.cpp)

We need to virtually create a single string with all the genomes. In order to do that, we define the following structure: A chromosome can be described as the concatenation of its alleles, then a genome can be described as the concatenation of its chromosomes (Named S). In that way, a genomic collection can be defined as a concatenation of each genome.   Given by the phrases generated by PARSE, we'll sort them by the following order: first by sample, then by chromosome, after that by allele and finally for position. This step is executed with an external memory sort algorithm provided by STXXL.  

### [BUILD](https://github.com/SanquirinoB/VCF_RLZ/blob/main/src/VCFParsingInterpreter.cpp)

Finally, with the sorted phrases, we'll construct LZ(S|R), the Lempel-Ziv factorization of S based on R. A minimal example is: for each pair of phrases we create one factor (p, l)  with p as the position within the reference, and l the length of the kept section (Green), then another factor (PE, LE) representing the edit. After that, given by the second phrase we create a fill factor representing the kept section between both edits (Yellow), and then again an edit factor, ending with a end factor (Red). With LZ(S|R), R and the meta data captured, we construct RLZ, getting our index, allowing us to search patterns, retrieving all occurrences characterized by sample, chromosome, allele and position of occurrence within the allele.

## Usage

Work in progress, , but at the sime time it has a restricted performance issue which doesn't allow you to index all that info in reasonable time.

I'll provide a new release which support a single chromosome instead of a collection. If you want it, you can send me a email to fsanchiricob@gmail.com :)

Any way you can use PARSE, SORT and partially BUILD for everything you want, the conversion of VCF into its Lempel-Ziv factorization is really useful for  another structures construction

### Build Index

```console
foo@bar:~$./build/fullbuild <DEST_PATH> <FASTA_PATH> -n <N_CHR> <VCF_PATH> ...
```

#### Parameters
1. <DEST_PATH>: Path to destination folder, where the binaries files will be placed.
2. <FASTA_PATH>: Path to reference file in FASTA format (.fa).
3. <N_CHR>: Amounf VCF files (autosomic chromosomes) to provide (.vcf).
4. <VCF_PATH>: Path to VCF files, must match with the N_CHR provided.

#### Output
Inside DEST_PATH youĺl finde the binary files that represent the genomic collection provided. These files can be used later to search.

## To Do
- Only support for autosomal chromosomes. Theorically X and Y have a different interpretation in VCF, open for assistance.
- Support to select a single chromosome instead of digest the full collection.
