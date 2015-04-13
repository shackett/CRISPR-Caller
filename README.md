CRISPR-Caller
=============

A method for confirming CRISPR knockouts of genes from sanger sequencing data. Inputs are:
- a gene name (either ensembl or hgnc), only tested for human
- an .ab1 sanger trace file

Two folders are used:
- ABIfiles
  - contains all .ab1 files of interest
  - gene names are associated with .ab1 files as the second field, seperated by underscore (e.g. xxxx_YFG_yyy.ab1) will
  analyze a gene YFG, while output will be tagged with xxx_YFK_yyy
- Output
  - for storing results

General approach is:
- call homozygous and heterozygous region of sanger sequence
- identify indels of the first allele
- subtract first allele from degenerate call
- identify indels of the second allele
