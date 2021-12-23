# Covid-Mutation-Pipieline

This repository contains scripts to extract mutations of interest from a vcf file.

**Requirements

tqdm>=4.62.2
pandas>=1.3.3
numpy>=1.20.3

**Usage
##To extract mutations of interest
usage: extract_mutations.py [-h] --mutFile MUTFILE --vcfFile VCFFILE
                            [--prefix PREFIX]

Script for extracting variants of interest

optional arguments:
  -h, --help            show this help message and exit
  --mutFile MUTFILE, -m MUTFILE
                        mutation list (csv)
  --vcfFile VCFFILE, -vcf VCFFILE
                        input vcf file
  --prefix PREFIX, -out PREFIX
                        output prefix, default is to data/ _mut and
                        haplist.txt
                        
## To merge GISAID metadata with mutations                        
usage: Merge_meta_and_mutrates.py [-h] --metaFile METAFILE [--prefix PREFIX]
                                  [--out OUT]

Merge metadata with mutation and hap output

optional arguments:
  -h, --help            show this help message and exit
  --metaFile METAFILE, -met METAFILE
                        Input metadata file (.json)
  --prefix PREFIX, -p PREFIX
                        Hapout prefix from previous script
  --out OUT, -o OUT     Output prefix
