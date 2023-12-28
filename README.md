# lolliprot
A script for drawing lollipop plots of proteins using genetic variant data.

This script draws lolliplots of genes showing sites of amino acid changes, given a VCF file and a gene symbol.
Plots are saved into a `lolliprot_output` folder.


![TET2 lolliprot](https://github.com/vaporised/lolliprot/blob/main/data/TET2_lolliprot.png)

## Pre-processing
All variants to be plotted should be in one VCF file. If variants across multiple VCFs should be plotted, merge them first.

## Usage
```
Usage:
   Rscript lolliprot.R vcf_path gene

Arguments:
   vcf_path       Path to the VCF containing the sequence variants
   gene           Gene symbol of the gene to be plotted. 
                  The program uses the protein isoform encoded by the canonical transcript. 
```

## Notes
Multinucleotide variants are removed before plotting. The canonical transcript is used for plotting. Protein domains are represented as grey boxes.
