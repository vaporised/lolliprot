# lolliprot
R package for drawing lollipop plots of proteins using genetic variant data.

Lolliplots show sites of amino acid changes, given a VCF file and a gene symbol.

![TET2 Lolliprot](https://github.com/vaporised/lolliprot/assets/86297267/fbd33bda-d2ef-4550-9cad-74d06d29b485)

## Pre-processing
All variants to be plotted should be in one VCF file. If variants across multiple VCFs should be plotted, merge them first. This can be done using the `./extdata/merged.sh` file.

## Usage
```R
plot_lolliprot(vcf_path, gene_symbol, remove_mnv = TRUE, to_pdf = FALSE)
```

## Arguments
- **vcf_path**: Path to the VCF file containing variants to be plotted.
- **gene_symbol**: Gene symbol of the gene to be plotted.
- **remove_mnv**: If multinucleotide variants should be omitted when plotting. Defaults to TRUE.
- **to_pdf**: If a PDF of the plot should be created instead in `./lolliprot_output`. Defaults to FALSE.

## Notes
The canonical transcript is used for plotting. Protein domains are represented as grey boxes. The GRCh38 assembly is used for reference.

## Examples
```R
example_path <- system.file("extdata", "example.vcf", package="lolliprot")

# Plot with the multinucleotide variant
plot_lolliprot(vcf_path = example_path, gene_symbol = "DNMT3A", remove_mnv = FALSE)

# Plot without the multinucleotide variant
plot_lolliprot(vcf_path = example_path, gene_symbol = "DNMT3A")

# Create a PDF
plot_lolliprot(vcf_path = example_path, gene_symbol = "DNMT3A", to_pdf = TRUE)
```
## Installation
```R
devtools::install_github("vaporised/lolliprot")
