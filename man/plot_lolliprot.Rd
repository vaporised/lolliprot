\name{plot_lolliprot}
\alias{plot_lolliprot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plot Lolliprot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Plots a lolliplot of a gene showing sites of amino acid changes, given a VCF file and a gene symbol.
}
\usage{
plot_lolliprot(vcf_path, gene_symbol, remove_mnv = TRUE, to_pdf = FALSE, width = 16, height = 6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcf_path}{
Path to the VCF file containing variants to be plotted
}
\item{gene_symbol}{
Gene symbol of the gene to be plotted
}
\item{remove_mnv}{
If multinucleotide variants should be omitted when plotting. Defaults to TRUE
}
\item{to_pdf}{
If a PDF of the plot should be created instead in \code{./lolliprot_output}. Defaults to FALSE
}
\item{width}{
Width of the PDF. Defaults to 16.
}
\item{height}{
Height of the PDF. Defaults to 6.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
All variants to be plotted should be in one VCF file. If variants across multiple VCFs should be plotted, merge them first.
The canonical transcript is used for plotting. Protein domains are represented as grey boxes. The GRCh38 assembly is used for reference.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
example_path <- system.file("extdata", "example.vcf", package="lolliprot")

# Plot with the multinucleotide variant
plot_lolliprot(vcf_path = example_path, gene_symbol = "DNMT3A", remove_mnv = F)

# Plot without the multinucleotide variant
plot_lolliprot(vcf_path = example_path, gene_symbol = "DNMT3A")

# Create a PDF
plot_lolliprot(vcf_path = example_path, gene_symbol = "DNMT3A", to_pdf = TRUE, width = 15, height = 5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
