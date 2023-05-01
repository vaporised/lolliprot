library("trackViewer")
library("VariantAnnotation")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("BSgenome.Hsapiens.UCSC.hg38")
library("Homo.sapiens")
library("dplyr")
library("plyranges")
library("drawProteins")
library("mygene")
library("AnnotationHub")
library("parallel")

# Return the abbreviation of a mutation
getAbbrev <- function(ref_aa, loc, type, new_aa) {
  start <- paste("p.", ref_aa, loc, sep="")
  end <- NULL
  
  if (type == "nonsynonymous") {
    end <- new_aa
  } else if (type == "nonsense") {
    end <- "X"
  } else if (type == "frameshift") {
    end <- "fs"
  } else if (type == "synonymous") {
    end <- "="
  } else if (type == "not translated") {
    end <- "?"
  }
  
  abbrev <- paste(start, end, sep = "")
  return(abbrev)
}

# Return the most severe consequence in a list
getConsequence <- function(consequence_list) {
  ordered_type <- c("frameshift", "nonsense", "nonsynonymous", "synonymous", "not translated")
  result <- NA
  consequence_list <- as.character(consequence_list)
  
  for (i in 1:length(ordered_type)) {
    if (ordered_type[i] %in% consequence_list) {
      result <- ordered_type[i]
      break
    }
  }
  
  return(result)
}

# If given condition is satisfied, exit program with error message
stopIfError <- function(message, condition) {
  if (condition) {
    message()
    message(message)
    options(show.error.messages = FALSE)
    stop() 
  }
}

loadVCFs <- function(path, subject, seqSource, cores) {
  chr_only <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
                "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
                "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 
                "chrX", "chrY")
  
  if (endsWith(path, ".rds")) {
    stopIfError("File does not exist", !file.exists(path))
    gr <- readRDS(path)
  } else {
    stopIfError("Directory does not exist", !dir.exists(path))
    setwd(path)
    
    vcf_list <- setdiff(list.files(pattern = "\\.vcf$"), list.dirs(recursive = FALSE, full.names = FALSE))
    stopIfError("No VCFs in directory", length(vcf_list) == 0)
    vcfs <- mclapply(vcf_list, mc.cores=cores, function(x) {
      vcf <- readVcf(x, "hg38")
      seqlevels(vcf) <- chr_only
      
      if (dim(vcf)[1] == 0) {
        return(NULL)
      }
      return(vcf)
    })
    names(vcfs) <- seq_along(vcfs)
    vcfs[sapply(vcfs, is.null)] <- NULL
    
    coding_list <- lapply(vcfs, function(x) predictCoding(x, txdb, Hsapiens))
    gr <- unlist(as(coding_list, "GRangesList"))
  }
  return(gr)
}

# Args (1. gene name 2. path to vcfs or rds path 3. num cores to use for mc functions)
args <- commandArgs(trailingOnly = TRUE)
stopIfError("Wrong number of arguments", (length(args) != 2) & (length(args) != 3))
gene <- args[1] # Gene symbol
path <- args[2] # Path to VCFs or RDS object
cores <- 1
if (length(args) == 3) {
  cores <- args[3]
  stopIfError("Number of cores should be a number greater than 1", (!all.equal(cores, as.integer(cores))) | (cores < 1 ) )
}
current_path <- getwd()

# Get gene data
ah <- AnnotationHub()
edb <- ah[["AH98047"]]
ebd_genes <- genes(edb)
gene_id <- ebd_genes[(ebd_genes$symbol == gene) & (ebd_genes$gene_biotype == "protein_coding") ]$gene_id
stopIfError(paste(gene, "does not exist"), identical(gene_id, character(0)))
transcripts <- transcripts(edb, filter = ~ tx_is_canonical == TRUE)
transcript <- transcripts[transcripts$gene_id == gene_id]$tx_id_version # Transcript ID

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
hsapiens <- BSgenome.Hsapiens.UCSC.hg38
t <- transcripts(txdb, filter=list(tx_name = c(transcript)))
tx_id <- t$tx_id[1] # tx_id
chr <- as.character(t@seqnames) # Chromosome
tx_df <- transcriptLengths(txdb, with.cds_len = TRUE)
length <- (tx_df[tx_df$tx_id == tx_id,]$cds_len / 3) - 1 # Length

gene_gr <- GRanges(chr, IRanges(t@ranges@start, t@ranges@start), names = c(gene))
protein_gr <- GRanges(chr, IRanges(1, length), names = c(gene))
output_filename <- paste(gene, "lolliprot.pdf", sep="_")

# Read vcfs or read existing R file with loaded vcfs
gr <- loadVCFs(path, txdb, Hsapiens, cores)

# Save RDS of gr that can be used later
dir.create(file.path(getwd(), "lolliprot_output"), showWarnings = FALSE)
setwd(file.path(getwd(), "lolliprot_output"))
saveRDS(gr, "lolliprot.rds") 

# Get variants for selected transcript only
txout_gr <- gr[gr$TXID == tx_id]

# Stop if no variants
stopIfError(paste("No mutations in", gene, sep = " "), length(txout_gr) == 0)

# Reshape coding grange
proteinloc <- txout_gr$PROTEINLOC
for (i in 1:length(proteinloc)) {
  if (length(proteinloc[[i]]) > 1) {
    proteinloc[i] <- proteinloc[[i]][1]
  }
}
proteinstart <- IRanges(start = unlist(proteinloc), width = 1, names=names(txout_gr@ranges))
txout_gr@ranges <- proteinstart

# Add nomenclature columns
mcols(txout_gr)$abbrev <- mapply(getAbbrev, ref_aa = txout_gr$REFAA, loc = txout_gr@ranges@start, type = txout_gr$CONSEQUENCE, new_aa = txout_gr$VARAA)
mcols(txout_gr)$abbrev_list <- mapply(function(abbrev) list(abbrev), abbrev = txout_gr$abbrev)

# Combine rows at the same site and variant type (DO NOT USE AA INFO AFTER THIS. 
# METADATA OF COMBINED ROWS ARE UNUSABLE)
site_merged_gr <- txout_gr
site_merged_gr$score <- c(1)
site_merged_gr$temp <- site_merged_gr@ranges@start

# Get rows that have sites with multiple mutations
dupes <- site_merged_gr %>%
  group_by(temp) %>%
  plyranges::filter(n() > 1) %>%
  ungroup()

# Remove these rows from the GRanges
site_merged_gr <- as.data.frame(site_merged_gr)
dupes <- as.data.frame(dupes)
site_merged_gr <- site_merged_gr[!(rownames(site_merged_gr) %in% rownames(dupes)),]
site_merged_gr <- as(site_merged_gr, "GRanges")
dupes <- as(dupes, "GRanges")

# Merge duplicate rows and put back into GRanges
merged <- dupes %>%
  group_by(temp) %>%
  mutate(abbrev_list = paste0(abbrev, collapse=", ")) %>%
  mutate(score = n()) %>% 
  mutate(CONSEQUENCE = getConsequence(CONSEQUENCE)) %>%
  ungroup()
merged <- unique(unname(plyranges::select(merged, temp, CONSEQUENCE, abbrev_list, score)))
site_merged_gr <- plyranges::select(site_merged_gr, temp, CONSEQUENCE, abbrev_list, score)
site_merged_gr <- append(site_merged_gr, merged)

# Change aesthetics of lollipop
final_gr <- site_merged_gr

## Legend
old <- c("synonymous", "nonsynonymous", "nonsense", "frameshift", "not translated")
new <- c("green", "purple", "red", "blue", "grey")
consequences <- as.character(final_gr$CONSEQUENCE)
final_gr$color <- consequences
final_gr$color[final_gr$color %in% old] <- new[match(final_gr$color, old, nomatch=0)]
xaxis <- c(1, length)

## Caterpillar layout
final_gr$SNPsideID <- consequences
new <- c("top", "top", "bottom", "bottom", "top")
final_gr$SNPsideID[final_gr$SNPsideID %in% old] <- new[match(final_gr$SNPsideID, old, nomatch=0)]

## Labels
names(final_gr) <- final_gr$abbrev_list
final_gr$label.parameter.rot <- 55
final_gr$label.parameter.gp <- gpar(fontsize=11)

## Domains
protein_ids <- queryMany(gene, scopes = "symbol", 
                         fields = c("name", "uniprot",  "ensemblgene"), 
                         species = "human", as_dataframe = "True")
uniprot_id <- protein_ids$uniprot.Swiss.Prot
gene_json <- drawProteins::get_features(uniprot_id)
gene_data <- drawProteins::feature_to_dataframe(gene_json)
domains <- gene_data[gene_data$type == "DOMAIN",]
domain_gr <- NULL

# No protein domains
if (nrow(domains) == 0) {
  domain_gr <- GRanges(chr, IRanges(c(1), width=c(1)))
  domain_gr$height <- 0
} else {
  domains$chr <- chr
  domain_gr <- unname(makeGRangesFromDataFrame(domains, start.field = "begin", end.field = "end"))
  domain_gr$fill <- "#a5a5a5"
  domain_gr$height <- 0.01
}

# Plot
pdf(output_filename, width = 16, height = 6)
lolliplot(final_gr, features=domain_gr, ranges = protein_gr, 
          legend="CONSEQUENCE", xaxis=xaxis, yaxis = c(0, max(final_gr$score)), ylab=transcript)
grid.text(gene, x=.5, y=.98, just="top", 
          gp=gpar(cex=1.5, fontface="bold"))
dev.off()
