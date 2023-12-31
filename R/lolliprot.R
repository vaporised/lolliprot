# ----------------------- Helper Functions -----------------------------#

# Return the abbreviation of a mutation
.getAbbrev <- function(ref_aa, loc, type, new_aa) {
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
.getConsequence <- function(consequence_list) {
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
.stopIfError <- function(message, condition) {
  if (condition) {
    message()
    message(message)
    stop(call. = FALSE)
  }
}


# ------------------------- Main -------------------------------#

plot_lolliprot <- function(vcf_path, gene_symbol, remove_mnv = T, to_pdf = F) {
  chr_only <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
                "chrX", "chrY")
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  hsapiens <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

  # Get gene_symbol data
  ah <- AnnotationHub::AnnotationHub()
  edb <- ah[["AH98047"]]
  ebd_genes <- GenomicFeatures::genes(edb)
  gene_id <- ebd_genes[(ebd_genes$symbol == gene_symbol) & (ebd_genes$gene_biotype == "protein_coding") ]$gene_id
  .stopIfError(paste(gene_symbol, "does not exist"), identical(gene_id, character(0)))
  transcripts <- GenomicFeatures::transcripts(edb, filter = ~ tx_is_canonical == TRUE)
  transcript <- transcripts[transcripts$gene_id == gene_id]$tx_id_version # Transcript ID

  t <- GenomicFeatures::transcripts(txdb, filter=list(tx_name = c(transcript)))
  tx_id <- t$tx_id[1] # tx_id
  chr <- as.character(t@seqnames) # Chromosome
  tx_df <- GenomicFeatures::transcriptLengths(txdb, with.cds_len = TRUE)
  length <- (tx_df[tx_df$tx_id == tx_id,]$cds_len / 3) - 1 # Length

  gene_gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(t@ranges@start, t@ranges@start), names = c(gene_symbol))
  protein_gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(1, length), names = c(gene_symbol))

  # Load VCF and predict amino acid changes
  .stopIfError("File given is not a VCF", !endsWith(vcf_path, ".vcf") & !endsWith(vcf_path, ".vcf.gz"))
  .stopIfError("VCF given does not exist", !file.exists(vcf_path))
  vcf <- VariantAnnotation::readVcf(vcf_path, "hg38")
  rownames(vcf) <- make.names(rownames(vcf), unique=TRUE)
  seqlevels(vcf) <- chr_only
  gr <- VariantAnnotation::predictCoding(vcf, txdb, BSgenome.Hsapiens.UCSC.hg38::Hsapiens)

  # Get variants for selected transcript only
  txout_gr <- gr[gr$TXID == tx_id]
  names(txout_gr) <- make.names(names(txout_gr), unique=TRUE)

  # Remove multinucleotide variants
  if (remove_mnv) {
    txout_gr <- txout_gr[txout_gr$REF@ranges@width <= 1 | txout_gr$ALT@unlistData@ranges@width <= 1,]
  }

  # Stop if no variants
  .stopIfError(paste("No mutations in", gene_symbol, sep = " "), length(txout_gr) == 0)

  # Reshape coding grange
  proteinloc <- txout_gr$PROTEINLOC
  for (i in 1:length(proteinloc)) {
    if (length(proteinloc[[i]]) > 1) {
      proteinloc[i] <- proteinloc[[i]][1]
    }
  }
  proteinstart <- IRanges::IRanges(start = unlist(proteinloc), width = 1, names=names(txout_gr@ranges))
  txout_gr@ranges <- proteinstart

  # Add nomenclature columns
  mcols(txout_gr)$abbrev <- mapply(.getAbbrev, ref_aa = txout_gr$REFAA, loc = txout_gr@ranges@start, type = txout_gr$CONSEQUENCE, new_aa = txout_gr$VARAA)
  mcols(txout_gr)$abbrev_list <- mapply(function(abbrev) list(abbrev), abbrev = txout_gr$abbrev)

  # Combine rows at the same site and variant type (DO NOT USE AA INFO AFTER THIS.
  # METADATA OF COMBINED ROWS ARE UNUSABLE)
  site_merged_gr <- txout_gr
  site_merged_gr$score <- c(1)
  site_merged_gr$temp <- site_merged_gr@ranges@start

  # Get rows that have sites with multiple mutations
  dupes <- site_merged_gr %>%
    plyranges::group_by(temp) %>%
    plyranges::filter(plyranges::n() > 1) %>%
    plyranges::ungroup()

  # Remove these rows from the gr, merge, and put back
  if (length(dupes) > 0) {
    site_merged_gr <- as.data.frame(site_merged_gr)
    dupes <- as.data.frame(dupes)
    site_merged_gr <- site_merged_gr[!(rownames(site_merged_gr) %in% rownames(dupes)),]
    site_merged_gr <- as(site_merged_gr, "GRanges")
    dupes <- as(dupes, "GRanges")
    merged <- dupes %>%
      plyranges::group_by(temp) %>%
      plyranges::mutate(abbrev_list = paste0(abbrev, collapse=", ")) %>%
      plyranges::mutate(score = plyranges::n()) %>%
      plyranges::mutate(CONSEQUENCE = .getConsequence(CONSEQUENCE)) %>%
      plyranges::ungroup()
    merged <- unique(unname(plyranges::select(merged, temp, CONSEQUENCE, abbrev_list, score)))
    site_merged_gr <- plyranges::select(site_merged_gr, temp, CONSEQUENCE, abbrev_list, score)
    site_merged_gr <- append(site_merged_gr, merged)
  }

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
  final_gr$label.parameter.gp <- grid::gpar(fontsize=11)

  ## Domains
  protein_ids <- mygene::queryMany(gene_symbol, scopes = "symbol",
                           fields = c("name", "uniprot",  "ensemblgene"),
                           species = "human", as_dataframe = "True")
  uniprot_id <- protein_ids$uniprot.Swiss.Prot
  gene_json <- drawProteins::get_features(uniprot_id)
  gene_data <- drawProteins::feature_to_dataframe(gene_json)
  domains <- gene_data[gene_data$type == "DOMAIN",]
  domain_gr <- NULL

  # If there are no protein domains
  if (nrow(domains) == 0) {
    domain_gr <- GenomicRanges::GenomicRanges(chr, IRanges::IRanges(c(1), width=c(1)))
    domain_gr$height <- 0
  } else {
    domains$chr <- chr
    domain_gr <- unname(GenomicRanges::makeGRangesFromDataFrame(domains, start.field = "begin", end.field = "end"))
    domain_gr$fill <- "#a5a5a5"
    domain_gr$height <- 0.01
  }

  # Plot
  if (to_pdf) {
    output_filename <- paste(gene_symbol, "lolliprot.pdf", sep="_")

    dir.create(file.path(getwd(), "lolliprot_output"), showWarnings = FALSE)
    setwd(file.path(getwd(), "lolliprot_output"))
    pdf(output_filename, width = 16, height = 6)
    trackViewer::lolliplot(final_gr, features=domain_gr, ranges = protein_gr,
                           legend="CONSEQUENCE", xaxis=xaxis, yaxis = c(0, max(final_gr$score)), ylab=transcript)
    grid::grid.text(gene_symbol, x=.5, y=.98, just="top", gp=grid::gpar(cex=1.5, fontface="bold"))
    dev.off()
  } else {
    trackViewer::lolliplot(final_gr, features=domain_gr, ranges = protein_gr,
                           legend="CONSEQUENCE", xaxis=xaxis, yaxis = c(0, max(final_gr$score)), ylab=transcript)
    grid::grid.text(gene_symbol, x=.5, y=.98, just="top", gp=grid::gpar(cex=1.5, fontface="bold"))
  }
}

