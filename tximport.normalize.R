# tximport.Combine.to.GCT
suppressMessages(suppressWarnings(library("devtools")))
suppressMessages(suppressWarnings(library("tools")))

suppressMessages(suppressWarnings(install.packages("foreign", repos = "https://cloud.r-project.org/",
 quiet = TRUE)))
suppressMessages(suppressWarnings(install.packages("Hmisc", repos = "https://cloud.r-project.org/",
 quiet = TRUE)))

suppressMessages(suppressWarnings(install.packages("getopt", repos = "https://cloud.r-project.org/",
 quiet = TRUE)))
suppressMessages(suppressWarnings(install.packages("optparse", repos = "https://cloud.r-project.org/",
 quiet = TRUE)))

if (!requireNamespace("BiocManager", quietly = TRUE)) suppressMessages(suppressWarnings(install.packages("BiocManager",
 repos = "https://cloud.r-project.org/", quiet = TRUE)))

suppressMessages(suppressWarnings(BiocManager::install("tximport", quiet = TRUE)))
suppressMessages(suppressWarnings(BiocManager::install("DESeq2", quiet = TRUE)))

suppressMessages(suppressWarnings(library("tximport")))
suppressMessages(suppressWarnings(library("DESeq2")))

suppressMessages(suppressWarnings(library("getopt")))
suppressMessages(suppressWarnings(library("optparse")))

arguments <- commandArgs(trailingOnly = TRUE)

option_list <- list(make_option("--quant", dest = "Quantifications"), make_option("--type",
 dest = "Quant.Type"), make_option("--txdb", dest = "Transcriptome.Database",
 default = NULL), make_option("--info", dest = "Sample.Info", default = NULL),
 make_option("--counts", dest = "Output.Normalized.Counts"), make_option("--tpm",
  dest = "Output.TPM"), make_option("--basename", dest = "output.file.base"),
 make_option("--annotate", dest = "annotate"), make_option("--reverse", dest = "reverse"),
 make_option("--split", dest = "Split.Identifiers"), make_option("--seed", dest = "random.seed",
  type = "integer", default = 779948241))

opt <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE,
 args = arguments)$options

quanttype = as.character(opt$Quant.Type)
counts = as.logical(opt$Output.Normalized.Counts)
TPM = as.logical(opt$Output.TPM)
outfile = as.character(opt$output.file.base)
split = as.logical(opt$Split.Identifiers)
seed = as.integer(opt$random.seed)
annotate = as.logical(opt$annotate)
reverse = as.logical(opt$reverse)

if (opt$Transcriptome.Database != "") {
 transcriptome = opt$Transcriptome.Database
 suppressMessages(suppressWarnings(BiocManager::install("GenomicFeatures", quiet = TRUE)))
 suppressMessages(suppressWarnings(library("GenomicFeatures")))
 TxDb <- suppressMessages(suppressWarnings(makeTxDbFromGFF(file = transcriptome)))
 k <- suppressMessages(suppressWarnings(keys(TxDb, keytype = "TXNAME")))
 tx2gene <- suppressMessages(suppressWarnings(select(TxDb, k, "GENEID", "TXNAME")))
}

set.seed(seed)

tximportfiles <- readLines(opt$Quantifications, warn = FALSE)

if (quanttype == "RSEM") {
 names(tximportfiles) <- gsub(".gz$", "", basename(tximportfiles))
 names(tximportfiles) <- gsub(".genes.results", "", names(tximportfiles))
 txi <- suppressMessages(suppressWarnings(tximport(tximportfiles, type = "rsem",
  txIn = FALSE, txOut = FALSE)))
 tpmmatrix <- txi$abundance
}

if (quanttype == "SALMON") {
 names(tximportfiles) <- gsub(".gz$", "", basename(tximportfiles))
 names(tximportfiles) <- gsub(".quant.sf", "", names(tximportfiles))
 txi <- suppressMessages(suppressWarnings(tximport(tximportfiles, type = "salmon",
  tx2gene = tx2gene, ignoreAfterBar = TRUE)))
 tpmmatrix <- txi$abundance
}

if (quanttype == "SAILFISH") {
 names(tximportfiles) <- gsub(".gz$", "", basename(tximportfiles))
 names(tximportfiles) <- gsub(".quant.sf", "", names(tximportfiles))
 txi <- suppressMessages(suppressWarnings(tximport(tximportfiles, type = "sailfish",
  tx2gene = tx2gene, ignoreAfterBar = TRUE)))
 tpmmatrix <- txi$abundance
}

if (quanttype == "KallistoH5") {
 BiocManager::install("rhdf5")
 library("rhdf5")
 names(tximportfiles) <- gsub(".gz$", "", basename(tximportfiles))
 names(tximportfiles) <- gsub(".abundance.h5", "", names(tximportfiles))
 txi <- suppressMessages(suppressWarnings(tximport(tximportfiles, type = "kallisto",
  tx2gene = tx2gene, ignoreAfterBar = TRUE)))
 tpmmatrix <- txi$abundance
}

if (quanttype == "KallistoTSV") {
 names(tximportfiles) <- gsub(".gz$", "", basename(tximportfiles))
 names(tximportfiles) <- gsub(".abundance.tsv", "", names(tximportfiles))
 txi <- suppressMessages(suppressWarnings(tximport(tximportfiles, type = "kallisto",
  tx2gene = tx2gene, ignoreAfterBar = TRUE)))
 tpmmatrix <- txi$abundance
}


# DESeq2 Normalize/Quant
txi$length[txi$length == 0] <- 1


coldata <- as.data.frame(colnames(txi$counts), stringsAsFactors = TRUE, header = FALSE)
colnames(coldata) <- c("EXPERIMENT")
rownames(coldata) <- colnames(txi$counts)
if (opt$Sample.Info == "") {
 dds <- suppressMessages(suppressWarnings(DESeqDataSetFromTximport(txi, colData = coldata,
  ~1)))
} else if (opt$Sample.Info != "") {
 info <- read.table(text = gsub(",", "\t", readLines(opt$Sample.Info, warn = FALSE)),
  stringsAsFactors = TRUE)
 if (dim(info)[1] == 1 + dim(coldata)[1]) {
  info <- read.table(text = gsub(",", "\t", readLines(opt$Sample.Info, warn = FALSE)),
   stringsAsFactors = TRUE, header = TRUE)
 }
 coldata <- merge(x = coldata, y = info, by.x = 1, by.y = 1)
 rownames(coldata) <- coldata[, 1]
 levels <- paste(rev(colnames(coldata[2:length(coldata)])), collapse = " + ")
 design = as.formula(paste("~", levels))

 txi$abundance <- txi$abundance[, rownames(coldata)]
 txi$counts <- txi$counts[, rownames(coldata)]
 txi$length <- txi$length[, rownames(coldata)]
 tpmmatrix <- txi$abundance
 if (txi$countsFromAbundance != "no") {
  txi$countsFromAbundance <- txi$countsFromAbundance[, rownames(coldata)]
 }

 dds <- suppressMessages(suppressWarnings(DESeqDataSetFromTximport(txi, colData = coldata,
  design = design)))
 dds <- suppressMessages(suppressWarnings(DESeq(dds)))
 if (reverse == FALSE) {
  # Compute Log2FC using the first factor vs. the second factor
  res <- results(dds, contrast = c(colnames(coldata)[2], levels(coldata[, c(2)])[1],
   levels(coldata[, c(2)])[2]))
 } else if (reverse == TRUE) {
  # Compute Log2FC using the second factor vs. the first factor
  res <- results(dds, contrast = c(colnames(coldata)[2], levels(coldata[, c(2)])[2],
   levels(coldata[, c(2)])[1]))
 }
 results <- as.data.frame(res, stringsAsFactors = FALSE)
 colnames(results) <- res@elementMetadata@listData$description

 if (split == TRUE) {
  splitids <- rownames(results)
  splitids <- sub("NM_", "NM.", splitids)
  mat <- do.call("rbind", strsplit(sub("_", ";", splitids), ";"))
  mat[, 1] <- sub("NM.", "NM_", mat[, 1])
  mat <- as.data.frame(mat, stringsAsFactors = FALSE)
  if (dim(mat)[2] == 1) {
   mat <- as.data.frame(cbind(mat, mat), stringsAsFactors = FALSE)
  }
  if (length(grep("^ENS", mat[, 1])) == length(mat[, 1])) {
   mat2 <- as.data.frame(do.call("rbind", strsplit(sub("[.]", ";", mat[,
    1]), ";")))
   mat[, 1] = mat2[, 1]
  }

  colnames(mat) <- c("NAME", "Description")

 } else if (split == FALSE) {
  mat <- as.data.frame(rownames(results), stringsAsFactors = FALSE)
  mat <- as.data.frame(cbind(mat, Description = "na"), stringsAsFactors = FALSE)

  if (length(grep("^ENS", mat[, 1])) == length(mat[, 1])) {
   mat2 <- as.data.frame(do.call("rbind", strsplit(sub("[.]", ";", mat[,
    1]), ";")))
   mat[, 2] = mat[, 1]
   mat[, 1] = mat2[, 1]
  }

  colnames(mat) <- c("NAME", "Description")
 }

 if (annotate == FALSE) {
  degmatrix <- as.data.frame(cbind(mat, results), stringsAsFactors = FALSE)
  colnames(degmatrix)[1] <- "Gene_ID"
  colnames(degmatrix)[2] <- "Details"
  degmatrix <- degmatrix[order(degmatrix$"BH adjusted p-values", decreasing = FALSE),
   ]
  suppressWarnings(write.table(degmatrix, paste0(outfile, ".Differential.Expression.txt"),
   sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE))
 } else if (annotate == TRUE) {
  if (opt$Transcriptome.Database != "") {
   gtf_db <- as.data.frame(rtracklayer::import(transcriptome), stringsAsFactors = FALSE)
   gtf_db2 <- unique(gtf_db[, c("gene_id", "gene_name")])
   degmatrix <- as.data.frame(cbind(mat, results), stringsAsFactors = FALSE)
   degmatrix_ann <- merge(x = gtf_db2, y = degmatrix, by.x = 1, by.y = 0,
    all.x = FALSE, all.y = TRUE)
   degmatrix_ann <- degmatrix_ann[, -c(1)]
   colnames(degmatrix_ann)[1] <- "Gene_Symbol"
   colnames(degmatrix_ann)[2] <- "Gene_ID"
   colnames(degmatrix_ann)[3] <- "Gene_Details"
   degmatrix_ann <- degmatrix_ann[order(degmatrix_ann$"BH adjusted p-values",
    decreasing = FALSE), ]
   suppressWarnings(write.table(degmatrix_ann, paste0(outfile, ".Differential.Expression.txt"),
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE))
  } else {
   message("No GTF provided. Unable to annotate results with gene symbols.")
   degmatrix <- as.data.frame(cbind(mat, results), stringsAsFactors = FALSE)
   colnames(degmatrix)[1] <- "Gene_ID"
   colnames(degmatrix)[2] <- "Details"
   degmatrix <- degmatrix[order(degmatrix$"BH adjusted p-values", decreasing = FALSE),
    ]
   suppressWarnings(write.table(degmatrix, paste0(outfile, ".Differential.Expression.txt"),
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE))
  }
 }

 write.table(paste(c(dim(coldata)[1], length(unique(coldata[, c(2)])), "1"), collapse = " "),
  paste0(outfile, ".cls"), sep = "\t", row.names = FALSE, col.names = FALSE,
  quote = FALSE)
 write.table(paste("#", paste(unique(coldata[, c(2)]), collapse = " "), collapse = " "),
  paste0(outfile, ".cls"), sep = "\t", row.names = FALSE, col.names = FALSE,
  quote = FALSE, append = TRUE)
 write.table(paste(coldata[, c(2)], collapse = " "), paste0(outfile, ".cls"),
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}

keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep, ]
dds <- suppressMessages(suppressWarnings(estimateSizeFactors(dds)))
normCounts <- suppressMessages(suppressWarnings(counts(dds, normalized = TRUE)))

if (split == TRUE) {
 splitids <- rownames(normCounts)
 splitids <- sub("NM_", "NM.", splitids)
 mat <- do.call("rbind", strsplit(sub("_", ";", splitids), ";"))
 mat[, 1] <- sub("NM.", "NM_", mat[, 1])
 mat <- as.data.frame(mat, stringsAsFactors = FALSE)
 if (dim(mat)[2] == 1) {
  mat <- as.data.frame(cbind(mat, mat), stringsAsFactors = FALSE)
 }
 if (length(grep("^ENS", mat[, 1])) == length(mat[, 1])) {
  mat2 <- as.data.frame(do.call("rbind", strsplit(sub("[.]", ";", mat[, 1]),
   ";")))
  mat[, 1] = mat2[, 1]
 }

 colnames(mat) <- c("NAME", "Description")

} else if (split == FALSE) {
 mat <- as.data.frame(rownames(normCounts), stringsAsFactors = FALSE)
 mat <- as.data.frame(cbind(mat, Description = "na"), stringsAsFactors = FALSE)

 if (length(grep("^ENS", mat[, 1])) == length(mat[, 1])) {
  mat2 <- as.data.frame(do.call("rbind", strsplit(sub("[.]", ";", mat[, 1]),
   ";")))
  mat[, 2] = mat[, 1]
  mat[, 1] = mat2[, 1]
 }

 colnames(mat) <- c("NAME", "Description")
}

if (counts == TRUE) {
 countmatrix <- as.data.frame(cbind(mat, normCounts), stringsAsFactors = FALSE)
 write.table("#1.2", paste0(outfile, ".Normalized.Counts.gct"), row.names = FALSE,
  col.names = FALSE, quote = FALSE)
 write.table(t(as.data.frame(dim(normCounts))), paste0(outfile, ".Normalized.Counts.gct"),
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
 suppressWarnings(write.table(countmatrix, paste0(outfile, ".Normalized.Counts.gct"),
  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE))
}

if (TPM == TRUE) {
 # TPM
 if (split == TRUE) {
  splitids <- rownames(tpmmatrix)
  splitids <- sub("NM_", "NM.", splitids)
  mat <- do.call("rbind", strsplit(sub("_", ";", splitids), ";"))
  mat[, 1] <- sub("NM.", "NM_", mat[, 1])
  mat <- as.data.frame(mat, stringsAsFactors = FALSE)
  if (dim(mat)[2] == 1) {
   mat <- as.data.frame(cbind(mat, mat), stringsAsFactors = FALSE)
  }
  if (length(grep("^ENS", mat[, 1])) == length(mat[, 1])) {
   mat2 <- as.data.frame(do.call("rbind", strsplit(sub("[.]", ";", mat[,
    1]), ";")))
   mat[, 1] = mat2[, 1]
  }

  colnames(mat) <- c("NAME", "Description")

 } else if (split == FALSE) {
  mat <- as.data.frame(rownames(tpmmatrix), stringsAsFactors = FALSE)
  mat <- as.data.frame(cbind(mat, Description = "na"), stringsAsFactors = FALSE)

  if (length(grep("^ENS", mat[, 1])) == length(mat[, 1])) {
   mat2 <- as.data.frame(do.call("rbind", strsplit(sub("[.]", ";", mat[,
    1]), ";")))
   mat[, 2] = mat[, 1]
   mat[, 1] = mat2[, 1]
  }

  colnames(mat) <- c("NAME", "Description")
 }
 tpmmatrix <- as.data.frame(cbind(mat, tpmmatrix), stringsAsFactors = FALSE)
 write.table("#1.2", paste0(outfile, ".TPM.gct"), row.names = FALSE, col.names = FALSE,
  quote = FALSE)
 write.table(t(as.data.frame(dim(txi$abundance))), paste0(outfile, ".TPM.gct"),
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
 suppressWarnings(write.table(tpmmatrix, paste0(outfile, ".TPM.gct"), sep = "\t",
  row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE))

}
