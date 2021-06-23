# ideally, this should be done before creating the Seurat object,
# but because we used already processed data (not raw fastqs), in most cases this wasn't possible
# this function filters dataset genes by removing genes expressed in few cells
# and replacing non accepted genes with their accepted version

library(Seurat)
library(tidyverse)

process_gene_names <- function(data, n_cells = NULL, remove_ensembl = TRUE) {

  # 1 - remove genes expressed in less than n cells
  # (default assay should be defined beforehand)
  # not recommended in most cases, as it can remove too many genes

  if (is.numeric(n_cells)) {
    genes_keep <- c()

    for (gene in rownames(data)) {
      # print(gene)
      if (sum(GetAssayData(object = data, slot = "counts")[gene, ] > n_cells)) {
        genes_keep <- append(genes_keep, gene)
      }
    }

    # keep those genes from the seurat object

    data <- subset(data, features = genes_keep)
  }


  # 2- convert ENSEMBL gene names to HGNC SYMBOL

  # check if there are ENSMBL genes, otherwise omit

  if (grepl(rownames(data), pattern = "^ENSG")) {
    require(org.Hs.eg.db)

    symbols_match <- mapIds(org.Hs.eg.db, keys = rownames(data), keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")

    symbols <- c()
    for (i in 1:length(symbols_match)) {
      correct_symbol <- if_else(is.na(symbols_match[i]), names(symbols_match)[i], symbols_match[i])
      symbols <- append(symbols, correct_symbol)
    }
  } else {
    symbols <- rownames(data)
  }


  # 3 - make sure the naming is homogeneous (synonym symbols, etc.)

  library(HGNChelper)

  symbols <- checkGeneSymbols(symbols,
    unmapped.as.na = FALSE, # if there is no mapping, return original gene name
    # map = getCurrentHumanMap(), # downloads latest version of aliases (needs internet connection)
    species = "human"
  )

  print(length(symbols$Suggested.Symbol))
  print(length(unique(symbols$Suggested.Symbol)))

  symbols <- make.unique(symbols$Suggested.Symbol) # keep only suggested symbols for renaming


  # rename

  if (nrow(data) == length(symbols)) {
    if (!all(is.na(data[[data@active.assay]]@counts))) {
      rownames(data[[data@active.assay]]@counts) <- symbols
    }
    if (!all(is.na(data[[data@active.assay]]@data))) {
      rownames(data[[data@active.assay]]@data) <- symbols
    }
  } else {
    "Different number of genes. Cannot rename"
  }


  # if after this step we still have (ENSG...) genes they probably are novel transcripts or other not-so-relevant genes, thus could be deleted
  if (remove_ensembl == TRUE) {
    print(paste(
      grep(rownames(data), pattern = "^ENSG", value = TRUE),
      "will be removed",
      "\\n",
      sep = " "
    ))
    data <- subset(data, features = grep(rownames(data), pattern = "^ENSG", value = TRUE, invert = TRUE))
  }

  # finally make sure the meta.features are correctly named
  data[["RNA"]]@meta.features <- data.frame(row.names = rownames(data[["RNA"]]))
  return(data)
}
