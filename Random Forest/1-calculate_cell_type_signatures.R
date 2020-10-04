# This script computes 25 cell type-specific signatures for each cell, based on the
# markers provided in the Supplementary Table 2. It also calculates signatures
# defined with random genes that will help us assess the significance of the classifier.
# Finally, it saves the dataframe that will be the input for the random forest.


# Load packages
library(Seurat)
library(readxl)
library(tidyverse)


# Load data
tica <- readRDS("path_to_TICA_seurat_object")
path_to_excel <- "data/supp_table_2.xlsx"
sheets <- excel_sheets(path_to_excel)
signatures_dfs <- purrr::map(sheets, ~ read_excel(path_to_excel, sheet = .x))
names(signatures_dfs) <- str_replace_all(sheets, " ", "_")
signatures_tica <- purrr::map(signatures_dfs, "gene")


# Create dataframe to run the random forest (signatures)
DefaultAssay(tica) <- "integrated"
tica@assays$RNA <- NULL
set.seed(123)
for (signature in names(signatures_tica)) {
  print(signature)
  start_time <- Sys.time()
  features <- signatures_tica[[signature]]
  features <- features[features %in% rownames(tica)]
  tica <- AddModuleScore(
    tica,
    features = list(features),
    name = str_c(signature, "signature", sep = "_")
  )
  random_features <- sample(
    rownames(tica),
    size = length(features),
    replace = FALSE
  )
  tica <- AddModuleScore(
    tica,
    features = list(random_features),
    name = str_c(signature, "random", sep = "_")
  )
  end_time <- Sys.time()
  print(end_time - start_time)
  print(head(tica@meta.data))
}
keep_cols_sign <- c(
  str_subset(colnames(tica@meta.data), "signature"),
  str_subset(colnames(tica@meta.data), "random"),
  "new_cell_types"
)
tica_df_signatures <- tica@meta.data[, keep_cols_sign]
tica_df_signatures$barcode <- colnames(tica)


# Save
dir.create("tmp", showWarnings = FALSE)
saveRDS(tica_df_signatures, "tmp/tica_input_dataframe_signatures_random_forest.rds")
