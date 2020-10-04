# This script creates the test sets for the 5-fold cross-validation.
# We ensure that each cell type is represented in each test set. Thus, we take
# 20% of cells for each cell type in each fold. Each cell in the dataset is found
# in 1 test set and 4 training sets.


# Load packages
library(tidyverse)


# Load data
tica_df_signatures <- readRDS("tmp/tica_input_dataframe_signatures_random_forest.rds")


# Create test sets in each fold
cells <- tica_df_signatures$barcode
cell_types <- levels(tica_df_signatures$new_cell_types)
test_cells_list <- list(
  fold1 = c(),
  fold2 = c(),
  fold3 = c(),
  fold4 = c(),
  fold5 = c()
)
for (cell_type in cell_types) {
  print(cell_type)
  cells <- tica_df_signatures[tica_df_signatures$new_cell_types == cell_type, "barcode"]
  steps_cross_valid <- round(
    seq(
      from = 1,
      to = length(cells),
      by = (length(cells) / 5)
    ),
    0
  )
  test_cells_list[["fold1"]] <- c(
    test_cells_list[["fold1"]],
    cells[steps_cross_valid[1]:steps_cross_valid[2]]
  )
  test_cells_list[["fold2"]] <- c(
    test_cells_list[["fold2"]],
    cells[(steps_cross_valid[2] + 1):steps_cross_valid[3]]
  )
  test_cells_list[["fold3"]] <- c(
    test_cells_list[["fold3"]],
    cells[(steps_cross_valid[3] + 1):steps_cross_valid[4]]
  )
  test_cells_list[["fold4"]] <- c(
    test_cells_list[["fold4"]],
    cells[(steps_cross_valid[4] + 1):steps_cross_valid[5]]
  )
  test_cells_list[["fold5"]] <- c(
    test_cells_list[["fold5"]],
    cells[(steps_cross_valid[5] + 1):length(cells)]
  )
}
saveRDS(test_cells_list, "tmp/test_set_per_fold.rds")

