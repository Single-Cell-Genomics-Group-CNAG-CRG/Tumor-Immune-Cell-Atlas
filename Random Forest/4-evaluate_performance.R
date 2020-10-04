# This script evaluates the performance of the random forest classifiers for all folds


# Load packages
library(caret)
library(e1071)
library(tidyverse)


# Load data
test_sets_cross_valid <- readRDS("tmp/test_set_per_fold.rds")
tica_df <- readRDS("tmp/tica_input_dataframe_signatures_random_forest.rds")
folds <- str_c("fold", 1:5, sep = "")
pred_class_sign_l <- purrr::map(folds, function(fold) {
  path_to_preds1 <- "tmp/random_forest_predicted_class_signatures_"
  path_to_preds <- str_c(path_to_preds1, fold, ".rds", sep = "")
  pred_class <- readRDS(path_to_preds)
  pred_class
})
names(pred_class_sign_l) <- folds
pred_class_rand_l <- purrr::map(folds, function(fold) {
  path_to_preds1 <- "tmp/random_forest_predicted_class_random_"
  path_to_preds <- str_c(path_to_preds1, fold, ".rds", sep = "")
  pred_class <- readRDS(path_to_preds)
  pred_class
})
names(pred_class_rand_l) <- folds

pred_probs_sign_l <- purrr::map(folds, function(fold) {
  path_to_probs1 <- "tmp/random_forest_predicted_probs_signatures_"
  path_to_probs <- str_c(path_to_probs1, fold, ".rds", sep = "")
  pred_probs <- readRDS(path_to_probs)
  pred_probs
})
names(pred_probs_sign_l) <- folds
pred_probs_rand_l <- purrr::map(folds, function(fold) {
  path_to_probs1 <- "tmp/random_forest_predicted_probs_random_"
  path_to_probs <- str_c(path_to_probs1, fold, ".rds", sep = "")
  pred_probs <- readRDS(path_to_probs)
  pred_probs
})
names(pred_probs_rand_l) <- folds


# Measure accuracy for each fold
accuracy_dfs <- purrr::map(folds, function(fold) {
  test_df <- tica_df[test_sets_cross_valid[[fold]], ]
  comparison_rownames_sign <- all(rownames(test_df) == names(pred_class_sign_l[[fold]]))
  comparison_rownames_rand <- all(rownames(test_df) == names(pred_class_rand_l[[fold]]))
  if (comparison_rownames_sign & comparison_rownames_rand) {
    conf_mat_sign <- caret::confusionMatrix(
      test_df$new_cell_types,
      pred_class_sign_l[[fold]]
    )
    conf_mat_rand <- caret::confusionMatrix(
      test_df$new_cell_types,
      pred_class_rand_l[[fold]]
    )
    df <- data.frame(
      accuracy = c(conf_mat_sign$overall["Accuracy"], conf_mat_rand$overall["Accuracy"]),
      kappa = c(conf_mat_sign$overall["Kappa"], conf_mat_rand$overall["Kappa"]),
      type = c("signatures", "random")
    )
    df

  } else {
    "Rownames do not match"
  }
})
names(accuracy_dfs) <- folds
accuracy_df <- bind_rows(accuracy_dfs, .id = "fold")
saveRDS(accuracy_df, "tmp/accuracy_random_forest_dataframe.rds")


# Calculate confusion matrix of probabilities for each fold
conf_mat_probs_l <- purrr::map(folds, function(fold) {
  test_df <- tica_df[test_sets_cross_valid[[fold]], ]
  if (all(rownames(test_df) == names(pred_probs_sign_l[[fold]]))) {
    cell_types <- levels(test_df$new_cell_types)
    average_probs_l <- purrr::map(cell_types, function(cell_type) {
      indices_cell_type <- rownames(test_df[test_df$new_cell_types == cell_type, ])
      mat <- pred_probs_sign_l[[fold]][indices_cell_type, ]
      average_prob <- colMeans(mat)
      average_prob
    })
    names(average_probs_l) <- cell_types
    average_probs_mat <- average_probs_l %>%
      bind_rows() %>%
      as.matrix()
    rownames(average_probs_mat) <- names(average_probs_l)
    average_probs_mat
  } else {
    "Rownames do not match"
  }
})
names(conf_mat_probs_l) <- folds
saveRDS(conf_mat_probs_l, "tmp/confusion_matrix_probabilities_list.rds")

