# This script trains and tests a random forest classifier on a specific fold
# of the cross validation specified in previous scripts. It saves both the trained model and
# the predicted labels.
start_time <- Sys.time()


# Load packages
library(Seurat)
library(tidyverse)
library(randomForest)
library(doParallel)
library(foreach)
library(readxl)


# Load data
tica_df <- readRDS("tmp/tica_input_dataframe_signatures_random_forest.rds")
test_cells_list <- readRDS("tmp/test_set_per_fold.rds")
print("Data loaded successfully!")


# Separate datframe in real and random signatures
signatures_cols <- c(
  str_subset(colnames(tica_df), "signature"),
  "new_cell_types"
)
random_cols <- c(
  str_subset(colnames(tica_df), "random"),
  "new_cell_types"
)
tica_df_signatures <- tica_df[, signatures_cols]
tica_df_random <- tica_df[, random_cols]
print("Separated real and random signatures")

  
# Process command-line arguments and subset to correct fold
args <- commandArgs(trailingOnly = TRUE)
fold <- args[[1]]
test_cells <- test_cells_list[[fold]]


# Divide cells in training and test sets
test_df_signatures <- tica_df_signatures[test_cells, ]
test_df_random <- tica_df_random[test_cells, ]
training_cells <- rownames(tica_df)[!(rownames(tica_df) %in% test_cells)]
training_df_signatures <- tica_df_signatures[training_cells, ]
training_df_random <- tica_df_random[training_cells, ]
print("Divided training and test sets")


# Train random forest
# ## Set parallelization backend
total_trees <- 600
parallel_trees <- 100
ncpu <- round(total_trees / parallel_trees, 0)
doParallel::registerDoParallel(ncpu)


########################################################################
############################## SIGNATURES ##############################
########################################################################

# Run n (ncpu) forest of ntrees (parallel_trees) each: signatures
rf_mod_signatures <- foreach::foreach(ntree = rep(parallel_trees, ncpu),
                           # Indicating we want to use combine from the randomForest to join all the parallelized trees
                           .combine = randomForest::combine,
                           # Indicate we need the randomForest within the ForeEach loop + dopar indicating parallelization
                           .packages="randomForest") %dopar% {
                             ## train RF model
                             model <- randomForest::randomForest(
                               new_cell_types ~ .,
                               data = training_df_signatures,
                               type = "classification",
                               # Keep importance in the RF object
                               importance = FALSE,
                               # How many trees to build
                               ntree = ntree,
                               # Return update every 10 trees
                               verbose = 10
                             )
                             return(model)
                           }
print("Finshed random forest signatures!")


# Test random forest
predicted_probs_signatures <- predict(
  object = rf_mod_signatures,
  newdata = test_df_signatures,
  type = "prob"
)
predicted_class_signatures <- predict(
  object = rf_mod_signatures,
  newdata = test_df_signatures,
  type = "class"
)


# Save model and predictions
saveRDS(rf_mod_signatures, str_c("tmp/random_forest_model_signatures", "_", fold, ".rds", sep = ""))
saveRDS(predicted_probs_signatures, str_c("tmp/random_forest_predicted_probs_signatures", "_", fold, ".rds"))
saveRDS(predicted_class_signatures, str_c("tmp/random_forest_predicted_class_signatures", "_", fold, ".rds"))
print("Saved results signatures!")


########################################################################
############################## RANDOM ##################################
########################################################################

# Run n (ncpu) forest of ntrees (parallel_trees) each: signatures
rf_mod_random <- foreach::foreach(ntree = rep(parallel_trees, ncpu),
                                      # Indicating we want to use combine from the randomForest to join all the parallelized trees
                                      .combine = randomForest::combine,
                                      # Indicate we need the randomForest within the ForeEach loop + dopar indicating parallelization
                                      .packages="randomForest") %dopar% {
                                        ## train RF model
                                        model <- randomForest::randomForest(
                                          new_cell_types ~ .,
                                          data = training_df_random,
                                          type = "classification",
                                          # Keep importance in the RF object
                                          importance = FALSE,
                                          # How many trees to build
                                          ntree = ntree,
                                          # Return update every 10 trees
                                          verbose = 10
                                        )
                                        return(model)
                                      }
print("Finshed random forest with random signatures!")


# Test random forest
predicted_probs_random <- predict(
  object = rf_mod_random,
  newdata = test_df_random,
  type = "prob"
)
predicted_class_random <- predict(
  object = rf_mod_random,
  newdata = test_df_random,
  type = "class"
)


# Save model and predictions
saveRDS(rf_mod_random, str_c("tmp/random_forest_model_random", "_", fold, ".rds", sep = ""))
saveRDS(predicted_probs_random, str_c("tmp/random_forest_predicted_probs_random", "_", fold, ".rds"))
saveRDS(predicted_class_random, str_c("tmp/random_forest_predicted_class_random", "_", fold, ".rds"))

print("Saved results random signatures!")

end_time <- Sys.time()
total_time <- end_time - start_time
print(str_c("Job took ", total_time, " to compute", sep = ""))
