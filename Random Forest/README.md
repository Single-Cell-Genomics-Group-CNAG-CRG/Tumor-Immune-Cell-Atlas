# Validation of TICA clusters using a random forest classifier


## Dependencies

* [R 3.6.0](https://cran.r-project.org/)
* [Seurat 3.2.0](https://cran.r-project.org/web/packages/Seurat/index.html)
* [readxl 1.3.1](https://cran.r-project.org/web/packages/readxl/index.html)
* [tidyverse 1.3.0](https://cran.r-project.org/web/packages/tidyverse/index.html)
* [randomForest 4.6.14](https://cran.r-project.org/web/packages/randomForest/index.html)
* [doParallel 1.0.15](https://cran.r-project.org/web/packages/doParallel/index.html)
* [foreach 1.5.0](https://cran.r-project.org/web/packages/foreach/index.html)
* [caret 6.0.86](https://cran.r-project.org/web/packages/caret/index.html)
* [e1071 1.7.3](https://cran.r-project.org/web/packages/e1071/index.html)
* [ggpubr 0.3.0](https://cran.r-project.org/web/packages/ggpubr/index.html)
* [pheatmap 1.0.12](https://cran.r-project.org/web/packages/pheatmap/index.html)


## Data

* TICA Seurat object: download it as specified in the publication
* Supplementary table 2: excel file with the markers for each cell type. For simplicity, we provide it in the ./data/ folder.


## Pipeline

The scripts should be run in the following order:


1. 1-calculate_cell_type_signatures.R: calculates 25 cell type-specific signatures + 25 random signatures for each cell in the dataset.
  * Input:
    * data/supp_table_2.xlsx: cell type specific markers (Supplementary Table 2 in the paper)
    * TICA Seurat object (path should be changed before running the script)
  * Output:
    * tmp/tica_input_dataframe_signatures_random_forest.rds: data frame with cells as observations, signatures as feautures and cell type as response variable.
2. 2-create_test_sets_for_cross_validation.R
  * Input: tmp/tica_input_dataframe_signatures_random_forest.rds
  * Output: tmp/test_set_per_fold.rds: list of 5 character vectors containing the barcodes of the cells that belong to each fold of the cross-validation.
3. 3-run_random_forest.R: train and test 2 RF (cell type-specific + random signatures).
  * Arguments: this script should have one argument specifying the fold to run the RF (ie "fold1"). The purpose of this is to parallelize the execution of the RF across folds.
  * Input:
    * tmp/tica_input_dataframe_signatures_random_forest.rds
    * tmp/test_set_per_fold.rds
  * Output: 
    * Cell type-specific:
      * tmp/random_forest_model_signatures_fold{1-5}.rds (RF models)
      * tmp/random_forest_predicted_probs_signatures_fold{1-5}.rds (RF probabilities)
      * tmp/random_forest_predicted_class_signatures_fold{1-5}.rds (RF predicted classes)
    * Random signatures:
      * tmp/random_forest_model_random_fold{1-5}.rds
      * tmp/random_forest_predicted_probs_random_fold{1-5}.rds
      * tmp/random_forest_predicted_class_random_fold{1-5}.rds
4. 4-evaluate_performance.R
  * Input: all outputs from previous script
  * Output:
  	* tmp/confusion_matrix_probabilities_list.rds: list with 5 confusion matrices (one per fold), where each CF is the probability that a cell from a cell type X is assigned to cell type Y.
    * tmp/accuracy_random_forest_dataframe.rds: accuracy and kappa statistics for each fold.
5. 5-plot_figures.R: 
  * Input:
    * tmp/accuracy_random_forest_dataframe.rds
    * tmp/confusion_matrix_probabilities_list.rds
  * Output: 
    * results/accuracy_kappa_random_forest.pdf (panels a and b)
    * results/confusion_matrix_probabilities.pdf (panel c)