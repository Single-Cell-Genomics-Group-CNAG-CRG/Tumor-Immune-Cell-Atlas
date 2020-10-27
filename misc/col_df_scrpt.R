# library(RColorBrewer)
library(pals)
library(dplyr)
# n <- 60
# qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
# Load color palette
col_pal_v <- readRDS(file = "data/immune_reference/cell_type_palette.rds")

# Save as df
col_df <- data.frame(ct_col = col_pal_v, stringsAsFactors = FALSE) %>%
  tibble::rownames_to_column("plt_name") %>%
  dplyr::mutate(
    ct_name =  gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                    x = plt_name, 
                    perl = TRUE),
    plt_name = dplyr::if_else(plt_name == "Central memory CD4 T cells", 
                               "Transitional Memory CD4 T cells",
                               plt_name),
    ct_name = dplyr::case_when(
      plt_name == "Recently activated CD4 T cells" ~ "Effector.precursor.CD4.T.cells",
      plt_name == "B cells" ~ "B.cell",
      plt_name == "Proliferative B cells" ~ "Proliferative.B.cell",
      plt_name == "Naive-memory CD4 T cells" ~ "Effector.memory.CD4.T.cells",
      plt_name == "Plasma B cells" ~ "Plasma.B",
      plt_name == "Proinflamatory TAMs" ~ "Proinflamatory.TAMs.and.neutrophils",
      # plt_name == "Central memory CD4 T cells" ~ "Transitional.Memory.CD4.T.cells",
      plt_name == "Mast cells" ~ "MAST",
      TRUE ~ ct_name)
    ) %>%
  dplyr::filter(plt_name != "Unknown")

# col_vector  <-  as.vector(pals::polychrome())
paula_order_old <- c("B cell", "Proliferative B cell", "Plasma B", "Naive T cells", "Regulatory T cells", "T helper cells", "Th17 cells", "Proliferative T cells", "Effector precursor CD4 T cells", "Effector memory CD4 T cells", "Central memory CD4 T cells", "Pre-exhausted CD8 T cells", "Cytotoxic CD8 T cells", "Effector memory CD8 T cells", "Terminally exhausted CD8 T cells", "NK", "SPP1 TAMs", "M2 TAMs", "Proinflamatory TAMs and neutrophils", "Proliferative monocytes and macrophages", "Monocytes", "cDC", "pDC", "mDC", "MAST")
paula_order <- c("B cells", "Proliferative B cells", "Plasma B cells", "Naive T cells", "Regulatory T cells", "T helper cells", "Th17 cells", "Proliferative T cells", "Recently activated CD4 T cells", "Naive-memory CD4 T cells", "Transitional Memory CD4 T cells", "Pre-exhausted CD8 T cells", "Cytotoxic CD8 T cells", "Effector memory CD8 T cells", "Terminally exhausted CD8 T cells", "NK", "SPP1 TAMs", "M2 TAMs", "Proinflamatory TAMs", "Proliferative monocytes and macrophages", "Monocytes", "cDC", "pDC", "mDC", "Mast cells")

# col_df <- data.frame(ct_name = gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
#                                     x = paula_order_old,
#                                     perl = TRUE),
#                      plt_name = paula_order_old,
#                      ct_col = col_vector[1:length(paula_order_old)],
#                      stringsAsFactors = FALSE)
# 
# # Change label names
# col_df <- col_df %>%
#   dplyr::mutate(plt_name = dplyr::if_else(plt_name == "Effector precursor CD4 T cells",
#                                           "Recently activated CD4 T cells", plt_name))
