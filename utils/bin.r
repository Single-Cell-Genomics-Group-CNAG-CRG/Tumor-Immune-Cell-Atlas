
# This function preprocesses the seurat object so that it can be used in QC_plots_fun and QC_UMAP_fun
se_preprocess_QC <- function(se_obj, vrs_names, mt_thresh, gene_thresh, count_thresh){
  
  # Relabel variables so that we can access them later
  se_obj@meta.data$nCount_RNA <- se_obj@meta.data[,vrs_names[1]]
  # Number of genes
  se_obj@meta.data$nFeature_RNA <- se_obj@meta.data[,vrs_names[2]]
  # Mito percentage
  se_obj@meta.data$percent.mt <- se_obj@meta.data[,vrs_names[3]]
  # Color cells by
  if(length(vrs_names) < 4){
    se_obj@meta.data$color_feat <- 'group1'
  } else {
    se_obj@meta.data$color_feat <- se_obj@meta.data[,vrs_names[4]]
  }
  
  # Label cells by good-poor quality
  se_obj$quality <- if_else(se_obj$percent.mt < mt_thresh &
                              se_obj$nFeature_RNA > gene_thresh &
                              se_obj$nCount_RNA > count_thresh, 'good', 'poor')
  
  return(se_obj)
}

####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

QC_seurat_hist <- function(se, assay, slot, nfeat, ncount, pctmt, pctrp) {
  # se: Seurat object
  # assay: assay to where to get the count matrix
  # slot: slot to where to get the count matrix
  require(Seurat)
  require(ggplot2)
  require(cowplot)
  
  p1 <- ggplot() +
    geom_histogram(data = se[[]],
                   aes_string(x = nfeat),
                   fill = "red",
                   alpha = 0.7,
                   color = "gray",
                   bins = 50) +
    ggplot2::theme_classic() +
    ggtitle("Genes per spot")
  
  p2 <- ggplot() +
    geom_histogram(data = se[[]],
                   aes_string(x = ncount),
                   fill = "red",
                   alpha = 0.7,
                   color = "gray",
                   bins = 50) +
    ggplot2::theme_classic() +
    ggtitle("Total counts per spots")
  
  p3 <- ggplot() +
    geom_histogram(data = se[[]],
                   aes_string(x = pctmt),
                   fill = "red",
                   alpha = 0.7,
                   color = "gray",
                   bins = 50) +
    ggplot2::theme_classic() +
    ggtitle("Mitochondrial % per spot")
  
  p4 <- ggplot() +
    geom_histogram(data = se[[]],
                   aes_string(x = pctrp),
                   fill = "red",
                   alpha = 0.7,
                   color = "gray",
                   bins = 50) +
    ggplot2::theme_classic() +
    ggtitle("Ribosomal % per spot")
  
  # Extract count_mtrx from seurat object
  count_mtrx <- GetAssayData(se,
                             slot = slot,
                             assay = assay)
  
  gene_attr <- data.frame(nUMI = Matrix::rowSums(count_mtrx),
                          nSpots = Matrix::rowSums(count_mtrx > 0))
  p5 <- ggplot() +
    geom_histogram(data = gene_attr,
                   aes(x = nUMI),
                   fill = "red",
                   alpha = 0.7,
                   color = "gray",
                   bins = 50) +
    ggplot2::theme_classic() +
    scale_x_log10() +
    ggtitle("Total counts per gene (log10 scale)")
  
  p6 <- ggplot() +
    geom_histogram(data = gene_attr,
                   aes(x = nSpots),
                   fill = "red",
                   alpha = 0.7,
                   color = "gray",
                   bins = 50) +
    ggplot2::theme_classic() +
    ggtitle("Total spots per gene")
  
  return(cowplot::plot_grid(p1, p2, p3, p4, p5, p6))
}

####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

# QC plots from Seurat object
QC_plots_fun <- function(se_obj, count_thresh, gene_thresh, mt_thresh, vrs_names){
  ###
  # This functions grabs a seurat object and returns basic exploratory plots to check the quality of the data
  # 
  # Parameters
  # se_object: seurat object 
  # count_thresh: threshold to filter out those cells with a smaller library size - count_thresh <- log10(750)
  # gene_thresh: minimum number of genes a cell needs to have detected - gene_thresh <- 200
  # mt_thresh: threshold to filter out to filter out those cells with a high MTpct - mt_thresh <- 25
  # vrs_names: pass a vector with the variable names in this order c(library_size, n_genes, MT_pct, color_feat)
  # 
  # Return
  # it returns a grid arranged plot with all the subplots
  ###
  
  se_obj <- se_preprocess_QC(se_obj = se_obj, vrs_names = vrs_names, mt_thresh, gene_thresh, count_thresh)
  
  ### Setting variable names ###
  # Library size
  # se_obj@meta.data$nCount_RNA <- se_obj@meta.data[,vrs_names[1]]
  # # Number of genes
  # se_obj@meta.data$nFeature_RNA <- se_obj@meta.data[,vrs_names[2]]
  # # Mito percentage
  # se_obj@meta.data$percent.mt <- se_obj@meta.data[,vrs_names[3]]
  # # Color cells by
  # if(length(vrs_names) < 4){
  #   se_obj@meta.data$color_feat <- 'group1'
  # } else {
  #   se_obj@meta.data$color_feat <- se_obj@meta.data[,vrs_names[4]]
  # }
  
  ####################
  ### Library size ###
  ####################
  ### We first filter out cells that have a library size (total number of RNA molecules) too small in comparison with other cells. Such cells are likely to have broken or failed to capture. To determine the threshold, we can visualize the library size distribution with a histogram. As there are outliers with a great deal of counts, we will plot the log distribution:
  
  x_titl <- expression("log"[10]*"(library size)")
  
  lib_size_log_qc <- data.frame(se_obj@meta.data) %>% 
    dplyr::mutate(exclude = ifelse(log10(nCount_RNA) < count_thresh, TRUE, FALSE)) %>%
    ggplot(aes(log10(nCount_RNA), fill = exclude, color = exclude)) + 
    geom_histogram(bins = 100, alpha = 0.65) +
    # geom_vline(xintercept = 2.85, color = "red", linetype = "dashed") +
    geom_vline(xintercept = count_thresh, color = "red", linetype = "dashed") +
    scale_x_continuous(x_titl) +
    scale_y_continuous(expand = c(0,0)) +
    scale_color_manual(values = c("black", "red2")) + 
    scale_fill_manual(values = c("black", "red2")) +
    theme_bw() +
    labs(y = 'Log transformed library size') +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = 'none')
  
  #####################
  ### Cell coverage ###
  #####################
  # We next filter by the cell coverage, which is the number of detected genes in each cell (i.e., number of genes with non-zero counts for a given cell). We want to ensure that the reads are distributed across the transcriptome. Thus, we rule out those cells that have an abnormally low number of detected genes.
  
  cell_coverage_hist <- data.frame(se_obj@meta.data) %>%
    mutate(exclude = ifelse(nFeature_RNA < gene_thresh, TRUE, FALSE)) %>%
    ggplot(aes(x = nFeature_RNA, fill = exclude, color = exclude)) + 
    geom_histogram(bins = 100, alpha = 0.65) +
    geom_vline(xintercept = gene_thresh, color = "red", linetype = "dashed") +
    scale_x_continuous("Number of detected genes") +
    scale_y_continuous(expand = c(0,0)) +
    scale_color_manual(values = c("black", "red2")) + 
    scale_fill_manual(values = c("black", "red2")) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = 'none')
  
  #######################
  ### Mitochondrial % ###
  #######################
  # We next filter by mitochondrial percentage  It is expected that poor-quality cells are enriched for the expression of mitochondrial genes, likely because cells underwent apoptosis and non-mitochondrial mRNA was released and lost
  
  mt_genes_qc <- data.frame(se_obj@meta.data) %>% 
    # Exclude samples with >=20% of counts mapped to mitochondrial DNA
    dplyr::mutate(exclude = ifelse(percent.mt >= mt_thresh, TRUE, FALSE)) %>%
    ggplot(aes(x = percent.mt, fill = exclude, color = exclude)) +
    geom_histogram(bins = 200, alpha = 0.65) +
    geom_vline(xintercept = mt_thresh, linetype = "dashed", color = "red") +
    scale_x_continuous("Mitochondrial proportion (%)") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_color_manual(values = c("black", "red2")) + 
    scale_fill_manual(values = c("black", "red2")) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = 'none')
  
  #############################################
  ### % Mitocondrial RNA vs number of genes ###
  #############################################
  # We wanna see if those cells with a high MT % have a low number of genes, it is also an indicator of low quality
  
  Mt.pct_vs_ngenes <- ggplot(se_obj@meta.data) +
    geom_point(aes(x = nFeature_RNA, y = percent.mt, color = quality), alpha = 0.65) +
    theme_classic() +
    labs(
      x = 'Number of detected genes',
      y = 'Mitochondrial proportion (%)') +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      legend.position = 'none'
    ) +
    scale_color_manual(values = c("black", "red2"))
  
  # top_row <- ggpubr::ggarrange(plotlist = list(lib_size_log_qc, cell_coverage_hist, mt_genes_qc), ncol = 3, align = 'hv')
  # bot_row <- ggpubr::ggarrange(plotlist = list(Mt.pct_vs_ngenes, UMAP_coord), ncol = 2, align = 'hv')
  # joint_plts <- ggpubr::ggarrange(plotlist = list(top_row, bot_row),
  #                   nrow = 2, align = 'hv')
  
  return(list(lib_size_log_qc, cell_coverage_hist, mt_genes_qc, Mt.pct_vs_ngenes))
  
}

####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

QC_UMAP_fun <- function(se_obj, vrs_names, mt_thresh, gene_thresh, count_thresh){
  # This function takes in a seurat objects and a list of variable names to carry out a UMAP labelling dots by good/poor quality.
  
  se_obj <- se_preprocess_QC(se_obj = se_obj, vrs_names = vrs_names, mt_thresh, gene_thresh, count_thresh)
  
  #########################################
  ### Clustering and color by good-poor ###
  #########################################
  
  # if(! "SCT" %in% names(se_obj@assays)) {
  #   # 1st scale and normalize the data
  #   se_obj <- SCTransform(se_obj, verbose = T)  
  # }
  
  if(! "pca" %in% names(se_obj@reductions)) {
    # 2nd perform PCA
    se_obj <- RunPCA(se_obj,
                     features = VariableFeatures(object = se_obj))
    
    # 3rd cluster the cells, we choose 1:20 from the PCs
    se_obj <- FindNeighbors(se_obj, dims = 1:40)
    se_obj <- FindClusters(se_obj, resolution = 0.6)
  }
  
  if(! "umap" %in% names(se_obj@reductions)) {
    # 4th Run and plot UMAP
    se_obj <- RunUMAP(se_obj, dims = 1:40)
  }
  
  UMAP_quality <- DimPlot(se_obj, reduction = "umap", group.by = 'quality')
  UMAP_quality2 <- UMAP_quality +
    scale_color_manual(values = c('black', 'red')) +
    geom_point(aes(shape = 'maturation_stage'))
  
  UMAP_coord <- data.frame(se_obj@reductions$umap@cell.embeddings) %>%
    rownames_to_column('id') %>%
    left_join(data.frame(se_obj@meta.data) %>% rownames_to_column('id'), by = 'id') %>%
    # mutate(maturation_stage = if_else(is.na(maturation_stage), 'NA', maturation_stage)) %>% 
    ggplot() +
    geom_point(aes(x=UMAP_1, y=UMAP_2, colour = quality), alpha = 0.65) +
    theme_classic() +
    scale_color_manual(values = c('black', 'red'))
  # UMAP_coord
  
  return(UMAP_coord)
  
}

####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

# test_spot_fun <- function(se_obj, n=10000, clust_vr, verbose=TRUE){
#   ######################################
#   # This functions takes in a seurat object and creates different mixtures resembling spots at different proportions
#   #
#   # Args
#   # se_obj: seurat object; The seurat object NEEDS to have the variable features calculated
#   # n: number of spots to generate
#   # clust_vr: name of the variables with the cluster assignment
#   # Returns:
#   # This function returnsa dataframe where each column is the proportion of each spot and
#   # a sparse matrix with the synthetic spots with the genenames as rownames and colnames being the spot names
#   #
#   ######################################
#   
#   library(DropletUtils) # For the downsampling
#   library(dtplyr) # To use dplyr commands with DT speed
#   
#   se_obj$seurat_clusters <- droplevels(se_obj@meta.data[,clust_vr])
#   
#   print('Generating synthetic test spots...')
#   start_gen <- Sys.time()
#   # create progress bar
#   pb <- txtProgressBar(min = 0, max = n, style = 3)
#   
#   # Defining a vector with the highly variable genes
#   # hvg <- VariableFeatures(se_obj)
#   
#   # Save count matrix
#   count_mtrx <- as.matrix(se_obj@assays$RNA@counts)
#   # Remove genes not expressed in at least 0.5% of the cells
#   min_cells <- 0.005*ncol(se_obj)
#   count_mtrx <- count_mtrx[rowSums(count_mtrx != 0) > min_cells,]
#   
#   
#   # Iterator to define the name of the spots generated
#   # id <- 1
#   
#   ds_spots <- lapply(1:n, function(i){
#     
#     # Select between 2 and 10 cells randomly from the count matrix
#     cell_pool <- sample(colnames(count_mtrx), sample(x = 2:10, size = 1))
#     
#     # Determine the weight each cell will have on the synthetic spot
#     # weigh <-runif(length(cell_pool))
#     # weigh <- weigh/sum(weigh)
#     
#     # We're not going to sum the reads as is bc spots are **enriched** 
#     # so we'll add up the counts and downsample to the ~depth of a typical spot.
#     
#     # Create a name for the spot with the necessary info to deconvolute it
#     pos <- which(colnames(count_mtrx) %in% cell_pool)
#     tmp_ds <- se_obj@meta.data[pos,] %>% mutate(weight = 1)
#     # tmp_ds[,'weight'] <- weigh
#     name_simp <- paste('spot_',i,sep='')
#     # name_long <- paste(with(tmp_ds, paste(rownames(tmp_ds), 'cluster', seurat_clusters, 'weight', weight, sep = '_')), collapse = '-')
#     
#     spot_ds <- tmp_ds %>%
#       select(seurat_clusters, weight) %>%
#       mutate(seurat_clusters = paste('clust_',seurat_clusters,sep = '')) %>% 
#       group_by(seurat_clusters) %>% 
#       summarise(sum_weights = sum(weight)) %>% 
#       ungroup() %>%
#       pivot_wider(names_from = seurat_clusters, values_from = sum_weights) %>% 
#       mutate(name = name_simp)
#     
#     # Generate synthetic spot
#     
#     ## Here we multiply each vector by its weight
#     # weighted <- lapply(1:length(cell_pool), function(ii) expr <- as.integer(round(weigh[[ii]]*count_mtrx[hvg,cell_pool[[ii]]],0)))
#     
#     ## Next step we add all the vectors by position
#     # syn_spot <- Reduce(`+`,weighted)
#     # ret_ds <- data.frame(gene=hvg, tmp=syn_spot)
#     
#     ## Here we add up the counts of each cell
#     syn_spot <- rowSums(as.matrix(count_mtrx[,cell_pool])); sum(syn_spot)
#     names_genes <- names(syn_spot)
#     ## Downsample
#     ### 25k is a bit above average 20k UMIs observed in spatial transcriptomics data then downsample to 20k
#     if(sum(syn_spot) > 25000){
#       syn_spot_sparse <- downsampleMatrix(Matrix::Matrix(syn_spot,sparse = T), prop = 20000/sum(syn_spot))
#     } else {
#       syn_spot_sparse <- Matrix::Matrix(syn_spot,sparse = T)
#     }
#     
#     rownames(syn_spot_sparse) <- names_genes
#     colnames(syn_spot_sparse) <- name_simp
#     # ret_ds <- data.frame(syn_spot) %>% mutate(gene=names_hvg)
#     
#     # colnames(ret_ds) <- c('gene', i)
#     # id <<- id+1
#     
#     # update progress bar
#     setTxtProgressBar(pb, i)
#     # Return the transpose so that documents are on the rows and genes are on the columns
#     return(list(syn_spot_sparse,spot_ds))
#   })
#   
#   # Generate sparse matrix of spots
#   # ds_syn_spots <- map(ds_spots, 1) %>%   
#   # Reduce(function(...) merge(..., by='gene', all.x=TRUE), .) %>%
#   # column_to_rownames(var="gene") %>%
#   # as.matrix() %>%
#   # Matrix(sparse = TRUE) %>% 
#   # t()
#   
#   ds_syn_spots <- map(ds_spots, 1) %>%
#     base::Reduce(function(m1,m2) cbind(unlist(m1),unlist(m2)), .)
#   
#   # Generate dataframe of spot characteristic
#   ds_spots_metadata <- map(ds_spots, 2) %>%
#     bind_rows() %>% 
#     lazy_dt(.) %>%  # Convert to lazy data so it can be avaluated as a data.table
#     # mutate_all(~replace_na(.,0)) %>%
#     data.frame()
#   
#   ds_spots_metadata[is.na(ds_spots_metadata)] <- 0
#   
#   # change column order so that its progressive
#   lev_mod <- gsub("[\\+|\\ ]", ".", levels(se_obj$seurat_clusters))
#   all_cn <- c(paste('clust_',lev_mod,sep = ''),'name')
#   if( sum(all_cn %in% colnames(ds_spots_metadata)) == (nlevels(se_obj$seurat_clusters)+1) ){
#     ds_spots_metadata <- ds_spots_metadata[,all_cn]
#   } else {
#     
#     # stringr::str_replace(levels(se_obj$seurat_clusters), pattern = '[\\+]', replacement = '.')
#     # lev_mod <- stringr::str_replace(levels(se_obj$seurat_clusters), pattern = '[\\+]', replacement = '.')
#     # lev_mod <- gsub(" ", ".", lev_mod, fixed = TRUE)
#     # lev_mod <- stringr::str_replace(lev_mod, pattern = '\\s', replacement = '.')
#     # lev_mod <- paste('clust',lev_mod,sep = '_')
#     
#     missing_cols <- all_cn[ which( ! all_cn %in% colnames(ds_spots_metadata) ) ]
#     ds_spots_metadata[missing_cols] <- 0
#     ds_spots_metadata <- ds_spots_metadata[,all_cn]
#   }
#   
#   # Close progress bar
#   close(pb)
#   
#   print(sprintf('Generation of %s test spots took %s mins', n, difftime(Sys.time(), start_gen, units = 'mins')))
#   print('output consists of a list with two dataframes, this first one has the weighted count matrix and the second has the metadata for each spot')
#   return(list(ds_syn_spots,ds_spots_metadata))
# }

####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

lda_performance_fun <- function(){
  results <- posterior(object=train_mod, newdata=test_dtm)
  
  train_mod@gamma[1:6,]
  train_mod@beta[1:6,1:6]
  train_mod@wordassignments$i
  results$topics[1,]
  
  syn_spts_metadata$parent$parent[1,]
  is()
} 

#############################################################################################################################################################################################
###########################################################################################################################################################################
####################################################################################################################################################################################

LDA_optimization_fun <- function(train_set,k,method='Gibbs',a=c(0.000001,0.00001,0.0001,0.001,0.01,0.1,1),verbose=TRUE){
  ######################################
  # This functions takes in a seurat object and optimizes for alpha by assessing which alpha is best to detect single cells and synthetic mixes
  #
  # Args
  # se_obj: seurat object
  # a: vectors of alpha parameter to tune, the higher the value the more homogeneous the dirichlet multinomial function will be (values towards the middle) and the lower the more polarized
  #
  # Returns:
  # This function returns a dataframe containing alpha, perplexity score and Jensen-Shannon divergence
  ######################################
  
  require(doParallel)
  
  # Set the number of cores
  det_cores <- parallel::detectCores()
  set_cores <- length(a)
  n_cores <- if_else(set_cores > det_cores, det_cores, set_cores)
  if (verbose) print(sprintf('The number of cores available are: %s, parallelization could use: %s. Calling to use: %s', det_cores, set_cores, n_cores)) 
  registerDoParallel(n_cores)
  
  # Parallelize the training
  r <- foreach(al=a) %dopar% {
    if (verbose) print(sprintf('Processing alpha: %s', al))
    # Setting control parameters
    control_LDA_Gibbs <- list(alpha = al, estimate.beta = TRUE,
                              verbose = 0, prefix = tempfile(), save = 0, keep = 1,
                              seed = as.integer(Sys.time()), nstart = 1, best = TRUE,
                              delta = 0.1, iter = 2000, burnin = 1000, thin = 10)
    
    # Run the model 10 times
    # multi_mod <- foreach(iter=1:10) %dopar% { 
    it <- 2000
    if (verbose) s_gibbs <- Sys.time()
    lda_mod_gibbs <- LDA(pbmc_lda_ready, k=k, method=method,control=control_LDA_Gibbs)
    if (verbose) print(sprintf('LDA Gibbs iteration: %s with alpha=%s took: %s', it, al, Sys.time() - s_gibbs)) # Takes ~10min
    # }
    return(lda_mod_gibbs)
  }
}


#############################################################################################################################################################################################
###########################################################################################################################################################################
####################################################################################################################################################################################

# lda_mod <- mod_ls[[3]]
get_top_topic_genes <- function(lda_mod,n_genes=100,perc,verbose=TRUE){
  ######################################
  # This functions takes in an topicmodels::lda model an returns a list with the most important genes for each topic by position
  # 
  # Args
  # lda_mod: topicmodels::lda model
  # perc: the top n% of genes we want to select
  # n: number of top genes to return for each topic, please note that if n is passed perc will be ignored
  #
  # Returns:
  # This function returns a list with the most important genes for each topic by position
  ######################################
  
  # Define if we take top n or top percentage
  if(!missing(n_genes)){ use = 'n'; if(!missing(perc)) warning('n_genes has been defines so perc is ignored')} else{  use = 'p'}
  
  # lda_mod@beta is a matrix of k*g; so k rows and g columns
  genes <- lda_mod@terms
  k <- paste('topic',1:lda_mod@k,sep = '')
  betas <- data.frame(t(lda_mod@beta))
  colnames(betas) <- k
  betas$gene <- genes
  
  top_gene_ls <- lapply(1:lda_mod@k, function(index){
    
    topic <- paste('topic',index,sep='')
    if (verbose) print(sprintf('Processing %s', topic))
    
    # Define if we take top n or top percentage
    if(use == 'n'){
      
      # Select genes above that threshold
      tmp_ds <- betas[order(betas[,index], decreasing = TRUE),][1:n_genes,c(index,lda_mod@k+1)]
      
    } else {
      
      # Defining likelyhood cutoof based on top perc
      thresh <- quantile(x = betas[,index],probs = 1-perc)
      # Select genes above that threshold
      tmp_ds <- betas[which(betas[,index] > thresh),c(topic,'gene')]
    }
    
    # return results to list to access faster downstream
    return(list(tmp_ds[,1],tmp_ds[,2]))
  })
  
  names(top_gene_ls) <- k
  return(top_gene_ls)
  
}


#############################################################################################################################################################################################
###########################################################################################################################################################################
####################################################################################################################################################################################

topic_assignment_maxlikelyhood <- function(lda_mod){
  ####
  # This function takes in a {topicmodels::lda} model and returns a topic assignment for each cell.
  #
  ####
  # Arguments
  # lda_mod: {topicmodels::lda} type object
  ####
  # Returns
  # returns a dataframe with 2 columns - cell barcode and topic assigned
  ####
  
  topic_assignment_ds <- lapply(1:nrow(lda_mod@gamma), function(r){
    
    # Select max likelyhood topic
    k <- which(max(lda_mod@gamma[r,]) == lda_mod@gamma[r,], lda_mod@gamma[r,])
    
    # Select cell barcode
    cell <- lda_mod@documents[r]
    
    # return dataframe with info
    return(data.frame('cell'=cell,topic=paste('topic',k,sep='')))
    
  }) %>% bind_rows() %>% data.frame()
  
  return(topic_assignment_ds)
}

#############################################################################################################################################################################################
###########################################################################################################################################################################
####################################################################################################################################################################################

convert_symb_entrez <- function(gene_vec){
  # This function converts a vector of gene symbols to its entrezID form so that it can be used to perform enrichment
  ###
  # gene_vec: vector of gene symbols to conver to entrezID
  ###
  # Returns a vecotr of entrezIDs
  ###
  
  library(org.Hs.eg.db)
  
  # Convert symbols to ENTREZID
  DE_entrezid <- mapIds(x = org.Hs.eg.db, keys = unlist(gene_vec), column = 'ENTREZID', keytype = 'SYMBOL')
  
  return(DE_entrezid)
  
}


#############################################################################################################################################################################################
###########################################################################################################################################################################
####################################################################################################################################################################################

gene_enrichment_GO <- function(gene_de, 
                               gene_universe,
                               pvalueCutoff = 0.05,
                               testDirection = 'over',
                               ontology = 'BP') {
  ####
  # From rmassonix github!
  # This function performs a hyper geometric test with a set of differentially expressed genes over a defined gene univers. 
  # It is set to be used for human genes only.
  # It requires the packages org.Hs.eg.db and GOstats to be installed.
  ####
  # gene_universe: vector with all the genes evaluated in order to set the gene universe. Must be ENTREZID names for the genes.
  # gene_de: vector with all the differentially expressed genes from the gene universe. Must be ENTREZID names for the genes.
  # pvalueCutoff: Cutoff to use for the Pvalue
  # testDirection: A string which can be either "over" or "under". This determines whether the test performed detects over or under represented GO terms.
  # gene_enrichment_GO: which ontology do we want to test BP (biological processes), CC (cellular component), or MF (molecular function)
  ####
  # Returns: hgOver with the enriched GO terms 
  ####
  # library(org.Hs.eg.db) # It is called in convert_symb_entrez
  library(GOstats)
  
  # Convert from symbol to entrezID
  target_enrich <- convert_symb_entrez(unique(gene_de))
  univers_g <- convert_symb_entrez(unique(gene_universe)) 
  
  # Create a GOHyperGParams instance
  params <- new("GOHyperGParams", geneIds=target_enrich, 
                universeGeneIds=univers_g,
                annotation="org.Hs.eg.db", ontology=ontology,
                pvalueCutoff=0.05, testDirection="over")
  
  # Carry out hyper geometric test
  hgOver <- hyperGTest(params)
  return(hgOver)
  
}

#############################################################################################################################################################################################
###########################################################################################################################################################################
####################################################################################################################################################################################

topic_profile_per_cluster <- function(lda_mod, se_obj, clust_vr){
  #####
  # This function takes in a seurat object and an LDA model generated from topic models and returns a list with a topic profile per cluster
  #
  ####
  # Params
  # lda_mod: LDA model from topicmodels package
  # se_obj: seurat object with the clustering performed
  # clust_vr: name of the cluster variable
  ####
  # Return
  # A dataframe of 
  ####
  library(tibble)
  library(dtplyr)
  
  se_obj$seurat_clusters <-  droplevels(se_obj@meta.data[,clust_vr])
  g_mtrx <- lda_mod@gamma # n of cells X n of topics
  colnames(g_mtrx) <- paste('topic_',1:ncol(g_mtrx),sep='')
  se_meta <- se_obj@meta.data
  
  clust_profiles <- cbind(se_meta,g_mtrx) %>%
    dtplyr::lazy_dt() %>%
    group_by(seurat_clusters) %>%
    select(seurat_clusters, paste('topic_',1:nlevels(se_obj$seurat_clusters),sep='')) %>%
    summarise_all(list(median)) %>%
    as.data.frame() %>% 
    column_to_rownames('seurat_clusters')
  
  colnames(clust_profiles) <- 1:ncol(clust_profiles)
  
  return(clust_profiles)
}

#############################################################################################################################################################################################
###########################################################################################################################################################################
####################################################################################################################################################################################

syn_spot_comb_topic <- function(clust_profiles, verbose = T){
  #####
  # This function takes in the cluster profiles and returns a sotchastic combination of all the possible ones 
  #
  ####
  # Params
  # clust_profiles: object returned from topic_profile_per_cluster
  # clust_vr: name of the variables with the cluster assignments
  ####
  # Return
  # Matrix with all possible combinations
  ####
  
  if(class(clust_profiles) != "matrix") clust_profiles <- as.matrix(clust_profiles)
  
  # If a cluster is 0 change it to 1
  if(sum(grepl(pattern = "0",rownames(clust_profiles))) != 0){
    rownames(clust_profiles) <- as.character(as.numeric(rownames(clust_profiles))+1)
  }
  
  # Compute all possible combinations up to grabbing round(nrow(comb)*0.66)
  # k_sub <- round(nrow(clust_profiles)*0.66)
  k_sub <- 8
  comb <- combinations(x = c(0:(nrow(clust_profiles))), k = k_sub, replace=TRUE)
  # comb <- combinations(x = rownames(clust_profiles), k = k_sub, replace=TRUE)
  # comb <- combinations(x = c(0:(ncol(clust_profiles))), k = 5, replace=TRUE)
  
  # Remove all those combinations that only include 1 or 2 cells
  comb <- comb[rowSums(comb != 0) > 2,]
  
  # Count 
  # comb_count <- apply(comb, 1, table)
  
  # Create all possible combinations
  ## Initialize matrix for increased speed so that it doesn't need to create indexes on the fly
  tmp_mtrx <- matrix(nrow = nrow(comb), ncol = ncol(clust_profiles))
  tmp_metadata <- matrix(nrow = nrow(comb), ncol = nrow(clust_profiles))
  
  if(verbose) print('Creating synthetic spots'); st_syn_spot <- Sys.time()
  if(verbose) pb_for <- txtProgressBar(min = 0, max = nrow(comb), style = 3) # Progress bar
  
  for (i in 1:nrow(comb)) {
    # Get how many cells of each type we have
    tt <- table(comb[i,][comb[i,]!=0])
    tmp_metadata[i,as.numeric(names(tt))] <- tt
    
    # Add all the profiles together
    row_i <- lapply(names(tt), function(nm){
      tmp_vec <- tt[[nm]]*clust_profiles[rownames(clust_profiles)[[as.numeric(nm)]],]
    }) %>% purrr::reduce(.,`+`)
    
    tmp_mtrx[i,] <- row_i/sum(tt)
    # update progress bar
    if(verbose) setTxtProgressBar(pb_for, i)
  }
  # rm(list(i,tt,row_i)) # For clean and good practice code, that there are no random tmp variables floating
  close(pb_for)
  if(verbose) print(sprintf('Creation of %s synthetic spot profiles took: %s minutes', 
                            nrow(comb), 
                            round(difftime(time1 = Sys.time(),time2 = st_syn_spot,units = 'mins'),2)))
  
  return(list(tmp_mtrx,tmp_metadata))
}

#############################################################################################################################################################################################
###########################################################################################################################################################################
####################################################################################################################################################################################

topic_viz <- function(lda_mod, k){
  
  require(wordcloud2)
  
  tmResult <- posterior(lda_mod)
  
  # visualize topics as word cloud
  topicToViz <- k # change for your own topic of interest
  
  # select to 40 most probable terms from the topic by sorting the term-topic-probability vector in decreasing order
  top40terms <- sort(tmResult$terms[topicToViz,], decreasing=TRUE)[1:40]
  words <- names(top40terms)
  
  # extract the probabilites of each of the 40 terms
  probabilities <- sort(tmResult$terms[topicToViz,], decreasing=TRUE)[1:40]
  
  # visualize the terms as wordcloud
  wc_plt <- wordcloud2(data.frame(words, probabilities), shuffle = FALSE, size = 0.8)  
  return(wc_plt)
}

#############################################################################################################################################################################################
###########################################################################################################################################################################
####################################################################################################################################################################################

top_n_predictions <- function(dist_mtrx,n){
  #####
  # This function gets a distance matrix and returns a list with the indices of the lowest (best) prediction indices
  # 
  #####
  # Params
  # dist_mtrx: distance matrix
  # n: number of best predictions wanted
  #####
  # Return
  # list of nrows(dist_mtrx) elements with ~n elements each
  #####
  
  #### Which is the 100th smalles value per row ####
  # We're subsetting it so that we can run JD on a subset of the data
  min_error <- Rfast::rownth(x = dist_mtrx, elems = rep(n, nrow(dist_mtrx)))
  
  # Get indices over which to calculate JD
  JD_indices <- lapply(1:length(min_error), function(i){
    which(dist_mtrx[i,] <= min_error[i])
  })
  return(JD_indices)
}

#############################################################################################################################################################################################
###########################################################################################################################################################################
####################################################################################################################################################################################

calculate_JSD_subset <- function(prediction, syn_spots_profiles, JSD_indices){
  #####
  # This function takes in 2 matrices of equal columns and calculates all pairwise JSD between them
  # 
  ####
  # Params
  # prediction: matrix with n rows and k columns where each column is a probabilistic array
  # syn_spots_profiles: matrix with m rows and k columns where each column is a probabilistic array
  # JSD_indices: list of m elements with with each element having row I indices of syn_spots_profiles
  ####
  # Return
  # this function returns a matrix with n rows and I columns the JSD values for the comparisons performed for each prediction
  ####
  
  #### Initialize JS matrix ####
  mtrx_JSD_full <- matrix(nrow = nrow(prediction), ncol = max(lengths(JSD_indices)))
  print('Calculating Jensen-Shannon Divergence')
  pb_JSD <- txtProgressBar(min = 0, max = nrow(prediction), style = 3)
  
  ##### Calculate Jensen-Shannon divergence of the subset of the data
  for (i in 1:nrow(prediction)) {
    # print(i)
    for (ii in 1:length(JSD_indices[[i]])) {
      # print(sprintf('nested:%s',ii))
      x <- rbind(prediction[i,], syn_spots_profiles[JSD_indices[[i]][ii],])
      mtrx_JSD_full[i,ii] <- JSD(x, unit = "log2")
    }
    # update progress bar
    setTxtProgressBar(pb_JSD,i)
  }
  close(pb_JSD)
  return(mtrx_JSD_full)
}

#############################################################################################################################################################################################
###########################################################################################################################################################################
####################################################################################################################################################################################

round_corners_grob <- function(plt) {
  library(ggplot2)
  library(grid)
  g <- ggplotGrob(plt)
  bg <- g$grobs[[1]]
  round_bg <- roundrectGrob(x = bg$x, 
                            y = bg$y, 
                            width = bg$width, 
                            height = bg$height,
                            r = unit(0.1, "snpc"),
                            just = bg$just, 
                            name = bg$name, 
                            gp = bg$gp, 
                            vp = bg$vp)
  g$grobs[[1]] <- round_bg
  return(g)
}
