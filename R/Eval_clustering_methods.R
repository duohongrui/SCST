#' Evaluate Clustering Methods Using Simulated Data
#'
#' @param SimulationResult A SCST object containing the simulated data.
#' @param eval_clustering Logical, whether to evaluate the clustering methods on the simulated data.
#' @param eval_cell_association Logical, whether to evaluate the methods for measuring the cell-cell association on the simulated data.
#' @param eval_annotation Logical, whether to evaluate the cell type annotation methods on the simulated data.
#' @param eval_DEGs Logical, whether to evaluate the methods for identifying DEGs on the simulated data.
#' @param eval_diff_abundance Logical, whether to evaluate the methods for testing differential abundance on the simulated data.
#' @param eval_batch_removal Logical, whether to evaluate the batch removal methods on the simulated data.
#' @param eval_trajectory Logical, whether to evaluate the methods of trajectory inference on the simulated data.
#' @param eval_spatial_domain Logical, whether to evaluate the methods for identifying spatial domains on the simulated spatial data.
#'
#' @export
EvaluateAppScenario <- function(
  SimulationResult,
  eval_clustering = FALSE,
  eval_cell_association = FALSE,
  eval_annotation = FALSE,
  eval_DEGs = FALSE,
  eval_diff_abundance = FALSE,
  eval_batch_removal = FALSE,
  eval_trajectory = FALSE,
  eval_spatial_domain = FALSE
){
  ############################## Clutering #####################################
  if(eval_clustering){
    eval_clustering_result <- eval_clustering(SimulationResult)
  }
  ############################# Cell-cell association###########################


  ############################## Annotation ####################################


  ############################## Detection of DEGs #############################


  ############################## Examine differential abundance ################


  ############################## Removal of batch effects ######################


  ############################## Trajectory ####################################


  ########################### Domain Detection #################################


}



#' Evaluation of Clustering Methods Based on Simulated Data
#'
#' @param SimulationResult A SCST object containing the simulated data.
#' @param min_cells Include genes detected in at least min_cells of cells. The genes detected below the number of min_cells are removed. Default no genes are removed.
#' @param min_genes Include cells where at least min_genes of genes are detected. The cells contain below the number of min_genes are removed. Default no genes are removed.
#' @param min_total_counts Include cells whose total counts are higher than min_total_counts.
#' @param min_counts_per_gene Include genes where more than min_counts_per_gene of counts are detected.
#' @param PCs The used number of principal components used. Default is 20.
#' @param n_cores The number of computational cores used in the clustering methods. Default is 1.
#' @param learning_rate Learning rate for training the artificial neural network. Default is 0.001.
#' @param batch_size Number of samples per training batch. Default is 64.
#' @param epochs Number of training epochs. Default is 200.
#' @param tsne_dims Output dimensionality of TSNE. Default is 3.
#' @param perplexity Perplexity parameter used for TSNE. Default is 30.
#' @param seed The random seed used for reproducibility.
#' @param verbose Whether to return messages or not. Default is TRUE.
#'
#' @export
#'
#' @importFrom stats prcomp kmeans
#' @importFrom tibble tibble
#' @importFrom Seurat AddMetaData
#'
EvalClusteringMethods <- function(SimulationResult,
                                  min_cells = 0,
                                  min_genes = 0,
                                  min_total_counts = 0,
                                  min_counts_per_gene = 1,
                                  PCs = 20,
                                  n_cores = 1,
                                  learning_rate = 0.001,
                                  batch_size = 64,
                                  epochs = 200,
                                  tsne_dims = 3,
                                  perplexity = 30,
                                  seed,
                                  verbose = TRUE){
  #-------- Check Necessary Information --------#
  if(length(SimulationResult@simulation_result) == 0){
    print_color_word("no No simulated data is found.")
    stop(error_output())
  }
  validated_methods <- check_info_pre_application(SimulationResult = SimulationResult, Group = TRUE)

  ### Iterate all simulation methods and perform clustering and sequent evaluation
  ClusteringResults <- purrr::map(validated_methods, .f = function(method){
    #-------- Filtering --------#
    data <- SimulationResult@simulation_result[[method]][["count_data"]]
    col_data <- SimulationResult@simulation_result[[method]][["col_meta"]]
    if(min_total_counts != 0){
      filter_index <- colSums(data) >= min_total_counts
      data <- data[, filter_index]
      col_data <- col_data[filter_index, ]
    }
    if(min_counts_per_gene != 0){
      data <- data[rowSums(data) >= min_counts_per_gene, ]
    }
    #-------- Methods --------#
    print_color_word(paste("------Perform Clustering On", method), color = "blue")
    ### 1. Seurat_Louvain
    seurat <- data %>%
      Seurat::CreateSeuratObject(min.cells = min_cells,
                                 min.features = min_genes,
                                 project = "simdata",
                                 meta.data = col_data) %>%
      Seurat::NormalizeData(verbose = FALSE) %>%
      Seurat::ScaleData(verbose = FALSE)
    ### Reset data and col_data
    if(packageVersion("Seurat") >= "5.0"){
      data <- SeuratObject::LayerData(seurat, layer = "counts") %>% as.matrix()
    }else{
      data <- methods::slot(Seurat::GetAssay(seurat), "counts") %>% as.matrix()
    }
    col_data <- seurat@meta.data
    ngroups <- length(unique(col_data$group))
    if(ncol(data) >= 8000){
      if(packageVersion("Seurat") >= "5.0"){
        seurat <- Seurat::FindVariableFeatures(seurat, verbose = FALSE)
        features <- SeuratObject::VariableFeatures(seurat)
      }else{
        seurat <- Seurat::FindVariableFeatures(seurat, verbose = FALSE)
        features <- Seurat::GetAssay(seurat)
        features <- features@var.features
      }
    }else{
      features <- rownames(data)
    }
    seurat <- seurat %>%
      Seurat::RunPCA(npcs = PCs,
                     features = features,
                     do.print = FALSE,
                     seed.use = seed,
                     verbose = FALSE) %>%
      Seurat::FindNeighbors(verbose = FALSE)
    resolution <- .find_resolution(seurat,
                                   groups = ngroups,
                                   algorithm = 1,
                                   seed = seed)
    seurat_louvain <- Seurat::FindClusters(object = seurat,
                                           algorithm = 1,
                                           resolution = resolution,
                                           random.seed = seed,
                                           verbose = FALSE)
    seurat_louvain_result <- as.numeric(seurat_louvain@meta.data$seurat_clusters)
    names(seurat_louvain_result) <- colnames(seurat)
    ### 2. Seurat_Leiden
    all_packages <- reticulate::py_list_packages()
    if("leidenalg" %in% all_packages$package){
      resolution <- .find_resolution(seurat,
                                     groups = ngroups,
                                     algorithm = 4,
                                     seed = seed)
      seurat_leiden <- Seurat::FindClusters(object = seurat,
                                            algorithm = 4,
                                            resolution = resolution,
                                            random.seed = seed,
                                            verbose = FALSE)
      seurat_leiden_result <- as.numeric(seurat_leiden@meta.data$seurat_clusters)
      names(seurat_leiden_result) <- colnames(seurat)
    }else{
      print_color_word("leidenalg python module has not been installed, so Seurat_Leiden will not be used for clustering", "yellow")
      seurat_leiden_result <- NULL
    }
    ### 3. SC3
    if(!requireNamespace("SC3")){
      message("SC3 is not installed on your device")
      message("Installing SC3...")
      BiocManager::install("SC3")
    }
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(
        counts = as.matrix(data),
        logcounts = log2(as.matrix(data) + 1)
      ),
      colData = col_data
    )
    SummarizedExperiment::rowData(sce)$feature_symbol <- rownames(sce)
    SC3_result <- SC3::sc3(sce,
                           pct_dropout_min = 5,
                           pct_dropout_max = 95,
                           biology = FALSE,
                           svm_max = 1e6,
                           k_estimator = FALSE,
                           ks = ngroups,
                           gene_filter = FALSE,
                           n_cores = n_cores,
                           rand_seed = seed)
    SC3_cluster <- as.numeric(SingleCellExperiment::colData(SC3_result)[, paste0("sc3_", ngroups, "_clusters")])
    names(SC3_cluster) <- colnames(sce)
    ### 4. SC3_SVM
    SC3_SVM_result <- SC3::sc3(sce,
                               pct_dropout_min = 5,
                               pct_dropout_max = 95,
                               biology = FALSE,
                               svm_max = 1,
                               svm_num_cells = round(ncol(sce) * 0.8),
                               k_estimator = FALSE,
                               ks = ngroups,
                               gene_filter = FALSE,
                               n_cores = n_cores,
                               rand_seed = seed)
    SC3_SVM_result <- SC3::sc3_run_svm(SC3_SVM_result, ks = ngroups)
    SC3_SVM_cluster <- as.numeric(SingleCellExperiment::colData(SC3_SVM_result)[, paste0("sc3_", ngroups, "_clusters")])
    names(SC3_SVM_cluster) <- colnames(sce)

    ### 5. scCCESS-kmeans
    if(!requireNamespace("scCCESS")){
      message("scCCESS is not installed on your device")
      message("Installing scCCESS...")
      devtools::install_github('PYangLab/scCCESS')
    }
    scCCESS_kmeans_result <- scCCESS::ensemble_cluster(
      data,
      seed = seed,
      cluster_func = function(x) {
        set.seed(seed)
        kmeans(x, centers = ngroups)
      },
      cores = n_cores,
      genes_as_rows = T,
      ensemble_sizes = 10,
      verbose = 0,
      scale = F,
      learning_rate = learning_rate,
      batch_size = batch_size,
      epochs = epochs
    )

    ### 6. scCCESS-SIMLR
    if(!requireNamespace("SIMLR")){
      message("SIMLR is not installed on your device")
      message("Installing SIMLR...")
      devtools::install_github("yulijia/SIMLR", ref = "master")
    }
    scCCESS_SIMLR_result <- scCCESS::ensemble_cluster(
      data,
      seed = seed,
      cluster_func = function(x) {
        set.seed(seed)
        SIMLR::SIMLR_Large_Scale(t(x), c = ngroups, kk = PCs)
      },
      cores = n_cores,
      genes_as_rows = T,
      ensemble_sizes = 10,
      verbose = 0,
      scale = F,
      learning_rate = learning_rate,
      batch_size = batch_size,
      epochs = epochs
    )
    names(scCCESS_SIMLR_result) <- colnames(data)

    ### 7. RtsneKmeans
    if(!requireNamespace("Rtsne")){
      message("Rtsne is not installed on your device")
      message("Installing Rtsne")
      utils::install.packages("Rtsne")
    }
    rtsne <- Rtsne::Rtsne(X = t(data),
                          dims = tsne_dims,
                          perplexity = perplexity,
                          pca = TRUE,
                          initial_dims = PCs,
                          check_duplicates = FALSE)
    RtsneKmeans_result <- stats::kmeans(rtsne$Y, centers = ngroups)$cluster
    names(RtsneKmeans_result) <- colnames(data)

    ### 8. PCAKmeans
    pca <- stats::prcomp(t(data), center = TRUE, scale. = FALSE, rank. = PCs)
    PCAKmeans_result <- kmeans(pca$x, centers = ngroups)$cluster

    ### 9. scLCA
    if(!requireNamespace("scLCA")){
      message("scLCA is not installed on your device")
      message("Installing scLCA")
      devtools::install_bitbucket("scLCA/single_cell_lca", ref = "master")
    }
    require(scLCA)
    scLCA_result <- scLCA::myscLCA(datmatrix = data,
                                   clust.max = ngroups,
                                   datBatch = NULL)
    scLCA_result <- scLCA_result[[1]]
    names(scLCA_result) <- colnames(data)

    ### 10. CIDR
    if(!requireNamespace("cidr")){
      message("CIDR is not installed on your device")
      message("Installing CIDR")
      devtools::install_github("VCCRI/CIDR")
    }
    sData <- data %>%
      cidr::scDataConstructor(tagType = "raw") %>%
      cidr::determineDropoutCandidates() %>%
      cidr::wThreshold() %>%
      cidr::scDissim(threads = n_cores) %>%
      cidr::scPCA(plotPC = FALSE) %>%
      cidr::nPC()
    ## Cluster with preset number of clusters
    sDataC <- cidr::scCluster(object = sData,
                              nCluster = ngroups,
                              nPC = sData@nPC,
                              cMethod = "ward.D2")
    CIDR_result <- sDataC@clusters
    names(CIDR_result) <- colnames(sDataC@tags)

    ### 11. Monocle3
    if(!requireNamespace("monocle3")){
      message("monocle3 is not installed on your device")
      message("Installing monocle3")
      devtools::install_github('cole-trapnell-lab/monocle3')
    }
    gene_metadata <- data.frame("gene_short_name" = rownames(data),
                                row.names = rownames(data))

    cds <- data %>%
      monocle3::new_cell_data_set(cell_metadata = col_data,
                                  gene_metadata = gene_metadata) %>%
      monocle3::preprocess_cds(method = "PCA",
                               num_dim = PCs,
                               verbose = verbose) %>%
      monocle3::reduce_dimension(reduction_method = "UMAP",
                                 preprocess_method = "PCA",
                                 verbose = verbose,
                                 cores = n_cores)
    resolution <- .find_resolution(cds,
                                   groups = ngroups,
                                   algorithm = 4,
                                   seed = seed)
    cds <- monocle3::cluster_cells(cds,
                                   reduction_method = "UMAP",
                                   k = 20,
                                   num_iter = 10,
                                   resolution = resolution,
                                   cluster_method = "leiden",
                                   verbose = FALSE)
    monocle3_result <- monocle3::clusters(cds) %>% as.numeric()
    names(monocle3_result) <- colnames(cds)

    ### collect all clustering methods
    all_clustering_results <- list(
      "seurat-louvain" = seurat_louvain_result,
      "seurat-leiden" = seurat_leiden_result,
      "SC3" = SC3_cluster,
      "SC3-SVM" = SC3_SVM_cluster,
      "scCCESS-kmeans" = scCCESS_kmeans_result,
      "scCCESS-SIMLR" = scCCESS_SIMLR_result,
      "RtsneKmeans" = RtsneKmeans_result,
      "PCAKmeans" = PCAKmeans_result,
      "scLCA" = scLCA_result,
      "CIDR" = CIDR_result,
      "monocle3" = monocle3_result
    )
    #-------- Evaluate Clustering Results --------#
    trueLabels <- col_data$"group"
    names(trueLabels) <- colnames(data)
    eval_clustering_table <- .CalculateClusteringMetrics(
      all_clustering_results,
      trueLabels,
      data,
      method)
    #-------- Add Clustering Results to Seurat For Visualization --------#
    all_clustering_results <- lapply(all_clustering_results, FUN = function(x){as.character(x)})
    seurat <- Seurat::AddMetaData(seurat, all_clustering_results)
    #-------- Outcome of one simulation method --------#
    list("seurat" = seurat,
         "eval_clustering_table" = eval_clustering_table)
  })
  names(ClusteringResults) <- validated_methods
  return(ClusteringResults)
}


.CalculateClusteringMetrics <- function(ClusteringResults, trueLabels, data, method){
  if(!requireNamespace("aricode")){
    message("aricode is not installed on your device")
    message("Installing aricode")
    utils::install.packages("aricode")
  }
  print_color_word(paste("------Clustering Evaluation For", method), color = "blue")
  cluster_eval_table <- purrr::map_dfr(1:length(ClusteringResults), .f = function(i){
    cluster_method <- names(ClusteringResults)[i]
    clustering_result <- ClusteringResults[[cluster_method]]
    ### 1. ARI
    ARI <- aricode::ARI(trueLabels, clustering_result)
    ### 2. NMI
    NMI <- aricode::NMI(trueLabels, clustering_result)
    ### 3. CDI
    if(!requireNamespace("CDI", quietly = TRUE)){
      message("Installing CDI package...")
      devtools::install_github('jichunxie/CDI')
    }
    ## cluster
    cell_label <- data.frame("cluster" = as.numeric(as.factor(clustering_result)))
    ## feature selection
    feature_gene_index <- CDI::feature_gene_selection(
      data,
      batch_label = NULL,
      method = "wds",
      nfeature = 1000
    )
    ## subset
    sub_data <- data[feature_gene_index, ]
    ## size factor
    size_factor <- CDI::size_factor(data)
    ## calculate CDI
    error <- try(
      CDI <- CDI::calculate_CDI(
        sub_data,
        cand_lab_df = cell_label,
        batch_label = NULL,
        cell_size_factor = size_factor
      ),
      silent = TRUE
    )
    if("try-error" %in% class(error)){
      print_color_word("no The calculation of CDI failed.")
      CDI <- NA
    }else{
      CDI <- min(CDI[1, 1], CDI[1, 2])
    }

    ### 4. Average Silhouette Width
    if(!requireNamespace("cluster", quietly = TRUE)){
      message("Installing cluster package...")
      utils::install.packages("cluster")
    }
    if(!requireNamespace("parallelDist", quietly = TRUE)){
      message("Installing parallelDist package...")
      utils::install.packages("parallelDist")
    }
    dist <- parallelDist::parDist(t(data))
    silhouette_width <- cluster::silhouette(x = as.numeric(as.factor(clustering_result)), dist)
    average_silhouette <- mean(silhouette_width[, 3])

    ### 5. Connectiity
    if(!requireNamespace("clValid", quietly = TRUE)){
      message("Installing clValid package...")
      utils::install.packages("clValid")
    }
    cluster_info <- as.numeric(as.factor(clustering_result))
    connectivity <- clValid::connectivity(distance = dist, clusters = cluster_info)

    ### 6. Dunn Index
    dunn <- clValid::dunn(distance = dist, clusters = cluster_info)

    #### Save results
    tibble::tibble("Simulation_Method" = method,
                   "Clustering_Method" = cluster_method,
                   "ARI" = ARI,
                   "NMI" = NMI,
                   "CDI" = CDI,
                   "ASW" = average_silhouette,
                   "Connectivity" = connectivity,
                   "Dunn_Index" = dunn)
  })
  ## normalize some values
  cluster_eval_table <- cluster_eval_table %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of("CDI"),
                    ~ replace(.x, .x == "Inf", values = NA))
    ) %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(c("CDI", "Connectivity", "Dunn_Index")),
                    ~ pnorm((.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
    )
  ### Subtract values by 1
  colume_name <- c("CDI", "Connectivity")
  cluster_eval_table <- cluster_eval_table %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(colume_name), ~ 1 - .x)
    )
  ### return
  return(cluster_eval_table)
}


