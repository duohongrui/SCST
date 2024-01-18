#' Evaluation of Batch-Removal Methods Based on Simulated Data
#'
#' @param conda_env The path of conda environment. Additionally, the packages of scvi/scanpy/scanorama/scgen are strongly recommended to be installed by users.
#' @param SimulationResult A SCST object containing the simulated data.
#' @param min_cells Include genes detected in at least min_cells of cells. The genes detected below the number of min_cells are removed. Default no genes are removed.
#' @param min_genes Include cells where at least min_genes of genes are detected. The cells contain below the number of min_genes are removed. Default no genes are removed.
#' @param min_total_counts Include cells whose total counts are higher than min_total_counts.
#' @param min_counts_per_gene Include genes where more than min_counts_per_gene of counts are detected.
#' @param nFeatures The number of highly variable genes. Default is 2000.
#' @param PCs The number of principal components used. Default is 20 and the value higher than 10 is recommended.
#' @param use_cuda Whether to use cuda to accelerate the training process or use cpu for computation. If true, make sure that you have installed the necessary packages with right versions in your Conda Environment.
#' @param nlayers Number of hidden layers used for encoder and decoder NNs. Default is 2.
#' @param nlatent The dimensionality of the latent space. Default is 20.
#' @param genelikelihood The distribution model for the data, which is one of nb (Negative binomial distribution), zinb (Zero-inflated negative binomial distribution) and poisson (Poisson distribution). Default is nb.
#' @param batch_size Number of samples per training batch. Default is 64.
#' @param epochs Number of training epochs. Default is NULL and it is determined automatically.
#' @param k_NNs An integer scalar specifying the number of nearest neighbors to consider when identifying MNNs.
#' @param seed The random seed used for reproducibility.
#' @param verbose Whether to return messages or not. Default is TRUE.
#'
#' @export
#'
#' @importFrom reticulate import use_condaenv py_install r_to_py py_to_r py_list_packages
#' @importFrom SeuratObject LayerData
#' @importFrom Seurat CreateDimReducObject DefaultAssay
#' @importFrom Matrix t
#'
EvalBatchRemovalMethods <- function(SimulationResult,
                                    conda_env,
                                    min_cells = 0,
                                    min_genes = 0,
                                    min_total_counts = 0,
                                    min_counts_per_gene = 1,
                                    nFeatures = 2000,
                                    PCs = 20,
                                    use_cuda = FALSE,
                                    nlayers = 2,
                                    nlatent = 20,
                                    genelikelihood = "nb",
                                    batch_size = 64,
                                    epochs = NULL,
                                    k_NNs = 20,
                                    seed,
                                    verbose = TRUE){
  #-------- Check Necessary Information --------#
  if(length(SimulationResult@simulation_result) == 0){
    print_color_word("no No simulated data is found.")
    stop(error_output())
  }
  validated_methods <- check_info_pre_application(SimulationResult = SimulationResult, Batch = TRUE)
  count_data <- get_count_data(SimulationResult)
  col_data <- get_cell_meta(SimulationResult)
  reticulate::use_condaenv(conda_env, required = TRUE)
  packages <- reticulate::py_list_packages(envname = conda_env)
  if(use_cuda){
    device <- "cuda"
    use_gpu <- reticulate::r_to_py(TRUE)
  }else{
    device <- "cpu"
    use_gpu <- reticulate::r_to_py(FALSE)
  }
  ### Iterate all simulation methods and perform clustering and sequent evaluation
  BatchRemovalResults <- purrr::map(validated_methods, .f = function(method){
    #-------- Preprocessing --------#
    print_color_word(paste("\u2192", "Preprocessing..."), color = "green")
    data <- count_data[[method]]
    cell_meta <- col_data[[method]]
    if(min_total_counts != 0){
      filter_index <- colSums(data) >= min_total_counts
      data <- data[, filter_index]
      cell_meta <- cell_meta[filter_index, ]
    }
    if(min_counts_per_gene != 0){
      data <- data[rowSums(data) >= min_counts_per_gene, ]
    }
    seurat <- data %>%
      methods::as("dgCMatrix") %>%
      Seurat::CreateSeuratObject(min.cells = min_cells,
                                 min.features = min_genes,
                                 project = "simdata",
                                 meta.data = cell_meta) %>%
      Seurat::NormalizeData(verbose = FALSE) %>%
      Seurat::ScaleData(verbose = FALSE)
    ### Reset data and col_data
    if(packageVersion("Seurat") >= "5.0"){
      seurat <- Seurat::FindVariableFeatures(seurat, verbose = FALSE, nfeatures = nFeatures)
      features <- SeuratObject::VariableFeatures(seurat)
    }else{
      seurat <- Seurat::FindVariableFeatures(seurat, verbose = FALSE, nfeatures = nFeatures)
      features <- Seurat::GetAssay(seurat)
      features <- features@var.features
    }
    seurat <- seurat %>%
      Seurat::RunPCA(npcs = PCs,
                     features = features,
                     do.print = FALSE,
                     seed.use = seed,
                     verbose = FALSE) %>%
      Seurat::FindNeighbors(verbose = FALSE) %>%
      Seurat::RunUMAP(dims = 1:PCs, seed.use = seed, verbose = FALSE) %>%
      Seurat::RunTSNE(dims = 1:PCs, seed.use = seed, verbose = FALSE)

    #### Python packages
    sc <-  reticulate::import("scanpy", convert = FALSE)
    scipy <-  reticulate::import("scipy", convert = FALSE)
    scvi <-  reticulate::import("scvi", convert = FALSE)
    if(packageVersion("Seurat") >= "5.0"){
      X <- SeuratObject::LayerData(seurat, layer = "counts")[features, ]
    }else{
      X <- methods::slot(Seurat::GetAssay(seurat), "counts")[features, ]
    }
    adata = sc$AnnData(
      X   = scipy$sparse$csr_matrix(Matrix::t(X)),
      obs = seurat@meta.data
    )
    sc$pp$normalize_total(adata)
    sc$pp$log1p(adata)
    sc$pp$highly_variable_genes(adata, n_top_genes = nFeatures)
    sc$tl$pca(adata, n_comps = as.integer(PCs))
    if(is.null(epochs)) {
      epochs <- reticulate::r_to_py(x = epochs)
    }else{
      epochs <- as.integer(x = epochs)
    }
    #-------- Methods --------#
    print_color_word(paste("------Perform Batch Removal On", method), color = "blue")

    ### 4. scGen
    scGen_moni <- peakRAM::peakRAM(
      seurat <- .scGen(adata = adata,
                       seurat = seurat,
                       conda_env = conda_env,
                       packages = packages,
                       use_gpu = use_gpu,
                       nlayers = nlayers,
                       nlatent = nlatent,
                       epochs = epochs,
                       batch_size = batch_size)
    )

    ### 1. scVI
    scVI_moni <- peakRAM::peakRAM(
      seurat <- .scVI(adata = adata,
                      seurat = seurat,
                      conda_env = conda_env,
                      packages = packages,
                      scvi = scvi,
                      use_gpu = use_gpu,
                      nlayers = nlayers,
                      nlatent = nlatent,
                      epochs = epochs,
                      batch_size = batch_size,
                      genelikelihood = genelikelihood)
    )
    model <- seurat[["model"]]
    seurat <- seurat[["seurat"]]

    ### 2. scANVI
    scANVI_moni <- peakRAM::peakRAM(
      seurat <- .scANVI(adata = adata,
                        seurat = seurat,
                        model = model,
                        scvi = scvi,
                        use_gpu = use_gpu,
                        epochs = epochs,
                        batch_size = batch_size)
    )

    ### 3. Scanorama
    Scanorama_moni <- peakRAM::peakRAM(
      seurat <- .Scanorama(adata = adata,
                           seurat = seurat,
                           conda_env = conda_env,
                           packages = packages,
                           sc = sc,
                           batch_size = batch_size)
    )

    ### 5. fastMNN
    fastMNN_moni <- peakRAM::peakRAM(
      seurat <- .fastMNN(seurat = seurat,
                         k_NNs = k_NNs,
                         PCs = PCs,
                         features = features)
    )

    ### 6. Harmony
    Harmony_moni <- peakRAM::peakRAM(
      seurat <- .Harmony(seurat = seurat,
                         PCs = PCs,
                         verbose = verbose,
                         seed = seed)
    )

    #-------- Evaluate Batch Removal Results --------#
    eval_batch_removal_table <- .CalculateBatchRemovalMetrics(
      seurat,
      method
    )
    #-------- Record Resource Occupation During Execution --------#
    resource_monitering <- tibble::tibble(
      "Simulation_Method" = method,
      "Batch_Removal_Method" = c("scVI", "scANVI", "Scanorama", "scGen", "fastMNN", "Harmony"),
      "Time" = c(scVI_moni[, 2],
                 scANVI_moni[, 2],
                 Scanorama_moni[, 2],
                 scGen_moni[, 2],
                 fastMNN_moni[, 2],
                 Harmony_moni[, 2]),
      "Memory" = c(scVI_moni[, 4],
                   scANVI_moni[, 4],
                   Scanorama_moni[, 4],
                   scGen_moni[, 4],
                   fastMNN_moni[, 4],
                   Harmony_moni[, 4]),
      "Device" = ifelse(use_cuda,
                        c("cuda", "cuda", "cpu", "cuda", "cpu", "cpu"),
                        "cpu")
    )
    #-------- Outcome of one simulation method --------#
    list("seurat" = seurat,
         "eval_batch_removal_table" = eval_batch_removal_table,
         "resource_monitering" = resource_monitering)
  })
  names(BatchRemovalResults) <- validated_methods
  return(BatchRemovalResults)
}


.CalculateBatchRemovalMetrics <- function(BatchRemovalResults, method){
  #-------- Check --------#
  assertthat::assert_that(methods::is(BatchRemovalResults, "Seurat"))
  print_color_word(paste("------Batch Removal Evaluation For", method), color = "blue")
  #-------- Preparing --------#
  col_data <- BatchRemovalResults@meta.data
  batch_info <- col_data[, "batch"]
  k <- sqrt(length(batch_info)) %>% round()
  reducedim_names <- Seurat::Reductions(BatchRemovalResults)
  reducedim_names <- reducedim_names[-grep("pca", reducedim_names)]
  reducedim_names <- reducedim_names[-grep("umap", reducedim_names)]
  reducedim_names <- reducedim_names[-grep("tsne", reducedim_names)]
  sce <- SingleCellExperiment::SingleCellExperiment(colData = col_data)
  #-------- Calculating Metrics --------#
  batch_removal_eval_table <- purrr::map_dfr(1:length(reducedim_names), .f = function(i){
    #### integration method
    integration_method <- reducedim_names[i]
    print_color_word(paste("\u2192", integration_method, "is being evaluated..."), color = "green")
    ##### Extract embeddings
    embeddings <- BatchRemovalResults@reductions[[integration_method]]@cell.embeddings
    ##### Perform PCA
    pca <- stats::prcomp(embeddings,
                         center = FALSE,
                         scale. = FALSE,
                         rank. = 10)
    ##### Add PCA reduction into sce
    SingleCellExperiment::reducedDim(sce, "PCA") <- pca$x %>% as.data.frame()

    ### 1. kBET
    if(!requireNamespace("kBET", quietly = TRUE)){
      message("Install kBET...")
      devtools::install_github('theislab/kBET')
    }
    batch_estimate <- kBET::kBET(df = embeddings,
                                 do.pca = FALSE,
                                 batch_info,
                                 plot = FALSE,
                                 k0 = 10)
    kBET <- batch_estimate$summary$kBET.observed %>% mean(na.rm = TRUE)

    ### 2. LISI
    if(!requireNamespace("lisi", quietly = TRUE)){
      message("Install lisi...")
      devtools::install_github("immunogenomics/lisi")
    }
    LISI <- lisi::compute_lisi(embeddings,
                               meta_data = col_data,
                               label_colnames = "batch",
                               perplexity = k)
    LISI <- mean(LISI$batch, na.rm = TRUE) / length(unique(batch_info))

    ### 3. ASW batch
    ASW_batch <- kBET::batch_sil(pca.data = pca,
                                 batch = as.factor(batch_info))
    ASW_batch <- 1 - abs(ASW_batch)

    ### 4. PCR
    batch_pcr <- kBET::pcRegression(pca.data = pca,
                                    batch = batch_info)
    PCR <- batch_pcr$R2Var

    ### 5. CMS
    if(!requireNamespace("CellMixS", quietly = TRUE)){
      message("Install CellMixS...")
      BiocManager::install('CellMixS')
    }
    results <- CellMixS::evalIntegration(metrics = c("cms",
                                                     "mixingMetric",
                                                     "entropy"),
                                         sce = sce,
                                         group = "batch",
                                         k = k)
    results <- as.data.frame(SingleCellExperiment::colData(results))
    CMS <- results$cms
    CMS <- mean(CMS, na.rm = TRUE)

    ### 6. Mixing Metric
    MM <- results$mm
    MM <- mean(MM, na.rm = TRUE)

    ### 7. Shannon Entropy
    Entropy <- results$entropy
    Entropy <- mean(Entropy, na.rm = TRUE)

    #### Save results
    tibble::tibble("Simulation_Method" = method,
                   "Integration_Method" = integration_method,
                   "kBET" = kBET,
                   "LISI" = LISI,
                   "PCR" = PCR,
                   "ASW_batch" = ASW_batch,
                   "CMS" = CMS,
                   "MM" = MM,
                   "Shannon_entropy" = Entropy)
  })

  ## normalize some values
  batch_removal_eval_table <- batch_removal_eval_table %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(c("CMS", "MM", "PCR")),
                    ~ pnorm((.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
    )
  ### Subtract values by 1
  colume_name <- c("kBET", "MM", "PCR")
  batch_removal_eval_table <- batch_removal_eval_table %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(colume_name), ~ 1 - .x)
    )
  ### return
  return(batch_removal_eval_table)
}


.scGen <- function(adata,
                   seurat,
                   conda_env,
                   packages,
                   use_gpu,
                   nlayers,
                   nlatent,
                   epochs,
                   batch_size){
  #### pip install git+https://github.com/theislab/scgen.git
  if(!"scgen" %in% packages$package){
    reticulate::py_install("scgen", envname = conda_env, pip = TRUE, ignore_installed = TRUE)
  }
  print_color_word(paste("\n \u2192", "scGen is running..."), color = "green")
  scgen <- reticulate::import("scgen")
  scgen$SCGEN$setup_anndata(adata, batch_key = "batch", labels_key = "group")
  scgen_model = scgen$SCGEN(adata, n_layers = as.integer(nlayers), n_latent = as.integer(nlatent))
  scgen_model$train(
    max_epochs = epochs,
    batch_size = as.integer(batch_size),
    use_gpu = use_gpu
  )
  scgen_result = scgen_model$batch_removal()
  scgen_result = scgen_result$obsm["corrected_latent"]
  scgen_result <- as.matrix(scgen_result)
  rownames(scgen_result) <- reticulate::py_to_r(adata$obs$index$values)
  colnames(scgen_result) <- paste0("scgen", "_", 1:ncol(scgen_result))
  reduction <- Seurat::CreateDimReducObject(
    embeddings = scgen_result,
    assay = Seurat::DefaultAssay(seurat),
    key = "scGen_"
  )
  seurat@reductions <- append(seurat@reductions, list("scGen" = reduction))
  return(seurat)
}


.scVI <- function(adata,
                  seurat,
                  conda_env,
                  packages,
                  scvi,
                  use_gpu,
                  nlayers,
                  nlatent,
                  epochs,
                  batch_size,
                  genelikelihood){
  if(!"scanpy" %in% packages$package){
    reticulate::py_install("scanpy", envname = conda_env, pip = TRUE)
  }
  if(!"scvi-tools" %in% packages$package){
    reticulate::py_install("scvi-tools", envname = conda_env, pip = TRUE)
  }
  print_color_word(paste("\n \u2192", "scVI is running..."), color = "green")
  scvi$model$SCVI$setup_anndata(adata, batch_key = "batch")
  model = scvi$model$SCVI(adata = adata,
                          n_latent = as.integer(x = nlatent),
                          n_layers = as.integer(x = nlayers),
                          gene_likelihood = genelikelihood)
  #### Training model
  model$train(max_epochs = epochs, batch_size = as.integer(batch_size), use_gpu = use_gpu)
  scvi_result = model$get_latent_representation()
  scvi_result <- as.matrix(scvi_result)
  rownames(scvi_result) <- reticulate::py_to_r(adata$obs$index$values)
  colnames(scvi_result) <- paste0("scvi", "_", 1:ncol(scvi_result))
  reduction <- Seurat::CreateDimReducObject(
    embeddings = scvi_result,
    assay = Seurat::DefaultAssay(seurat),
    key = "scVI_"
  )
  seurat@reductions <- append(seurat@reductions, list("scVI" = reduction))
  return(list(seurat = seurat,
              model = model))
}


.scANVI <- function(adata,
                    seurat,
                    model,
                    scvi,
                    use_gpu,
                    epochs,
                    batch_size){
  print_color_word(paste("\n \u2192", "scANVI is running..."), color = "green")
  scanvi_model = scvi$model$SCANVI$from_scvi_model(
    model,
    adata = adata,
    unlabeled_category = "Unknown",
    labels_key = "group"
  )
  scanvi_model$train(max_epochs = epochs,
                     batch_size = as.integer(batch_size),
                     use_gpu = use_gpu)
  scANVI_result = scanvi_model$get_latent_representation()
  scANVI_result <- reticulate::py_to_r(scANVI_result)
  rownames(scANVI_result) <- reticulate::py_to_r(adata$obs$index$values)
  colnames(scANVI_result) <- paste0("scANVI", "_", 1:ncol(scANVI_result))
  scANVI_result <- as.matrix(scANVI_result)
  reduction <- Seurat::CreateDimReducObject(
    embeddings = scANVI_result,
    assay = Seurat::DefaultAssay(seurat),
    key = "scANVI_"
  )
  seurat@reductions <- append(seurat@reductions, list("scANVI" = reduction))
  return(seurat)
}


.Scanorama <- function(adata,
                       seurat,
                       conda_env,
                       packages,
                       sc,
                       batch_size){
  if(!"scanorama" %in% packages$package){
    reticulate::py_install("scanorama", envname = conda_env, pip = TRUE)
  }
  print_color_word(paste("\n \u2192", "Scanorama is running..."), color = "green")
  #### Train models
  sc$external$pp$scanorama_integrate(adata, key = "batch", batch_size = as.integer(batch_size))
  Scanorama_result <- as.matrix(adata$obsm["X_scanorama"])
  rownames(Scanorama_result) <- reticulate::py_to_r(adata$obs$index$values)
  colnames(Scanorama_result) <- paste0("scanorama", "_", 1:ncol(Scanorama_result))
  reduction <- Seurat::CreateDimReducObject(
    embeddings = Scanorama_result,
    assay = Seurat::DefaultAssay(seurat),
    key = "Scanorama_"
  )
  seurat@reductions <- append(seurat@reductions, list("Scanorama" = reduction))
  return(seurat)
}


.fastMNN <- function(seurat,
                     k_NNs,
                     PCs,
                     features){
  if(!requireNamespace("batchelor")){
    message("batchelor is not installed on your device")
    message("Installing batchelor...")
    BiocManager::install("batchelor")
  }
  if(packageVersion("Seurat") >= "5.0"){
    data <- SeuratObject::LayerData(seurat, layer = "counts") %>% as.matrix()
  }else{
    data <- methods::slot(Seurat::GetAssay(seurat), "counts") %>% as.matrix()
  }
  print_color_word(paste("\n \u2192", "fastMNN is running..."), color = "green")
  cell_meta <- seurat@meta.data
  batch_label <- unique(cell_meta[, "batch"])
  fastMNN_input <- lapply(batch_label, FUN = function(x){
    index <- grep(x, cell_meta[, "batch"])
    tmp_data <- data[features, index]
    tmp_data
  })
  names(fastMNN_input) <- batch_label
  suppressWarnings(fastMNN <- batchelor::fastMNN(fastMNN_input, k = k_NNs, d = PCs))
  fastMNN_result <- SingleCellExperiment::reducedDim(fastMNN)
  rownames(fastMNN_result) <- colnames(fastMNN)
  colnames(fastMNN_result) <- paste0("fastMNN", "_", 1:ncol(fastMNN_result))
  reduction <- Seurat::CreateDimReducObject(
    embeddings = fastMNN_result,
    assay = Seurat::DefaultAssay(seurat),
    key = "fastMNN_"
  )
  seurat@reductions <- append(seurat@reductions, list("fastMNN" = reduction))
  return(seurat)
}


.Harmony <- function(seurat,
                     PCs,
                     verbose,
                     seed){
  if(!requireNamespace("harmony")){
    message("harmony is not installed on your device")
    message("Installing harmony...")
    devtools::install_github("immunogenomics/harmony")
  }
  print_color_word(paste("\n \u2192", "Harmony is running..."), color = "green")
  set.seed(seed)
  seurat <- harmony::RunHarmony(seurat,
                                group.by.vars = "batch",
                                reduction = "pca",
                                dims.use = 1:PCs,
                                verbose = verbose,
                                reduction.save = "Harmony",
                                max.iter.harmony = 50)
  return(seurat)
}
