#' Evaluation of Domain Detection Methods Based on Simulated Spatial Data
#'
#' @param SimulationResult A SCST object containing the simulated spatial data.
#' @param conda_env The path of conda environment. Additionally, the packages or the environment (Python packages: DeepST, STAGATE, GraphST and SCGDL) should be installed by users previously. The detail configuration in conda is demonstrated below.
#' @param min_cells Include genes detected in at least min_cells of cells. The genes detected below the number of min_cells are removed. Default no genes are removed.
#' @param min_genes Include cells where at least min_genes of genes are detected. The cells contain below the number of min_genes are removed. Default no genes are removed.
#' @param min_total_counts Include cells whose total counts are higher than min_total_counts.
#' @param min_counts_per_gene Include genes where more than min_counts_per_gene of counts are detected.
#' @param nFeatures The number of highly variable genes. Default is 2000.
#' @param PCs The number of principal components used. Default is 20.
#' @param use_cuda Whether to use cuda to accelerate the training process or use cpu for computation. If true, make sure that you have installed the necessary packages with right versions in your Conda Environment.
#' @param learning_rate Learning rate for training the artificial neural network. Default is 0.001.
#' @param epochs Number of training epochs. Default is 1000.
#' @param seed The random seed used for reproducibility.
#' @param verbose Whether to return messages or not. Default is TRUE.
#'
#' @export

EvalDomainDetectionMethods <- function(SimulationResult,
                                       conda_env,
                                       min_cells = 0,
                                       min_genes = 0,
                                       min_total_counts = 0,
                                       min_counts_per_gene = 1,
                                       nFeatures = 2000,
                                       PCs = 20,
                                       use_cuda = FALSE,
                                       learning_rate = 0.001,
                                       epochs = 1000,
                                       seed,
                                       verbose = TRUE){
  #-------- Check Necessary Information --------#
  if(length(SimulationResult@simulation_result) == 0){
    print_color_word("no No simulated data is found.")
    stop(error_output())
  }
  validated_methods <- check_info_pre_application(SimulationResult = SimulationResult,
                                                  Group = TRUE,
                                                  Spatial = TRUE)
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
  DomainDetectionResults <- purrr::map(validated_methods, .f = function(method){
    #-------- Preprocessing --------#
    print_color_word(paste("\n \u2192", "Preprocessing..."), color = "green")
    count <- count_data[[method]]
    cell_meta <- col_data[[method]]
    n_domains <- length(unique(cell_meta[, "group"])) %>% as.integer()
    location <- cell_meta %>%
      dplyr::select("spatial_x", "spatial_y") %>%
      dplyr::rename("imagecol" = "spatial_x",
                    "imagerow" = "spatial_y")
    cell_meta <- cell_meta %>%
      dplyr::rename("col" = "spatial_x",
                    "row" = "spatial_y")
    #-------- Filtering --------#
    seurat <- Seurat::CreateSeuratObject(counts = count %>% methods::as("dgCMatrix"),
                                         project = "SlideSeq",
                                         assay = "Spatial",
                                         min.cells = min_cells,
                                         min.features = min_genes,
                                         meta.data = cell_meta)
    seurat[['image']] = methods::new(Class = "SlideSeq",
                                     assay = "Spatial",
                                     coordinates = location)
    seurat <- seurat %>%
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
      Seurat::RunTSNE(dims = 1:PCs, seed.use = seed, verbose = FALSE, check_duplicates = FALSE)

    #-------- Python Object --------#
    sc <-  reticulate::import("scanpy", convert = FALSE)
    scipy <-  reticulate::import("scipy", convert = FALSE)
    np <- reticulate::import("numpy", convert = FALSE)
    if(packageVersion("Seurat") >= "5.0"){
      X <- SeuratObject::LayerData(seurat, layer = "counts")[features, ]
    }else{
      X <- methods::slot(Seurat::GetAssay(seurat), "counts")[features, ]
    }
    ann_cell_meta <- seurat@meta.data %>%
      dplyr::rename("arraycol" = "col",
                    "arrayrow" = "row")
    spatial_coordiantes <- reticulate::r_to_py(list("spatial" = ann_cell_meta %>%
                                                      dplyr::select("arraycol", "arrayrow") %>%
                                                      as.matrix()))
    adata = sc$AnnData(
      X   = scipy$sparse$csr_matrix(Matrix::t(X)),
      obs = ann_cell_meta,
      obsm = spatial_coordiantes,
      dtype = np$float32
    )
    sc$pp$normalize_total(adata)
    sc$pp$log1p(adata)
    sc$pp$highly_variable_genes(adata, n_top_genes = as.integer(nFeatures))
    sc$tl$pca(adata, n_comps = as.integer(PCs))

    #-------- Domain Detection Methods --------#
    ### 1. STAGATE (STAGATE_path)
    STAGATE_moni <- peakRAM::peakRAM(
      STAGATE_result <- .STAGATE(adata = adata,
                                 seurat = seurat,
                                 epochs = epochs,
                                 learning_rate = learning_rate,
                                 device = device,
                                 n_domains = n_domains,
                                 seed = seed,
                                 verbose = verbose)
    )

    ### 2. DeepST
    DeepST_moni <- peakRAM::peakRAM(
      DeepST_result <- .DeepST(conda_env = conda_env,
                               packages = packages,
                               adata = adata,
                               seurat = seurat,
                               epochs = epochs,
                               use_gpu = use_gpu,
                               n_domains = n_domains)
    )

    ### 3. Seurat_Louvain
    seurat_louvain_moni <- peakRAM::peakRAM(
      seurat_louvain_result <- .seurat_louvain(seurat = seurat,
                                               n_domains = n_domains,
                                               seed = seed)
    )

    ### 4. Seurat_Leiden
    seurat_leiden_moni <- peakRAM::peakRAM(
      seurat_leiden_result <- .seurat_leiden(seurat = seurat,
                                             packages = packages,
                                             conda_env = conda_env,
                                             n_domains = n_domains,
                                             seed = seed)
    )

    ### 5. BASS
    BASS_moni <- peakRAM::peakRAM(
      BASS_result <- .BASS(X = X,
                           location = location,
                           n_domains = n_domains,
                           nFeatures = nFeatures,
                           PCs = PCs,
                           seed = seed)
    )

    ### 6. SCGDL
    SCGDL_moni <- peakRAM::peakRAM(
      SCGDL_result <- .SCGDL(adata = adata,
                             seurat = seurat,
                             epochs = epochs,
                             learning_rate = learning_rate,
                             n_domains = n_domains,
                             seed = seed)
    )

    ### 7. GraphST
    GraphST_moni <- peakRAM::peakRAM(
      GraphST_result <- .GraphST(conda_env = conda_env,
                                 packages = packages,
                                 adata = adata,
                                 seurat = seurat,
                                 epochs = epochs,
                                 learning_rate = learning_rate,
                                 n_domains = n_domains,
                                 seed = seed,
                                 device = device)
    )

    ### 8. DR-SC
    DR.SC_moni <- peakRAM::peakRAM(
      DR.SC_result <- .DR.SC(seurat = seurat,
                             n_domains = n_domains,
                             verbose = verbose)
    )

    ### collect all clustering methods
    all_domain_detection_results <- list(
      "DeepST" = DeepST_result,
      "seurat-leiden" = seurat_leiden_result,
      "seurat-louvain" = seurat_louvain_result,
      "BASS" = BASS_result,
      "GraphST" = GraphST_result,
      "STAGATE" = STAGATE_result,
      "SCGDL" = SCGDL_result,
      "DR.SC" = DR.SC_result
    )
    #-------- Evaluate Clustering Results --------#
    trueLabels <- seurat$group
    names(trueLabels) <- colnames(seurat)
    eval_domain_detection_table <- .CalculateDomainDetectionMetrics(
      all_domain_detection_results,
      trueLabels,
      method)
    #-------- Record Resource Occupation During Execution --------#
    resource_monitering <- tibble::tibble(
      "Simulation_Method" = method,
      "Clustering_Method" = names(all_domain_detection_results),
      "Time" = c(DeepST_moni[, 2],
                 seurat_leiden_moni[, 2],
                 seurat_louvain_moni[, 2],
                 BASS_moni[, 2],
                 GraphST_moni[, 2],
                 STAGATE_moni[, 2],
                 SCGDL_moni[, 2],
                 DR.SC_moni[, 2]),
      "Memory" = c(DeepST_moni[, 4],
                   seurat_leiden_moni[, 4],
                   seurat_louvain_moni[, 4],
                   BASS_moni[, 4],
                   GraphST_moni[, 4],
                   STAGATE_moni[, 4],
                   SCGDL_moni[, 4],
                   DR.SC_moni[, 4]),
      "Device" = ifelse(use_cuda,
                        c("cuda", "cpu", "cpu", "cpu", "cuda", "cuda", "cuda", "cpu"),
                        "cpu")
    )
    #-------- Add Clustering Results to Seurat For Visualization --------#
    all_domain_detection_results <- lapply(all_domain_detection_results, FUN = function(x){as.character(x)})
    seurat <- Seurat::AddMetaData(seurat, all_domain_detection_results)
    #-------- Outcome of one simulation method --------#
    list("seurat" = seurat,
         "eval_domain_detection_table" = eval_domain_detection_table,
         "resource_monitering" = resource_monitering)
  })
  names(DomainDetectionResults) <- validated_methods
  return(DomainDetectionResults)
}


.CalculateDomainDetectionMetrics <- function(DomainDetectionResults, trueLabels, method){
  if(!requireNamespace("aricode")){
    message("aricode is not installed on your device")
    message("Installing aricode")
    utils::install.packages("aricode")
  }
  print_color_word(paste("------Evaluation of Domain Detection For", method), color = "blue")
  domain_detection_eval_table <- purrr::map_dfr(1:length(DomainDetectionResults), .f = function(i){
    domain_detection_method <- names(DomainDetectionResults)[i]
    domain_detection_result <- DomainDetectionResults[[domain_detection_method]]
    ### 1. ARI
    ARI <- aricode::ARI(trueLabels, domain_detection_result)
    ### 2. NMI
    NMI <- aricode::NMI(trueLabels, domain_detection_result)
    #### Save results
    tibble::tibble("Simulation_Method" = method,
                   "Clustering_Method" = domain_detection_method,
                   "ARI" = ARI,
                   "NMI" = NMI)
  })
  ### return
  return(domain_detection_eval_table)
}


.STAGATE <- function(adata,
                     seurat,
                     epochs,
                     learning_rate,
                     device,
                     n_domains,
                     seed,
                     verbose){
  if(!requireNamespace("mclust")){
    message("mclust is not installed on your device")
    message("Installing mclust...")
    utils::install.packages("mclust")
  }
  print_color_word(paste("\n \u2192", "STAGATE is running..."), color = "green")
  STAGATE_path <- system.file("STAGATE", package = "SCST")
  stagate <- reticulate::import_from_path("STAGATE_pyG", path = STAGATE_path, convert = FALSE)
  stagate$Cal_Spatial_Net(adata, verbose = reticulate::r_to_py(verbose), model = "KNN", k_cutoff = as.integer(20))
  stagate$Stats_Spatial_Net(adata)
  adata = stagate$train_STAGATE(adata,
                                n_epochs = as.integer(epochs),
                                lr = learning_rate,
                                verbose = reticulate::r_to_py(verbose),
                                random_seed = as.integer(seed),
                                device = device)
  # sc$pp$neighbors(adata, use_rep = 'STAGATE', key_added = "STAGATE_neighbors")
  # adata = stagate$mclust_R(adata, used_obsm = 'STAGATE', num_cluster = n_domains)
  STAGATE_obsm <- reticulate::py_to_r(adata$obsm["STAGATE"])
  require("mclust")
  model_names <- c("EEE", "EEI", "EII", "EVI", "VEI", "VII", "VVI")
  for(i in model_names){
    res <- mclust::Mclust(data = STAGATE_obsm, G = n_domains, modelNames = i)
    if(!is.null(res)){
      print(i)
      break
    }
  }
  STAGATE_result <- res$classification
  print(STAGATE_result)
  names(STAGATE_result) <- colnames(seurat)
  return(STAGATE_result)
}


.DeepST <- function(conda_env,
                    packages,
                    adata,
                    seurat,
                    epochs,
                    use_gpu,
                    n_domains){
  if(!"louvain" %in% packages$package){
    reticulate::py_install("louvain", envname = conda_env, pip = TRUE, ignore_installed = TRUE)
  }
  print_color_word(paste("\n \u2192", "DeepST is running..."), color = "green")
  # new_path <- file.path(DeepST_path, "DeepST.py")
  # reticulate::source_python(new_path, convert = FALSE)
  DeepST_path <- system.file("DeepST", package = "SCST")
  DeepST_function <- reticulate::import_from_path("DeepST", path = DeepST_path, convert = FALSE)
  deepen = DeepST_function$run(task = "Identify_Domain",
                               pre_epochs = as.integer(epochs),
                               epochs = as.integer(epochs),
                               use_gpu = use_gpu)
  adata = deepen$"_get_augment"(adata, use_morphological = reticulate::r_to_py(FALSE))
  graph_dict = deepen$"_get_graph"(adata$obsm["spatial"], distType = "BallTree")
  ##### Enhanced data preprocessing
  processed_data = deepen$"_data_process"(adata, pca_n_comps = as.integer(100))
  ##### Training models
  DeepST = deepen$"_fit"(processed_data, graph_dict)
  ##### DeepST outputs
  adata$obsm$setdefault(key = "DeepST_embed", default = DeepST)
  # adata$obsm["DeepST_embed"] = DeepST
  ##### Define the number of space domains, and the model can also be customized.
  adata = deepen$"_get_cluster_data"(adata, n_domains = n_domains, priori = reticulate::r_to_py(TRUE))
  DeepST_result <- reticulate::py_to_r(adata$obs["DeepST_refine_domain"])
  DeepST_result <- as.character(DeepST_result)
  names(DeepST_result) <- colnames(seurat)
  return(DeepST_result)
}


.seurat_louvain <- function(seurat,
                            n_domains,
                            seed){
  print_color_word(paste("\n \u2192", "Seurat_Louvain is running..."), color = "green")
  resolution <- .find_resolution(seurat,
                                 groups = n_domains,
                                 algorithm = 1,
                                 seed = seed)
  seurat_louvain <- Seurat::FindClusters(object = seurat,
                                         algorithm = 1,
                                         resolution = resolution,
                                         random.seed = seed,
                                         verbose = FALSE)
  seurat_louvain_result <- as.numeric(seurat_louvain@meta.data$seurat_clusters)
  names(seurat_louvain_result) <- colnames(seurat)
  return(seurat_louvain_result)
}


.seurat_leiden <- function(seurat,
                           packages,
                           conda_env,
                           n_domains,
                           seed){
  if(!"leidenalg" %in% packages$package){
    reticulate::py_install("leidenalg", envname = conda_env, pip = TRUE, ignore_installed = TRUE)
  }
  print_color_word(paste("\n \u2192", "Seurat_Leiden is running..."), color = "green")
  resolution <- .find_resolution(seurat,
                                 groups = n_domains,
                                 algorithm = 4,
                                 seed = seed)
  seurat_leiden <- Seurat::FindClusters(object = seurat,
                                        algorithm = 4,
                                        resolution = resolution,
                                        random.seed = seed,
                                        verbose = FALSE)
  seurat_leiden_result <- as.numeric(seurat_leiden@meta.data$seurat_clusters)
  names(seurat_leiden_result) <- colnames(seurat)
  return(seurat_leiden_result)
}


.BASS <- function(X,
                  location,
                  n_domains,
                  nFeatures,
                  PCs,
                  seed){
  if(!requireNamespace("BASS")){
    message("BASS is not installed on your device")
    message("Installing BASS...")
    devtools::install_github("zhengli09/BASS")
  }
  print_color_word(paste("\n \u2192", "BASS is running..."), color = "green")
  #### Create a BASS object
  set.seed(seed)
  require(BASS)
  BASS <- BASS::createBASSObject(list("A" = log2(as.matrix(X) + 1)),
                                 list("A" = as.matrix(location)),
                                 C = n_domains,
                                 R = n_domains,
                                 beta_method = "SW")
  BASS <- BASS::BASS.preprocess(BASS,
                                doLogNormalize = FALSE,
                                nSE = nFeatures,
                                geneSelect = "sparkx",
                                doPCA = TRUE,
                                scaleFeature = FALSE,
                                nPC = PCs,
                                doBatchCorrect = FALSE)
  BASS <- BASS::BASS.run(BASS)
  BASS <- BASS::BASS.postprocess(BASS)
  BASS_result <- BASS@results[["z"]][[1]]
  names(BASS_result) <- colnames(X)
  return(BASS_result)
}


.SCGDL <- function(adata,
                   seurat,
                   epochs,
                   learning_rate,
                   n_domains,
                   seed){
  print_color_word(paste("\n \u2192", "SCGDL is running..."), color = "green")
  SCGDL_auxiliary_path <- system.file("SCGDL", package = "SCST")
  SCGDL_auxiliary <- reticulate::import_from_path("SCGDL_auxiliary", path = SCGDL_auxiliary_path)
  # SCGDL <- reticulate::import_from_path("SCGDL", path = "/Users/duohongrui/Downloads/SCGDL")
  SCGDL <- reticulate::import_from_path("SCGDL", path = SCGDL_auxiliary_path)
  reticulate::source_python(paste0(SCGDL_auxiliary_path, "/SCGDL_Train.py"), convert = FALSE)
  SCGDL_auxiliary$Spatial_Dis_Cal(adata, knn_dis = as.integer(20), model = "KNN")
  adata = SCGDL_Train(adata,
                      lr = learning_rate,
                      random_seed = as.integer(seed),
                      num_epochs = as.integer(epochs))
  knowledge = SCGDL_auxiliary$BayesianGaussianMixture(n_components = n_domains,
                                                      weight_concentration_prior_type = 'dirichlet_process',
                                                      weight_concentration_prior = 50)$fit(adata$obsm["SCGDL"])
  SCGDL_result = knowledge$predict(adata$obsm["SCGDL"]) + 1
  SCGDL_result <- as.character(SCGDL_result)
  # SCGDL_result <- reticulate::py_to_r(SCGDL_result)
  names(SCGDL_result) <- colnames(seurat)
  return(SCGDL_result)
}


.GraphST <- function(conda_env,
                     packages,
                     adata,
                     seurat,
                     n_domains,
                     epochs,
                     learning_rate,
                     seed,
                     device){
  if(!"graphst" %in% packages$package){
    reticulate::py_install("GraphST", envname = conda_env, pip = TRUE, ignore_installed = TRUE)
  }
  if(!"pot" %in% packages$package){
    reticulate::py_install("pot", envname = conda_env, pip = TRUE, ignore_installed = TRUE)
  }
  if(!"scikit-misc" %in% packages$package){
    reticulate::py_install("scikit-misc", envname = conda_env, pip = TRUE, ignore_installed = TRUE)
  }
  # reticulate::source_python(file.path(conda_env, "lib/python3.9/site-packages/GraphST/preprocess.py"), convert = FALSE)
  # reticulate::source_python(file.path(conda_env, "lib/python3.9/site-packages/GraphST/GraphST.py"), convert = FALSE)
  print_color_word(paste("\n \u2192", "GraphST is running..."), color = "green")
  GraphST <- reticulate::import("GraphST", convert = FALSE)
  model = GraphST$GraphST$GraphST(adata,
                                  device = device,
                                  epochs = as.integer(epochs),
                                  learning_rate = as.integer(learning_rate),
                                  random_seed = as.integer(seed))
  adata = model$train()
  GraphST$utils$clustering(adata,
                           n_clusters = n_domains,
                           radius = as.integer(20),
                           method = "louvain",
                           start = 0.01,
                           end = 2.0,
                           increment = 0.01,
                           refinement = reticulate::r_to_py(FALSE))
  GraphST_result = adata$obs["domain"]
  GraphST_result = GraphST_result$astype("object")
  GraphST_result <- reticulate::py_to_r(GraphST_result) %>% as.character()
  names(GraphST_result) <- colnames(seurat)
  return(GraphST_result)
}


.DR.SC <- function(seurat,
                   n_domains,
                   verbose){
  if(!requireNamespace("DR.SC")){
    message("DR.SC is not installed on your device")
    message("Installing DR.SC...")
    utils::install.packages("DR.SC")
  }
  print_color_word(paste("\n \u2192", "DR.SC is running..."), color = "green")
  DR.SC_result <- DR.SC::DR.SC(seurat,
                               K = n_domains,
                               platform = "Visium",
                               verbose = verbose)
  DR.SC_result <- DR.SC_result$spatial.drsc.cluster
  return(DR.SC_result)
}


#### Environment Configuration
#' conda create -n spatial python=3.9
#' conda activate spatial
#' git clone https://github.com/spatial-Transcriptomics/DeepST.git
#' cd DeepST
#'
#' #### CPU
#' pip install torch==1.13.0 torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cpu
#' pip torch_scatter torch_sparse torch_cluster torch_spline_conv torch_geometric -f https://data.pyg.org/whl/torch-1.13.0+cpu.html
#' pip install -r requirements.txt
#' pip install louvain
#'
#' #### CUDA
#'
#'
#' https://data.pyg.org/whl/index.html
#' Driver https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html
#'
#' https://www.cnblogs.com/Wanggcong/p/12625540/html
#'
#' torch_clusters/torch_scatter/torch_sparse
#'
#'
#' pip install pot
#' pip install scikit-misc
#' pip install GraphST
#'
#' git clone https://github.com/QIFEIDKN/STAGATE_pyG.git
#' cd STAGATE_pyG
#' python setup.py build
#' python setup.py install
#'
