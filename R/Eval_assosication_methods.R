#' Evaluation of Clustering Methods Based on Simulated Data
#' Evaluate Cell Association Measurements Using Simulated Data
#'
#' @param SimulationResult A SCST object containing the simulated data.
#' @param n_cores The number of computational cores used in the clustering methods. Default is 1.
#' @param verbose Whether to return messages or not. Default is TRUE.
#'
#' @importFrom stats hclust cutree as.dist
#'
#' @export
#'
#'
#' @references Skinnider M A, Squair J W, Foster L J. Evaluating measures of association for single-cell transcriptomics[J]. Nature methods, 2019, 16(5): 381-386.
#'
EvalAssociationMethods <- function(SimulationResult,
                                   n_cores = 1,
                                   verbose = TRUE){
  #-------- Check Necessary Information --------#
  if(length(SimulationResult@simulation_result) == 0){
    print_color_word("no No simulated data is found.")
    stop(error_output())
  }
  validated_methods <- check_info_pre_application(SimulationResult = SimulationResult, Group = TRUE)
  count_data <- get_count_data(SimulationResult)
  col_data <- get_cell_meta(SimulationResult)
  ### Iterate all simulation methods and perform clustering and sequent evaluation
  AssociationResults <- purrr::map_dfr(validated_methods, .f = function(method){
    #-------- Methods --------#
    print_color_word(paste("------Measuring Cell Association On", method), color = "blue")
    data <- count_data[[method]]
    cell_meta <- col_data[[method]]
    ngroup <- length(unique(cell_meta[, "group"]))
    if(!is.matrix(data)){
      data <- as.matrix(data)
    }
    all_association_results <- list()
    ### 1. Pearson Correlation
    if(!requireNamespace("WGCNA")){
      message("WGCNA is not installed on your device")
      message("Installing WGCNA...")
      BiocManager::install("WGCNA")
    }
    print_color_word(paste("\u2192", "Calculating Pearson Correlation..."), color = "green")
    PearsonCorrelation_moni <- peakRAM::peakRAM(
      all_association_results[["PearsonCorrelation"]] <- WGCNA::cor(data, method = "pearson", nThreads = n_cores)
    )

    ### 2. Spearman correlation
    print_color_word(paste("\u2192", "Calculating Spearman Correlation..."), color = "green")
    SpearmanCorrelation_moni <- peakRAM::peakRAM(
      all_association_results[["SpearmanCorrelation"]] <- WGCNA::cor(data, method = "spearman", nThreads = n_cores)
    )

    ### 3. Kendall correlation
    if(!requireNamespace("pcaPP")){
      message("pcaPP is not installed on your device")
      message("Installing pcaPP...")
      utils::install.packages("pcaPP")
    }
    print_color_word(paste("\u2192", "Calculating Kendall Correlation..."), color = "green")
    KendallCorrelation_moni <- peakRAM::peakRAM(
      all_association_results[["KendallCorrelation"]] <- pcaPP::cor.fk(data)
    )

    ### 4. Weighted Rank Correlation
    print_color_word(paste("\u2192", "Calculating Weighted Rank Correlation..."), color = "green")
    WeightedRankCorrelation_moni <- peakRAM::peakRAM(
      all_association_results <- .WeightedRankCorrelation(all_association_results = all_association_results,
                                                          data = data,
                                                          n_cores = n_cores)
    )

    ### 5. Mutual information
    print_color_word(paste("\u2192", "Calculating Mutual Information..."), color = "green")
    MutualInfo_moni <- peakRAM::peakRAM(
      all_association_results[["MutualInfo"]] <- WGCNA::mutualInfoAdjacency(data)$AdjacencyUniversalVersion1
    )

    ### 6. Cosine distance
    print_color_word(paste("\u2192", "Calculating Cosine Distance..."), color = "green")
    CosineCorrelation_moni <- peakRAM::peakRAM(
      all_association_results[["CosineCorrelation"]] <- WGCNA::cor(data, method = "pearson", cosine = TRUE, nThreads = n_cores)
    )

    ### 7. Canberra distance
    if(!requireNamespace("parallelDist")){
      message("parallelDist is not installed on your device")
      message("Installing parallelDist...")
      utils::install.packages("parallelDist")
    }
    print_color_word(paste("\u2192", "Calculating Canberra Distance..."), color = "green")
    CanberraDistance_moni <- peakRAM::peakRAM(
      all_association_results[["CanberraDistance"]] <- -1.0 * parallelDist::parallelDist(t(data),
                                                                                         method = "canberra",
                                                                                         threads = n_cores) %>%
        as.matrix()
    )

    ### 8. Euclidean distance
    print_color_word(paste("\u2192", "Calculating Euclidean Distance..."), color = "green")
    EuclideanDistance_moni <- peakRAM::peakRAM(
      all_association_results[["EuclideanDistance"]] <- -1.0 * parallelDist::parallelDist(t(data),
                                                                                          method = "euclidean",
                                                                                          threads = n_cores) %>%
        as.matrix()
    )

    ### 9. Manhattan distance
    print_color_word(paste("\u2192", "Calculating Manhattan Distance..."), color = "green")
    ManhattanDistance_moni <- peakRAM::peakRAM(
      all_association_results[["ManhattanDistance"]] <- -1.0 * parallelDist::parallelDist(t(data),
                                                                                          method = "manhattan",
                                                                                          threads = n_cores) %>%
        as.matrix()
    )

    ### 10. Hamming distance
    print_color_word(paste("\u2192", "Calculating Hamming Distance..."), color = "green")
    HammingDistance_moni <- peakRAM::peakRAM(
      all_association_results[["HammingDistance"]] <- -1.0 * parallelDist::parallelDist(t(data),
                                                                                        method = "hamming",
                                                                                        threads = n_cores) %>%
        as.matrix()
    )

    ### 11. Phi S
    if(!requireNamespace("propr")){
      message("propr is not installed on your device")
      message("Installing propr...")
      devtools::install_github("tpq/propr")
    }
    print_color_word(paste("\u2192", "Calculating Phi..."), color = "green")
    Phi_S_moni <- peakRAM::peakRAM(
      all_association_results[["Phi_S"]] <- -1.0 * propr::propr(counts = data, metric = "phs")@matrix
    )

    ### 12. Rho
    print_color_word(paste("\u2192", "Calculating Rho..."), color = "green")
    Rho_moni <- peakRAM::peakRAM(
      all_association_results[["Rho"]] <- propr::propr(counts = data, metric = "rho")@matrix
    )

    #-------- Evaluate Cell Association Results --------#
    trueLabels <- cell_meta[, "group"] %>% as.character()
    eval_association_table <- .CalculateAssociationMetrics(
      all_association_results,
      method,
      trueLabels,
      ngroup
    )
    #-------- Record Resource Occupation During Execution --------#
    resource_monitering <- tibble::tibble(
      "Simulation_Method" = method,
      "Cell_Association_Method" = names(all_association_results),
      "Time" = c(PearsonCorrelation_moni[, 2],
                 SpearmanCorrelation_moni[, 2],
                 KendallCorrelation_moni[, 2],
                 WeightedRankCorrelation_moni[, 2],
                 MutualInfo_moni[, 2],
                 CosineCorrelation_moni[, 2],
                 CanberraDistance_moni[, 2],
                 EuclideanDistance_moni[, 2],
                 ManhattanDistance_moni[, 2],
                 HammingDistance_moni[, 2],
                 Phi_S_moni[, 2],
                 Rho_moni[, 2]),
      "Memory" = c(PearsonCorrelation_moni[, 4],
                   SpearmanCorrelation_moni[, 4],
                   KendallCorrelation_moni[, 4],
                   WeightedRankCorrelation_moni[, 4],
                   MutualInfo_moni[, 4],
                   CosineCorrelation_moni[, 4],
                   CanberraDistance_moni[, 4],
                   EuclideanDistance_moni[, 4],
                   ManhattanDistance_moni[, 4],
                   HammingDistance_moni[, 4],
                   Phi_S_moni[, 4],
                   Rho_moni[, 4]),
      "Device" = "cpu"
    )
    eval_association_table <- eval_association_table %>%
      dplyr::full_join(resource_monitering, by = c("Simulation_Method", "Cell_Association_Method"))
    eval_association_table
  })
  return(AssociationResults)
}


.CalculateAssociationMetrics <- function(
    all_association_results,
    method,
    trueLabels,
    ngroup
){
  if(!requireNamespace("scales")){
    message("scales is not installed on your device")
    message("Installing scales...")
    utils::install.packages("scales")
  }
  print_color_word(paste("------Evaluation of Domain Detection For", method), color = "blue")
  association_methods <- names(all_association_results)
  cell_association_eval_table <- purrr::map_dfr(association_methods, .f = function(associate_method){
    cor <- all_association_results[[associate_method]]
    if(associate_method == "MutualInfo" |
       associate_method == "CosineCorrelation" |
       associate_method == "CanberraDistance" |
       associate_method == "EuclideanDistance" |
       associate_method == "ManhattanDistance" |
       associate_method == "HammingDistance" |
       associate_method == "Phi_S"){
      cor <- scales::rescale(cor)
    }
    clust <- stats::hclust(stats::as.dist(-cor))
    clusters <- stats::cutree(clust, k = ngroup)
    ARI <- aricode::ARI(clusters, trueLabels)
    NMI <- aricode::NMI(clusters, trueLabels)
    tibble::tibble("Simulation_Method" = method,
                   "Cell_Association_Method" = associate_method,
                   "ARI" = ARI,
                   "NMI" = NMI)
  })
  return(cell_association_eval_table)
}


.WeightedRankCorrelation <- function(all_association_results,
                                     data,
                                     n_cores){
  ranks <- apply(data, 2, rank, ties = "average")
  # weight the ranks
  # calculate the savage scores
  n <- nrow(data)
  reciprocals <- 1 / seq_len(n)
  savage <- sapply(seq_len(n), function(i){sum(reciprocals[i:n])})
  # replace each rank with the savage score
  savages <- ranks
  savages[] <- savage[ranks]
  # calculate pearson correlation
  all_association_results[["WeightedRankCorrelation"]] <- WGCNA::cor(savages, method = "pearson", nThreads = n_cores)
  return(all_association_results)
}

