#' #' Evaluation of Clustering Methods Based on Simulated Data
#' #' Evaluate Cell Association Measurements Using Simulated Data
#' #'
#' #' @param SimulationResult A SCST object containing the simulated data.
#' #' @param n_cores The number of computational cores used in the clustering methods. Default is 1.
#' #' @param verbose Whether to return messages or not. Default is TRUE.
#' #'
#' #' @export
#' #'
#' #'
#' #' @references Skinnider M A, Squair J W, Foster L J. Evaluating measures of association for single-cell transcriptomics[J]. Nature methods, 2019, 16(5): 381-386.
#' #'
#' EvalAssociationMethods <- function(SimulationResult,
#'                                    n_cores = 1,
#'                                    verbose = TRUE){
#'   #-------- Check Necessary Information --------#
#'   if(S4Vectors::isEmpty(SimulationResult@simulation_result)){
#'     print_color_word("no No simulated data is found.")
#'     stop(error_output())
#'   }
#'   validated_methods <- check_info_pre_application(SimulationResult = SimulationResult, Group = TRUE)
#'
#'   ### Iterate all simulation methods and perform clustering and sequent evaluation
#'   methods_names <- names(SimulationResult@simulation_result)
#'   AssociationResults <- purrr::map(methods_names, .f = function(methods_names){
#'     #-------- Methods --------#
#'     print_color_word(paste("------Measuring Cell Association On", methods_names), color = "blue")
#'     data <- SimulationResult@simulation_result[[methods_names]][["count_data"]]
#'     col_data <- SimulationResult@simulation_result[[methods_names]][["col_meta"]]
#'     ngroup <- length(unique(col_data[, "group"]))
#'     all_association_results <- list()
#'     ### 1. Pearson Correlation
#'     if(!requireNamespace("WGCNA")){
#'       message("WGCNA is not installed on your device")
#'       message("Installing WGCNA...")
#'       BiocManager::install("WGCNA")
#'     }
#'     all_association_results[["PearsonCorrelation"]] <- WGCNA::cor(data, method = "pearson", nThreads = n_cores)
#'
#'     ### 2. Spearman correlation
#'     all_association_results[["SpearmanCorrelation"]] <- WGCNA::cor(data, method = "spearman", nThreads = n_cores)
#'
#'     ### 3. Kendall correlation
#'     if(!requireNamespace("pcaPP")){
#'       message("pcaPP is not installed on your device")
#'       message("Installing pcaPP...")
#'       utils::install.packages("pcaPP")
#'     }
#'     all_association_results[["KendallCorrelation"]] <- pcaPP::cor.fk(data)
#'
#'     ### 4. Weighted Rank Correlation
#'     ranks <- apply(data, 2, rank, ties = "average")
#'     # weight the ranks
#'     # calculate the savage scores
#'     n <- nrow(data)
#'     reciprocals <- 1 / seq_len(n)
#'     savage <- sapply(seq_len(n), function(i){sum(reciprocals[i:n])})
#'     # replace each rank with the savage score
#'     savages <- ranks
#'     savages[] <- savage[ranks]
#'     # calculate pearson correlation
#'     all_association_results[["WeightedRankCorrelation"]] <- WGCNA::cor(savages, method = "pearson", nThreads = n_cores)
#'
#'     ### 5. Mutual information
#'     all_association_results[["MutualInfo"]] <- WGCNA::mutualInfoAdjacency(data)$AdjacencyUniversalVersion1
#'
#'     ### 6. Cosine distance
#'     all_association_results[["CosineCorrelation"]] <- WGCNA::cor(data, method = "pearson", cosine = TRUE, nThreads = n_cores)
#'
#'     ### 7. Canberra distance
#'     if(!requireNamespace("parallelDist")){
#'       message("parallelDist is not installed on your device")
#'       message("Installing parallelDist...")
#'       utils::install.packages("parallelDist")
#'     }
#'     all_association_results[["CanberraDistance"]] <- -1.0 * parallelDist::parallelDist(t(data),
#'                                                                                        method = "canberra",
#'                                                                                        threads = n_cores) %>%
#'       as.matrix()
#'
#'
#'     ### 8. Euclidean distance
#'     all_association_results[["EuclideanDistance"]] <- -1.0 * parallelDist::parallelDist(t(data),
#'                                                                                         method = "euclidean",
#'                                                                                         threads = n_cores) %>%
#'       as.matrix()
#'
#'
#'     ### 9. Manhattan distance
#'     all_association_results[["ManhattanDistance"]] <- -1.0 * parallelDist::parallelDist(t(data),
#'                                                                                         method = "manhattan",
#'                                                                                         threads = n_cores) %>%
#'       as.matrix()
#'
#'     ### 10. Hamming distance
#'     all_association_results[["HammingDistance"]] <- -1.0 * parallelDist::parallelDist(t(data),
#'                                                                                       method = "hamming",
#'                                                                                       threads = n_cores) %>%
#'       as.matrix()
#'
#'     ### 11. Phi S
#'     if(!requireNamespace("propr")){
#'       message("propr is not installed on your device")
#'       message("Installing propr...")
#'       devtools::install_github("tpq/propr")
#'     }
#'     all_association_results[["Phi_S"]] <- -1.0 * propr::propr(counts = data, metric = "phs")@matrix
#'
#'     ### 12. Rho
#'     all_association_results[["Rho"]] <- propr::propr(counts = data, metric = "rho")@matrix
#'
#'     #-------- Evaluate Clustering Results --------#
#'     eval_association_table <- .CalculateAssociationMetrics(
#'       all_association_results,
#'       methods_names
#'     )
#'     #-------- Add Clustering Results to Seurat For Visualization --------#
#'     all_clustering_results <- lapply(all_clustering_results, FUN = function(x){as.character(x)})
#'     seurat <- Seurat::AddMetaData(seurat, all_clustering_results)
#'     #-------- Outcome of one simulation method --------#
#'     list("seurat" = seurat,
#'          "eval_clustering_table" = eval_clustering_table)
#'   })
#'   names(ClusteringResults) <- methods_names
#' }
#'
#'
#' .CalculateAssociationMetrics <- function(
#'     all_association_results,
#'     methods_names,
#'     ngroup
#' ){
#'   association_methods <- names(all_association_results)
#'
#'   hclust_result <- lapply(all_association_results, FUN = function(per_ass_result){
#'     clust <- hclust(as.dist(per_ass_result-1))
#'     clusters <- cutree(clust, k = ngroup)
#'     aricode::ARI(clusters, as.character(col_data[, "group"]))
#'   })
#'
#' }
#'
#'
#'
#'
#'
