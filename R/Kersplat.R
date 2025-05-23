#' Estimation Parameters Using Kersplat
#'
#' @param SCST_Object A SCST object generated by [SCST::Create_SCST_Object()]
#' @param verbose Whether to return messages or not
#'
#' @export
#'

Kersplat_estimation <- function(SCST_Object,
                                verbose = FALSE
){
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  ref_data <- as.matrix(SCST_Object@reference)
  seed <- SCST_Object@seed
  ##############################################################################
  ####                            Estimation                                 ###
  ##############################################################################
  if(verbose){
    message("Estimating parameters using Kersplat")
  }
  # Seed
  set.seed(seed)
  # Estimation
  estimate_detection <- peakRAM::peakRAM(
    estimate_result <- splatter::kersplatEstimate(ref_data, verbose = FALSE)
  )
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  slot(SCST_Object, "estimation_time_memory") <- list("Kersplat" = list("estimation_time" = estimate_detection$Elapsed_Time_sec,
                                                                        "estimation_memory" = estimate_detection$Peak_RAM_Used_MiB))
  slot(SCST_Object, "estimation_result") <- list("Kersplat" = estimate_result)
  slot(SCST_Object, "estimation_parameters") <- list("Kersplat" = NULL)
  return(SCST_Object)
}


#' Simulate ScRNA-seq Data Using Kersplat
#'
#' @param estimated_result The SCST object after estimating parameters using [SCST::Kersplat_estimation()]
#' @param cell_num The expected number of simulated cells. If NULL, the original cell number in the reference data will be adopted.
#' @param gene_num The expected number of simulated genes. If NULL, the original gene number in the reference data will be adopted.
#' @param verbose Whether to return messages or not
#' @param seed Random seed
#' @param return_format The format of returned simulation data. Choices: list, SingleCellExperiment and Seurat.
#' @param ... Other parameters represented in Kersplat, see [splatter::KersplatParams()]
#'
#' @export
#'
Kersplat_simulation <- function(estimated_result,
                                cell_num = NULL,
                                gene_num = NULL,
                                return_format,
                                verbose = FALSE,
                                seed,
                                ...
){
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  parameters <- estimated_result@estimation_result$Kersplat
  if(!is.null(cell_num)){
    parameters <- splatter::setParam(parameters, "nCells", cell_num)
  }
  if(!is.null(gene_num)){
    parameters <- splatter::setParam(parameters, "nGenes", gene_num)
  }
  # Seed
  parameters <- splatter::setParam(parameters, name = "seed", value = seed)
  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    message("Simulating datasets using Kersplat")
  }
  # Simulation
  simulate_detection <- peakRAM::peakRAM(
    simulate_result <- splatter::kersplatSimulate(parameters,
                                                  verbose = verbose)
  )
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  ## counts
  counts <- as.matrix(SingleCellExperiment::counts(simulate_result))
  ## col_data
  col_data <- data.frame("cell_name" = colnames(counts),
                         row.names = colnames(counts))
  ## row_data
  row_data <- as.data.frame(SingleCellExperiment::rowData(simulate_result)[, 1])
  rownames(row_data) <- row_data[, 1]
  colnames(row_data) <- "gene_name"
  # Establish SingleCellExperiment
  simulate_result <- SingleCellExperiment::SingleCellExperiment(list(counts = counts),
                                                                colData = col_data,
                                                                rowData = row_data)
  simulate_result <- data_conversion(SCE_object = simulate_result, return_format = return_format)
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  slot(estimated_result, "simulation_time_memory") <- list("Kersplat" = list("estimation_time" = simulate_detection$Elapsed_Time_sec,
                                                                             "estimation_memory" = simulate_detection$Peak_RAM_Used_MiB))
  slot(estimated_result, "simulation_result") <- list("Kersplat" = simulate_result)
  slot(estimated_result, "simulation_parameters") <- list("Kersplat" = parameters)
  return(estimated_result)
}
