#' Estimation Parameters Using Lun
#'
#' @param SCST_Object A SCST object generated by [SCST::Create_SCST_Object()]
#' @param verbose Whether to return messages or not
#'
#' @export
#'
Lun_estimation <- function(SCST_Object,
                           verbose = FALSE
){
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  ref_data <- SCST_Object@reference %>% as.matrix()
  seed <- SCST_Object@seed
  ##############################################################################
  ####                            Estimation                                 ###
  ##############################################################################
  if(verbose){
    message("Estimating parameters using Lun")
  }
  # Seed
  set.seed(seed)
  # Estimation
  estimate_detection <- peakRAM::peakRAM(
    estimate_result <- splatter::lunEstimate(ref_data)
  )
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  slot(SCST_Object, "estimation_time_memory") <- list("Lun" = list("estimation_time" = estimate_detection$Elapsed_Time_sec,
                                                                   "estimation_memory" = estimate_detection$Peak_RAM_Used_MiB))
  slot(SCST_Object, "estimation_result") <- list("Lun" = estimate_result)
  slot(SCST_Object, "estimation_parameters") <- list("Lun" = NULL)
  return(SCST_Object)
}


#' Simulate ScRNA-seq Data Using Lun
#'
#' @param estimated_result The SCST object after estimating parameters using [SCST::Lun_estimation()]
#' @param cell_num The expected number of simulated cells. If NULL, the original cell number in the reference data will be adopted.
#' @param nGroups The expected number of groups in the simulated data. Default is 1.
#' @param prop_group The proportion of cells in each group. Default is 1.
#' @param de_prop The proportion of DEGs over all genes between different simulated cell groups. Default is 0.2.
#' @param fc_group The fold change of DEGs between cell groups. Default is 2.
#' @param verbose Whether to return messages or not
#' @param seed Random seed
#' @param return_format The format of returned simulation data. Choices: list, SingleCellExperiment and Seurat.
#' @param ... Other parameters represented in Lun, see [splatter::LunParams()]
#'
#' @export
#'
Lun_simulation <- function(estimated_result,
                           cell_num = NULL,
                           nGroups = 1,
                           prop_group = 1,
                           de_prop = 0.2,
                           fc_group = 2,
                           return_format,
                           verbose = FALSE,
                           seed,
                           ...
){
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  parameters <- estimated_result@estimation_result$Lun
  if(is.null(cell_num)){
    cell_num <- estimated_result@cell_num
  }
  # prop_group
  if(nGroups > 1){
    if(nGroups == length(prop_group)){
      cell_num <- proportionate(
        number = cell_num,
        result_sum_strict = cell_num,
        prop = prop_group,
        prop_sum_strict = 1,
        digits = 0)
    }else{
      cell_num <- proportionate(
        number = cell_num,
        result_sum_strict = cell_num,
        prop = rep(1/nGroups, nGroups),
        digits = 0)
    }
    parameters <- splatter::setParam(parameters,
                                     name = "groupCells",
                                     value = cell_num)
  }
  # de_prop
  parameters <- splatter::setParam(parameters,
                                   name = "de.nGenes",
                                   value = round(de_prop * estimated_result@gene_num/nGroups))
  # fc.up.group
  parameters <- splatter::setParam(parameters,
                                   name = "de.upFC",
                                   value = fc_group)
  # fc.down.group
  parameters <- splatter::setParam(parameters,
                                   name = "de.upFC",
                                   value = 1/fc_group)
  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    message("Simulating datasets using Lun")
  }
  # Seed
  parameters <- splatter::setParam(parameters, name = "seed", value = seed)

  # Simulation
  simulate_detection <- peakRAM::peakRAM(
    simulate_result <- splatter::lunSimulate(parameters,
                                             verbose = verbose)
  )
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  # counts
  counts <- as.matrix(SingleCellExperiment::counts(simulate_result))
  # col_data
  col_data <- as.data.frame(SummarizedExperiment::colData(simulate_result))
  if(nGroups == 1){
    col_data[, 2] <- rep("Group1", ncol(col_data))
    colnames(col_data) <- c("cell_name", "group")
  }else{
    col_data <- col_data[, c("Cell", "Group")]
    colnames(col_data) <- c("cell_name", "group")
  }

  # row_data
  row_data <- as.data.frame(SummarizedExperiment::rowData(simulate_result))
  if(nGroups == 1){
    row_data <- data.frame("gene_name" = row_data$Gene,
                           row.names = rownames(counts))
  }else{
    group_fac <- row_data[, grep(colnames(row_data), pattern = "^DEFac")]
    total_sum <- rowSums(group_fac)
    de_gene <- ifelse(total_sum == nGroups, "no", "yes")
    row_data[, 2] <- de_gene
    row_data <- row_data[, 1:2]
    row_data <- BiocGenerics::cbind(row_data, group_fac)
    colnames(row_data) <- c("gene_name", "de_gene", colnames(group_fac))
  }
  # Establish SingleCellExperiment
  simulate_result <- SingleCellExperiment::SingleCellExperiment(list(counts = counts),
                                                                colData = col_data,
                                                                rowData = row_data)
  simulate_result <- data_conversion(SCE_object = simulate_result, return_format = return_format)
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  slot(estimated_result, "simulation_time_memory") <- list("Lun" = list("estimation_time" = simulate_detection$Elapsed_Time_sec,
                                                                        "estimation_memory" = simulate_detection$Peak_RAM_Used_MiB))
  slot(estimated_result, "simulation_result") <- list("Lun" = simulate_result)
  slot(estimated_result, "simulation_parameters") <- list("Lun" = parameters)
  return(estimated_result)
}
