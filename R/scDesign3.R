#' Estimation Parameters Using scDesign3
#'
#' @param SCST_Object A SCST object generated by [SCST::Create_SCST_Object()]
#' @param family_use A string or a vector of strings of the marginal distribution. Must be one of "poisson", "nb", "zip", "zinb" or "gaussian". Default is "zinb".
#' @param n_cores An integer. The number of cores to use.
#' @param copula A string of the copula choice. Must be one of "gaussian" or "vine". Default is "gaussian".
#' @param path A Boolean symbol (TRUE or FLASE) which determines whether to simulate datasets with trajectory along the path. Default is false.
#' @param verbose Whether to return messages or not
#' @param ... Other parameters represented in scDesign3, see [scDesign3::fit_marginal()], [scDesign3::fit_copula()], [scDesign3::extract_para()]
#'
#' @export
#'
#' @references Song, D., Wang, Q., Yan, G. et al. scDesign3 generates realistic in silico data for multimodal single-cell and spatial omics. Nat Biotechnol (2023). https://doi.org/10.1038/s41587-023-01772-1
#'
scDesign3_estimation <- function(SCST_Object,
                                 family_use = "zinb",
                                 n_cores = 1,
                                 copula = "gaussian",
                                 path = FALSE,
                                 verbose = FALSE,
                                 ...
){
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  if(!requireNamespace("scDesign3", quietly = TRUE)){
    message("scDesign3 is not installed on your device...")
    message("Installing scDesign3...")
    devtools::install_github("SONGDONGYUAN1994/scDesign3")
  }
  ref_data <- as.matrix(SCST_Object@reference)
  seed <- SCST_Object@seed
  ### group
  if(!S4Vectors::isEmpty(SCST_Object@group_label)){
    col_data <- data.frame("group_label" = paste0("Group", as.numeric(factor(SCST_Object@group_label))))
    print("Group...")
  }else{
    col_data <- data.frame("group_label" = rep("A", ncol(ref_data)))
  }
  rownames(col_data) <- colnames(ref_data)
  col_data$"cell_name" <- colnames(ref_data)
  ### batch
  if(!S4Vectors::isEmpty(SCST_Object@batch_label)){
    col_data$batch_label <- paste0("Batch", as.numeric(factor(SCST_Object@batch_label)))
    print("Batch...")
  }
  if(!S4Vectors::isEmpty(SCST_Object@spatial_x) & !S4Vectors::isEmpty(SCST_Object@spatial_y)){
    col_data$spatial.x <- SCST_Object@spatial_x
    col_data$spatial.y <- SCST_Object@spatial_y
    spatial <- c("spatial.x", "spatial.y")
  }else{
    spatial <- NULL
  }
  ### trajectory
  if(path){
    if(!requireNamespace("dyndimred", quietly = TRUE)){
      devtools::install_github("dynverse/dyndimred")
    }
    cat("Constructing lineages for the data...\n")
    traj_info <- .pseudotime_info(SCST_Object = SCST_Object,
                                  dynwrap_data = NULL,
                                  start_cell_id = NULL,
                                  col_data = col_data,
                                  seed = seed)
    col_data <- traj_info$col_data
    mu_formula <- traj_info$mu_formula
    pseudotime <- traj_info$pseudotime
  }else{
    pseudotime <- NULL
  }

  # Establish SingleCellExperiment
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = ref_data), colData = col_data)
  ### step by step
  #### other_covariates
  other_covariates <- NULL
  if(!S4Vectors::isEmpty(SCST_Object@batch_label)){
    other_covariates <- "batch_label"
  }
  scDesign3_data <- scDesign3::construct_data(
    sce = sce,
    assay_use = "counts",
    celltype = "group_label",
    pseudotime = pseudotime,
    spatial = spatial,
    other_covariates = other_covariates,
    corr_by = "group_label",
    ncell = SCST_Object@cell_num
  )

  ### mu_formula
  if(!exists("mu_formula")){
    mu_formula <- NULL
    for(i in c("group_label", "batch_label")){
      if(!S4Vectors::isEmpty(methods::slot(SCST_Object, i))){
        if(is.null(mu_formula)){
          mu_formula <- append(mu_formula, i)
        }else{
          mu_formula <- paste0(mu_formula, " + ", i)
        }
      }
    }
    if(!S4Vectors::isEmpty(SCST_Object@spatial_x) & !S4Vectors::isEmpty(SCST_Object@spatial_y)){
      mu_formula <- "s(spatial.x, spatial.y, bs = 'gp')"
    }
    if(path){
      mu_formula <- traj_info$mu_formula
    }
  }
  scDesign3_est <- function(){
    scDesign3_marginal <- scDesign3::fit_marginal(
      data = scDesign3_data,
      predictor = ifelse(exists("predictor"), get("predictor"), "gene"),
      mu_formula = ifelse(!is.null(mu_formula), mu_formula, "1"),
      sigma_formula = ifelse(exists("sigma_formula"), get("sigma_formula"), "1"),
      family_use = family_use,
      n_cores = n_cores,
      usebam = FALSE
    )
    set.seed(seed)
    scDesign3_copula <- scDesign3::fit_copula(
      sce = sce,
      assay_use = "counts",
      marginal_list = scDesign3_marginal,
      family_use = family_use,
      copula = copula,
      n_cores = n_cores,
      input_data = scDesign3_data$dat
    )
    scDesign3_para <- scDesign3::extract_para(
      sce = sce,
      marginal_list = scDesign3_marginal,
      n_cores = 1,
      family_use = family_use,
      new_covariate = NULL,
      data = scDesign3_data$dat
    )
    estimate_result <- list(scDesign3_data = scDesign3_data,
                            scDesign3_marginal = scDesign3_marginal,
                            scDesign3_copula = scDesign3_copula,
                            scDesign3_para = scDesign3_para)
    return(estimate_result)

  }
  ##############################################################################
  ####                            Estimation                                 ###
  ##############################################################################
  if(verbose){
    message("Estimating parameters using scDesign3")
  }
  # Seed
  set.seed(seed)
  # Estimation
  estimate_detection <- peakRAM::peakRAM(
    estimate_result <- scDesign3_est()
  )
  estimate_result <- list(sce = sce,
                          scDesign3_data = estimate_result$scDesign3_data,
                          scDesign3_marginal = estimate_result$scDesign3_marginal,
                          scDesign3_copula = estimate_result$scDesign3_copula,
                          scDesign3_para = estimate_result$scDesign3_para)
  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  slot(SCST_Object, "estimation_time_memory") <- list("scDesign3" = list("estimation_time" = estimate_detection$Elapsed_Time_sec,
                                                                         "estimation_memory" = estimate_detection$Peak_RAM_Used_MiB))
  slot(SCST_Object, "estimation_result") <- list("scDesign3" = estimate_result)
  slot(SCST_Object, "estimation_parameters") <- list("scDesign3" = NULL)
  return(SCST_Object)
}



#' Simulate ScRNA-seq Data Using scDesign3
#'
#' @param estimated_result The SCST object after estimating parameters using [scDesign3::fit_marginal()], [scDesign3::fit_copula()], [scDesign3::extract_para()]
#' @param n_cores An integer. The number of cores to use.
#' @param family_use A string or a vector of strings of the marginal distribution. Must be one of "poisson", "nb", "zip", "zinb" or "gaussian". Default is "zinb".
#' @param verbose Whether to return messages or not
#' @param seed Random seed
#' @param return_format The format of returned simulation data. Choices: list, SingleCellExperiment and Seurat.
#' @param ... Other parameters represented in scDesign3, see [scDesign3::simu_new()]
#'
#' @export
#'
#' @references Song, D., Wang, Q., Yan, G. et al. scDesign3 generates realistic in silico data for multimodal single-cell and spatial omics. Nat Biotechnol (2023). https://doi.org/10.1038/s41587-023-01772-1
#'
scDesign3_simulation <- function(estimated_result,
                                 n_cores = 1,
                                 family_use = "zinb",
                                 return_format,
                                 verbose = FALSE,
                                 seed,
                                 ...
){
  ##############################################################################
  ####                            Environment                                ###
  ##############################################################################
  if(!requireNamespace("scDesign3", quietly = TRUE)){
    message("scDesign3 is not installed on your device...")
    message("Installing scDesign3...")
    devtools::install_github("SONGDONGYUAN1994/scDesign3")
  }
  ##############################################################################
  ####                               Check                                   ###
  ##############################################################################
  parameters <- estimated_result@estimation_result$scDesign3
  existed_params <- list(sce = parameters[["sce"]],
                         mean_mat =  parameters[["scDesign3_para"]][["mean_mat"]],
                         sigma_mat = parameters[["scDesign3_para"]][["sigma_mat"]],
                         zero_mat = parameters[["scDesign3_para"]][["zero_mat"]],
                         copula_list = parameters[["scDesign3_copula"]][["copula_list"]],
                         n_cores = n_cores,
                         family_use = family_use,
                         input_data = parameters[["scDesign3_data"]][["dat"]],
                         new_covariate = parameters[["scDesign3_data"]][["newCovariate"]],
                         important_feature = parameters[["scDesign3_copula"]][["important_feature"]],
                         filtered_gene = NULL)
  extra_params <- list(...)
  existed_params <- append(existed_params, extra_params)
  argument <- formals(scDesign3::simu_new)
  used_argument <- change_parameters(existed_params, argument)
  used_argument <- append(used_argument, list(filtered_gene = NULL))
  if(!"new_covariate" %in% names(used_argument)){
    used_argument <- append(used_argument, list(new_covariate = NULL))
  }
  ##############################################################################
  ####                            Simulation                                 ###
  ##############################################################################
  if(verbose){
    message("Simulating datasets using scDesign3")
  }
  # Seed
  set.seed(seed)
  # Estimation
  simulate_detection <- peakRAM::peakRAM(
    simulate_result <- do.call(scDesign3::simu_new, used_argument)
  )
  ##############################################################################
  ####                        Format Conversion                              ###
  ##############################################################################
  counts <- as.matrix(simulate_result)
  col_data <- data.frame("cell_name" = colnames(counts))
  ### spatial information
  if(!is.null(existed_params[["new_covariate"]][["spatial.x"]])){
    col_data$"spatial.x" <- existed_params[["new_covariate"]][["spatial.x"]]
    col_data$"spatial.y" <- existed_params[["new_covariate"]][["spatial.y"]]
  }
  ### group
  nGroups <- ifelse(S4Vectors::isEmpty(estimated_result@group_label), 1, length(unique(estimated_result@group_label)))
  if(nGroups != 1){
    if(!is.null(existed_params[["new_covariate"]])){
      if(is.null(existed_params[["new_covariate"]][["corr_group"]])){
        group <- existed_params[["input_data"]][["corr_group"]]
        col_data$group <- group
      }else{
        group <- as.character(existed_params[["new_covariate"]][["corr_group"]])
        col_data$group <- group
      }
    }else{
      group <- existed_params[["input_data"]][["corr_group"]]
      col_data$group <- group
    }
  }
  ### batch
  nBatches <- ifelse(S4Vectors::isEmpty(estimated_result@batch_label), 1, length(unique(estimated_result@batch_label)))
  if(nBatches != 1){
    if(!is.null(existed_params[["new_covariate"]])){
      if(is.null(existed_params[["new_covariate"]][["batch_label"]])){
        batch <- existed_params[["input_data"]][["batch_label"]]
        col_data$batch <- batch
      }else{
        batch <- existed_params[["new_covariate"]][["batch_label"]]
        col_data$batch <- batch
      }
    }else{
      batch <- existed_params[["input_data"]][["batch_label"]]
      col_data$batch <- batch
    }
  }
  ### row_data
  row_data <- data.frame("gene_name" = rownames(counts))
  # Establish SingleCellExperiment
  simulate_result <- SingleCellExperiment::SingleCellExperiment(list(counts = counts),
                                                                colData = col_data,
                                                                rowData = row_data)
  simulate_result <- data_conversion(SCE_object = simulate_result, return_format = return_format)

  ##############################################################################
  ####                           Ouput                                       ###
  ##############################################################################
  slot(estimated_result, "simulation_time_memory") <- list("scDesign3" = list("estimation_time" = simulate_detection$Elapsed_Time_sec,
                                                                              "estimation_memory" = simulate_detection$Peak_RAM_Used_MiB))
  slot(estimated_result, "simulation_result") <- list("scDesign3" = simulate_result)
  slot(estimated_result, "simulation_parameters") <- list("scDesign3" = used_argument)
  return(estimated_result)
}
