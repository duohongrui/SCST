change_parameters <- function(input, argument){
  if(methods::is(input, "Params")){
    alternative_methods <- c("SCRIP",
                             "Splat",
                             "Lun",
                             "Lun2")
    method <- gsub("Params", "", methods::is(input)[1])
    if(method %in% alternative_methods){
      # match names in input and in parameters
      index <- methods::slotNames(input)[stats::na.omit(match(names(argument), methods::slotNames(input)))]
      # filter
      argument <- argument[index]
      # reset customed paramters
      parameters <- splatter::setParams(input, argument)
      # return
      return(parameters)
    }else{
      print_color_word(paste0("no You can not specify the parameters in ", method))
      stop(error_output())
    }
  }else{
    if(!is.null(input)){
      input_names <- names(input)
      for(i in input_names){
        if(i %in% names(argument)){
          argument[[i]] <- input[[i]]
        }
      }
    }
    param_name <- names(argument)
    argument <- methods::as(argument, "list")
    names(argument) <- param_name
    return(argument)
  }
}

random_seed <- function(){
  max_num <- .Machine$integer.max
  seed <- sample(1:max_num, 1)
  return(seed)
}


#' @importFrom SingleCellExperiment SingleCellExperiment colData rowData
#' @importFrom Seurat as.Seurat
#' @importFrom dplyr lst
#'
data_conversion <- function(
    SCE_object,
    return_format = "list"
){
  simulate_result <- SCE_object
  if(return_format == "SingleCellExperiment"){
    simulate_result <- simulate_result
  }
  if(return_format == "list"){
    count_data <- SingleCellExperiment::counts(simulate_result)
    col_meta <- as.data.frame(SingleCellExperiment::colData(simulate_result))
    row_meta <- as.data.frame(SingleCellExperiment::rowData(simulate_result))
    simulate_result <- dplyr::lst(count_data,
                                  col_meta,
                                  row_meta)
  }
  if(return_format == "Seurat"){
    simulate_result <- Seurat::as.Seurat(simulate_result,
                                         counts = "counts",
                                         data = NULL)
  }
  return(simulate_result)
}


#' A Number Multiplies Proportions Under Limitations
#'
#' @param number A number
#' @param result_sum_strict The limitation of the sum of results
#' @param prop The proportions to be multiplied by the number
#' @param prop_sum_strict The limitation of the sum of proportions
#' @param digits The number of decimal places
#'
#' @importFrom assertthat assert_that
#'
#' @export
#'
#' @examples
#' ## Save 1 decimal place
#' a <- proportionate(number = 355,
#'                    prop = c(0.2, 0.6, 0.15, 0.36),
#'                    digits = 1)
#' ## The sum of the proportions is 1
#' b <- proportionate(number = 355,
#'                    prop = c(0.2, 0.6, 0.15, 0.05),
#'                    prop_sum_strict = 1,
#'                    digits = 1)
#' ## Save 0 decimal place
#' c <- proportionate(number = 355,
#'                    prop = c(0.2, 0.6, 0.15, 0.05),
#'                    prop_sum_strict = 1,
#'                    digits = 0)
#' ## The sum of the results is 355
#' d <- proportionate(number = 355,
#'                    result_sum_strict = 355,
#'                    prop = c(0.2, 0.6, 0.15, 0.05),
#'                    prop_sum_strict = 1,
#'                    digits = 0)
proportionate <- function(
    number,
    result_sum_strict = NULL,
    prop,
    prop_sum_strict = NULL,
    digits = 0
){
  # Check
  if(!is.null(prop_sum_strict)){
    if(sum(prop) != prop_sum_strict){
      prop[length(prop)] <- prop_sum_strict - sum(prop[1:(length(prop) - 1)])
    }
  }
  # Assign
  len_prop <- length(prop)
  assign_num <- number * prop
  assign_num <- round(assign_num, digits = digits)
  if(!is.null(result_sum_strict)){
    if(sum(assign_num) != result_sum_strict){
      assign_num <- c(assign_num[1:(len_prop-1)],
                      result_sum_strict - sum(assign_num[1:(len_prop-1)]))
      assertthat::assert_that(sum(assign_num) == result_sum_strict)
    }
  }
  return(assign_num)
}

# Necessary prior information for methods is determined
chech_prior_info <- function(SCST_Object, method){
  #### group label
  if(method == "Lun2"){
    if(S4Vectors::isEmpty(SCST_Object@batch_label)){
      print_color_word(paste0("no Batch labels (slot of batch_label in SCST object) of cells are required for ", method))
      stop(error_output(), call. = FALSE)
    }
  }
  #### spatial coordinates
  if(method == "SRTsim"){
    if(any(S4Vectors::isEmpty(SCST_Object@spatial_x), S4Vectors::isEmpty(SCST_Object@spatial_y))){
      print_color_word(paste0("no Spatial coordinates of cells/spots are required for ", method))
      stop(error_output(), call. = FALSE)
    }
  }
  #### ERCC
  if(method == "SPARSim"){
    if(!identical(S4Vectors::isEmpty(SCST_Object@dilution_factor), S4Vectors::isEmpty(SCST_Object@volume))){
      print_color_word(paste0("no The dilution factor of ERCC spike-in RNA mix and volume are both required for ", method))
      stop(error_output(), call. = FALSE)
    }
  }
}


# Colorful notes inform users the outcome of checkout
print_color_word <- function(text, color = "white"){
  code <- switch(color,
                 red = 31,
                 green = 32,
                 yellow = 33,
                 blue = 34,
                 white = 37)
  if(startsWith(text, "yes")){
    text <- gsub(pattern = "^yes", replacement = "\u2714", x = text)
    code <- 34
  }
  if(startsWith(text, "no")){
    text <- gsub(pattern = "^no", replacement = "\u2716", x = text)
    code <- 31
  }
  cat(paste0("\033[", code, "m", text, "\033[0m","\n"))
}

# Output error when the message has been raised
error_output <- function(){return("procedure failed. See above information and run again.")}


# combine all results generated by multiple methods
integrate_multi_result <- function(
    SCST_Object,
    result,
    step = "estimation"
){
  method_name <- names(slot(result, paste0(step, "_result")))

  slot(SCST_Object, paste0(step, "_result")) <- append(slot(SCST_Object, paste0(step, "_result")),
                                                       slot(result, paste0(step, "_result")))
  slot(SCST_Object, paste0(step, "_time_memory")) <- append(slot(SCST_Object, paste0(step, "_time_memory")),
                                                            slot(result, paste0(step, "_time_memory")))
  slot(SCST_Object, paste0(step, "_parameters")) <- append(slot(SCST_Object, paste0(step, "_parameters")),
                                                           slot(result, paste0(step, "_parameters")))

  # if(is.null(slot(result, paste0(step, "_result")))){
  #   ### When the method does not contain the estimation step
  #
  # }else{
  #   ### When the estimation result is available
  # }
  return(SCST_Object)
}



#' Set Optional Parameters For Simulation
#'
#' @param SCST_Object The estimation result generated by [SCST::Estimate_parameters] function or the estimation function for an individual method.
#' @param return_format The format of returned simulation data. Choices: list, SingleCellExperiment and Seurat. Default is list.
#' @param cell_num The expected number of cells to be simulated.
#' @param gene_num The expected number of genes to be simulated.
#' @param nGroups The expected number of cell groups to be simulated. Default is 1.
#' @param prop_group The proportion of cells in each group. Default is 1.
#' @param de_prop The proportion of DEGs over all genes between different simulated cell groups. Default is 0.2.
#' @param fc_group The fold change of the generated DEGs. Default is 2.
#' @param nBatches The expected number of cell batches to be simulated. Default is 1.
#' @param prop_batch The proportion of cells in each batch. Default is 1.
#' @param path A Boolean (TRUE or FLASE) which determines whether to simulate datasets with trajectory along the path. Default is false.
#'
#' @return A SCST object
#' @export
#'
Set_customed_parameters <- function(
  SCST_Object,
  return_format = "list",
  cell_num = NULL,
  gene_num = NULL,
  nGroups = NULL,
  prop_group = NULL,
  de_prop = NULL,
  fc_group = NULL,
  nBatches = NULL,
  prop_batch = NULL,
  path = NULL
){
  ### Change default parameter values
  if(!is.null(cell_num)) {SCST_Object@customed_setting$cell_num <- cell_num}
  if(!is.null(gene_num)) {SCST_Object@customed_setting$gene_num <- gene_num}
  if(!is.null(nGroups)) {SCST_Object@customed_setting$nGroups <- nGroups}
  if(!is.null(prop_group)) {SCST_Object@customed_setting$prop_group <- prop_group}
  if(!is.null(de_prop)) {SCST_Object@customed_setting$de_prop <- de_prop}
  if(!is.null(fc_group)) {SCST_Object@customed_setting$fc_group <- fc_group}
  if(!is.null(nBatches)) {SCST_Object@customed_setting$nBatches <- nBatches}
  if(!is.null(prop_batch)) {SCST_Object@customed_setting$prop_batch <- prop_batch}
  if(!is.null(path)) {SCST_Object@customed_setting$path <- path}

  ### Determine useful parameters for each method
  method_names <- methods::slot(SCST_Object, "methods")
  env <- asNamespace("SCST")
  params <- SCST_Object@customed_setting
  for(i in method_names){
    right_method <- paste0(i, "_simulation")
    assign(right_method, get(right_method, envir = env))
    method_params <- change_parameters(params, formals(get(right_method)))
    method_params <- method_params[-1]
    method_params <- method_params[-grep("verbose", names(method_params))]
    method_params[["return_format"]] <- return_format
    tmp_list <- list()
    tmp_list[[i]] <- method_params
    methods::slot(SCST_Object, "customed_setting") <- append(methods::slot(SCST_Object, "customed_setting"),
                                                             tmp_list)
  }
  return(SCST_Object)
}


