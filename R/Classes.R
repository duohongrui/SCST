#' The SCST Class
#'
#' S4 class that contains the settings and information used for estimating important parameters from real reference data
#'
#' @slot reference A reference gene expression matrix
#' @slot methods The method(s) used for simulation
#' @slot cell_num The cell/spots/bin number of the reference matrix
#' @slot gene_num The gene number of the reference matrix
#' @slot group_label The group labels of cells/spots/bins
#' @slot batch_label The batch labels of cells/spots/bins
#' @slot dilution_factor If ERCC spike-in RNAs were added into the single-cell solution, the dilution factor can be used for estimating parameters from real using certain methods (such as SPARSim)
#' @slot volume The volume of ERCC spike-in mix added into the single-cell solution
#' @slot spatial_x The coordinates of the spots/bins in the spatial transcriptome data on the x-axis
#' @slot spatial_y The coordinates of the spots/bins in the spatial transcriptome data on the y-axis
#' @slot seed The random seed used for reproducibility
#' @slot estimation_result The outcome of parameter estimation step for simulating new datasets
#' @slot estimation_time_memory Time consumption and memory usage in the parameter estimation step
#' @slot estimation_parameters The input parameters used in the estimation step
#' @slot simulation_result The outcome of simulated data
#' @slot simulation_time_memory Time consumption and memory usage in the data simulation step
#' @slot simulation_parameters The input parameters used in the simulation step
#' @slot customed_setting The prior settings of different simulation scenarios are specified for generating single-cell or ST datasets with different sizes, groups, DEGs, batches and trajectories.
#'
#' @name SCST
#' @aliases SCST-class
#' @exportClass SCST
#'
#' @importClassesFrom Matrix dgCMatrix
#'
setClass(Class = "SCST",
         slots = c(reference = "dgCMatrix",
                   methods = "character",
                   cell_num = "numeric",
                   gene_num = "numeric",
                   group_label = "character",
                   batch_label = "character",
                   dilution_factor = "numeric",
                   volume = "numeric",
                   spatial_x = "numeric",
                   spatial_y = "numeric",
                   seed = "numeric",
                   estimation_result = "list",
                   estimation_time_memory = "list",
                   estimation_parameters = "list",
                   simulation_result = "list",
                   simulation_time_memory = "list",
                   simulation_parameters = "list",
                   customed_setting = "list"),
         prototype = prototype(customed_setting = list(cell_num = 1000,
                                                       gene_num = 1000,
                                                       mode = "GP-trendedBCV",
                                                       nGroups = 1,
                                                       prop_group = 1,
                                                       de_prop = 0.2,
                                                       fc_group = 2,
                                                       de_facLoc = 0.1,
                                                       de_facScale = 0.4,
                                                       nBatches = 1,
                                                       prop_batch = 1,
                                                       batch_facLoc = 0.1,
                                                       batch_facScale = 0.1,
                                                       path = FALSE,
                                                       seed = 111)))


#' Create SCST Object
#'
#' @param reference The reference gene expression matrix used for estimating parameters
#' @param methods The simulation methods used for estimating parameters and data simulation. Users can specify more than 1 methods to use and the procedure will check the prior information which are indispensable for each method.
#' @param group_label The group labels of cells/spots/bins
#' @param batch_label The batch labels of cells/spots/bins
#' @param batch_label The batch labels of cells/spots/bins
#' @param dilution_factor If ERCC spike-in RNAs were added into the single-cell solution, the dilution factor can be used for estimating parameters from real using certain methods (such as SPARSim)
#' @param volume The volume of ERCC spike-in mix added into the single-cell solution
#' @param spatial_x The coordinates of the spots/bins on the x-axis for the spatial transcriptome data
#' @param spatial_y The coordinates of the spots/bins on the y-axis for the spatial transcriptome data
#' @param seed The random seed used for reproducibility
#'
#' @importFrom methods new as
#' @export
Create_SCST_Object <- function(
    reference,
    methods,
    group_label = NULL,
    batch_label = NULL,
    dilution_factor = NULL,
    volume = NULL,
    spatial_x = NULL,
    spatial_y = NULL,
    seed = random_seed()
){
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Input Reference Data%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SCST <- methods::new("SCST")
  if(!is.data.frame(reference) & !is.matrix(reference)){
    stop("Please input a gene expression matrix (or dataframe)!")
  }
  if(is.data.frame(reference)){
    reference <- as.matrix(reference)
  }
  SCST@reference <- methods::as(reference, "dgCMatrix")
  SCST@methods <- methods
  SCST@cell_num <- ncol(SCST@reference)
  SCST@gene_num <- nrow(SCST@reference)
  SCST@customed_setting$cell_num <- ncol(SCST@reference)
  SCST@customed_setting$gene_num <- nrow(SCST@reference)
  SCST@customed_setting$seed <- random_seed()
  if(!is.null(group_label)){SCST@group_label <- group_label}
  if(!is.null(batch_label)){SCST@batch_label <- batch_label}
  if(!is.null(dilution_factor)){SCST@dilution_factor <- dilution_factor}
  if(!is.null(volume)){SCST@volume <- volume}
  if(!is.null(spatial_x)){SCST@spatial_x <- spatial_x}
  if(!is.null(spatial_y)){SCST@spatial_y <- spatial_y}
  if(!is.null(seed)){SCST@seed <- seed}
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%Check Necessary Prior Info%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### check method name
  env <- asNamespace("SCST")
  all_functions <- ls("package:SCST", envir = env, pattern = "(simulation)$")
  all_methods <- gsub(pattern = "_.*", replacement = "", x = all_functions)
  ### check prior info
  message("Check Method Name & Prior Information")
  for(i in methods){
    if(i %in% all_methods){
      chech_prior_info(SCST, i)
      print_color_word(paste0("yes ", i))
    }else{
      print_color_word(paste0("no ", i, " is not included in the package, please check out the spelling and input the right name!"))
      stop(error_output())
    }
  }
  return(SCST)
}
