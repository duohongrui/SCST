
#' get_cell_meta <- function(SCSTObject){
#'   UseMethod("get_cell_meta")
#' }
#'

#' get_cell_meta.SCST <- function(SCSTObject){
#'   col_data <- SingleCellExperiment::colData(SCSTObject) %>% as.data.frame()
#'   return(col_data)
#' }
#'

#' get_cell_meta.SingleCellExperiment <- function(SCSTObject){
#'   col_data <- SingleCellExperiment::colData(SCSTObject) %>% as.data.frame()
#'   return(col_data)
#' }
#'

#' get_cell_meta.Seurat <- function(SCSTObject){
#'   col_data <- SCSTObject$meta.data %>% as.data.frame()
#'   return(col_data)
#' }
#'
#'

#' get_cell_meta.list <- function(SCSTObject){
#'   col_data <- SCSTObject$col_meta %>% as.data.frame()
#'   return(col_data)
#' }
#'

#' get_cell_meta.default <- function(SCSTObject){
#'   col_data <- SCSTObject$col_meta %>% as.data.frame()
#'   return(col_data)
#' }
