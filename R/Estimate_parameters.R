#' Estimate Parameters From Reference Data
#'
#' @param SCST_Object A SCST object for downstream procedure, which is generated by [SCST::Create_SCST_Object]
#' @param verbose Whether to return messages or not
#' @param ... Other customed parameters represented in the function of methods
#'
#' @export
#'
Estimate_parameters <- function(
    SCST_Object,
    verbose = TRUE,
    ...
){
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Parameter Estimation%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  message("(1)---Parameter estimation")
  extra_params <- list(...)
  final_estimation_result <- SCST_Object
  methods <- methods::slot(SCST_Object, "methods")
  env <- asNamespace("SCST")
  for(i in 1:length(methods)){
    used_method <- methods[i]
    print_color_word(paste("\u2192", used_method, "is running..."), color = "green")
    right_method <- paste0(used_method, "_estimation")
    determine_error <- try(expr = assign(right_method, get(right_method, envir = env)),
                           silent = TRUE)
    if(is(determine_error, "try-error")){
      print_color_word(paste("Note:", used_method, "does not contain individual function for parameter estimation"), color = "yellow")
      estimation_result <- SCST_Object
      estimation_time_memory <- list()
      estimation_time_memory[[used_method]] <- list("estimation_time" = NULL,
                                                    "estimation_memory" = NULL)
      methods::slot(estimation_result, "estimation_time_memory") <- estimation_time_memory
      estimated_result <- list()
      estimated_result[[used_method]] <- list()
      methods::slot(estimation_result, "estimation_result") <- estimated_result
      estimation_parameters <- list()
      estimation_parameters[[used_method]] <- list()
      methods::slot(estimation_result, "estimation_parameters") <- estimation_parameters
      estimation_error <- NULL
    }else{
      used_params <- change_parameters(extra_params, formals(right_method))
      used_params[["SCST_Object"]] <- SCST_Object
      used_params[["verbose"]] <- TRUE
      if("..." %in% names(formals(right_method))){
        used_params[["..."]] <- "..."
      }
      estimation_error <- try(
        expr = estimation_result <- do.call(right_method, used_params),
        silent = TRUE
      )
    }
    ## If error in estimation
    if(methods::is(estimation_error, "try-error")){
      print_color_word(paste("no The estimation process failed using", used_method))
      #### return error messages
      message(estimation_error[1], appendLF = FALSE)
      ## print messages for those methods that successfully complete the estimation
      if(i > 1){
        print_color_word(paste("However, methods of",
                               paste(methods[1:(i-1)], collapse = "/"),
                               "have completed the estimation and results will be returned."),
                         color = "blue")
        return(final_estimation_result)
      }else{
        #### If the first method encountered error, stop the procedure
        stop(estimation_error[1])
      }
    }
    ## If not error in estimation
    ### When the candidate methods are more than 1
    if(length(methods) > 1){
      final_estimation_result <- integrate_multi_result(final_estimation_result, estimation_result, "estimation")
    }else{
      ### When the only one method is specified
      final_estimation_result <- estimation_result
    }
  }
  message("\n(2)---All estimations are done..")
  return(final_estimation_result)
}
