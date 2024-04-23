#' Simulate Datasets
#'
#' @param SCST_Object A SCST object genaated by [SCST::Estimate_parameters] and [SCST::Set_customed_parameters] functions
#' @param verbose Whether to return messages or not
#'
#' @return A SCST object with all results of parameter estimation, simulation and used parameters.
#' @export
#'
Simulate_datasets <- function(
  SCST_Object,
  verbose = TRUE
){
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%Check Simulation Methods%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  message("(1)---Check methods and estimation results")
  ## Simulation_datasets is the core function which needs formal programming object to prevent errors
  methods <- methods::slot(SCST_Object, "methods")
  for(i in methods){
    if(S4Vectors::isEmpty(grep(i, names(SCST_Object@customed_setting)))){
      print_color_word("no Please run `Set_customed_parameters` function before simulating datasets")
      stop(error_output())
    }
    if(!i %in% names(SCST_Object@estimation_result)){
      print_color_word(paste("no No estimation results detected for", i, "method, please check and try again"))
      stop(error_output())
    }
    print_color_word(paste0("yes ", i))
  }
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Simulate New Datasets%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  message("(2)---Simulating datasets")
  env <- asNamespace("SCST")
  for(i in 1:length(methods)){
    used_method <- methods[i]
    print_color_word(paste("\u2192", used_method, "is running..."), color = "green")
    right_method <- paste0(used_method, "_simulation")
    assign(right_method, get(right_method, envir = env))
    ### get params for method
    sim_params <- SCST_Object@customed_setting[[used_method]]
    sim_params[["estimated_result"]] <- SCST_Object
    sim_params[["verbose"]] <- verbose
    if("..." %in% names(sim_params)){
      sim_params <- sim_params[-grep(pattern = "[...]", x = names(sim_params))]
    }
    if("SCST_Object" %in% names(formals(right_method))){
      sim_params[["SCST_Object"]] <- SCST_Object
    }
    ### simulate
    errors <- try(
      expr = sim_result <- do.call(right_method, sim_params),
      silent = TRUE
    )
    #### If error
    if(is(errors, "try-error")){
      print_color_word(paste("no The simulation process failed using", used_method))
      #### return error messages
      message(errors[1], appendLF = FALSE)
      ## print messages for those methods that successfully complete the estimation
      if(i > 1){
        print_color_word(paste("However, methods of",
                               paste(methods[1:(i-1)], collapse = "/"),
                               "have completed the simulation task and results will be returned."),
                         color = "blue")
        return(SCST_Object)
      }else{
        #### If the first method encountered error, stop the procedure
        stop(error_output(), call. = FALSE)
      }
    }else{
      #### If not error
      SCST_Object@simulation_result <- append(SCST_Object@simulation_result,
                                              sim_result@simulation_result)
      SCST_Object@simulation_parameters <- append(SCST_Object@simulation_parameters,
                                                  sim_result@simulation_parameters)
      SCST_Object@simulation_time_memory <- append(SCST_Object@simulation_time_memory,
                                                   sim_result@simulation_time_memory)
    }
  }
  message("\n(3)---All simulations complete")
  return(SCST_Object)
}

