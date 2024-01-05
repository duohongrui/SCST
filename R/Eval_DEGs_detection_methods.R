#' Perform Differential Expression Analysis
#'
#' @param SimulationResult A SCST object containing the simulated results.
#' @param verbose Whether the process massages are returned.
#' @importFrom stats model.matrix p.adjust t.test var predict
#' @importFrom methods new
#' @importFrom utils combn
#' @references
#' Soneson C, Robinson M D. Bias, robustness and scalability in single-cell differential expression analysis. Nature methods, 2018, 15(4): 255-261. <https://doi.org/10.1038/nmeth.4612>
perform_DEA <- function(
    SimulationResult,
    verbose = TRUE
){
  #-------- Check Necessary Information --------#
  validated_methods <- check_info_pre_application(SimulationResult,
                                                  Group = TRUE,
                                                  DEGs = TRUE)
  #-------- Iterate All Validated Methods --------#
  DEA_results <- purrr::map(1:length(validated_methods), .f = function(x){
    used_method <- validated_methods[x]
    print_color_word(paste("------Perform Defferential Expression Analysis On", used_method), color = "blue")
    data <- SimulationResult@simulation_result[[used_method]][["count_data"]]
    col_data <- SimulationResult@simulation_result[[used_method]][["col_meta"]]
    group <- col_data[, "group"]
    ## Filter
    filter_index <- rowSums(data) > 0
    if(verbose){
      message(paste0(nrow(data) - sum(filter_index), " genes are removed when filtering"))
    }
    data <- data[filter_index, ]
    ## paired groups
    group_unique <- unique(group)
    group_paired <- utils::combn(group_unique, 2)
    ## blank result list
    result_list <- list()
    for(i in 1:ncol(group_paired)){
      group_candidate <- group_paired[, i]
      index <- group %in% group_candidate
      sub_group <- group[index]
      sub_data <- data[, index]
      sublist_name <- paste0(group_candidate, collapse = "vs")
      result_list[[sublist_name]] <- list()
      if(verbose){
        message("--------------------------------------------------")
        message(paste("Performing DEA with the", i, "/", ncol(group_paired), "paired of groups."))
      }
      ### 1. edgeRQLF
      if(!requireNamespace("edgeR")){
        message("edgeR is not installed on your device")
        message("Installing edgeR...")
        BiocManager::install("edgeR")
      }
      if(verbose){
        message("Performing DEA using edgeR QLF model")
      }
      dge <- edgeR::DGEList(sub_data, group = sub_group)
      dge <- edgeR::calcNormFactors(dge)
      design <- stats::model.matrix(~ sub_group)
      dge <- edgeR::estimateDisp(dge, design = design)
      fit <- edgeR::glmQLFit(dge, design = design)
      qlf <- edgeR::glmQLFTest(fit)
      tt <- edgeR::topTags(qlf, n = Inf)
      result_list[[sublist_name]][["edgeRQLF"]] <- tt$table

      ### 2. edgeRQLFDetRate
      if(verbose){
        message("Performing DEA using edgeR QLF model including the cellular detection rate")
      }
      cdr <- scale(colMeans(sub_data > 0))
      design <- stats::model.matrix(~ cdr + sub_group)
      dge <- edgeR::estimateDisp(dge, design = design)
      fit <- edgeR::glmQLFit(dge, design = design)
      qlf <- edgeR::glmQLFTest(fit)
      tt <- edgeR::topTags(qlf, n = Inf)
      result_list[[sublist_name]][["edgeRQLFDetRate"]] <- tt$table

      ### 3. MASTcpmDetRate
      if(!requireNamespace("MAST", quietly = TRUE)){
        message("Installing MAST...")
        BiocManager::install("MAST")
      }
      if(verbose){
        message("Performing DEA using MAST including the cellular detection rate with CPM data")
      }
      names(sub_group) <- colnames(sub_data)
      dge <- edgeR::DGEList(counts = sub_data)
      dge <- edgeR::calcNormFactors(dge)
      cpms <- edgeR::cpm(dge)
      sca <- MAST::FromMatrix(exprsArray = log2(cpms + 1),
                              cData = data.frame(wellKey = names(sub_group),
                                                 group = sub_group, cdr = cdr))
      zlmdata <- MAST::zlm(~ cdr + group, sca)
      mast <- MAST::lrTest(zlmdata, "group")
      df <- data.frame(PValue = mast[, "hurdle", "Pr(>Chisq)"],
                       FDR = stats::p.adjust(mast[, "hurdle", "Pr(>Chisq)"], method = "BH"),
                       row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
      result_list[[sublist_name]][["MASTcpmDetRate"]] <- df

      ### 4. MASTcpm
      if(verbose){
        message("Performing DEA using MAST with CPM data")
      }
      sca <- MAST::FromMatrix(exprsArray = log2(cpms + 1),
                              cData = data.frame(wellKey = names(sub_group),
                                                 group = sub_group))
      zlmdata <- MAST::zlm(~ group, sca)
      mast <- MAST::lrTest(zlmdata, "group")
      df <- data.frame(PValue = mast[, "hurdle", "Pr(>Chisq)"],
                       FDR = stats::p.adjust(mast[, "hurdle", "Pr(>Chisq)"], method = "BH"),
                       row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
      result_list[[sublist_name]][["MASTcpm"]] <- df

      ### 5. limmatrend
      if(!requireNamespace("limma", quietly = TRUE)){
        message("Installing limma...")
        BiocManager::install("limma")
      }
      if(verbose){
        message("Performing DEA using limma-trend")
      }
      design <- stats::model.matrix(~ sub_group)
      y <- methods::new("EList")
      y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
      fit <- limma::lmFit(y, design = design)
      fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
      tt <- limma::topTable(fit, n = Inf, adjust.method = "BH")
      colnames(tt)[c(4, 5)] <- c("PValue", "FDR")
      result_list[[sublist_name]][["limmatrend"]] <- tt

      ### 6. limmavoom
      if(verbose){
        message("Performing DEA using limma-voom")
      }
      design <- model.matrix(~ sub_group)
      vm <- limma::voom(dge, design = design, plot = FALSE)
      fit <- limma::lmFit(vm, design = design)
      fit <- limma::eBayes(fit)
      tt <- limma::topTable(fit, n = Inf, adjust.method = "BH")
      colnames(tt)[c(4, 5)] <- c("PValue", "FDR")
      result_list[[sublist_name]][["limmavoom"]] <- tt

      ### 7. ttest
      if(verbose){
        message("Performing DEA using t-test")
      }
      timing <- system.time({
        logcpms <- log2(cpms + 1)
        idx <- seq_len(nrow(logcpms))
        names(idx) <- rownames(logcpms)
        ttest_p <- sapply(idx, function(i) {
          stats::t.test(logcpms[i, ] ~ sub_group)$p.value
        })
      })
      df <- data.frame(PValue = ttest_p,
                       FDR = stats::p.adjust(ttest_p, method = "BH"),
                       row.names = names(ttest_p))
      result_list[[sublist_name]][["ttest"]] <- df

      ### 8. wilcox
      if(verbose){
        message("Performing DEA using wilcox-test")
      }
      wilcox_p <- sapply(idx, function(i) {
        stats::wilcox.test(logcpms[i, ] ~ sub_group)$p.value
      })
      df <- data.frame(PValue = wilcox_p,
                       FDR = stats::p.adjust(wilcox_p, method = "BH"),
                       row.names = names(wilcox_p))
      result_list[[sublist_name]][["wilcox"]] <- df
    }
    result_list
  })
  names(DEA_results) <- validated_methods

  #-------- Evaluated The Validity of DEGs --------#
  Eval_DEGs_validity <- .CalculateDEGDetectionMetrics(
    SimulationResult = SimulationResult,
    DEAResult = DEA_results
  )
  #-------- Return Results --------#
  return(Eval_DEGs_validity)
}


.CalculateDEGDetectionMetrics <- function(SimulationResult, DEAResult){
  #-------- Check Necessary Information --------#
  validated_methods <- check_info_pre_application(SimulationResult,
                                                  Group = TRUE,
                                                  DEGs = TRUE)
  #-------- Evaluating Reliability of DEGs --------#
  counts_data <- get_count_data(SimulationResult)
  col_data <- get_cell_meta(SimulationResult)
  row_data <- get_gene_meta(SimulationResult)

  DEGs_evaluation_results <- purrr::map_dfr(.x = validated_methods, .f = function(used_method){
    print_color_word(paste("------Perform Evaluation of DEGs On", used_method), color = "blue")
    predicted_result <- DEAResult[[used_method]]
    method_count <- counts_data[[used_method]]
    method_col_data <- col_data[[used_method]]
    method_row_data <- row_data[[used_method]]

    ### Iterate all group pairs and calculate the weighted metric values
    group_pairs <- length(predicted_result)
    true_DEGs_prop <- c()
    P_distribution_score <- c()
    print_color_word("1 & 2---True DEGs proportion and The Distribution of P Values...", color = "green")
    for(i in 1:group_pairs){
      ### 1. True proportion of DEGs
      paired_name <- names(predicted_result)[i]
      group1 <- stringr::str_split(paired_name, pattern = "vs", simplify = TRUE)[1]
      group2 <- stringr::str_split(paired_name, pattern = "vs", simplify = TRUE)[2]
      ### subset data
      index1 <- grep(group1, method_col_data[, "group"])
      index2 <- grep(group2, method_col_data[, "group"])
      sub_data <- method_count[, c(index1, index2)]
      sub_group <- c(rep(group1, length(index1)), rep(group2, length(index2)))
      true_DEGs <- get_true_DEGs(row_meta = method_row_data,
                                 group_pairs = group_pairs,
                                 group1_name = group1,
                                 group2_name = group2,
                                 method_name = used_method)
      if(is.null(true_DEGs)){
        weighted_true_DEGs_prop <- tibble::tibble("Simulation_Method" = used_method,
                                                  "DEA_Method" = names(predicted_result[[paired_name]]),
                                                  "prop_true_DEGs" = NA)
        weighted_P_distribution_score <- tibble::tibble("Simulation_Method" = used_method,
                                                        "DEA_Method" = names(predicted_result[[paired_name]]),
                                                        "PValue_distribution" = NA)
        next
      }
      prop <- lapply(predicted_result[[paired_name]], FUN = function(each_DEGs_result){
        each_DEGs_result <- each_DEGs_result[each_DEGs_result$"PValue" < 0.05, ]
        each_DEGs_result <- rownames(each_DEGs_result)
        intertect_genes <- BiocGenerics::intersect(true_DEGs, each_DEGs_result)
        union_genes <- BiocGenerics::union(true_DEGs, each_DEGs_result)
        # t_prop <- length(intertect_genes) / length(each_DEGs_result)
        # tt_prop <- length(intertect_genes) / length(true_DEGs)
        # t_prop * tt_prop
        length(intertect_genes) / length(union_genes)
      })
      prop <- dplyr::bind_cols(prop) %>% t()
      true_DEGs_prop <- cbind(true_DEGs_prop, prop)
      ### If the iterations on the group pairs are complete, calculate the weighted values
      if(i == group_pairs){
        weighted_true_DEGs_prop <- apply(true_DEGs_prop, MARGIN = 1, FUN = mean, na.omit = TRUE)
        weighted_true_DEGs_prop <- tibble::tibble("Simulation_Method" = used_method,
                                                  "DEA_Method" = names(weighted_true_DEGs_prop),
                                                  "prop_true_DEGs" = unname(weighted_true_DEGs_prop))
      }

      ### 2. Distribution of P Values
      distribution_p <- lapply(predicted_result[[paired_name]], FUN = function(each_DEGs_result){
        if(!requireNamespace("spgs", quietly = TRUE)){
          message("spgs is not installed on your device...")
          message("Installing spgs...")
          utils::install.packages("spgs")
        }
        #### Extract DEGs-removal data
        each_DEGs_result <- each_DEGs_result[each_DEGs_result$"PValue" < 0.05,, ]
        each_DEGs_result <- rownames(each_DEGs_result)
        index_DEGs <- which(rownames(sub_data) %in% each_DEGs_result)
        remove_DEGs_count <- sub_data[-index_DEGs, ]
        #### Perform DEA for DEGs-removal data
        ## Filter cells
        filter_cells <- colSums(remove_DEGs_count) > 0
        remove_DEGs_count <- remove_DEGs_count[, filter_cells]
        filtered_group <- sub_group[filter_cells]
        ## DEA
        dge <- edgeR::DGEList(remove_DEGs_count, group = filtered_group)
        dge <- edgeR::calcNormFactors(dge)
        cdr <- scale(colMeans(remove_DEGs_count > 0))
        design <- stats::model.matrix(~ cdr + filtered_group)
        dge <- edgeR::estimateDisp(dge, design = design)
        fit <- edgeR::glmQLFit(dge, design = design)
        qlf <- edgeR::glmQLFTest(fit)
        tt <- edgeR::topTags(qlf, n = Inf)
        p <- tt$table$PValue
        if(length(p) < 30){
          score <- NA
        }else{
          ### Chi-squared test
          result <- spgs::chisq.unif.test(p, bins = 20)
          if(result$p.value < 0.05){
            score <- 0
          }else{
            score <- 1
          }
        }
        score
      })
      distribution_p <- dplyr::bind_cols(distribution_p) %>% t()
      P_distribution_score <- cbind(P_distribution_score, distribution_p)
      ### If the iterations on the group pairs are complete, calculate the weighted values
      if(i == group_pairs){
        weighted_P_distribution_score <- apply(P_distribution_score, MARGIN = 1, FUN = mean, na.omit = TRUE)
        weighted_P_distribution_score <- tibble::tibble("Simulation_Method" = used_method,
                                                        "DEA_Method" = names(weighted_P_distribution_score),
                                                        "PValue_distribution" = unname(weighted_P_distribution_score))
      }
    }
    ### 3. Capality of identifying cell groups
    #### filter genes with constant value across cells
    gene_var <- apply(method_count, 1, stats::var)
    model_data <- method_count[gene_var != 0, ]
    ## scale
    model_data <- scale(model_data, center = FALSE)
    ## split data into train and test subsets
    set.seed(1111)
    train_index <- sample(1:ncol(model_data), ncol(model_data) * 0.8, replace = FALSE)
    model_data <- as.data.frame(t(model_data))
    group <- method_col_data[, "group"]
    group <- as.factor(group)
    ## SVM
    print_color_word("3---SVM...", color = "green")
    SVM_result <- lapply(predicted_result[[paired_name]], FUN = function(each_DEGs_result){
      each_DEGs_result <- each_DEGs_result[each_DEGs_result$"PValue" < 0.05, ]
      each_DEGs_result <- rownames(each_DEGs_result)
      gene_for_model <- intersect(each_DEGs_result, names(gene_var)[gene_var != 0])
      train_data <- model_data[train_index, gene_for_model]
      test_data <- model_data[-train_index, gene_for_model]
      ## group information
      train_group <- group[train_index]
      test_group <- group[-train_index]
      ## SVM
      error <- try(
        SVM_result <- model_predict(train_data = train_data,
                                    train_group = train_group,
                                    test_data = test_data,
                                    test_group = test_group),
        silent = TRUE)
      if(methods::is(error, "try-error")){
        print_color_word("The model traning or predicting failed and the result is NA", color = "yellow")
        SVM_result <- NULL
        AUC <- NA
        Accuracy <- NA
        Precision <- NA
        Recall <- NA
        F1 <- NA
      }else{
        AUC <- as.numeric(SVM_result$roc$auc)
        Accuracy <- unname(SVM_result[["conf_matrix"]][["overall"]][1])
        if(length(unique(group)) == 2){
          Precision <- unname(SVM_result[["conf_matrix"]][["byClass"]]["Precision"])
          Recall <- unname(SVM_result[["conf_matrix"]][["byClass"]]["Recall"])
          F1 <- unname(SVM_result[["conf_matrix"]][["byClass"]]["F1"])
        }else{
          Precision <- mean(SVM_result[["conf_matrix"]][["byClass"]][, "Precision"], na.rm = TRUE)
          Recall <- mean(SVM_result[["conf_matrix"]][["byClass"]][, "Recall"], na.rm = TRUE)
          F1 <- mean(SVM_result[["conf_matrix"]][["byClass"]][, "F1"], na.rm = TRUE)
        }
      }
      list(AUC = AUC,
           Accuracy = Accuracy,
           Precision = Precision,
           Recall = Recall,
           F1 = F1)
    })
    SVM_result <- purrr::map_df(SVM_result, dplyr::bind_rows)
    prediction_result <- tibble::tibble("Simulation_Method" = used_method,
                                        "DEA_Method" = weighted_P_distribution_score$DEA_Method,
                                        "AUC" = SVM_result$AUC,
                                        "Accuracy" = SVM_result$Accuracy,
                                        "Precision" = SVM_result$Precision,
                                        "Recall" = SVM_result$Recall,
                                        "F1" = SVM_result$F1)
    all_results <- weighted_true_DEGs_prop %>%
      dplyr::full_join(weighted_P_distribution_score, by = c("Simulation_Method", "DEA_Method")) %>%
      dplyr::full_join(prediction_result, by = c("Simulation_Method", "DEA_Method"))
    ## return results from a simulation method
    all_results
  })
  return(DEGs_evaluation_results)
}
