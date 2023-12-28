.pseudotime_info <- function(SCST_Object,
                             dynwrap_data = NULL,
                             start_cell_id = NULL,
                             col_data,
                             seed){
  if(!requireNamespace("parallelDist", quietly = TRUE)){
    utils::install.packages("parallelDist")
  }
  if(!requireNamespace("dyndimred", quietly = TRUE)){
    devtools::install_github("dynverse/dyndimred")
  }
  ### dimension reduction
  dimred <- dyndimred::dimred_umap(t(as.matrix(SCST_Object@reference)))
  ### start cell
  if(is.null(start_cell_id)){
    set.seed(seed)
    start_cell_id <- colnames(SCST_Object@reference)[sample(1:ncol(SCST_Object@reference), size = 1)]
  }
  ### If dynwrap data is NULL
  if(is.null(dynwrap_data)){
    if(!requireNamespace("dynwrap", quietly = TRUE)){
      devtools::install_github("dynverse/dynwrap")
    }
    if(!requireNamespace("tislingshot", quietly = TRUE)){
      devtools::install_github("dynverse/ti_slingshot/package")
    }
    if(!requireNamespace("tidyr", quietly = TRUE)){
      utils::install.packages("tidyr")
    }
    dynwrap_data <- dynwrap::wrap_expression(expression = log2(t(as.matrix(SCST_Object@reference)) + 1),
                                             counts = t(as.matrix(SCST_Object@reference)))
    if(!S4Vectors::isEmpty(SCST_Object@group_label)){
      dynwrap_data <- dynwrap::add_grouping(dynwrap_data, SCST_Object@group_label)
    }
    dynwrap_data <- dynwrap::add_dimred(dynwrap_data, dimred)
    dynwrap_data <- dynwrap::infer_trajectory(dynwrap_data, tislingshot::ti_slingshot())
  }
  ### dynwrap_expression data
  traj_type <- dynwrap_data$trajectory_type
  milestone_network <- dynwrap_data$milestone_network
  start_milestone <- dynwrap_data$root_milestone_id
  if(is.null(start_milestone)){
    start_milestone <- milestone_network$from[1]
  }
  if(is.null(dynwrap_data$grouping)){
    start_cell <- col_data$cell_name[which(SCST_Object@group_label %in% start_milestone)]
    if(S4Vectors::isEmpty(start_cell)){
      start_cell <- col_data$cell_name[which(SCST_Object@group_label %in% unique(SCST_Object@group_label)[1])]
    }
    groups <- SCST_Object@group_label
  }else{
    start_cell <- names(dynwrap_data$grouping)[which(dynwrap_data$grouping %in% start_milestone)]
    groups <- dynwrap_data$grouping
  }
  set.seed(seed)
  start_cell_id <- start_cell[sample(1:length(start_cell), size = 1)]
  if(traj_type == "bifurcation"){
    inter_name <- unique(milestone_network$to)[1]
    li1_name <- c(start_milestone, milestone_network$from[-grep(start_milestone, milestone_network$from)][1],
                  inter_name)
    li2_name <- c(start_milestone, milestone_network$from[-grep(start_milestone, milestone_network$from)][2],
                  inter_name)
    dist <- parallelDist::parallelDist(dimred) %>% as.matrix()
    dist1 <- parallelDist::parallelDist(dimred[which(groups %in% li1_name), ]) %>% as.matrix()
    dist2 <- parallelDist::parallelDist(dimred[which(groups %in% li2_name), ]) %>% as.matrix()
    pseudotime1 <- data.frame("cell_name" = names(dist[which(groups %in% li1_name), start_cell_id]),
                              "pseudotime1" = S4Vectors::unname(dist[which(groups %in% li1_name), start_cell_id]))
    pseudotime2 <- data.frame("cell_name" = names(dist[which(groups %in% li2_name), start_cell_id]),
                              "pseudotime2" = S4Vectors::unname(dist[which(groups %in% li2_name), start_cell_id]))
    col_data <- col_data %>%
      dplyr::full_join(pseudotime1, by = "cell_name") %>%
      dplyr::full_join(pseudotime2, by = "cell_name") %>%
      dplyr::mutate(
        dplyr::across(c("pseudotime1", "pseudotime2"), ~ tidyr::replace_na(.x, -1)),
        l1 = dplyr::case_when(
          pseudotime1 >= 0 & pseudotime2 < 0 ~ TRUE,
          pseudotime1 >= 0 & pseudotime2 >= 0 ~ TRUE,
          TRUE ~ FALSE
        ),
        l2 = dplyr::case_when(
          pseudotime1 < 0 & pseudotime2 >= 0 ~ TRUE,
          pseudotime1 >= 0 & pseudotime2 >= 0 ~ TRUE,
          TRUE ~ FALSE
        )
      )
    col_data$l1 <- factor(col_data$l1)
    col_data$l2 <- factor(col_data$l2)
    mu_formula = "s(pseudotime1, k = 10, by = l1, bs = 'cr') + s(pseudotime2, k = 10, by = l2, bs = 'cr')"
    pseudotime <- c("pseudotime1", "pseudotime2", "l1", "l2")
  }
  if(traj_type == "linear" | traj_type == 'cycle'){
    dist <- parallelDist::parallelDist(dimred) %>% as.matrix()
    pseudotime <- data.frame("cell_name" = names(dist[, start_cell_id]),
                             "pseudotime" = unname(dist[, start_cell_id]))
    col_data <- col_data %>%
      dplyr::full_join(pseudotime, by = "cell_name")
    mu_formula = "s(pseudotime, k = 10, bs = 'cr')"
    pseudotime <- c("pseudotime")
  }
  if(traj_type == "multifurcation" | traj_type == "tree"){
    dist <- parallelDist::parallelDist(dimred) %>% as.matrix()
    pseudotime <- data.frame("cell_name" = colnames(dist))
    for(i in 1:nrow(milestone_network)){
      from_cell <- milestone_network[i, ] %>% dplyr::pull("from")
      to_cell <- milestone_network[i, ] %>% dplyr::pull("to")
      pseudotime <- data.frame("cell_name" = names(dist[grep(from_cell, dynwrap_data$grouping), grep(start_cell_id, colnames(dist))]),
                               "pseudotime" = S4Vectors::unname(dist[grep(from_cell, dynwrap_data$grouping), grep(start_cell_id, colnames(dist))]))
      col_data <- col_data %>%
        dplyr::full_join(pseudotime, by = "cell_name")
    }
    col_data <- col_data %>%
      dplyr::mutate(
        dplyr::across(3:ncol(col_data), ~ tidyr::replace_na(.x, -1))
      )
    colnames(col_data)[3:ncol(col_data)] <- paste0("pseudotime", 1:nrow(milestone_network))
    pseudotime <- paste0("pseudotime", 1:nrow(milestone_network))
    mu_formula <- paste0(paste0("s(", pseudotime, ", k = 10, bs = 'cr')"), collapse = " + ")
  }
  return(
    dplyr::lst(col_data, mu_formula, pseudotime)
  )
}
