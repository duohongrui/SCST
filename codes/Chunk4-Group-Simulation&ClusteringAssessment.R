library(Seurat)
library(SCST)
library(dplyr)
group_datasets <- list.files("../SCST-article/Data/Group/")

#### Prior information of DEGs in real reference datasets
for(data_id in group_datasets){
  print(data_id)
  data <- readRDS(file.path("../SCST-article/Data/Group", data_id))
  group <- data$cell_type

  ## paired groups
  group_unique <- unique(group)
  group_paired <- utils::combn(group_unique, 2)
  ## blank result list
  result_list <- list()
  for(group_pair in 1:ncol(group_paired)){
    message("--------------------------------------------------")
    message(paste("Performing DEA with the", group_pair, "/", ncol(group_paired), "paired of groups."))
    group_candidate <- group_paired[, group_pair]
    index <- group %in% group_candidate
    sub_group <- group[index]
    sub_data <- data[, index]
    sublist_name <- paste0(group_candidate, collapse = "/")
    result_list[[sublist_name]] <- list()
    result_list <- .edgeRQLFDetRate(sub_data = as.data.frame(GetAssayData(sub_data, layer = "data")),
                                    sub_group = sub_group,
                                    sublist_name = sublist_name,
                                    result_list = result_list)
  }
  saveRDS(result_list, paste0("../SCST-article/Data/DEGs/", data_id))
}


#### Simulation Pipeline
for(data_id in group_datasets){
  print(data_id)
  data <- readRDS(file.path("../SCST-article/Data/Group", data_id))

  ### prior information of real batches
  nGroups <- length(unique(data$cell_type))
  prop_group <- as.numeric(table(data$cell_type))/length(data$cell_type)
  ### DEGs
  real_DEGs <- readRDS(file.path("../SCST-article/Data/DEGs", data_id))
  DEGs_number <- purrr::map(real_DEGs, .f = function(x){
    x[["edgeRQLFDetRate"]] %>%
      filter(FDR < 0.05) %>%
      rownames()
  })
  DEGs_number <- unique(BiocGenerics::Reduce(x = DEGs_number, f = union))
  de_prob <- length(DEGs_number)/nrow(data)

  ### Splat, SCRIP, SPARSim, Lun, scDesign, muscat, scDesign3, scDesign2
  SCSTObject <- Create_SCST_Object(reference = as.matrix(SeuratObject::GetAssayData(data, layer = "counts")),
                                   methods = c("scDesign3", "scDesign", "scDesign2", "muscat", "Splat", "SPARSim", "Lun"),
                                   group_label = as.character(data$cell_type),
                                   seed = 888)
  # estimation <- SPARSim_estimation(SCSTObject, verbose = TRUE)
  estimation <- Estimate_parameters(SCSTObject, verbose = TRUE)
  estimation <- Set_customed_parameters(estimation,
                                        return_format = "Seurat",
                                        nGroups = nGroups,
                                        prop_group = prop_group,
                                        de_facLoc = 1.2,
                                        de_facScale = 0.4,
                                        de_prop = de_prob,
                                        fc_group = 2,
                                        nBatches = 1)
  simulation <- Simulate_datasets(estimation, verbose = TRUE)
  saveRDS(simulation,
          file = paste0("../simulation_results/Group_results/", data_id))
}

