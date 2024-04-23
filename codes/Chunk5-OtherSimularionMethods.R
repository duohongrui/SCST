library(Seurat)
library(SCST)
library(dplyr)
group_datasets <- list.files("../SCST-article/Data/Group/")

#### Simulation Pipeline
for(data_id in group_datasets){
  print(data_id)
  data <- readRDS(file.path("../SCST-article/Data/Group", data_id))
  ### ZINB-WaVE and dropsim
  SCSTObject <- Create_SCST_Object(reference = as.matrix(SeuratObject::GetAssayData(data, layer = "counts")),
                                   methods = c("zinbwave", "dropsim"),
                                   group_label = as.character(data$cell_type),
                                   seed = 888)
  estimation <- Estimate_parameters(SCSTObject, verbose = TRUE)
  estimation <- Set_customed_parameters(estimation,
                                        return_format = "Seurat")
  simulation <- Simulate_datasets(estimation, verbose = TRUE)
  saveRDS(simulation,
          file = paste0("../simulation_results/Other_results/", data_id))
}
