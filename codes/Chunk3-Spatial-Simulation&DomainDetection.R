library(Seurat)
library(SCST)
library(dplyr)
spatial_datasets <- list.files("../SCST-article/Data/Spatial/")

for(data_id in spatial_datasets){
  print(data_id)
  data <- readRDS(file.path("../SCST-article/Data/Spatial", data_id))
  ### SRTsim
  SCSTObject <- Create_SCST_Object(reference = as.matrix(SeuratObject::GetAssayData(data, layer = "counts")),
                                   methods = c("SRTsim"),
                                   group_label = as.character(data$group),
                                   spatial_x = data$imagerow,
                                   spatial_y = data$imagecol,
                                   seed = 888)
  # estimation <- SPARSim_estimation(SCSTObject, verbose = TRUE)
  estimation <- Estimate_parameters(SCSTObject, verbose = TRUE)
  estimation <- Set_customed_parameters(estimation,
                                        return_format = "Seurat")
  simulation <- Simulate_datasets(estimation, verbose = TRUE)
  simulation@simulation_result$SRTsim[['image']] <- new(Class = "SlideSeq",
                                                        assay = "RNA",
                                                        coordinates = data@meta.data %>% select(c("imagerow", "imagecol")))
  saveRDS(simulation@simulation_result$SRTsim,
          file = paste0("../simulation_results/Spatial_results/", "SRTsim_", data_id))


  ### scDesign3
  SCSTObject <- Create_SCST_Object(reference = as.matrix(SeuratObject::GetAssayData(data, layer = "counts")),
                                   methods = c("scDesign3"),
                                   group_label = as.character(data$group),
                                   spatial_x = data$imagerow,
                                   spatial_y = data$imagecol,
                                   seed = 888)
  # estimation <- SPARSim_estimation(SCSTObject, verbose = TRUE)
  estimation <- Estimate_parameters(SCSTObject, verbose = TRUE)
  estimation <- Set_customed_parameters(estimation,
                                        return_format = "Seurat")
  simulation <- Simulate_datasets(estimation, verbose = TRUE)
  simulation@simulation_result$scDesign3[['image']] <- new(Class = "SlideSeq",
                                                           assay = "RNA",
                                                           coordinates = data@meta.data %>% select(c("imagerow", "imagecol")))
  saveRDS(simulation@simulation_result$scDesign3,
          file = paste0("../simulation_results/Spatial_results/", "scDesign3_", data_id))

}

