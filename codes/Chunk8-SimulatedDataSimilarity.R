library(dplyr)
library(Seurat)
ref_data_list <- c(list.files("../SCST-article/Data/Batch", full.names = TRUE, pattern = "batch"),
                   list.files("../SCST-article/Data/Group", full.names = TRUE),
                   list.files("../SCST-article/Data/Spatial", full.names = TRUE))

for(i in ref_data_list){
  print(i)
  data <- readRDS(i)
  if("batch" %in% i){
    data <- JoinLayers(data)
  }
  data_name <- strsplit(i, "/") %>% sapply(FUN = function(x){x[5]})
  ### cell-level properties
  cell_properties <- .cell_properties(GetAssayData(data, layer = "counts") %>% as.matrix(), verbose = TRUE, nCores = 4)
  ### gene-level properties
  gene_properties <- .gene_properties(GetAssayData(data, layer = "counts") %>% as.matrix(), cpm_norm = TRUE, verbose = TRUE)
  saveRDS(list("cell_properties" = cell_properties,
               "gene_properties" = gene_properties),
          file = paste0("../SCST-article/Data/ref_data_properties/", data_name))
}





