library(Seurat)
library(SCST)
library(dplyr)
batch_datasets <- list.files("../SCST-article/Data/Batch/", pattern = "batch")

for(data_id in batch_datasets){
  data <- readRDS(file.path("../SCST-article/Data/Batch", data_id))
  data <- JoinLayers(data)

  ### prior information of real batches
  nBatches <- length(unique(data$orig.ident))
  prop_batch <- table(data$orig.ident)/length(data$orig.ident)

  ### scDesign3 and Splat
  SCSTObject <- Create_SCST_Object(reference = as.matrix(SeuratObject::GetAssayData(data)),
                                   methods = c("scDesign3", "Splat"),
                                   batch_label = as.character(data$orig.ident),
                                   seed = 888)
  # estimation <- SPARSim_estimation(SCSTObject, verbose = TRUE)
  estimation <- Estimate_parameters(SCSTObject, verbose = TRUE)
  estimation <- Set_customed_parameters(estimation,
                                        return_format = "Seurat",
                                        cell_num = estimation@customed_setting$cell_num,
                                        gene_num = estimation@customed_setting$gene_num,
                                        nBatches = nBatches,
                                        prop_batch = prop_batch,
                                        batch_facLoc = 0.4,
                                        batch_facScale = 0.4)
  simulation <- Simulate_datasets(estimation, verbose = TRUE)
  saveRDS(simulation@simulation_result$Splat,
          file = paste0("../simulation_results/Batch_results/", "Splat_", data_id))
  saveRDS(simulation@simulation_result$scDesign3,
          file = paste0("../simulation_results/Batch_results/", "scDesign3_", data_id))


  ### SPARSim
  param_A <- rep(1, nBatches)
  param_B <- rep(0.3, nBatches)
  estimation <- SPARSim_estimation(SCSTObject, verbose = TRUE)
  simulation <- SPARSim_simulation(estimation,
                                   nBatches = nBatches,
                                   prop_batch = prop_batch,
                                   param_A = param_A,
                                   param_B = param_B,
                                   verbose = TRUE,
                                   seed = 888,
                                   return_format = "Seurat")
  saveRDS(simulation@simulation_result$SPARSim,
          file = paste0("../simulation_results/Batch_results/", "SPARSim_", data_id))


  ### SCRIP
  modes <- c("GP-trendedBCV", "GP-trendedBCV", "GP-commonBCV", "BGP-commonBCV", "BP", "BGP-trendedBCV")
  SCSTObject <- Create_SCST_Object(reference = as.matrix(SeuratObject::GetAssayData(data)),
                                   methods = c("SCRIP"),
                                   batch_label = as.character(data$orig.ident),
                                   seed = 888)
  estimation <- Estimate_parameters(SCSTObject, verbose = TRUE)
  for(mode in modes){
    estimation <- Set_customed_parameters(estimation,
                                          mode = mode,
                                          return_format = "Seurat",
                                          cell_num = estimation@customed_setting$cell_num,
                                          gene_num = estimation@customed_setting$gene_num,
                                          nBatches = nBatches,
                                          prop_batch = prop_batch,
                                          batch_facLoc = 0.2,
                                          batch_facScale = 0.2)
    simulation <- Simulate_datasets(estimation, verbose = TRUE)
    saveRDS(simulation@simulation_result$SCRIP,
            file = paste0("../simulation_results/Batch_results/", "SCRIP-", mode, "_", data_id))
  }
}


# SCSTObject <- Create_SCST_Object(reference = as.matrix(SeuratObject::GetAssayData(data)),
#                                  methods = c("Lun2"),
#                                  batch_label = as.character(data$orig.ident),
#                                  seed = 888)
# estimation <- Estimate_parameters(SCSTObject, verbose = TRUE)
# estimation <- Set_customed_parameters(estimation, return_format = "Seurat")
# estimation@estimation_result$Lun2@plate.var <- 10
# simulation <- Simulate_datasets(estimation, verbose = TRUE)

# sim1 <- simulation1@simulation_result$Splat
# sim1 <- simulation1@simulation_result$Lun2
# sim1 <- simulation1@simulation_result$scDesign3
# sim1 <- simulation@simulation_result$Lun2
# sim1 <- fastMNN(sim1, k_NNs = 30, PCs = 20, features = rownames(sim1))
# sim1[["RNA"]] <- split(sim1[["RNA"]], f = sim1$batch)
# sim1 <- sim1 %>%
#   NormalizeData() %>%
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#   ScaleData(verbose = FALSE) %>%
#   RunPCA(npcs = 20, verbose = FALSE) %>%
#   RunUMAP(reduction = "pca", dims = 1:20, seed.use = 111)
# DimPlot(sim1, group.by = "batch", reduction = "umap")
#
#
# sim1 <- IntegrateLayers(
#   object = sim1,
#   method = HarmonyIntegration,
#   orig.reduction = "pca",
#   new.reduction = "harmony",
#   verbose = TRUE,
#   max.iter.harmony = 50
# )
# sim1 <- sim1 %>%
#   RunUMAP(reduction = "harmony", dims = 1:20, seed.use = 111)
# DimPlot(sim1, group.by = "batch")

# conda_env <- "/Users/duohongrui/opt/miniconda3/envs/scvi-env"
# Batch_removal_result <- EvalBatchRemovalMethods(simulation1,
#                                                 conda_env = conda_env,
#                                                 PCs = 20,
#                                                 nlatent = 20,
#                                                 epochs = 100,
#                                                 batch_size = 32,
#                                                 seed = 111,
#                                                 verbose = TRUE)
