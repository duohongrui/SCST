library(Seurat)
library(dplyr)
library(ggplot2)
preprocessing_function <- function(matrix,
                                   meta_data = NULL,
                                   project = NULL,
                                   species = "human"){
  seurat <- CreateSeuratObject(counts = matrix,
                               project = project,
                               min.cells = 5,
                               min.features = 300,
                               meta.data = meta_data)
  seurat$orig.ident <- factor(project, levels = project)
  if(species == "human"){
    seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "mt")
  }else{
    seurat <- PercentageFeatureSet(seurat, pattern = "^mt", col.name = "mt")
  }
  print(VlnPlot(seurat, features = c("nCount_RNA", "nFeature_RNA", "mt"), group.by = "orig.ident"))
  if(project == "Batch_10X_V2"){
    seurat <- subset(seurat, nFeature_RNA > 200 & nFeature_RNA < 3000 & mt < 10)
  }
  if(project == "Batch_10X_V3"){
    seurat <- subset(seurat, nFeature_RNA > 200 & nFeature_RNA < 3000 & mt < 15)
  }
  if(project %in% c("Batch_10X_bladder", "Group_Smart-seq2_trachea")){
    seurat <- subset(seurat, nFeature_RNA > 500 & nFeature_RNA < 5000)
  }
  if(project == "Batch_Smart-seq2_bladder"){
    seurat <- subset(seurat, nFeature_RNA > 2000 & nFeature_RNA < 7000)
  }
  if(project %in% c("Batch_10X_spleen", "Batch_Smart-seq2_spleen", "Batch_10X_liver")){
    seurat <- subset(seurat, nFeature_RNA > 200 & nFeature_RNA < 3000)
  }
  if(project == "Batch_Smart-seq2_liver"){
    seurat <- subset(seurat, nFeature_RNA < 7500)
  }
  if(project == "Batch_CEL-seq_pancreas"){
    seurat <- subset(seurat, nFeature_RNA > 400 & nFeature_RNA < 6000)
  }
  if(project == "Batch_CEL-seq2_pancreas"){
    seurat <- subset(seurat, nFeature_RNA > 500 & nFeature_RNA < 7500)
  }
  if(project == "Batch_Smart-seq2_pancreas"){
    seurat <- subset(seurat, nFeature_RNA > 500 & nFeature_RNA < 10000)
  }
  if(project == "Batch_C1_pancreas"){
    seurat <- subset(seurat, mt < 20)
  }
  if(project == "Batch_inDrop_pancreas"){
    seurat <- subset(seurat, nFeature_RNA < 3500)
  }
  if(project == "Group_10X_lung"){
    seurat <- subset(seurat, nFeature_RNA < 3500)
  }
  if(project == "Group_10X_marrow"){
    seurat <- subset(seurat, nFeature_RNA < 3000)
  }
  if(project == "Group_Smart-seq2_tongue"){
    seurat <- subset(seurat, nFeature_RNA > 2500 & nFeature_RNA < 6500)
  }
  seurat <- seurat %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 20, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:20, seed.use = 111) %>%
    FindNeighbors(reduction = "pca", dims = 1:20) %>%
    FindClusters(resolution = 0.6, random.seed = 111)
  DimPlot(seurat, group.by = "group")
  return(seurat)
}



###################------------------------------------------###################
##################               Batch_1_data
###################------------------------------------------###################
Batch_10X_V2 <- Read10X("../SCST-article/Data/Batch_10X_V2/raw/")
Batch_10X_V2 <- preprocessing_function(Batch_10X_V2, project = "Batch_10X_V2")

Batch_10X_V3 <- Read10X("../SCST-article/Data/Batch_10X_V3/raw/")
Batch_10X_V3 <- preprocessing_function(Batch_10X_V3, project = "Batch_10X_V3")


DimPlot(Batch_10X_V2, label = TRUE)
DotPlot(Batch_10X_V2, features = c("CD3D", ### T cells
                                   "CD4", ### CD4 T
                                   "CD8A", ### CD8 T
                                   "GNLY", "NKG7", "KLRF1", "NCR1", "NCAM1", ### NK
                                   "MS4A1", "CD79A", ### B cells
                                   "CD27", ### naive B (MS4A1+/CD27−) and memory B cells (MS4A1+/CD27+)
                                   "CD14", "FCGR3A", ### Monocyte
                                   "FCER1A", "CST3", ### dendritic cells
                                   "PPBP"), scale.by = "size") + ### Platelet
  coord_flip()
Batch_10X_V3_cluster <- c("CD4 T","CD4 T","CD4 T","Monocyte","CD8 T","Naive B",
                          "Memory B","CD4 T","NK","NK","CD8 T","DC","Platelet")
Batch_10X_V2_cluster <- c("Monocyte","CD4 T","CD4 T","Naive B","CD8 T","CD8 T",
                          "CD8 T","NK", "Memory B","Monocyte", "DC", "DC","Platelet")
names(Batch_10X_V2_cluster) <- levels(Batch_10X_V2)
Batch_10X_V2 <- RenameIdents(Batch_10X_V2, Batch_10X_V2_cluster)
Batch_10X_V2$"cell_type" <- Idents(Batch_10X_V2)
DimPlot(Batch_10X_V2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(Batch_10X_V2, file = "../SCST-article/Data/Batch/Batch_10X_V2.rds")


#### Read both datasets
Batch_10X_V3 <- readRDS("../SCST-article/Data/Batch/Batch_10X_V3.rds")
Batch_10X_V2 <- readRDS("../SCST-article/Data/Batch/Batch_10X_V2.rds")
batch1_data <- merge(Batch_10X_V2, Batch_10X_V3)
batch1_data <- batch1_data %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 20, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:20, seed.use = 111) %>%
  RunTSNE(reduction = "pca", dims = 1:20, seed.use = 111) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.6, random.seed = 111)
DimPlot(batch1_data, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "orig.ident") + NoLegend()
DimPlot(batch1_data, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "cell_type") + NoLegend()
saveRDS(batch1_data, file = "../SCST-article/Data/Batch/batch1_data.rds")

# batch1_data <- batch1_data %>%
#   IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony")
# batch1_data <- batch1_data %>%
#   FindNeighbors(reduction = "harmony", dims = 1:20) %>%
#   FindClusters(resolution = 0.6, random.seed = 111) %>%
#   RunUMAP(reduction = "harmony", dims = 1:20, seed.use = 111) %>%
#   RunTSNE(reduction = "harmony", dims = 1:20, seed.use = 111)
# batch1_data <- JoinLayers(batch1_data)


###################------------------------------------------###################
##################               Batch_2_data
###################------------------------------------------###################
### 10X
bladder_10X_file <- file.path("../SCST-article/Data/Batch_10X_Bladder",
                              list.files("../SCST-article/Data/Batch_10X_Bladder/"))
names(bladder_10X_file) <- c("10X_P4_3", "10X_P4_4", "10X_P7_7")
batch_10X_bladder <- Read10X(data.dir = bladder_10X_file, strip.suffix = TRUE)
annotation_10X <- read.csv("../SCST-article/Data/Annotation/annotations_droplet.csv")
intersect_index <- intersect(colnames(batch_10X_bladder), annotation_10X$cell)
bladder_10X_metadata <- annotation_10X %>%
  filter(cell %in% intersect_index) %>%
  select(cell, free_annotation) %>%
  rename(cell_type = free_annotation) %>%
  tibble::column_to_rownames("cell")
batch_10X_bladder <- preprocessing_function(matrix = batch_10X_bladder,
                                            meta_data = bladder_10X_metadata,
                                            project = "Batch_10X_bladder",
                                            species = "mouse")
batch_10X_bladder$orig.ident <- factor("Batch_10X_bladder", levels = "Batch_10X_bladder")
saveRDS(batch_10X_bladder, file = "../SCST-article/Data/Batch/Batch_10X_bladder.rds")

### Smart-seq2
batch_Smartseq2_bladder <- read.csv("../SCST-article/Data/Batch_Smart-seq2_Bladder/Bladder-counts.csv", row.names = 1)
annotation_Smartseq2 <- read.csv("../SCST-article/Data/Annotation/annotations_FACS.csv") %>%
  filter(tissue == "Bladder")
intersect_index <- intersect(colnames(batch_Smartseq2_bladder), annotation_Smartseq2$cell)
batch_Smartseq2_bladder <- batch_Smartseq2_bladder[, intersect_index]
bladder_Smartseq2_metadata <- annotation_Smartseq2 %>%
  filter(cell %in% intersect_index) %>%
  select(cell, free_annotation) %>%
  rename(cell_type = free_annotation) %>%
  tibble::column_to_rownames("cell")
batch_Smartseq2_bladder <- preprocessing_function(matrix = batch_Smartseq2_bladder,
                                                  meta_data = bladder_Smartseq2_metadata,
                                                  project = "Batch_Smart-seq2_bladder",
                                                  species = "mouse")
saveRDS(batch_Smartseq2_bladder, file = "../SCST-article/Data/Batch/Batch_Smart-seq2_bladder.rds")


### Merge
batch_10X_bladder <- readRDS("../SCST-article/Data/Batch/Batch_10X_bladder.rds")
batch_Smartseq2_bladder <- readRDS("../SCST-article/Data/Batch/Batch_Smart-seq2_bladder.rds")
batch2_data <- merge(batch_10X_bladder, batch_Smartseq2_bladder)
batch2_data <- batch2_data %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 20, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:20, seed.use = 111) %>%
  RunTSNE(reduction = "pca", dims = 1:20, seed.use = 111) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.6, random.seed = 111)
DimPlot(batch2_data, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "orig.ident") + NoLegend()
DimPlot(batch2_data, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "cell_type") + NoLegend()
saveRDS(batch2_data, file = "../SCST-article/Data/Batch/batch2_data.rds")



###################------------------------------------------###################
##################               Batch_3_data
###################------------------------------------------###################
### 10X
spleen_10X_file <- file.path("../SCST-article/Data/Batch_10X_Spleen",
                             list.files("../SCST-article/Data/Batch_10X_Spleen/"))
names(spleen_10X_file) <- c("10X_P4_7", "10X_P7_6")
batch_10X_spleen <- Read10X(data.dir = spleen_10X_file, strip.suffix = TRUE)
annotation_10X <- read.csv("../SCST-article/Data/Annotation/annotations_droplet.csv")
intersect_index <- intersect(colnames(batch_10X_spleen), annotation_10X$cell)
batch_10X_spleen <- batch_10X_spleen[, intersect_index]
spleen_10X_metadata <- annotation_10X %>%
  filter(cell %in% intersect_index) %>%
  select(cell, cell_ontology_class) %>%
  rename(cell_type = cell_ontology_class) %>%
  tibble::column_to_rownames("cell")
batch_10X_spleen <- preprocessing_function(matrix = batch_10X_spleen,
                                           meta_data = spleen_10X_metadata,
                                           project = "Batch_10X_spleen",
                                           species = "mouse")
batch_10X_spleen$orig.ident <- factor("Batch_10X_spleen", levels = "Batch_10X_spleen")
saveRDS(batch_10X_spleen, file = "../SCST-article/Data/Batch/Batch_10X_spleen.rds")

### Smart-seq2
batch_Smartseq2_spleen <- read.csv("../SCST-article/Data/Batch_Smart-seq2_Spleen/Spleen-counts.csv", row.names = 1)
annotation_Smartseq2 <- read.csv("../SCST-article/Data/Annotation/annotations_FACS.csv") %>%
  filter(tissue == "Spleen")
intersect_index <- intersect(colnames(batch_Smartseq2_spleen), annotation_Smartseq2$cell)
batch_Smartseq2_spleen <- batch_Smartseq2_spleen[, intersect_index]
spleen_Smartseq2_metadata <- annotation_Smartseq2 %>%
  filter(cell %in% intersect_index) %>%
  select(cell, cell_ontology_class) %>%
  rename(cell_type = cell_ontology_class) %>%
  tibble::column_to_rownames("cell")
batch_Smartseq2_spleen <- preprocessing_function(matrix = batch_Smartseq2_spleen,
                                                 meta_data = spleen_Smartseq2_metadata,
                                                 project = "Batch_Smart-seq2_spleen",
                                                 species = "mouse")
saveRDS(batch_Smartseq2_spleen, file = "../SCST-article/Data/Batch/Batch_Smart-seq2_spleen.rds")


### Merge
batch_10X_spleen <- readRDS("../SCST-article/Data/Batch/Batch_10X_spleen.rds")
batch_Smartseq2_spleen <- readRDS("../SCST-article/Data/Batch/Batch_Smart-seq2_spleen.rds")
batch3_data <- merge(batch_10X_spleen, batch_Smartseq2_spleen)
batch3_data <- batch3_data %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 20, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:20, seed.use = 111) %>%
  RunTSNE(reduction = "pca", dims = 1:20, seed.use = 111) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.6, random.seed = 111)
DimPlot(batch3_data, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "orig.ident") + NoLegend()
DimPlot(batch3_data, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "cell_type")
saveRDS(batch3_data, file = "../SCST-article/Data/Batch/batch3_data.rds")



###################------------------------------------------###################
##################               Batch_4_data
###################------------------------------------------###################
### 10X
liver_10X_file <- file.path("../SCST-article/Data/Batch_10X_liver",
                            list.files("../SCST-article/Data/Batch_10X_liver/"))
names(liver_10X_file) <- c("10X_P4_2", "10X_P7_0", "10X_P7_1")
batch_10X_liver <- Read10X(data.dir = liver_10X_file, strip.suffix = TRUE)
annotation_10X <- read.csv("../SCST-article/Data/Annotation/annotations_droplet.csv")
intersect_index <- intersect(colnames(batch_10X_liver), annotation_10X$cell)
batch_10X_liver <- batch_10X_liver[, intersect_index]
liver_10X_metadata <- annotation_10X %>%
  filter(cell %in% intersect_index) %>%
  select(cell, cell_ontology_class) %>%
  rename(cell_type = cell_ontology_class) %>%
  tibble::column_to_rownames("cell")
batch_10X_liver <- preprocessing_function(matrix = batch_10X_liver,
                                          meta_data = liver_10X_metadata,
                                          project = "Batch_10X_liver",
                                          species = "mouse")
batch_10X_liver$orig.ident <- factor("Batch_10X_liver", levels = "Batch_10X_liver")
saveRDS(batch_10X_liver, file = "../SCST-article/Data/Batch/Batch_10X_liver.rds")

### Smart-seq2
batch_Smartseq2_liver <- read.csv("../SCST-article/Data/Batch_Smart-seq2_liver/Liver-counts.csv", row.names = 1)
annotation_Smartseq2 <- read.csv("../SCST-article/Data/Annotation/annotations_FACS.csv") %>%
  filter(tissue == "Liver")
intersect_index <- intersect(colnames(batch_Smartseq2_liver), annotation_Smartseq2$cell)
batch_Smartseq2_liver <- batch_Smartseq2_liver[, intersect_index]
liver_Smartseq2_metadata <- annotation_Smartseq2 %>%
  filter(cell %in% intersect_index) %>%
  select(cell, free_annotation) %>%
  rename(cell_type = free_annotation) %>%
  tibble::column_to_rownames("cell")
batch_Smartseq2_liver <- preprocessing_function(matrix = batch_Smartseq2_liver,
                                                meta_data = liver_Smartseq2_metadata,
                                                project = "Batch_Smart-seq2_liver",
                                                species = "mouse")
saveRDS(batch_Smartseq2_liver, file = "../SCST-article/Data/Batch/Batch_Smart-seq2_liver.rds")

### Merge
batch4_data <- merge(batch_10X_liver, batch_Smartseq2_liver)
batch4_data <- batch4_data %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 20, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:20, seed.use = 111) %>%
  RunTSNE(reduction = "pca", dims = 1:20, seed.use = 111) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.6, random.seed = 111)
DimPlot(batch4_data, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "orig.ident") + NoLegend()
DimPlot(batch4_data, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "cell_type")
saveRDS(batch4_data, file = "../SCST-article/Data/Batch/batch4_data.rds")


###################------------------------------------------###################
##################               Batch_5_data
###################------------------------------------------###################
meta_data <- read.table("../SCST-article/Data/Annotation/metadata.txt", header = TRUE, sep = ",") %>%
  mutate(
    batch = case_when(
      dataset == "celseq" ~ "Batch_CEL-seq_Pancreas",
      dataset == "celseq2" ~ "Batch_CEL-seq2_Pancreas",
      dataset == "smartseq" ~ "Batch_Smart-seq2_Pancreas",
      dataset == "c1" ~ "Batch_C1_Pancreas",
      dataset == "indrop" ~ "Batch_inDrop_Pancreas"
    ),
    cell_id = stringr::str_split(cell_id, "_", n = 2, simplify = TRUE)[, 2]
  ) %>%
  select(cell_id, batch, cell_type)


pancreas_markers <- c("GCG", ### alpha
                      "MAFA", ### beta
                      "PPY", ### gamma
                      "SST", ### delta
                      "PRSS1", ### acinar
                      "KRT19", ### ductal
                      "CDH5", ### endothelial
                      "COL1A2", ### stellate
                      "PTPRC") ### immune

### Batch_CEL-seq_Pancreas
batch_celseq_pancreas <- readRDS("../SCST-article/Data/Batch_CEL-seq_Pancreas/Batch_CEL-seq_Pancreas.rds")
batch_celseq_pancreas <- preprocessing_function(batch_celseq_pancreas, project = "Batch_CEL-seq_pancreas")
DimPlot(batch_celseq_pancreas, label = TRUE) + NoLegend()
DotPlot(batch_celseq_pancreas, features = pancreas_markers)
new_idents <- c("ductal", "acinar", "alpha", "beta", "alpha", "acinar", "mix", "acinar", "stellate")
names(new_idents) <- levels(batch_celseq_pancreas)
batch_celseq_pancreas <- RenameIdents(batch_celseq_pancreas, new_idents)
batch_celseq_pancreas$"cell_type" <- Idents(batch_celseq_pancreas)

### extract cluster6
cluster6 <- subset(batch_celseq_pancreas, seurat_clusters == "6")
cluster6 <- preprocessing_function(GetAssayData(cluster6), project = "cluster6")
DimPlot(cluster6, label = TRUE) + NoLegend()
DotPlot(cluster6, features = pancreas_markers)
cluster6_idents <- c("delta", "gamma")
names(cluster6_idents) <- levels(cluster6)
cluster6 <- RenameIdents(cluster6, cluster6_idents)
cluster6$"cell_type" <- Idents(cluster6)

### Add cluster6 annotation
meta_data <- batch_celseq_pancreas@meta.data
meta_data$cell_type <- as.character(meta_data$cell_type)
meta_data[names(Idents(cluster6)), "cell_type"] <- as.character(Idents(cluster6))
meta_data$cell_type <- factor(meta_data$cell_type,
                              levels = c("ductal", "acinar", "alpha", "beta", "delta", "gamma", "stellate"))
batch_celseq_pancreas@meta.data <- meta_data
Idents(batch_celseq_pancreas) <- batch_celseq_pancreas$cell_type
saveRDS(batch_celseq_pancreas, file = "../SCST-article/Data/Batch/Batch_CEL-seq_pancreas.rds")

### Batch_CEL-seq2_Pancreas
batch_celseq2_pancreas <- readRDS("../SCST-article/Data/Batch_CEL-seq2_Pancreas/Batch_CEL-seq2_Pancreas.rds")
batch_celseq2_pancreas <- preprocessing_function(batch_celseq2_pancreas, project = "Batch_CEL-seq2_pancreas")
DimPlot(batch_celseq2_pancreas, label = TRUE) + NoLegend()
DotPlot(batch_celseq2_pancreas, features = pancreas_markers)
new_idents <- c("alpha", "alpha", "beta", "ductal", "delta", "acinar",
                "gamma", "immune", "stellate", "delta", "acinar", "endothelial")
names(new_idents) <- levels(batch_celseq2_pancreas)
batch_celseq2_pancreas <- RenameIdents(batch_celseq2_pancreas, new_idents)
batch_celseq2_pancreas$"cell_type" <- Idents(batch_celseq2_pancreas)
saveRDS(batch_celseq2_pancreas, file = "../SCST-article/Data/Batch/Batch_CEL-seq2_pancreas.rds")

### Batch_Smart-seq2_Pancreas
batch_smartseq2_pancreas <- readRDS("../SCST-article/Data/Batch_Smart-seq2_Pancreas/Batch_Smart-seq2_Pancreas.rds")
batch_smartseq2_pancreas <- preprocessing_function(batch_smartseq2_pancreas, project = "Batch_Smart-seq2_pancreas")
DimPlot(batch_smartseq2_pancreas, label = TRUE) + NoLegend()
DotPlot(batch_smartseq2_pancreas, features = pancreas_markers)
new_idents <- c("alpha", "ductal", "alpha", "alpha", "beta", "gamma",
                "acinar", "alpha", "delta", "stellate", "alpha", "endothelial", "alpha")
names(new_idents) <- levels(batch_smartseq2_pancreas)
batch_smartseq2_pancreas <- RenameIdents(batch_smartseq2_pancreas, new_idents)
batch_smartseq2_pancreas$"cell_type" <- Idents(batch_smartseq2_pancreas)
saveRDS(batch_smartseq2_pancreas, file = "../SCST-article/Data/Batch/Batch_Smart-seq2_pancreas.rds")


### Batch_C1_Pancreas
batch_C1_pancreas <- readRDS("../SCST-article/Data/Batch_C1_Pancreas/Batch_C1_Pancreas.rds")
batch_C1_pancreas <- preprocessing_function(batch_C1_pancreas, project = "Batch_C1_pancreas")
cell_type <- meta_data %>%
  filter(cell_id %in% colnames(batch_C1_pancreas)) %>%
  pull(cell_type)
batch_C1_pancreas$"cell_type" <- cell_type
batch_C1_pancreas$cell_type[batch_C1_pancreas$cell_type == ""] <- "Unknown"
Idents(batch_C1_pancreas) <- batch_C1_pancreas$cell_type
DimPlot(batch_C1_pancreas, label = TRUE) + NoLegend()
saveRDS(batch_C1_pancreas, file = "../SCST-article/Data/Batch/Batch_C1_pancreas.rds")

### Batch_inDrop_Pancreas
batch_inDrop_pancreas <- readRDS("../SCST-article/Data/Batch_inDrop_Pancreas/Batch_inDrop_Pancreas.rds")
batch_inDrop_pancreas <- preprocessing_function(batch_inDrop_pancreas, project = "Batch_inDrop_pancreas")
cell_type <- meta_data %>%
  filter(cell_id %in% colnames(batch_inDrop_pancreas)) %>%
  pull(cell_type)
batch_inDrop_pancreas$"cell_type" <- cell_type
Idents(batch_inDrop_pancreas) <- batch_inDrop_pancreas$cell_type
DimPlot(batch_inDrop_pancreas, label = TRUE) + NoLegend()
saveRDS(batch_inDrop_pancreas, file = "../SCST-article/Data/Batch/Batch_inDrop_pancreas.rds")

### Merge
batch5_data <- merge(batch_celseq_pancreas, list(batch_celseq2_pancreas,
                                                 batch_smartseq2_pancreas,
                                                 batch_C1_pancreas,
                                                 batch_inDrop_pancreas))
batch5_data <- batch5_data %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 20, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:20, seed.use = 111) %>%
  RunTSNE(reduction = "pca", dims = 1:20, seed.use = 111) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.6, random.seed = 111)
DimPlot(batch5_data, label = TRUE, group.by = "orig.ident") + NoLegend()
DimPlot(batch5_data, label = TRUE, group.by = "cell_type") + NoLegend()
saveRDS(batch5_data, file = "../SCST-article/Data/Batch/batch5_data.rds")

# batch5_data <- batch5_data %>%
#   IntegrateLayers(orig.reduction = "pca", method = HarmonyIntegration, new.reduction = "harmony")
# batch5_data <- batch5_data %>%
#   FindNeighbors(reduction = "harmony", dims = 1:20) %>%
#   FindClusters(resolution = 0.6, random.seed = 111) %>%
#   RunUMAP(reduction = "harmony", dims = 1:20, seed.use = 111) %>%
#   RunTSNE(reduction = "harmony", dims = 1:20, seed.use = 111)


###################------------------------------------------###################
##################               Spatial Data
###################------------------------------------------###################
spatial_data_preprocessing <- function(spatial_data){
  spatial_data <- spatial_data %>%
    SCTransform(assay = "Spatial", variable.features.n = 2000) %>%
    RunPCA(npcs = 20, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:20, seed.use = 111)
  return(spatial_data)
}

### Spatial_DPLFC
DPLFC_files <- list.dirs("../SCST-article/Data", recursive = FALSE)
DPLFC_files <- DPLFC_files[grep(pattern = "DLPFC", x = DPLFC_files)]
for(i in DPLFC_files){
  print(i)
  spatial_data <- Load10X_Spatial(data.dir = i)
  data_name <- stringr::str_split(i, pattern = "\\/", simplify = TRUE)[4]
  #### Add annotated layers
  metadata <- read.table(paste0(i, "/metadata.tsv"))
  metadata <- metadata[colnames(spatial_data), ]
  spatial_data$"group" <- metadata$layer_guess
  print(SpatialDimPlot(spatial_data, group.by = "group"))
  #### Add coordinates to metadata
  spatial_coordinates <- GetTissueCoordinates(spatial_data)
  spatial_data@meta.data <- spatial_data@meta.data %>%
    mutate(imagerow = spatial_coordinates$imagerow,
           imagecol = spatial_coordinates$imagecol)
  #### Filtering
  genes_keep <- rowSums(GetAssayData(spatial_data, layer = "counts") != 0)
  genes_keep <- rownames(spatial_data)[genes_keep >= 20]
  spatial_data <- subset(spatial_data, features = genes_keep)
  spatial_data <- spatial_data_preprocessing(spatial_data)
  print(spatial_data)
  saveRDS(spatial_data, file = paste0("../SCST-article/Data/Spatial/", data_name, ".rds"))
}

### Slide-seq
library(SpatialPCA)
library(peakRAM)
load("../SCST-article/Data/Spatial_Slide-seq_cerebellum/slideseq.rds")

mem_sparse1 <- peakRAM({
  start_time <- Sys.time()
  slideseq = CreateSpatialPCAObject(counts = sp_count,
                                    location = location,
                                    project = "SpatialPCA",
                                    gene.type ="spatial",
                                    sparkversion = "sparkx",
                                    numCores_spark = 5,
                                    customGenelist = NULL,
                                    min.loctions = 20,
                                    min.features = 20)
  end_time <- Sys.time()
  T_sparse1 = end_time - start_time
})

mem_sparse1 <- peakRAM({
  start_time <- Sys.time()
  slideseq <- SpatialPCA_buildKernel(slideseq,
                                     kerneltype = "gaussian",
                                     bandwidthtype = "Silverman",
                                     bandwidth.set.by.user = NULL,
                                     sparseKernel = TRUE,
                                     sparseKernel_tol = 1e-20,
                                     sparseKernel_ncore = 5)
  slideseq <- SpatialPCA_EstimateLoading(slideseq,
                                         fast = TRUE,
                                         SpatialPCnum = 20)
  slideseq <- SpatialPCA_SpatialPCs(slideseq, fast = TRUE)
  end_time <- Sys.time()
  T_sparse1 <- end_time - start_time
})
clusterlabel <- louvain_clustering(clusternum = 8,
                                   latent_dat = as.matrix(slideseq@SpatialPCs),
                                   knearest = round(sqrt(dim(slideseq@SpatialPCs)[2])))
cbp_spatialpca <- c("lightyellow2",
                    "coral",
                    "lightcyan2",
                    "#66C2A5",
                    "cornflowerblue",
                    "#FFD92F",
                    "#E78AC3",
                    "skyblue1")
p <- plot_cluster(legend = "right",
                  location = slideseq@location,
                  clusterlabel,
                  pointsize = 1,
                  text_size = 20,
                  title_in = paste0("SpatialPCA"),
                  color_in = cbp_spatialpca)

meta_data <- tibble("spot" = colnames(slideseq@counts[, rownames(slideseq@location)]),
                    "cluster" = clusterlabel) %>%
  mutate(
    group = case_when(
      cluster == "1" ~ "Choroid plexus",
      cluster == "2" ~ "White matter",
      cluster == "3" ~ "Granule middle sublayer",
      cluster == "4" ~ "Granule inner sublayer",
      cluster == "5" ~ "Cerebellar nucleus",
      cluster == "6" ~ "Molecular layer",
      cluster == "7" ~ "Granule outer sublayer",
      cluster == "8" ~ "Purkinje layer"
    ),
    imagerow = slideseq@location[, 1],
    imagecol = slideseq@location[, 2]
  ) %>%
  tibble::column_to_rownames("spot") %>%
  select(-1)
spatial_data <- CreateSeuratObject(counts = sp_count,
                                   project = "SlideSeq",
                                   assay = "Spatial",
                                   min.cells = 3,
                                   min.features = 0,
                                   meta.data = meta_data)
spatial_data[['image']] = new(Class = "SlideSeq",
                              assay = "Spatial",
                              coordinates = meta_data %>% select(c("imagerow", "imagecol")))
SpatialDimPlot(spatial_data, group.by = "group", cols = "white") + scale_fill_manual(values = cbp_spatialpca)
saveRDS(spatial_data, file = "../SCST-article/Data/Spatial/Spatial_Slide-seq_cerebellum.rds")


### Slide-seqv2
sceasy::convertFormat("../SCST-article/Data/Spatial_Slide-seqv2_OB/GSM5173943_OB1_Slide19.h5ad",
                      from = "anndata", to = "seurat",
                      outFile = "../SCST-article/Data/Spatial_Slide-seqv2_OB/GSM5173943_OB1_Slide19.rds")
spatial_data <- readRDS("../SCST-article/Data/Spatial_Slide-seqv2_OB/GSM5173943_OB1_Slide19.rds")
spatial_data <- subset(spatial_data, layer != "unknown")
coordinates <- spatial_data@reductions$spatial@cell.embeddings %>%
  as.data.frame() %>%
  rename(imagerow = SPATIAL_1, imagecol = SPATIAL_2)
meta_data <- spatial_data@meta.data %>%
  select(layer) %>%
  rename(group = layer) %>%
  cbind(coordinates)
spatial_data <- CreateSeuratObject(counts = spatial_data@assays$RNA$counts,
                                   project = "SlideSeq",
                                   assay = "Spatial",
                                   min.cells = 20,
                                   min.features = 0,
                                   meta.data = meta_data)
spatial_data[['image']] = new(Class = "SlideSeq",
                              assay = "Spatial",
                              coordinates = coordinates)
SpatialDimPlot(spatial_data, group.by = "group", stroke = 0)
spatial_data <- spatial_data_preprocessing(spatial_data)
saveRDS(spatial_data, file = "../SCST-article/Data/Spatial/Spatial_Slide-seqv2_OB.rds")


# load("../SCST-article/Data/Spatial_Slide-seqv2_hippocampus/SlideseqV2_SpatialPCA_result.rdata")
# clusterlabel <- louvain_clustering(14,
#                                    latent_dat = as.matrix(SpatialPCA_result$SpatialPCs),
#                                    310)
# # spatial domain cluster label for each location
# cbp_spatialpca <- c("#FD7446",
#                     "#709AE1",
#                     "#31A354",
#                     "#9EDAE5",
#                     "#DE9ED6",
#                     "#BCBD22",
#                     "#CE6DBD",
#                     "#DADAEB",
#                     "yellow",
#                     "#FF9896",
#                     "#91D1C2",
#                     "#C7E9C0",
#                     "#6B6ECF",
#                     "#7B4173")
# cluster <- as.character(clusterlabel)
#
# load("../SCST-article/Data/Spatial_Slide-seqv2_hippocampus/Puck_200115_08_count_location.rdata")
# counts <- countmat[, rownames(SpatialPCA_result$location)]
# rm(countmat)
# meta_data <- tibble("spot" = colnames(counts),
#                     "cluster" = clusterlabel) %>%
#   mutate(
#     group = case_when(
#       cluster == "1" ~ "Layer 5",
#       cluster == "2" ~ "Hippocampus(slm)",
#       cluster == "3" ~ "Layer 6",
#       cluster == "4" ~ "Dentate gyrus",
#       cluster == "5" ~ "Layer 4",
#       cluster == "6" ~ "Hippocampus(so/sr)",
#       cluster == "7" ~ "Thalamus subregion1",
#       cluster == "8" ~ "Thalamus subregion3",
#       cluster == "9" ~ "CA3",
#       cluster == "10" ~ "CA1",
#       cluster == "11" ~ "Third ventricle",
#       cluster == "12" ~ "Thalamus subregion2",
#       cluster == "13" ~ "Hippocampus(so)",
#       cluster == "14" ~ "Corpus callosum"
#     ),
#     imagerow = SpatialPCA_result$location[, 1],
#     imagecol = SpatialPCA_result$location[, 2]
#   ) %>%
#   tibble::column_to_rownames("spot") %>%
#   select(-1)
# spatial_data <- CreateSeuratObject(counts = counts,
#                                    project = "SlideSeq",
#                                    assay = "Spatial",
#                                    min.cells = 3,
#                                    min.features = 0,
#                                    meta.data = meta_data)
# spatial_data[['image']] = new(Class = "SlideSeq",
#                               assay = "Spatial",
#                               coordinates = meta_data %>% select(c("imagerow", "imagecol")))
# spatial_data <- subset(spatial_data, group == "Dentate gyrus" |
#                          group == "CA3" |
#                          group == "CA1" |
#                          group == "Hippocampus(so/sr)" |
#                          group == "Hippocampus(slm)" |
#                          group == "Hippocampus(so)")
# SpatialDimPlot(spatial_data, group.by = "group", stroke = 0, pt.size.factor = 1) + scale_fill_manual(values = cbp_spatialpca)
# saveRDS(spatial_data, file = "../SCST-article/Data/Spatial/Spatial_Slide-seqv2_hippocampus.rds")
#
#
# fov <- SeuratObject::CreateCentroids(meta_data %>% filter(group == "Dentate gyrus" |
#                                                                 group == "CA3" |
#                                                                 group == "CA1") %>% select(c(2:3)))
# fov <- SeuratObject::CreateFOV(fov)
# spatial_data@images <- list(slideseq = fov)
# ImageDimPlot(spatial_data, fov = "slideseq", group.by = "group")
# ImageFeaturePlot(spatial_data, fov = "slideseq", features = c("Myl4", "Map4"))
#
# spatial_data <- spatial_data %>%
#   NormalizeData(normalization.method = "CLR", margin = 2) %>%
#   ScaleData()


### osmFISH
library(SeuratDisk)
loom_data <- Connect(filename = "../SCST-article/Data/Spatial_osmFISH_cortex/osmFISH_SScortex_mouse_all_cells.loom", mode = "r")
data <- as(t(loom_data[["matrix"]][,]), "dgCMatrix")
cell_name <- loom_data[["/col_attrs/CellID"]][]
gene_name <- loom_data[["/row_attrs/Gene"]][]
colnames(data) <- cell_name
rownames(data) <- gene_name

## meta_data
meta_data <- data.frame("imagerow" = max(loom_data[["/col_attrs/Y"]][]) - loom_data[["/col_attrs/Y"]][],
                        "imagecol" = loom_data[["/col_attrs/X"]][],
                        "group" = loom_data[["/col_attrs/Region"]][],
                        "cell_type" = loom_data[["/col_attrs/ClusterName"]][],
                        row.names = loom_data[["/col_attrs/CellID"]][])
meta_data <- meta_data[colnames(data), ]
spatial_data <- CreateSeuratObject(counts = data,
                                   project = "SlideSeq",
                                   assay = "Spatial",
                                   min.cells = 20,
                                   min.features = 0,
                                   meta.data = meta_data)
spatial_data[['image']] = new(Class = "SlideSeq",
                              assay = "Spatial",
                              coordinates = meta_data %>% select(c("imagerow", "imagecol")))
spatial_data <- subset(spatial_data, group != "Excluded" &
                         group != "Hippocampus" &
                         group != "Ventricle" &
                         group != "White matter" &
                         group != "Internal Capsule Caudoputamen")
spatial_data <- spatial_data_preprocessing(spatial_data)
SpatialDimPlot(spatial_data, group.by = "group")
saveRDS(spatial_data, file = "../SCST-article/Data/Spatial/Spatial_osmFISH_cortex.rds")


#### STARmap
sceasy::convertFormat("../SCST-article/Data/Spatial_STARmap_cortex/STARmap_20180505_BY3_1k.h5ad",
                      from = "anndata", to = "seurat",
                      outFile = "../SCST-article/Data/Spatial_STARmap_cortex/Spatial_STARmap_cortex.rds")
spatial_data <- readRDS("../SCST-article/Data/Spatial_STARmap_cortex/Spatial_STARmap_cortex.rds")
coordinates <- spatial_data@reductions$spatial@cell.embeddings %>%
  as.data.frame() %>%
  rename(imagerow = SPATIAL_1, imagecol = SPATIAL_2)
meta_data <- spatial_data@meta.data %>%
  select(label) %>%
  rename(group = label) %>%
  cbind(coordinates)
spatial_data <- CreateSeuratObject(counts = spatial_data@assays$RNA$counts,
                                   project = "SlideSeq",
                                   assay = "Spatial",
                                   min.cells = 20,
                                   min.features = 0,
                                   meta.data = meta_data)
spatial_data[['image']] = new(Class = "SlideSeq",
                              assay = "Spatial",
                              coordinates = coordinates)
spatial_data <- spatial_data_preprocessing(spatial_data)
SpatialDimPlot(spatial_data, group.by = "group", stroke = 0)
saveRDS(spatial_data, file = "../SCST-article/Data/Spatial/Spatial_STARmap_cortex.rds")



#### Stereo-seq
sceasy::convertFormat("../SCST-article/Data/Spatial_Stereo-seq_OB/Spatial_Stereo-seq_OB.h5ad",
                      from = "anndata", to = "seurat",
                      outFile = "../SCST-article/Data/Spatial_Stereo-seq_OB/Spatial_Stereo-seq_OB.rds")
spatial_data <- readRDS("../SCST-article/Data/Spatial_Stereo-seq_OB/Spatial_Stereo-seq_OB.rds")
coordinates <- spatial_data@reductions$spatial@cell.embeddings %>%
  as.data.frame() %>%
  rename(imagerow = SPATIAL_1, imagecol = SPATIAL_2)
spatial_data <- CreateSeuratObject(counts = spatial_data@assays$RNA$counts,
                                   project = "SlideSeq",
                                   assay = "Spatial",
                                   min.cells = 10,
                                   min.features = 200,
                                   meta.data = coordinates)
spatial_data[['image']] = new(Class = "SlideSeq",
                              assay = "Spatial",
                              coordinates = coordinates[colnames(spatial_data), ])
spatial_data <- spatial_data %>%
  SCTransform(assay = "Spatial", variable.features.n = 5000) %>%
  RunPCA(npcs = 20, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:20, seed.use = 111) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.8, random.seed = 111)
SpatialDimPlot(spatial_data, stroke = 0)
DotPlot(spatial_data, features = c("Mbp", "Nrgn", "Pcp4", "Gabra1", "Spp1", "Slc6a11", "Cck", "Apod", "Ptn"))
new_ident_ids <- c("EPL", "GCL", "LPL", "ONL", "MCL", "GL", "RMS", "LPL")
names(new_ident_ids) <- levels(spatial_data)
spatial_data <- RenameIdents(spatial_data, new_ident_ids)
SpatialDimPlot(spatial_data, stroke = 0)
spatial_data$group <- unname(S4Vectors::unfactor(Idents(spatial_data)))
saveRDS(spatial_data, file = "../SCST-article/Data/Spatial/Spatial_Stereo-seq_OB.rds")

#
# conda_env <- "/Users/duohongrui/opt/miniconda3/envs/spatial"
# reticulate::use_condaenv(conda_env, required = TRUE)
# packages <- reticulate::py_list_packages(envname = conda_env)
# if(use_cuda){
#   device <- "cuda"
#   use_gpu <- reticulate::r_to_py(TRUE)
# }else{
#   device <- "cpu"
#   use_gpu <- reticulate::r_to_py(FALSE)
# }
#
# sc <-  reticulate::import("scanpy", convert = FALSE)
# scipy <-  reticulate::import("scipy", convert = FALSE)
# np <- reticulate::import("numpy", convert = FALSE)
# if(packageVersion("Seurat") >= "5.0"){
#   X <- SeuratObject::LayerData(seurat, layer = "counts")
# }else{
#   X <- methods::slot(Seurat::GetAssay(seurat), "counts")
# }
# ann_cell_meta <- seurat@meta.data %>%
#   dplyr::rename("arraycol" = "imagecol",
#                 "arrayrow" = "imagerow")
# spatial_coordiantes <- reticulate::r_to_py(list("spatial" = ann_cell_meta %>%
#                                                   dplyr::select("arraycol", "arrayrow") %>%
#                                                   as.matrix()))
# adata = sc$AnnData(
#   X   = scipy$sparse$csr_matrix(Matrix::t(X)),
#   obs = ann_cell_meta,
#   obsm = spatial_coordiantes,
#   dtype = np$float32
# )
# sc$pp$highly_variable_genes(adata, n_top_genes = as.integer(3000))
# sc$pp$normalize_total(adata, target_sum = 1e4)
# sc$pp$log1p(adata)
# STAGATE_path <- system.file("STAGATE", package = "SCST")
# stagate <- reticulate::import_from_path("STAGATE_pyG", path = STAGATE_path, convert = FALSE)
# stagate$Cal_Spatial_Net(adata, verbose = reticulate::r_to_py(TRUE), model = "KNN", k_cutoff = as.integer(50))
# stagate$Stats_Spatial_Net(adata)
# adata = stagate$train_STAGATE(adata,
#                               n_epochs = as.integer(1000),
#                               lr = 0.001,
#                               verbose = reticulate::r_to_py(TRUE),
#                               random_seed = as.integer(888),
#                               device = device)
# adata$write_h5ad("../STAGATE.h5ad")
# sc$pp$neighbors(adata, use_rep = 'STAGATE')
# sc$tl$umap(adata)
# sc$tl$louvain(adata, resolution = 0.8)
# sc$pl$embedding(adata, basis="spatial", color="louvain", title='STAGATE')
# # sc$pp$neighbors(adata, use_rep = 'STAGATE', key_added = "STAGATE_neighbors")
# # adata = stagate$mclust_R(adata, used_obsm = 'STAGATE', num_cluster = n_domains)
# STAGATE_obsm <- reticulate::py_to_r(adata$obsm["STAGATE"])
# require("mclust")
# model_names <- c("EEE", "EEI", "EII", "EVI", "VEI", "VII", "VVI")
# for(i in model_names){
#   res <- mclust::Mclust(data = STAGATE_obsm, G = 7, modelNames = i)
#   if(!is.null(res)){
#     print(i)
#     break
#   }
# }
#
# rownames(STAGATE_obsm) <- reticulate::py_to_r(adata$obs$index$values)
# colnames(STAGATE_obsm) <- paste0("scANVI", "_", 1:ncol(STAGATE_obsm))
# reduction <- Seurat::CreateDimReducObject(
#   embeddings = STAGATE_obsm,
#   assay = "Spatial",
#   key = "STAGATE_"
# )
# spatial_data@reductions <- append(spatial_data@reductions, list("STAGATE" = reduction))
# spatial_data <- spatial_data %>%
#   SCTransform(assay = "Spatial", variable.features.n = 3000) %>%
#   RunPCA(npcs = 20, verbose = FALSE) %>%
#   RunUMAP(reduction = "STAGATE", dims = 1:20, seed.use = 111) %>%
#   FindNeighbors(reduction = "STAGATE", dims = 1:20, graph.name = c("NN", "SNN")) %>%
#   FindClusters(resolution = 0.4, random.seed = 111, graph.name = "NN")
#
# # STAGATE_result <- res$classification
# STAGATE_result <- adata$obs["louvain"]$astype('str')$tolist() %>% reticulate::py_to_r()
# spatial_data$"STAGATE_result" <- STAGATE_result
# Idents(spatial_data) <- STAGATE_result
# SpatialDimPlot(spatial_data, stroke = 0)
# DotPlot(spatial_data, features = c("Mbp", "Nrgn", "Pcp4", "Gabra1", "Spp1", "Slc6a11", "Cck", "Apod", "Ptn"))
# new_ident_ids <- c("EPL", "GCL", "LPL", "ONL", "MCL", "GL", "RMS", "LPL")
# names(new_ident_ids) <- levels(spatial_data)
# spatial_data <- RenameIdents(spatial_data, new_ident_ids)
# SpatialDimPlot(spatial_data, stroke = 0)


###################------------------------------------------###################
##################         Data with multiple groups/clusters
###################------------------------------------------###################

### Group_10X_Lung
file.list <- file.path("../SCST-article/Data/Group_10X_Lung",
                       list.files("../SCST-article/Data/Group_10X_Lung/"))
names(file.list) <- c("10X_P7_0", "10X_P7_1", "10X_P7_8", "10X_P7_9")
data <- Read10X(data.dir = file.list, strip.suffix = TRUE)
annotation <- read.csv("../SCST-article/Data/Annotation/annotations_droplet.csv")
intersect_index <- intersect(colnames(data), annotation$cell)
metadata <- annotation %>%
  filter(cell %in% intersect_index) %>%
  select(cell, cell_ontology_class) %>%
  rename(cell_type = cell_ontology_class) %>%
  tibble::column_to_rownames("cell") %>%
  filter(cell_type != "")
reserved_cell_type <- names(table(metadata$cell_type))[table(metadata$cell_type) > 50]
metadata <- metadata %>%
  filter(cell_type %in% reserved_cell_type)
data <- data[, rownames(metadata)]
data <- preprocessing_function(matrix = data,
                               meta_data = metadata,
                               project = "Group_10X_lung",
                               species = "mouse")
saveRDS(data, file = "../SCST-article/Data/Group/Group_10X_lung.rds")



### Group_10X_Marrow
file.list <- file.path("../SCST-article/Data/Group_10X_Marrow",
                       list.files("../SCST-article/Data/Group_10X_Marrow/"))
names(file.list) <- c("10X_P7_2", "10X_P7_3")
data <- Read10X(data.dir = file.list, strip.suffix = TRUE)
annotation <- read.csv("../SCST-article/Data/Annotation/annotations_droplet.csv")
intersect_index <- intersect(colnames(data), annotation$cell)
metadata <- annotation %>%
  filter(cell %in% intersect_index) %>%
  select(cell, cell_ontology_class) %>%
  rename(cell_type = cell_ontology_class) %>%
  tibble::column_to_rownames("cell") %>%
  filter(cell_type != "")
data <- data[, rownames(metadata)]
data <- preprocessing_function(matrix = data,
                               meta_data = metadata,
                               project = "Group_10X_marrow",
                               species = "mouse")
saveRDS(data, file = "../SCST-article/Data/Group/Group_10X_marrow.rds")


### Group_Smart-seq2_Trachea
data <- read.csv("../SCST-article/Data/Group_Smart-seq2_Trachea/Trachea-counts.csv", row.names = 1)
annotation <- read.csv("../SCST-article/Data/Annotation/annotations_FACS.csv") %>%
  filter(tissue == "Trachea")
intersect_index <- intersect(colnames(data), annotation$cell)
metadata <- annotation %>%
  filter(cell %in% intersect_index) %>%
  select(cell, cell_ontology_class) %>%
  rename(cell_type = cell_ontology_class) %>%
  tibble::column_to_rownames("cell") %>%
  filter(cell_type != "")
data <- data[, rownames(metadata)]
data <- preprocessing_function(matrix = data,
                               meta_data = metadata,
                               project = "Group_Smart-seq2_trachea",
                               species = "mouse")
saveRDS(data, file = "../SCST-article/Data/Group/Group_Smart-seq2_trachea.rds")


### Group_Smart-seq2_Tongue
data <- read.csv("../SCST-article/Data/Group_Smart-seq2_Tongue/Tongue-counts.csv", row.names = 1)
annotation <- read.csv("../SCST-article/Data/Annotation/annotations_FACS.csv") %>%
  filter(tissue == "Tongue")
intersect_index <- intersect(colnames(data), annotation$cell)
metadata <- annotation %>%
  filter(cell %in% intersect_index) %>%
  select(cell, cell_ontology_class) %>%
  rename(cell_type = cell_ontology_class) %>%
  tibble::column_to_rownames("cell") %>%
  filter(cell_type != "")
data <- data[, rownames(metadata)]
data <- preprocessing_function(matrix = data,
                               meta_data = metadata,
                               project = "Group_Smart-seq2_tongue",
                               species = "mouse")
saveRDS(data, file = "../SCST-article/Data/Group/Group_Smart-seq2_tongue.rds")


### Group_Smart-seq
data <- readRDS("../SCST-article/Data/Group_Smart-seq/data79_aging-hsc-old_kowalczyk.rds")
counts <- as(t(data$counts), "dgCMatrix")
metadata <- data.frame("cell_type" = data[["grouping"]],
                       row.names = data[["cell_ids"]])
data <- preprocessing_function(matrix = counts,
                               meta_data = metadata,
                               project = "Group_Smart-seq",
                               species = "mouse")
saveRDS(data, file = "../SCST-article/Data/Group/Group_Smart-seq.rds")



### Group_CEL-seq2
data <- readRDS("../SCST-article/Data/Group_CEL-seq2/data82_cellbench-SC1_luyitian.rds")
counts <- as(t(data$counts), "dgCMatrix")
metadata <- data.frame("cell_type" = data[["grouping"]],
                       row.names = data[["cell_ids"]])
data <- preprocessing_function(matrix = counts,
                               meta_data = metadata,
                               project = "Group_CEL-seq2",
                               species = "mouse")
saveRDS(data, file = "../SCST-article/Data/Group/Group_CEL-seq2.rds")

