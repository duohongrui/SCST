.make_trees <- function(ref_data,
                        group = NULL,
                        is_Newick = TRUE,
                        is_parenthetic = FALSE,
                        return_group = FALSE){
  if(!requireNamespace("ctc", quietly = TRUE)){
    message("Install ctc...")
    BiocManager::install("ctc")
  }
  if(!requireNamespace("ape", quietly = TRUE)){
    message("Install ape...")
    utils::install.packages('ape')
  }
  data <- Seurat::CreateSeuratObject(counts = ref_data)
  data <- Seurat::NormalizeData(data,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000,
                                verbose=FALSE)
  data <- Seurat::FindVariableFeatures(data,
                                       selection.method = "vst",
                                       nfeatures = 2000,
                                       verbose = FALSE)
  all.genes <- rownames(data)
  data <- Seurat::ScaleData(data, features = all.genes, verbose=FALSE)
  if(is.null(group)){
    data <- Seurat::RunPCA(data, features = Seurat::VariableFeatures(object = data), verbose = FALSE)
    data <- Seurat::FindNeighbors(data, dims = 1:10)
    data <- Seurat::FindClusters(data, resolution = 0.5, verbose = FALSE)
    group_num <- length(unique(data@meta.data[["seurat_clusters"]]))
    cat(paste0("Your data has ", group_num, " groups \n"))
    group <- paste0("group", as.numeric(data@meta.data[["seurat_clusters"]]))
  }
  data@meta.data$'group' <- stringr::str_remove_all(group, pattern = "[(].+[)]")
  #Get state tree by hierachical clustering on the state means
  if(utils::packageVersion("Seurat") >= "5.0"){
    exp_data <- Seurat::AggregateExpression(data, group.by = 'group')
  }else{
    exp_data <- Seurat::AverageExpression(data, slot = 'data', group.by = 'group')
  }
  clu <- stats::hclust(stats::dist(t(as.matrix(exp_data$RNA))), method = 'ward.D')
  for(i in 1:length(clu[["labels"]])){
    clu[["labels"]][i] <- stringr::str_replace_all(clu[["labels"]][i], ",", "_")
  }
  phyla <- ctc::hc2Newick(clu, flat=TRUE)
  if(is_Newick){
    if(return_group) return(list(phyla = phyla, group = group)) else return(phyla)
  }
  phyla <- ape::read.tree(text = phyla)
  if(is_parenthetic){
    if(return_group) return(list(phyla = phyla <- list(phyla), group = group)) else return(phyla <- list(phyla))
  }
  phyla$edge.length <- ceiling(phyla$edge.length)
  phyla$edge.length[phyla$edge.length == 0] <- 1
  if(return_group) return(list(phyla = phyla, group = group)) else return(phyla)
}
