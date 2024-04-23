# VisualizeEvalReults <- function(seurat,
#                                 aspect = "cluster",
#                                 method = "all",
#                                 reduction = "umap",
#                                 alpha = 0.8,
#                                 colors = NULL,
#                                 ncol = 3,
#                                 legend_position = NULL){
#
# }
#
# seurat <- Batch_removal_result$Splat$seurat
#
#
# reductions <- Seurat::Reductions(seurat)
# reductions <- reductions[-grep("pca", reductions)]
#
#
# vis_list <- lapply(reductions[-c(2)], FUN = function(x){
#   Seurat::DimPlot(seurat, reduction = x, group.by = "batch", alpha = 0.7) +
#     ggplot2::theme_test() +
#     ggplot2::theme(
#       axis.text = ggplot2::element_blank(),
#       axis.ticks = ggplot2::element_blank(),
#       panel.border = ggplot2::element_rect(color = "black"),
#       plot.title = ggplot2::element_text(hjust = 0.5),
#       legend.position = "right"
#     ) +
#     ggplot2::xlab("UMAP 1") +
#     ggplot2::ylab("UMAP 2") +
#     ggplot2::ggtitle(ifelse(x == "umap" | x == "tsne", "Raw", x)) +
#     ggplot2::scale_color_manual(values = plot_colors())
# })
#
# vis_list <- patchwork::wrap_plots(vis_list, ncol = 3) +
#   patchwork::guide_area() +
#   patchwork::plot_layout(guides = "collect")
# vis_list
#
# Seurat::DimPlot(seurat, reduction = "scVI", group.by = "batch")
# Seurat::DimPlot(seurat, reduction = "umap", group.by = "batch")
#
#
# plot_colors <- function(){
#   colors <- c("#ec671e", "#2e8d3a", "#1d5e9b", "#774d9c", "#29a8b9", "#9cb5db", "#f4a965",
#               "#ee7e80", "#b38880", "#f5c290", "#838384", "#cdce79", "#be1a22", "#f4c33f",
#               "#37b380", "#4c98c8", "#36546d", "#b98519", "#dc4530", "#c5e6eb", "#5d5098",
#               "#edec73", "#f18a1c", "#c27cad", "#ee7634", "#661115", "#bcd531", "#dbbc80",
#               "#f6bbbd")
#   return(colors)
# }
