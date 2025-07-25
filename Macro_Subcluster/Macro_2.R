library(openxlsx)
library(openxlsx) 
library(Seurat)
library(tidyverse)
library(patchwork)
library(sctransform)
library(ggplot2)
message("libs loaded")

# Process the data

MacroCluster <- RunUMAP(MacroCluster, reduction = "pca", dims = 1:11)

MacroCluster <- FindNeighbors(MacroCluster, reduction = "pca", dims = 1:11)

#Resolutions 
MacroCluster <- FindClusters(MacroCluster, resolution = 0.5)
MacroCluster <- FindClusters(MacroCluster, resolution = 0.1)
MacroCluster <- FindClusters(MacroCluster, resolution = 0.2)
MacroCluster <- FindClusters(MacroCluster, resolution = 0.3)
MacroCluster <- FindClusters(MacroCluster, resolution = 0.05)


# Umaps
dim_plots <- list(
  DimPlot_clusters0.5 = DimPlot(MacroCluster, reduction = "umap", group.by = "integrated_snn_res.0.5") + ggtitle("Resolution 0.5"),
  DimPlot_clusters0.1 = DimPlot(MacroCluster, reduction = "umap", group.by = "integrated_snn_res.0.1") + ggtitle("Resolution 0.1"),
  DimPlot_clusters0.05 = DimPlot(MacroCluster, reduction = "umap", group.by = "integrated_snn_res.0.05") + ggtitle("Resolution 0.05"),
  DimPlot_clusters0.2 = DimPlot(MacroCluster, reduction = "umap", group.by = "integrated_snn_res.0.2") + ggtitle("Resolution 0.2"),
  DimPlot_clusters0.3 = DimPlot(MacroCluster, reduction = "umap", group.by = "integrated_snn_res.0.3") + ggtitle("Resolution 0.3")
)

DimPlot_clusters0.5 = DimPlot(MacroCluster, reduction = "umap", group.by = "integrated_snn_res.0.5") + ggtitle("Resolution 0.5")
DimPlot_clusters0.1 = DimPlot(MacroCluster, reduction = "umap", group.by = "integrated_snn_res.0.1") + ggtitle("Resolution 0.1")
DimPlot_clusters0.05 = DimPlot(MacroCluster, reduction = "umap", group.by = "integrated_snn_res.0.05") + ggtitle("Resolution 0.05")
DimPlot_clusters0.2 = DimPlot(MacroCluster, reduction = "umap", group.by = "integrated_snn_res.0.2") + ggtitle("Resolution 0.2")
DimPlot_clusters0.3 = DimPlot(MacroCluster, reduction = "umap", group.by = "integrated_snn_res.0.3") + ggtitle("Resolution 0.3")


DimPlot_clusters0.3
DimPlot_clusters0.5
DimPlot_clusters0.1
DimPlot_clusters0.05
message("UMAPs created!")

#Saving the UMAPs

#PDF
for (plot_name in names(dim_plots)) {
  ggsave(paste0(plot_name, ".pdf"), 
         plot = dim_plots[[plot_name]], 
         width = 8, height = 6)}

#PNG
for (plot_name in names(dim_plots)) {
  ggsave(paste0(plot_name, ".png"), 
         plot = dim_plots[[plot_name]], 
         width = 8, height = 6)}

message("UMAPs saved")

MacroCluster <- PrepSCTFindMarkers(MacroCluster)

#FindAllMarkers
Markers0.5 <- FindAllMarkers(MacroCluster, assay = "SCT", min.pct = 0.25, logfc.threshold = 0.25, group.by = "integrated_snn_res.0.5")
write.xlsx(Markers0.5,"/home/gabriel.batzli/jamie_project/Markers_0.5.xlsx")

Markers0.1 <- FindAllMarkers(MacroCluster, assay = "SCT", min.pct = 0.25, logfc.threshold = 0.25, group.by = "integrated_snn_res.0.1")
write.xlsx(Markers0.1,"/home/gabriel.batzli/jamie_project/Markers_0.1.xlsx")

Markers0.05 <- FindAllMarkers(MacroCluster, assay = "SCT", min.pct = 0.25, logfc.threshold = 0.25, group.by = "integrated_snn_res.0.05")
write.xlsx(Markers0.05,"/home/gabriel.batzli/jamie_project/Markers_0.05.xlsx")
message("Find All Markers Finished and Saved")


Markers0.3 <- FindAllMarkers(MacroCluster, assay = "SCT", min.pct = 0.25, logfc.threshold = 0.25, group.by = "integrated_snn_res.0.3", only.pos = TRUE)
write.xlsx(Markers0.3,"MacroMarkers_0.3.xlsx")


#Save Object
saveRDS(MacroCluster,"MacroCluster.rds")

message("obj saved")

























