library(openxlsx) 
library(Seurat)
library(tidyverse)
library(patchwork)
library(sctransform)
library(ggplot2)
message("libs loaded")

Merged_Donors <- readRDS("/home/gabriel.batzli/jamie_project/Merged_Donors_Object.rds")

# Process the data

Merged_Donors <- RunUMAP(Merged_Donors, reduction = "pca", dims = 1:11)

Merged_Donors <- FindNeighbors(Merged_Donors, reduction = "pca", dims = 1:11)

#Resolutions 
Merged_Donors <- FindClusters(Merged_Donors, resolution = 0.5)
Merged_Donors <- FindClusters(Merged_Donors, resolution = 0.1)
Merged_Donors <- FindClusters(Merged_Donors, resolution = 0.05)


# Umaps
dim_plots <- list(
  DimPlot_clusters0.5 = DimPlot(Merged_Donors, reduction = "umap", group.by = "integrated_snn_res.0.5") + ggtitle("Resolution 0.5"),
  DimPlot_clusters0.1 = DimPlot(Merged_Donors, reduction = "umap", group.by = "integrated_snn_res.0.1") + ggtitle("Resolution 0.1"),
  DimPlot_clusters0.05 = DimPlot(Merged_Donors, reduction = "umap", group.by = "integrated_snn_res.0.05") + ggtitle("Resolution 0.05")
)

message("UMAPs created!")

#Saving the UMAPs

#PDF
for (plot_name in names(dim_plots)) {
  ggsave(paste0("/home/gabriel.batzli/jamie_project/", plot_name, ".pdf"), 
         plot = dim_plots[[plot_name]], 
         width = 8, height = 6)}

#PNG
for (plot_name in names(dim_plots)) {
  ggsave(paste0("/home/gabriel.batzli/jamie_project/", plot_name, ".png"), 
         plot = dim_plots[[plot_name]], 
         width = 8, height = 6)}

message("UMAPs saved")

Merged_Donors <- PrepSCTFindMarkers(Merged_Donors)

#FindAllMarkers
Markers0.5 <- FindAllMarkers(Merged_Donors, assay = "SCT", min.pct = 0.25, logfc.threshold = 0.25, group.by = "integrated_snn_res.0.5")
write.xlsx(Markers0.5,"/home/gabriel.batzli/jamie_project/Markers_0.5.xlsx")

Markers0.1 <- FindAllMarkers(Merged_Donors, assay = "SCT", min.pct = 0.25, logfc.threshold = 0.25, group.by = "integrated_snn_res.0.1")
write.xlsx(Markers0.1,"/home/gabriel.batzli/jamie_project/Markers_0.1.xlsx")

Markers0.05 <- FindAllMarkers(Merged_Donors, assay = "SCT", min.pct = 0.25, logfc.threshold = 0.25, group.by = "integrated_snn_res.0.05")
write.xlsx(Markers0.05,"/home/gabriel.batzli/jamie_project/Markers_0.05.xlsx")
message("Find All Markers Finished and Saved")


#Save Object
saveRDS(Merged_Donors, "/home/gabriel.batzli/jamie_project/Merged_Donors_Object_2.rds")

message("obj saved")
