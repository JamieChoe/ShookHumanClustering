

install.packages("openxlsx", dependencies = TRUE)
library(openxlsx)
library(writexl)

MyeloidCluster <- readRDS("/home/gabriel.batzli/jamie_project/Myeloid_Cluster_Object.rds")


# Process the data

MyeloidIntegrated <- RunUMAP(MyeloidIntegrated, reduction = "pca", dims = 1:11)

MyeloidIntegrated <- FindNeighbors(MyeloidIntegrated, reduction = "pca", dims = 1:11)


#Resolutions 
MyeloidIntegrated <- FindClusters(MyeloidIntegrated, resolution = 0.5)
MyeloidIntegrated <- FindClusters(MyeloidIntegrated, resolution = 0.4)
MyeloidIntegrated <- FindClusters(MyeloidIntegrated, resolution = 0.3)
MyeloidIntegrated <- FindClusters(MyeloidIntegrated, resolution = 0.2)
MyeloidIntegrated <- FindClusters(MyeloidIntegrated, resolution = 0.1)
MyeloidIntegrated <- FindClusters(MyeloidIntegrated, resolution = 0.05)
message("FindClusters complete. Beginning UMAPs")


# Umaps
dim_plots <- list(
  DimPlot_clusters0.5 = DimPlot(MyeloidIntegrated, reduction = "umap", group.by = "integrated_snn_res.0.5") + ggtitle("Resolution 0.5"),
  DimPlot_clusters0.1 = DimPlot(MyeloidIntegrated, reduction = "umap", group.by = "integrated_snn_res.0.1") + ggtitle("Resolution 0.1"),
  DimPlot_clusters0.05 = DimPlot(MyeloidIntegrated, reduction = "umap", group.by = "integrated_snn_res.0.05") + ggtitle("Resolution 0.05"),
  DimPlot_clusters0.3 = DimPlot(MyeloidIntegrated, reduction = "umap", group.by = "integrated_snn_res.0.3") + ggtitle("Resolution 0.3"),
  DimPlot_clusters0.2 = DimPlot(MyeloidIntegrated, reduction = "umap", group.by = "integrated_snn_res.0.2") + ggtitle("Resolution 0.2")
)

DimPlot_clusters0.3 = DimPlot(MyeloidIntegrated, reduction = "umap", group.by = "integrated_snn_res.0.3") + ggtitle("Resolution 0.3")
DimPlot_clusters0.4 = DimPlot(MyeloidIntegrated, reduction = "umap", group.by = "integrated_snn_res.0.4") + ggtitle("Resolution 0.4")

DimPlot_clusters0.5 = DimPlot(MyeloidIntegrated, reduction = "umap", group.by = "integrated_snn_res.0.5") + ggtitle("Resolution 0.5")
DimPlot_clusters0.1 = DimPlot(MyeloidIntegrated, reduction = "umap", group.by = "integrated_snn_res.0.1") + ggtitle("Resolution 0.1")
DimPlot_clusters0.05 = DimPlot(MyeloidIntegrated, reduction = "umap", group.by = "integrated_snn_res.0.05") + ggtitle("Resolution 0.05")
DimPlot_clusters0.2 = DimPlot(MyeloidIntegrated, reduction = "umap", group.by = "integrated_snn_res.0.2") + ggtitle("Resolution 0.2")



DimPlot_clusters0.4
DimPlot_clusters0.5 
DimPlot_clusters0.1 
DimPlot_clusters0.05 
DimPlot_clusters0.3 
DimPlot_clusters0.2 

DimPlot_clusters0.2 = DimPlot(MyeloidIntegrated, reduction = "umap", group.by = "integrated_snn_res.0.2") + ggtitle("Resolution 0.2")


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

MyeloidIntegrated <- PrepSCTFindMarkers(MyeloidIntegrated)

#FindAllMarkers
Markers0.5 <- FindAllMarkers(MyeloidIntegrated, assay = "SCT", min.pct = 0.25, logfc.threshold = 0.25, group.by = "integrated_snn_res.0.5")
write_xlsx(Markers0.5,"/home/gabriel.batzli/jamie_project/Markers_0.5.xlsx",col_names = TRUE,format_headers = TRUE)

Markers0.1 <- FindAllMarkers(MyeloidIntegrated, assay = "SCT", min.pct = 0.25, logfc.threshold = 0.25, group.by = "integrated_snn_res.0.1")
write_xlsx(Markers0.1,"/home/gabriel.batzli/jamie_project/Markers_0.1.xlsx",col_names = TRUE,format_headers = TRUE)

Markers0.05 <- FindAllMarkers(MyeloidIntegrated, assay = "SCT", min.pct = 0.25, logfc.threshold = 0.25, group.by = "integrated_snn_res.0.05")
write_xlsx(Markers0.05,"Markers0.05Macro.xlsx",col_names = TRUE,format_headers = TRUE)
message("Find All Markers Finished and Saved")


Markers0.2 <- FindAllMarkers(MyeloidIntegrated, assay = "SCT", min.pct = 0.25, logfc.threshold = 0.25, group.by = "integrated_snn_res.0.2")

write_xlsx(Markers0.2,"Markers_0.2.xlsx",col_names = TRUE,format_headers = TRUE)

saveRDS(MyeloidIntegrated, "Myeloid_Integrated.rds")


Markers0.4 <- FindAllMarkers(MyeloidIntegrated, assay = "SCT", min.pct = 0.25, logfc.threshold = 0.25, group.by = "integrated_snn_res.0.4")
write_xlsx(Markers0.4,"Markers_0.4.xlsx",col_names = TRUE,format_headers = TRUE)

SaveSeuratRds(MyeloidIntegrated, "MyeloidIntegrated.rds")




















