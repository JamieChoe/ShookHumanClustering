library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(SeuratObject)

MacroCluster <- subset(x = MyeloidIntegrated, idents = c("Mac 1","Mac 2"))

MacroCluster <- SplitObject(MacroCluster, split.by = "day")

#Using SCTransform on each object
MacroCluster <- lapply(MacroCluster, SCTransform, vst.flavor = "v2")
message("SCTransform Finished!")

#Integration
features <- SelectIntegrationFeatures(MacroCluster, nfeatures = 3000)
MacroCluster <- PrepSCTIntegration(MacroCluster, anchor.features = features)
message("Integration prepped, finding anchors now.")

#Anchors
DonorAnchors <- FindIntegrationAnchors(MacroCluster, normalization.method = "SCT", anchor.features = features)
#saveRDS(DonorAnchors, "/directory")
message("Anchors made and saved")

MacroCluster <- IntegrateData(anchorset = DonorAnchors, normalization.method = "SCT")
message("Integration complete.")

#Run PCA to find PCs
MacroCluster <- RunPCA(MacroCluster)
#saveRDS(MacroCluster, "/directory")
message("Finished running PCA")

#Creating Elbowplot
ElbowPlot(MacroCluster)
#ggsave("/directory", plot = EP, width = 8, height = 6)
message("Elbowplot created and saved")

SaveSeuratRds(MacroCluster)
































