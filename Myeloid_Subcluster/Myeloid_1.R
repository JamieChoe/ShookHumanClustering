library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(SeuratObject)

#Merged_Donors <- readRDS("/home/gabriel.batzli/jamie_project/Merged_Donors_Object.rds")

#Subsetting from Merged Donors to isolate Myeloid cluster
MyeloidCluster <- subset(x = Merged_Donors, idents = c("Myeloid"))

MyeloidCluster <- SplitObject(MyeloidCluster, split.by = "day")

#Using SCTransform
MyeloidCluster <- lapply(MyeloidCluster, SCTransform, vst.flavor = "v2")
message("SCTransform Finished!")

#Integration
features <- SelectIntegrationFeatures(object.list = MyeloidCluster, nfeatures = 3000)
MyeloidCluster <- PrepSCTIntegration(object.list = MyeloidCluster, anchor.features = features)
message("Integration prepped, finding anchors now.")

#Anchors
DonorAnchors <- FindIntegrationAnchors(object.list = MyeloidCluster, normalization.method = "SCT", anchor.features = features)
#saveRDS(DonorAnchors, "/home/gabriel.batzli/jamie_project/Myeloid_Anchors.rds")
message("Anchors made and saved")

MyeloidCluster <- IntegrateData(anchorset = DonorAnchors, normalization.method = "SCT")
message("Integration complete.")


#Run PCA to find PCs
MyeloidIntegrated <- RunPCA(MyeloidCluster)
#saveRDS(MyeloidCluster, "/home/gabriel.batzli/jamie_project/Myeloid_Cluster_Object.rds")
message("Finished running PCA")

#Creating Elbowplot
ElbowPlot(MyeloidIntegrated)



































