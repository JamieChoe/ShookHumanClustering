library(Seurat)
library(tidyverse)
library(patchwork)
library(sctransform)
library(ggplot2)
message("libs loaded")


options(future.globals.maxSize = 16 * 1024^3)  # 16 GB

# Load the PBMC datasets
Donor5_D0.data <- Read10X(data.dir = "/home/gabriel.batzli/jamie_project/GSE241132_RAW/PWH26D0")
Donor5_D1.data <- Read10X(data.dir = "/home/gabriel.batzli/jamie_project/GSE241132_RAW/PWH26D1")
Donor4_D0.data <- Read10X(data.dir = "/home/gabriel.batzli/jamie_project/GSE241132_RAW/PWH27D0")
Donor4_D1.data <- Read10X(data.dir = "/home/gabriel.batzli/jamie_project/GSE241132_RAW/PWH27D1")
Donor3_D0.data <- Read10X(data.dir = "/home/gabriel.batzli/jamie_project/GSE241132_RAW/PWH28D0")
Donor3_D1.data <- Read10X(data.dir = "/home/gabriel.batzli/jamie_project/GSE241132_RAW/PWH28D1")

# Create objects. Filter out cells expressing <500 genes, and genes expressed in < 10 cells. 
Donor5_D0 <- CreateSeuratObject(counts = Donor5_D0.data, project = "Donor5_Day0", min.cells = 10, min.features = 500)
Donor5_D1 <- CreateSeuratObject(counts = Donor5_D1.data, project = "Donor5_Day1",min.cells = 10, min.features = 500)
Donor4_D0 <- CreateSeuratObject(counts = Donor4_D0.data, project = "Donor4_Day0",min.cells = 10, min.features = 500)
Donor4_D1 <- CreateSeuratObject(counts = Donor4_D1.data, project = "Donor4_Day1",min.cells = 10, min.features = 500)
Donor3_D0 <- CreateSeuratObject(counts = Donor3_D0.data, project = "Donor3_Day0",min.cells = 10, min.features = 500)
Donor3_D1 <- CreateSeuratObject(counts = Donor3_D1.data, project = "Donor3_Day1",min.cells = 10, min.features = 500)

Donor5_D0[["percent.mt"]] <- PercentageFeatureSet(Donor5_D0, pattern = "^MT-")
Donor5_D1[["percent.mt"]] <- PercentageFeatureSet(Donor5_D1, pattern = "^MT-")
Donor4_D0[["percent.mt"]] <- PercentageFeatureSet(Donor4_D0, pattern = "^MT-")
Donor4_D1[["percent.mt"]] <- PercentageFeatureSet(Donor4_D1, pattern = "^MT-")
Donor3_D0[["percent.mt"]] <- PercentageFeatureSet(Donor3_D0, pattern = "^MT-")
Donor3_D1[["percent.mt"]] <- PercentageFeatureSet(Donor3_D1, pattern = "^MT-")


#Create list of objects
DonorList = list(Donor5_D0,Donor5_D1,Donor4_D0,Donor4_D1,Donor3_D0,Donor3_D1)
message("List created")

#Filtering out each object
DonorList <- lapply(DonorList, function(obj) {
  subset(obj, subset = nCount_RNA > 1000 & percent.mt < 20)})
message("Filering each object finished!")

#Using SCTransform on each object
Merged_Donors <- lapply(X = DonorList, FUN = SCTransform, vst.flavor = "v2", vars.to.regress = "percent.mt")
message("SCTransform Finished!")

#Integration
features <- SelectIntegrationFeatures(object.list = Merged_Donors, nfeatures = 3000)
Merged_Donors <- PrepSCTIntegration(object.list = Merged_Donors, anchor.features = features)
message("Integration prepped, finding anchors now.")

#Anchors
DonorAnchors <- FindIntegrationAnchors(object.list = Merged_Donors, normalization.method = "SCT", anchor.features = features)
saveRDS(DonorAnchors, "/home/gabriel.batzli/jamie_project/Merged_1_Anchors.rds")
message("Anchors made and saved")

Merged_Donors <- IntegrateData(anchorset = DonorAnchors, normalization.method = "SCT")
message("Integration complete.")

#Run PCA to find PCs
Merged_Donors <- RunPCA(Merged_Donors, verbose = FALSE)
saveRDS(Merged_Donors, "/home/gabriel.batzli/jamie_project/Merged_Donors_Object.rds")
message("Finished running PCA")

#Creating Elbowplot
EP <- ElbowPlot(Merged_Donors)
ggsave("/home/gabriel.batzli/jamie_project/Merged_1_Elbow", plot = EP, width = 8, height = 6)
message("Elbowplot created and saved")

































