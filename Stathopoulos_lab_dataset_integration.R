##Codes for Integrating three wildtype datasets
# The datasets include: wildtype#1, wildtype#2, fixed_wildtype

library(Seurat)
library(monocle)

#Load Cassie's WT#1 data
Doe_dm_wt1<-readRDS("R_input/76476877_220901Sta.rds")

#Load WT#2 data
Doe_dm_wt2<-readRDS("R_input/Doe_dm_06022023Sta.rds")

#Load WT_fixed data
Doe_dm_fixed<-readRDS("R_input/Doe_dm_fixed_06022023Sta.rds")

#Changing the orginal identities of each datasets
Doe_dm_wt1$orig.ident = "wt1"
Doe_dm_wt2$orig.ident = "wt2"
Doe_dm_fixed$orig.ident = "wtfixed"


# Combine dataset
Doe_dm.list <- list(Doe_dm_wt1, Doe_dm_wt2,Doe_dm_fixed)

# normalize and identify variable features for each dataset independently
Doe_dm.list <- lapply(X = Doe_dm.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Doe_dm.list)

#Perform integration
Doe_dm.anchors <- FindIntegrationAnchors(object.list = Doe_dm.list, anchor.features = features)

# this command creates an 'integrated' data assay
wt.combined <- IntegrateData(anchorset = Doe_dm.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(wt.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
wt.combined <- ScaleData(wt.combined, verbose = FALSE)
wt.combined <- RunPCA(wt.combined, npcs = 30, verbose = FALSE)
ElbowPlot(wt.combined, ndims = 60, reduction = "pca")
wt.combined <- RunUMAP(wt.combined, reduction = "pca", dims = 1:30)
wt.combined <- FindNeighbors(wt.combined, reduction = "pca", dims = 1:30)
wt.combined <- FindClusters(wt.combined, resolution = 0.5)

saveRDS(wt.combined, file ="R_output/wt_combinedPC20_081723.rds")

#Visualization
p1 <- DimPlot(wt.combined, reduction = "umap",group.by="orig.ident")

p2 <- DimPlot(wt.combined, reduction = "umap", split.by = "orig.ident")

p1
p2

#Write markers
library(dplyr)

markers <- FindAllMarkers(wt.combined, only.pos = T)

write.table(markers, file="wt_all_combined_markers_0817.txt", quote=F, sep="\t", col.names=T)

###############################################################################################################

##Codes for integrating stage 12 dataset  by Doe et al. with Wildtype #1 dataset
##The wildtype dataset#1 was previously filtered and analyzed, and it is named Fang dataset for the following analysis.
## Stage 12 dataset(GSE202987: GSM6145582) is downloaded from (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202987)

library(Seurat)
library(monocle)

dm.data<-Read10X("GSM6145582_stg12_aggregated_barcodes_filtered_matrix", unique.features = TRUE)
dm<-CreateSeuratObject(dm.data, project = "cellranger")
dm[["percent.rRNA"]] <- PercentageFeatureSet(dm, pattern = "rRNA")
VlnPlot(dm, features = c("nFeature_RNA", "nCount_RNA", "percent.rRNA"), ncol = 3, pt.size = 0)
dm <- subset(dm, subset = nFeature_RNA < 7000 & nCount_RNA > 250 & nCount_RNA < 100000 & percent.rRNA<25)
dm <- NormalizeData(dm, normalization.method = "LogNormalize", scale.factor = 10000)
dm <- FindVariableFeatures(dm, selection.method = "vst", nfeatures = 2000)
dm <- ScaleData(dm, verbose = FALSE, vars.to.regress = c("nCount_RNA"))
dm <- RunPCA(dm, npcs = 60, ndims.print = 1:5, nfeatures.print = 5)
# t-SNE and Clustering
dm <- FindNeighbors(dm, reduction = "pca", dims = 1:30)
dm <- RunTSNE(dm, dims = 1:30, check_duplicates = FALSE)
dm <- RunUMAP(dm, dims = 1:30, check_duplicates = FALSE)
DimPlot(dm)
saveRDS(dm,file="DOE_WTSta.rds",compress=F)

##Load wildtype dataset and integrate with the filtered stage 12 dataset 
#Load Fang's data
dm_1<-readRDS("76476877_220901Sta.rds")

#Load Doe's data
dm_2<-readRDS("DOE_WTSta.rds")

#Changing the orginal identities of each datasets
dm_1$orig.ident = "Fang"
dm_2$orig.ident = "Doe"

# Combine dataset
dm.list <- list(dm_1, dm_2)

# normalize and identify variable features for each dataset independently
dm.list <- lapply(X = dm.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = dm.list)

#Perform integration
dm.anchors <- FindIntegrationAnchors(object.list = dm.list, anchor.features = features)

# this command creates an 'integrated' data assay
wt.combined <- IntegrateData(anchorset = dm.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(wt.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
wt.combined <- ScaleData(wt.combined, verbose = FALSE)
wt.combined <- RunPCA(wt.combined, npcs = 20, verbose = FALSE)
wt.combined <- RunUMAP(wt.combined, reduction = "pca", dims = 1:20)
wt.combined <- FindNeighbors(wt.combined, reduction = "pca", dims = 1:20)
wt.combined <- FindClusters(wt.combined, resolution = 0.5)

saveRDS(wt.combined, file ="FD_combined_081723.rds")

#Visualization
p1 <- DimPlot(wt.combined, reduction = "umap",group.by="orig.ident")

p2 <- DimPlot(wt.combined, reduction = "umap", split.by = "orig.ident")

p3 <- DimPlot(dm_1, reduction = "umap")
p4 <- DimPlot(dm_2, reduction = "umap")

p1
p2
p3 + p4 

#Write markers
library(dplyr)

markers <- FindAllMarkers(wt.combined, only.pos = T)

write.table(markers, file="FDwtcombined_markers_081723.txt", quote=F, sep="\t", col.names=T)
