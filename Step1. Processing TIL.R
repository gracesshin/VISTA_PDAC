############################################################################
# File 1. Processing TIL data ----------------------------------------------
############################################################################
### Step 1
### Extract TILs from whole tumor cells (Pan02) ###
# (Note)
# The UMAP coordinates provided include all cells used during the initial embedding computation. 
# Cells labeled as 'masked' were not used in any downstream analysis and do not appear in the manuscript figures or results. 
# These cells have been retained in the coordinate file to ensure that the UMAP structure is reproducible.

### Step 1-1. Create Seurat Object
# wt: Wild-type # vko: Vsir-/- #
Directories  <- list.dirs(path = Path.DataDir, full.names = FALSE, recursive = FALSE)
for(EachSet in Directories){
  CountMatrix = ReadMtx(mtx      = file.path(Path.DataDir,EachSet,"matrix.mtx.gz"),
                        features = file.path(Path.DataDir,EachSet,"features.tsv.gz"),
                        cells    = file.path(Path.DataDir,EachSet,"barcodes.tsv.gz"))
  assign(x = gsub(pattern = "_filtered_feature_bc_matrix", replacement = "", EachSet),
         value = CreateSeuratObject(counts = CountMatrix, min.cells = 3, min.features = 200))
}

### Step 1-2. Merge Seurat Objects
Pan02_Total <- merge(wt, y = c(vko, masked), add.cell.ids = c("wt", "vko", "masked"), project = "Pan02_Total")

### Step 1-3. Create a column titled "sample"
Pan02_Total %<>% AddMetaData(metadata = gsub("\\_.*","",colnames(Pan02_Total)), col.name = "sample")

### Step 1-4. Prepare Quality Control for each cell data
Pan02_Total %<>% AddMetaData(metadata = PercentageFeatureSet(Pan02_Total, pattern = "^mt-"), col.name = "mitoPercent")
Pan02_Total <- subset(Pan02_Total, subset = nFeature_RNA > 200 & nFeature_RNA < 9500 & nCount_RNA < 100000 & mitoPercent <10)

### Step 1-5. Normalize counts using SCTransform (Embedded in Seurat) and Execute PCA, UMAP sequentially
Pan02_SCTrans <- SCTransform(Pan02_Total, vst.flavor = "v2", verbose = FALSE)
Pan02_SCTrans %<>% RunPCA(npcs = 20, verbose = FALSE)
Pan02_SCTrans %<>% RunUMAP(reduction = "pca", dims = 1:12, verbose = FALSE)
Pan02_SCTrans %<>% FindNeighbors(dims = 1:12, verbose = FALSE)
Pan02_SCTrans %<>% FindClusters(resolution = 0.8, verbose = FALSE)

### Visualize Step 1. Cells from whole tumor
DimPlot(Pan02_SCTrans, label = TRUE)


### Step 2 
### Isolate immune cells from whole tumor cells ###
# Strategy : a-e
# a. Count cell numbers and Ptptrc distributions for each cluster
# b. In case where a cluster exhibited a bimodal distribution of Ptprc expression, split the cluster and the sub-cluster with higher Ptprc expression was selected.
# c. Excluding clusters counted less than 100
# d. Exclude cluster co-expressing  both Cd14 and Cd3g
# e. As mentioned, wild-type and Vista KO cells are going to be used in downstream analysis.

### Step 2-1. Collecting variables
# Strategy a. Count cell numbers and Ptptrc distributions for each cluster
table(Pan02_SCTrans$seurat_clusters)
VlnPlot(Pan02_SCTrans, "Ptprc", pt.size = 0) # Ptprc = Cd45

# Strategy b. Cluster 14 showed the bimodal distribution of Ptprc expression
# Note : Because cluster 14 have small cell numbers, reduced dimention and resolution parmeters
Pan02_cl14 <- subset(Pan02_SCTrans, subset = seurat_clusters == "14")
Pan02_cl14 %<>% RunUMAP(reduction = "pca", dims = 1:5, verbose = FALSE)
Pan02_cl14 %<>% FindNeighbors(dims = 1:5, verbose = FALSE)
Pan02_cl14 %<>% FindClusters(resolution = 0.5, verbose = FALSE)
Pan02_cl14_Ptprc.Low = subset(Pan02_cl14, subset = seurat_clusters %in% c(2,4), invert=TRUE)
Pan02_cl14_Excluding = which(colnames(Pan02_SCTrans) %in% colnames(Pan02_cl14_Ptprc.Low)) # Strategy b.

# Variable to apply to isolate TILs
Pan02_Low.Ptprc  <- c(1,4,9,12,19,18,27) # Strategy a.
Pan02_Bimodal    <- Pan02_cl14_Excluding # Strategy b. [not cluster, but cells]
Pan02_Low.Number <- c(26,27,28)  # Strategy c.
Pan02_Cd3.Cd14   <- c(24) # Strategy d.
Pan02_Groups     <- c("wt", "vko") # Strategy e. [Exclude masked group]

### Step 2-2. Apply variables
# Isolating TILs
Pan02_SCTrans <- Pan02_SCTrans[,-Pan02_Bimodal] # Apply b.
Pan02_TILs_WT_and_KO <- subset(Pan02_SCTrans, 
                               subset = seurat_clusters %in% c(Pan02_Low.Ptprc,    # Apply a.
                                                               Pan02_Low.Number,   # Apply c.
                                                               Pan02_Cd3.Cd14),    # Apply d,
                               invert=TRUE) %>% subset(., subset = sample %in% Pan02_Groups) # Apply e.

# Save the Seurat object into RDS file
saveRDS(Pan02_TILs_WT_and_KO, "Pan02_TILs_WT_and_KO.rds")

### Step 3
### Annotate Seurat clusters in TILs
Pan02_TILs_WT_and_KO %<>% AddMetaData(metadata = "", col.name = "Annotation")
Pan02_TILs_WT_and_KO$Annotation[Pan02_TILs_WT_and_KO$seurat_clusters %in% c(10,2,3,17)]    = "TAN"
Pan02_TILs_WT_and_KO$Annotation[Pan02_TILs_WT_and_KO$seurat_clusters %in% c(8,0,25,14)] = "MoMx"
Pan02_TILs_WT_and_KO$Annotation[Pan02_TILs_WT_and_KO$seurat_clusters %in% c(22)] = "Prolif.MoMx"
Pan02_TILs_WT_and_KO$Annotation[Pan02_TILs_WT_and_KO$seurat_clusters %in% c(23,13)] = "DC"
Pan02_TILs_WT_and_KO$Annotation[Pan02_TILs_WT_and_KO$seurat_clusters %in% c(5,7,20)] = "CD8"
Pan02_TILs_WT_and_KO$Annotation[Pan02_TILs_WT_and_KO$seurat_clusters %in% c(16)] = "Prolif.CD8"
Pan02_TILs_WT_and_KO$Annotation[Pan02_TILs_WT_and_KO$seurat_clusters %in% c(11)] = "NK"
Pan02_TILs_WT_and_KO$Annotation[Pan02_TILs_WT_and_KO$seurat_clusters %in% c(6)] = "CD4"
Pan02_TILs_WT_and_KO$Annotation[Pan02_TILs_WT_and_KO$seurat_clusters %in% c(15)] = "B"
Pan02_TILs_WT_and_KO$Annotation[Pan02_TILs_WT_and_KO$seurat_clusters %in% c(21)] = "Mast"

# Suppl. Figure S8a
DimPlot(Pan02_TILs_WT_and_KO, label = TRUE)+xlim(-15,10)+ylim(-15,10)+NoLegend()

# Suppl. Figure S8b
DimPlot(Pan02_TILs_WT_and_KO, group.by = "Annotation")+xlim(-15,10)+ylim(-15,10)+ggtitle(NULL)+NoLegend()

# Suppl. Figure S8c
TIL.ReLevels <- c( 10,2,3,17,8,0,25,14,22,23,13,5,7,16,20,11,6,15,21)
Idents(Pan02_TILs_WT_and_KO) <- factor(Pan02_TILs_WT_and_KO$seurat_clusters,levels = TIL.ReLevels)
TIL.Markers <- c("Mki67","Fcer1a",'Cd34',"Ms4a1",'Cd19','Cd4','Foxp3',"Itga2", "Klrb1c", "Ncr1", 'Xcl1',
                 'Cd247',"Tox",'Nkg7','Cd8b1',"Cd8a",'Cd3g','Xcr1',"Marco",'Adgre1','Cd86','Ccr2',"Cd14",
                 'Cd80',"Neat1", 'S100a9') 
DotPlot(Pan02_TILs_WT_and_KO, TIL.Markers, col.min = 0) + coord_flip()

# Save the Seurat object into RDS file
saveRDS(Pan02_TILs_WT_and_KO, "Pan02_TILs_WT_and_KO.RDS")


# Configure MoMx (Monocyte/Macrophage) -----------------------------------
# Because the variable name 'Pan02_MoMx_WT_and_KO' is long, we will use 'MoMx' as a shorthand from this point onward.
# Figure 3a: Highlight Monocyte/Macrophages in TIL
Pan02_TILs_WT_and_KO %<>% AddMetaData(metadata = ifelse(Pan02_TILs_WT_and_KO$Annotation %in% "MoMx", 1, 0), col.name = "MoMx")
FeaturePlot(Pan02_TILs_WT_and_KO, "MoMx", cols = c("gray70","red"))+xlim(-15,10)+ylim(-15,10)+NoLegend()+ ggtitle(NULL)

MoMx <- subset(Pan02_TILs_WT_and_KO, subset = Annotation == "MoMx")
MoMx %<>% RunUMAP(reduction = "pca", dims = 1:15, verbose = FALSE)
MoMx %<>% FindNeighbors(dims = 1:15, verbose = FALSE)
MoMx %<>% FindClusters(resolution = 0.5, verbose = FALSE)
saveRDS(MoMx, "Pan02_MoMx_WT_and_KO.RDS")



# Configure CD8T (CD8+ T cells) -------------------------------------------
# Because the variable name 'Pan02_CD8T_WT_and_KO' is long, we will use 'CD8T' as a shorthand from this point onward.
CD8T <- subset(Pan02_TILs_WT_and_KO, subset = Annotation == "CD8")
CD8T %<>% RunUMAP(reduction = "pca", dims = 1:15, verbose = FALSE)
CD8T %<>% FindNeighbors(dims = 1:15, verbose = FALSE)
CD8T %<>% FindClusters(resolution = 0.5, verbose = FALSE)
saveRDS(CD8T, "Pan02_CD8T_WT_and_KO.RDS")