############################################################################
# File 4. Validation using human Dataset -----------------------------------
############################################################################
# To validate our findings, we obtained all relevant data from the HTAN Pancreas NOS cohort. 
# In the final step, we compared cell population frequencies between the VSIR^High and VSIR^Low groups, using the sample IDs provided below.
#  HT064P1, HT115P1, HT121P1, HT185P1, HT123P1, HT190P1, HT191P1, HT168P1, HT125P1"

# Creating Seurat Object
reRun = FALSE # reRun = TRUE if creating Seurat Object at the first time
if(reRun){
  library(foreach)
  Directories <- list.dirs()[-1] %>% gsub("\\.\\/","",.)
  cluster <- parallel::makeCluster(40)
  doParallel::registerDoParallel(cluster)
  Human.Seurat_LIST <- foreach(s = 1:length(Directories), 
                               .packages = c("dplyr","data.table","Seurat","magrittr")) %dopar% {
                                 
                                 x = Directories[s]
                                 Read = ReadMtx(features = file.path(x,paste(x, "features.tsv.gz", sep = "-")),
                                                mtx =      file.path(x,paste(x, "matrix.mtx.gz", sep = "-")),
                                                cells =    file.path(x,paste(x, "barcodes.tsv.gz", sep = "-")))
                                 Human.Seurat_tmp = CreateSeuratObject(counts = Read)
                                 Human.Seurat_tmp %<>% AddMetaData(metadata = gsub("\\_1$","",x) %>% gsub(".*\\_","",.), col.name = "Donor")
                                 Human.Seurat_tmp = Human.Seurat_tmp[,which(Human.Seurat_tmp@assays$RNA$counts["PTPRC",] > 0)]
                                 return(Human.Seurat_tmp)
                               }
  parallel::stopCluster(cluster); gc()
  
  Human.Seurat_prep <- merge(x = Human.Seurat_LIST[[1]], y = c(Human.Seurat_LIST[[2]],  Human.Seurat_LIST[[3]],  Human.Seurat_LIST[[4]],
                                                               Human.Seurat_LIST[[5]],  Human.Seurat_LIST[[6]],  Human.Seurat_LIST[[7]],
                                                               Human.Seurat_LIST[[8]],  Human.Seurat_LIST[[9]],  Human.Seurat_LIST[[10]],
                                                               Human.Seurat_LIST[[11]], Human.Seurat_LIST[[12]], Human.Seurat_LIST[[13]],
                                                               Human.Seurat_LIST[[14]], Human.Seurat_LIST[[15]], Human.Seurat_LIST[[16]],
                                                               Human.Seurat_LIST[[17]], Human.Seurat_LIST[[18]], Human.Seurat_LIST[[19]],
                                                               Human.Seurat_LIST[[20]], Human.Seurat_LIST[[21]], Human.Seurat_LIST[[22]],
                                                               Human.Seurat_LIST[[23]], Human.Seurat_LIST[[24]], Human.Seurat_LIST[[25]],
                                                               Human.Seurat_LIST[[26]], Human.Seurat_LIST[[27]], Human.Seurat_LIST[[28]],
                                                               Human.Seurat_LIST[[29]], Human.Seurat_LIST[[30]], Human.Seurat_LIST[[31]],
                                                               Human.Seurat_LIST[[32]], Human.Seurat_LIST[[33]], Human.Seurat_LIST[[34]],
                                                               Human.Seurat_LIST[[35]], Human.Seurat_LIST[[36]], Human.Seurat_LIST[[37]],
                                                               Human.Seurat_LIST[[38]], Human.Seurat_LIST[[39]], Human.Seurat_LIST[[40]],
                                                               Human.Seurat_LIST[[41]], Human.Seurat_LIST[[42]], Human.Seurat_LIST[[43]],
                                                               Human.Seurat_LIST[[44]], Human.Seurat_LIST[[45]], Human.Seurat_LIST[[46]],
                                                               Human.Seurat_LIST[[47]], Human.Seurat_LIST[[48]], Human.Seurat_LIST[[49]],
                                                               Human.Seurat_LIST[[50]], Human.Seurat_LIST[[51]], Human.Seurat_LIST[[52]],
                                                               Human.Seurat_LIST[[53]], Human.Seurat_LIST[[54]], Human.Seurat_LIST[[55]],
                                                               Human.Seurat_LIST[[56]], Human.Seurat_LIST[[57]], Human.Seurat_LIST[[58]],
                                                               Human.Seurat_LIST[[59]], Human.Seurat_LIST[[60]], Human.Seurat_LIST[[61]],
                                                               Human.Seurat_LIST[[62]], Human.Seurat_LIST[[63]], Human.Seurat_LIST[[64]],
                                                               Human.Seurat_LIST[[65]], Human.Seurat_LIST[[66]], Human.Seurat_LIST[[67]],
                                                               Human.Seurat_LIST[[68]], Human.Seurat_LIST[[69]], Human.Seurat_LIST[[70]],
                                                               Human.Seurat_LIST[[71]], Human.Seurat_LIST[[72]], Human.Seurat_LIST[[73]],
                                                               Human.Seurat_LIST[[74]], Human.Seurat_LIST[[75]]),
                             add.cell.ids = c(gsub(".*\\_","",Directories)))
  
  Human.Seurat_prep <- JoinLayers(Human.Seurat_prep)
  Human.Seurat_prep@meta.data$orig.ident = rownames(Human.Seurat_prep@meta.data$Donor) %>% gsub("\\_.*","",.)
  
  Idents(Human.Seurat_prep) = Human.Seurat_prep$orig.ident
  Human.Seurat_prep[["percent.mt"]] <- PercentageFeatureSet(Human.Seurat_prep, pattern = "^MT-")
  saveRDS(object = Human.Seurat_prep, file = paste0("Human.Seurat_prep.rds"))
  
  Human.Seurat <- Seurat::SCTransform(subset(Human.Seurat_prep, subset = (percent.mt < 50) & (nFeature_RNA > 100) & (nFeature_RNA < 3000) & (nCount_RNA > 100) & (nCount_RNA   < 10000)),
                                      vst.flavor = "v2", ncells = 3000,
                                      conserve.memory = TRUE, do.correct.umi = TRUE,
                                      vars.to.regress = c("nCount_RNA", "nFeature_RNA"),
                                      return.only.var.genes = TRUE)
  
  Human.Seurat %<>% RunPCA(verbose = FALSE, npcs = 30)
  Human.Seurat %<>% RunUMAP(reduction = "pca", dims = 1:15)
  Human.Seurat %<>% FindNeighbors(dims = 1:15, reduction = "pca")
  Human.Seurat %<>% FindClusters(resolution = 0.8)
  
  saveRDS(object = Human.Seurat, file = paste0("Human.Seurat.rds"))
  Human.Seurat = readRDS("Human.Seurat.rds")
}else{
  Human.Seurat = readRDS("Human.Seurat.rds")
}

Human.Seurat %<>% AddMetaData(metadata = Human.Seurat$orig.ident %>% gsub("\\-.*","",.), col.name = "Case")
Human.Seurat <- subset(Human.Seurat, subset = Case != "TWCE")

# Explore seurat_clusters in Human dataset
DimPlot(Human.Seurat, label = TRUE)
VlnPlot(Human.Seurat, features = c("PTPRC","TRAC","CD3G","CD4","CD8B","CD19", "XCL1","CD163","ADGRE1","VSIR"), pt.size = 0)


# Macophages from PDAC TIL ------------------------------------------------
# Processing Macrophage cells from human PDAC TILs
MacrophageSeurat <- subset(Human.Seurat, subset = (seurat_clusters %in% "8") & (CD3E < 0.01) & (CD1C < 0.01) & (FCER1A < 0.01))
MacrophageSeurat %<>% RunPCA(verbose = FALSE, npcs = 30)
MacrophageSeurat %<>% RunUMAP(reduction = "pca", dims = 1:30)
MacrophageSeurat %<>% FindNeighbors(dims = 1:30, reduction = "pca")
MacrophageSeurat %<>% FindClusters(resolution = 0.3)


# Because Cluster 6 showed low CD83 and Cd86, we decided to exclude cluster 6
VlnPlot(MacrophageSeurat, c("CD83","CD80","CD86"))
MacrophageSeurat = subset(MacrophageSeurat, subset = (seurat_clusters %in% "6"), invert = TRUE)

MacrophageSeurat %<>% RunPCA(verbose = FALSE, npcs = 30)
MacrophageSeurat %<>% RunUMAP(reduction = "pca", dims = 1:30)
MacrophageSeurat %<>% FindNeighbors(dims = 1:30, reduction = "pca")
MacrophageSeurat %<>% FindClusters(resolution = 0.5)

# Suppl. Figure S10a -------------------------------------------------------
DimPlot(MacrophageSeurat, label=TRUE)+NoLegend()


# Suppl. Figure S10b -------------------------------------------------------
library(viridis)
Sorting.cluster <- c(1,2,0,4,3,6,5,7)
MacrophageSeurat$seurat_clusters <- factor(MacrophageSeurat$seurat_clusters, levels = Sorting.cluster)
SFig10b.Markers <- c("S100A8", "S100A9", "S100A12", "SELL", "RETN", "TIMP1", "AREG", "ZFAND2A", "RNASE1", 
                     "CCL2", "MAF", "F13A1", "SPP1", "MARCO", "FABP5","CSTB","HLA-DQA1","HLA-DPA1","HLA-DPB1", 
                     "C1QB", "APOE", "APOC1", "GPNMB", "PLA2G7", "IL1B", "TNF", "CXCL2", "CCL20", "CDKN1C",
                     "ADGRE1", "CD79B", "RAP1GAP2") %>% rev
DotPlot(object = MacrophageSeurat, features = SFig10b.Markers, group.by = "seurat_clusters") + scale_colour_viridis(option = "viridis") + coord_flip()


# Suppl. Figure S10c -------------------------------------------------------
# Export count cell numbers by Seurat cluster and Donor
Donors <- unique(MacrophageSeurat@meta.data$Case)
MetaData <- MacrophageSeurat@meta.data 
Number <- lapply(MacrophageSeurat, function(Case){
  Cell   = MetaData$seurat_clusters[(MetaData$Case == Case)]
  Counts = table(Cell)
  Export = data.table(Cluster = levels(Cell), Number = Counts[match(levels(Cell), names(Counts))])
  return(Export)
}) %>% bind_cols %>% .[,grep("Number.N", colnames(.)),with=FALSE]
colnames(Number) <- Donors 
ReportCounts <- data.table(Cluster = levels(MetaData$seurat_clusters), Number)
write.table(ReportCounts,"Macrophage.ReportCounts.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


# CD8+ T cells from PDAC TILs ---------------------------------------------

# Because cluster 7 is mixed with gdT cells, we first defined cell barcoded of gdT cells in cluster 7
Human.Seurat_Cluster7 <- subset(Human.Seurat, subset = seurat_clusters == "7")
Human.Seurat_Cluster7 %<>% RunPCA(verbose = FALSE, npcs = 30)
Human.Seurat_Cluster7 %<>% RunUMAP(reduction = "pca", dims = 1:30)
Human.Seurat_Cluster7 %<>% FindNeighbors(dims = 1:30, reduction = "pca")
Human.Seurat_Cluster7 %<>% FindClusters(resolution = 0.8)
DimPlot(Human.Seurat_Cluster7, label = TRUE)

# gdT cellsL: TRDC positive cluster 
VlnPlot(Human.Seurat_Cluster7, "TRDC", pt.size = 0)
TRDC_positive_cluster <- c("4","5","7","9","10","11","12","13")
Barcodes_gdTcells  <- colnames(Human.Seurat_Cluster7)[Human.Seurat_Cluster7$seurat_clusters %in%  TRDC_positive_cluster]

# Select CD8 T cells from TILs & exclude gdT cells
CD8TSuerat <- subset(Human.Seurat, subset = seurat_clusters %in% c(2,4,7,15))
CD8TSuerat <- CD8TSuerat[,-which(colnames(CD8TSuerat) %in% Barcodes_gdTcells)]

# Re analysis of CD8T cell clusters 
CD8TSuerat %<>% RunPCA(verbose = FALSE, npcs = 30)
CD8TSuerat %<>% RunUMAP(reduction = "pca", dims = 1:30)
CD8TSuerat %<>% FindNeighbors(dims = 1:30, reduction = "pca")
CD8TSuerat %<>% FindClusters(resolution = 0.5)
DimPlot(CD8TSuerat, label = TRUE)

# Because cluster 9 is CD8B- exclude the cluster
CD8TSuerat.Fin <- subset(CD8TSuerat, subset = seurat_clusters %in% c("9"), invert = TRUE)
CD8TSuerat.Fin %<>% RunPCA(verbose = FALSE, npcs = 30)
CD8TSuerat.Fin %<>% RunUMAP(reduction = "pca", dims = 1:30)
CD8TSuerat.Fin %<>% FindNeighbors(dims = 1:30, reduction = "pca")
CD8TSuerat.Fin %<>% FindClusters(resolution = 0.5) # 0.3 prev
DimPlot(CD8TSuerat.Fin, label = TRUE)


# Suppl. Figure S17a -------------------------------------------------------
DimPlot(CD8TSuerat.Fin, label = TRUE) + xlim(c(-10,6))+ylim(c(-8,8))+NoLegend()


# Suppl. Figure S17b -------------------------------------------------------
library(viridis)
Sorting.cluster <- c(0,1,6,9,11,5,12,2,3,10,13,7,4,8)
CD8TSuerat.Fin$seurat_clusters <- factor(CD8TSuerat.Fin$seurat_clusters, levels = Sorting.cluster)
SFig17b.Markers <- c("GZMK","GZMB","KLRG1","GPR183","TRAT1","CCL4","TNF","CCL4L2","IFNG","CCL3","XCL1","XCL2","PPP1R2C",
                     "KCNQ1OT1","KLRC1","IFI6","TOX","MKI67","CXCL13","CTLA4","GNLY","KLRF1") %>% rev
DotPlot(object = CD8TSuerat.Fin, features = SFig17b.Markers, group.by = "seurat_clusters") + scale_colour_viridis(option = "viridis") + coord_flip()


# Suppl. Figure S17c -------------------------------------------------------
# Export count cell numbers by Seurat cluster and Donor
Donors <- unique(CD8TSuerat.Fin@meta.data$Case)
MetaData <- CD8TSuerat.Fin@meta.data 
Number <- lapply(CD8TSuerat.Fin, function(Case){
  Cell   = MetaData$seurat_clusters[(MetaData$Case == Case)]
  Counts = table(Cell)
  Export = data.table(Cluster = levels(Cell), Number = Counts[match(levels(Cell), names(Counts))])
  return(Export)
}) %>% bind_cols %>% .[,grep("Number.N", colnames(.)),with=FALSE]
colnames(Number) <- Donors 
ReportCounts <- data.table(Cluster = levels(MetaData$seurat_clusters), Number)
write.table(ReportCounts,"CD8.ReportCounts.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
