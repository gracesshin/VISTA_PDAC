############################################################################
# File 2. Processing Monocyte/Macrophage -----------------------------------
############################################################################

############################################################################
# Monocyte/Macrophage (MoMx) -----------------------------------------------
# Because the variable name 'Pan02_MoMx_WT_and_KO' is long, we will use 'MoMx' as a shorthand from this point onward.
############################################################################

MoMx <- readRDS("Pan02_MoMx_WT_and_KO.RDS.RDS")
# Figure 3b ---------------------------------------------------------------
DimPlot(MoMx, group.by = "seurat_clusters")+xlim(-8,8)+ylim(-6,6)+NoLegend()+ ggtitle(NULL)


# Figure 3c-1 -------------------------------------------------------------
DimPlot(MoMx, group.by = "sample", cols = c("#F04635", "#368AF0"))+xlim(-8,8)+ylim(-6,6)+NoLegend()+ ggtitle(NULL)


# Figure 3c-2 -------------------------------------------------------------
# Create trajectory lines dependent to variable 'sample'
MoMx_WT  <- subset(MoMx, subset = sample %in% "wt")
MoMx_VKO <- subset(MoMx, subset = sample %in% "vko")

# Because cluster 6 is not expressing Cx3cr1 or Ccr2 draw trajectory lines without cluster 6
MoMx_WT_Seurat  <- subset(MoMx_WT,  idents = "6", invert = TRUE)
MoMx_VKO_Seurat <- subset(MoMx_VKO, idents = "6", invert = TRUE)


# Downsampling samples based on the genotypes
# Calculate trajectory line in WT and VKO  
MoMx_Sizes         <- c(ncol(MoMx_WT_Seurat), ncol(MoMx_VKO_Seurat))
MoMx_Minimum_Sizes <- min(MoMx_Sizes)

## Wild type
set.seed(39)
Index.MoMx_WT <- sample(1:ncol(MoMx_WT_Seurat),  MoMx_Minimum_Sizes, replace = FALSE) %>% sort
dn.MoMx_WT_Seurat <- MoMx_WT_Seurat[,Index.MoMx_WT]
dn.MoMx_WT_CDS <- as.cell_data_set(dn.MoMx_WT_Seurat)
dn.MoMx_WT_CDS %<>% cluster_cells(resolution = 8e-4)
dn.MoMx_WT_CDS %<>% learn_graph() %>% order_cells()
dn.MoMx_WT_Pseudotime <- dn.MoMx_WT_CDS@principal_graph_aux$UMAP$pseudotime
dn.MoMx_WT_Seurat %<>% AddMetaData(metadata = dn.MoMx_WT_Pseudotime, col.name = "pseudotime")

## Vsir-/-
set.seed(39)
Index.MoMx_VKO <- sample(1:ncol(MoMx_VKO_Seurat), MoMx_Minimum_Sizes, replace = FALSE) %>% sort
dn.MoMx_VKO_Seurat <- MoMx_VKO_Seurat[,Index.MoMx_VKO]
dn.MoMx_VKO_CDS <- as.cell_data_set(dn.MoMx_VKO_Seurat)
dn.MoMx_VKO_CDS %<>% cluster_cells(resolution = 8e-4)
dn.MoMx_VKO_CDS %<>% learn_graph() %>% order_cells()
dn.MoMx_VKO_Pseudotime = dn.MoMx_VKO_CDS@principal_graph_aux$UMAP$pseudotime
dn.MoMx_VKO_Seurat %<>% AddMetaData(metadata = dn.MoMx_VKO_Pseudotime, col.name = "pseudotime")


# Figure 3d - WT ----------------------------------------------------------
Genexpression_timepoints(seuratobjt = dn.MoMx_WT_Seurat, genesymbol = "Spp1",   divider = 0.7)$mean %>% write.table(., "Spp1_WT_pseudotime.txt",  col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
Genexpression_timepoints(seuratobjt = dn.MoMx_WT_Seurat, genesymbol = "Cxcl9",  divider = 0.7)$mean %>% write.table(., "Cxcl9_WT_pseudotime.txt", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
plot(dn.MoMx_WT_Seurat$pseudotime, as.numeric(dn.MoMx_WT_Seurat$seurat_clusters),  xlim = c(0, max(dn.MoMx_WT_Pseudotime)))


# Figure 3d - VKO ---------------------------------------------------------
Genexpression_timepoints(seuratobjt = dn.MoMx_VKO_Seurat, genesymbol = "Spp1",   divider = 0.7)$mean %>% write.table(., "Spp1_VKO_pseudotime.txt",  col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
Genexpression_timepoints(seuratobjt = dn.MoMx_VKO_Seurat, genesymbol = "Cxcl9",  divider = 0.7)$mean %>% write.table(., "Cxcl9_VKO_pseudotime.txt", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
plot(dn.MoMx_VKO_Seurat$pseudotime,as.numeric(dn.MoMx_VKO_Seurat$seurat_clusters), xlim = c(0, max(dn.MoMx_VKO_Pseudotime)))


# Figure 3e ---------------------------------------------------------------
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = MoMx_WT,  Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode ="Freq")
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = MoMx_VKO, Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode ="Freq")


# Figure 3f - 1st column --------------------------------------------------
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "4"), Mode = "Freq")
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "0"), Mode = "Freq")
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "6"), Mode = "Freq")
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "1"), Mode = "Freq")
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "3"), Mode = "Freq")
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "7"), Mode = "Freq")
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "5"), Mode = "Freq")
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "2"), Mode = "Freq")


# Figure 3f - 2nd column --------------------------------------------------
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "4"), Mode = "Number")
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "0"), Mode = "Number")
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "6"), Mode = "Number")
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "1"), Mode = "Number")
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "3"), Mode = "Number")
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "7"), Mode = "Number")
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "5"), Mode = "Number")
Calculate_TwoGenotype(Seurat = subset(MoMx, subset = seurat_clusters == "2"), Mode = "Number")


# Figure 3f - 3rd column --------------------------------------------------
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "4"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "Number")
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "0"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "Number")
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "6"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "Number")
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "1"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "Number")
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "3"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "Number")
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "7"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "Number")
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "5"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "Number")
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "2"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "Number")


# Figure 3f - 4th column --------------------------------------------------
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "4"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "GeneRatio")
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "0"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "GeneRatio")
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "6"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "GeneRatio")
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "1"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "GeneRatio")
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "3"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "GeneRatio")
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "7"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "GeneRatio")
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "5"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "GeneRatio")
Calculate_TwoGene_Expressing(Original = MoMx, Seurat = subset(MoMx, subset = seurat_clusters == "2"), Gene1 ="Cxcl9", Gene2="Spp1", method = "Mean", Mode = "GeneRatio")


# Suppl. Figure S9a -------------------------------------------------------
Idents(MoMx) <- factor(Idents(MoMx), levels = c(4,0,6,1,3,7,5,2))
VlnPlot(MoMx, features = c("Vsir","Adgre1"), pt.size = 0)


# Suppl. Figure S9b -------------------------------------------------------
Idents(MoMx) <- factor(Idents(MoMx), levels = c(4,0,6,1,3,7,5,2))
VlnPlot(MoMx, features = c("Apoe","Trem2","Abca1","Mertk"), pt.size = 0)


# Suppl. Figure S9c -------------------------------------------------------
library(viridis)
Idents(MoMx) = factor(Idents(MoMx), levels = c(4,0,6,1,3,7,5,2))
SFig9.markers = c("Stat1", "Stat2", "Stat4", "Jak1", "Jak2", "Irf1", "Irf7", "Irf8", "Irf9", "Socs1", "Socs3", "Pim1", "Cxcl9", "Cxcl10", "Vcam1", "Ccl9", "Ccl7", "Ccl5", "Il18bp", "Itgax", "Slamf7")
DoHeatmap(MoMx, features = SFig9.markers, draw.line = T, label = F, assay = "SCT", disp.min = -2, disp.max = 2, slot = "scale.data", lines.width = 10) + 
  scale_fill_viridis(na.value = "white")


# Suppl. Figure S9d (Cxcl9 vs Spp1) -------------------------------------------------------
ProcessingTable <- FetchData(MoMx, vars = c("Cxcl9", "Spp1", "seurat_clusters"))
set.seed(42)
# Because NA value do not displayed in contour plot, replace NA value to a small random number ranged from 0 to 0.01
ProcessingTable[,1] <- ifelse(is.na(ProcessingTable[,1]),0,ProcessingTable[,1]) + runif(n = nrow(ProcessingTable), min = 0, max = 0.01)
ProcessingTable[,2] <- ifelse(is.na(ProcessingTable[,2]),0,ProcessingTable[,2]) + runif(n = nrow(ProcessingTable), min = 0, max = 0.01)
colnames(ProcessingTable) <- c("Cxcl9", "Spp1", "seurat_clusters")

Group_lists = data.table(Group = paste(ifelse(ProcessingTable[,"Cxcl9"] > mean_cutoff(Gene = "Cxcl9", Seurat = MoMx),"Cxcl9+","Cxcl9-"),
                                       ifelse(ProcessingTable[,"Spp1"]  > mean_cutoff(Gene = "Spp1",  Seurat = MoMx),"Spp1+","Spp1-"), sep = ""),
                         Cluster = ProcessingTable[,"seurat_clusters"])

Group_numbers <- lapply(Clusters_to_Select, function(EachCluster){
  Cluster_Annotated = Group_lists[Cluster == EachCluster]$Group %>% sort %>% table %>% as.data.table 
  colnames(Cluster_Annotated)[2] = paste0("Cluster.",EachCluster)
  return(Cluster_Annotated)
}) %>% bind_cols()
colnames(Group_numbers)[1] = "Groups"
Group_numbers <- Group_numbers[,-grep("\\.\\.\\.", colnames(Group_numbers)), with=FALSE]

# Cxcl9 vs Spp1 Contour plot
ggplot(ProcessingTable, aes(x = Cxcl9, y = Spp1)) + geom_density_2d(aes(color = seurat_clusters), size = 0.5, contour_var = "ndensity") +
  geom_vline(xintercept= mean_cutoff(Gene = "Cxcl9", Seurat = MoMx), linetype="dashed", color = "red", size=0.5)+
  geom_hline(yintercept= mean_cutoff(Gene = "Spp1",  Seurat = MoMx), linetype="dashed", color = "red", size=0.5)+
  facet_wrap(~seurat_clusters) + labs(x = "Cxcl9", y = "Spp1") + scale_color_brewer(palette = "Dark2") + theme_minimal()+ NoLegend()


# Figure 4a ---------------------------------------------------------------
# Differential genes between MoMx3vs1 and MoMx2vs1
MoMx <- PrepSCTFindMarkers(MoMx)
pct.cutoff = 0.25

# Cluster 3vs1
MoMx_Compare_3vs1 <- FindMarkers(object= MoMx, ident.1 = "3", ident.2 = "1", assay = "SCT", group_by = 'seurat_clusters', recorrect_umi = FALSE)
MoMx_3vs1.Raw <- data.table(symbol = rownames(MoMx_Compare_3vs1),MoMx_Compare_3vs1)
MoMx_3vs1.DEG <-  rbindlist(list(MoMx_3vs1.Raw[avg_log2FC > 0] %>% .[pct.1 > pct.cutoff], MoMx_3vs1.Raw[avg_log2FC < 0] %>% .[pct.2 > pct.cutoff])) 
write.table(MoMx_3vs1.DEG, "MoMx_3vs1.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Cluster 2vs1
MoMx_Compare_2vs1 <- FindMarkers(object= MoMx, ident.1 = "2", ident.2 = "1", assay = "SCT", group_by = 'seurat_clusters', recorrect_umi = FALSE)
MoMx_2vs1.Raw <- data.table(symbol = rownames(MoMx_Compare_2vs1),MoMx_Compare_2vs1)
MoMx_2vs1.DEG <-  rbindlist(list(MoMx_2vs1.Raw[avg_log2FC > 0] %>% .[pct.1 > pct.cutoff], MoMx_2vs1.Raw[avg_log2FC < 0] %>% .[pct.2 > pct.cutoff])) 
write.table(MoMx_2vs1.DEG, "MoMx_2vs1.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Merge Cluster 3vs1 and 2vs1
Cluster3vs1.2vs1 <- merge.data.table(x= MoMx_3vs1.DEG, y = MoMx_2vs1.DEG, by.x = "symbol", by.y = "symbol", all.x = TRUE, all.y = TRUE)
Cluster3vs1.2vs1$avg_log2FC.x[is.na(Cluster3vs1.2vs1$avg_log2FC.x)] = 0
Cluster3vs1.2vs1$avg_log2FC.y[is.na(Cluster3vs1.2vs1$avg_log2FC.y)] = 0
Cluster3vs1.2vs1$p_val_adj.x[is.na(Cluster3vs1.2vs1$p_val_adj.x)] = 1
Cluster3vs1.2vs1$p_val_adj.y[is.na(Cluster3vs1.2vs1$p_val_adj.y)] = 1

# Figure 4a (gt glance)
plot(Cluster3vs1.2vs1$avg_log2FC.x, Cluster3vs1.2vs1$avg_log2FC.y)
write.table(Cluster3vs1.2vs1, "Cluster3vs1.2vs1.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


# Suppl. Figure S12a ------------------------------------------------------
# PCA based on the gene expression level (SCT)
Avg.Expression <- (AverageExpression(MoMx, assays = "SCT")$SCT+0.1) %>% data.matrix 
Avg.Expression <- Avg.Expression[apply(Avg.Expression,1,mean) > 0.1,]
Avg.Expression <- Avg.Expression[apply(Avg.Expression,1,mean) < 100,]
Avg.Expression_PCA <- prcomp(t(Avg.Expression))
plot(Avg.Expression_PCA$x[,1:2])
text(Avg.Expression_PCA$x[,1:2], label = rownames(Avg.Expression_PCA$x))
summary(Avg.Expression_PCA)



# Suppl. Figure S12b ------------------------------------------------------
library(viridis)
Idents(MoMx) <- factor(Idents(MoMx), levels = c(4,0,6,1,3,7,5,2))
MoMx$seurat_clusters
SF12.Markers <- c("Cxcl9", "Spp1", "Cxcl10", "Vcam1", "Ccl9", "Ccl7", "Ccl5", "Il18bp", "Itgax", "Slamf7", "Slamf9")
DoHeatmap(MoMx, features = SF12.Markers, draw.line = T, label = F, assay = "SCT", disp.min = -2, disp.max = 2, slot = "scale.data", lines.width = 10) + 
  scale_fill_viridis(na.value = "white")


# Suppl. Figure S12c ------------------------------------------------------
library(ggridges)
KEGG_APC_Geneset <- list(c("Lgmn","Tap2","Nfya","Creb1","Nfyc","Ctsl","Psme3","Ctsb","Pdia3","Ctss","Canx","Ifi30","Hspa4","Psme1","Tap1","Hsp90ab1","Ciita","Hsp90aa1","Cd74","Psme2","Hspa5","Calr","Tapbp","Hspa8"))
MoMx <- AddModuleScore(MoMx, features = KEGG_APC_Geneset, nbin = 20, ctrl = 5, assay = "SCT", name = "APCs_")
MoMx_APCs_Data.Table = data.table(Clusters = MoMx$seurat_clusters, APC_score = MoMx$APCs_1)
MoMx_APCs_Data.Table$Clusters = factor(MoMx_APCs_Data.Table$Clusters, levels = Idents(MoMx) %>% levels %>% rev )
ggplot(MoMx_APCs_Data.Table, aes(x = APC_score, y = Clusters, fill = Clusters)) + geom_density_ridges() + theme_ridges() + xlim(-0.4,1.4)+ theme(legend.position = "none")
# P-values
p.value_MoMX_1vs2 <- wilcox.test(MoMx_APCs_Data.Table[Clusters == 1]$APC_score, MoMx_APCs_Data.Table[Clusters == 2]$APC_score)$p.value
p.value_MoMX_1vs3 <- wilcox.test(MoMx_APCs_Data.Table[Clusters == 1]$APC_score, MoMx_APCs_Data.Table[Clusters == 3]$APC_score)$p.value
p.value_MoMX_2vs3 <- wilcox.test(MoMx_APCs_Data.Table[Clusters == 2]$APC_score, MoMx_APCs_Data.Table[Clusters == 3]$APC_score)$p.value
