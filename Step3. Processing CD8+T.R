############################################################################
# File 3. Processing CD8+ T cells ------------------------------------------
############################################################################

############################################################################
# CD8 ---------------------------------------------------------------------
# Because the variable name 'Pan02_CD8T_WT_and_KO' is long, we will use 'CD8T' as a shorthand from this point onward.
############################################################################
CD8T <- readRDS("Pan02_CD8T_WT_and_KO.RDS")

# Prepare Contour Plot in CD8T (using Four clusters) -------------------------------------------------------
Clusters_to_Select <- c(0,1,2,3)
CD8T_CellNumbers  <- CD8T_FourClusters$seurat_clusters %>% table %>%  .[names(.) %in% Clusters_to_Select]
CD8T_MinNumbers   <- min(CD8T_CellNumbers)
CD8T_FourClusters <- subset(CD8T, subset = seurat_clusters %in% Clusters_to_Select)
CD8T_FourClusters <- CD8T_FourClusters[,lapply(Clusters_to_Select, function(EachCluster){ sample(which(CD8T_FourClusters$seurat_clusters %in% EachCluster), size = CD8T_MinNumbers, replace = FALSE) }) %>% unlist %>% sort]


# Figure 5a ---------------------------------------------------------------
DimPlot(CD8T, label = TRUE)+xlim(-6,10)+ylim(-6,8)+NoLegend()+ggtitle(NULL)


# Figure 5b ---------------------------------------------------------------
Idents(CD8T) <- factor(Idents(CD8T),c(4,5,0,1,2,3))
F5b.Markers <- c("Fos","Gab2","Maf","Il2ra","Tox","Havcr2","Pdcd1","Igf2r","Cx3cr1","Klrd1","Gzmk","Slamf7","Eomes","Klrc1","Cd101","Cd38","Cxcr3")
DotPlot(CD8T, F5b.Markers, col.min = 0)


# Figure 5c ---------------------------------------------------------------
DimPlot(CD8T, label = FALSE, cols = c("#F04635", "#368AF0"), group.by = "sample")+xlim(-6,10)+ylim(-6,8)+NoLegend()+ggtitle(NULL)
CD8T_WT  <- subset(CD8T, subset = sample == "wt")
CD8T_VKO <- subset(CD8T, subset = sample == "vko")


# Figure 5d ---------------------------------------------------------------
CD8T_WT_Freq  <- table(CD8T_WT$seurat_clusters)/ncol(CD8T_WT)
CD8T_VKO_Freq <- table(CD8T_VKO$seurat_clusters)/ncol(CD8T_VKO)


# Figure 5e (left; Eomes vs Cd38) -----------------------------------------
ProcessingTable   <- FetchData(CD8T_FourClusters, vars = c("Eomes", "Cd38", "seurat_clusters"))
set.seed(42)
# Because NA value do not displayed in contour plot, replace NA value to a small random number ranged from 0 to 0.01
ProcessingTable[,1] <- ifelse(is.na(ProcessingTable[,1]),0,ProcessingTable[,1]) + runif(n = nrow(ProcessingTable), min = 0, max = 0.01)
ProcessingTable[,2] <- ifelse(is.na(ProcessingTable[,2]),0,ProcessingTable[,2]) + runif(n = nrow(ProcessingTable), min = 0, max = 0.01)
colnames(ProcessingTable) <- c("Eomes", "Cd38", "seurat_clusters")

Group_lists <- data.table(Group = paste(ifelse(ProcessingTable[,"Eomes"] > mean_cutoff(Gene = "Eomes", Seurat = CD8T_FourClusters),"Eomes+","Eomes-"),
                                        ifelse(ProcessingTable[,"Cd38"]  > mean_cutoff(Gene = "Cd38",  Seurat = CD8T_FourClusters),"Cd38+","Cd38-"), sep = ""),
                          Cluster = ProcessingTable[,"seurat_clusters"])

Group_numbers <- lapply(Clusters_to_Select, function(EachCluster){
  Cluster_Annotated = Group_lists[Cluster == EachCluster]$Group %>% sort %>% table %>% as.data.table 
  colnames(Cluster_Annotated)[2] = paste0("Cluster.",EachCluster)
  return(Cluster_Annotated)
}) %>% bind_cols()
colnames(Group_numbers)[1] = "Groups"
Group_numbers = Group_numbers[,-grep("\\.\\.\\.", colnames(Group_numbers)), with=FALSE]

# Eomes vs Cd38 Contour
ggplot(ProcessingTable, aes(x = Eomes, y = Cd38)) + geom_density_2d(aes(color = seurat_clusters), size = 0.5, contour_var = "ndensity") +
  geom_vline(xintercept= mean_cutoff(Gene = "Eomes", Seurat = CD8T_FourClusters), linetype="dashed", color = "red", size=0.5)+
  geom_hline(yintercept= mean_cutoff(Gene = "Cd38",  Seurat = CD8T_FourClusters), linetype="dashed", color = "red", size=0.5)+
  facet_wrap(~seurat_clusters) + labs(x = "Eomes", y = "Cd38") + scale_color_brewer(palette = "Dark2") + theme_minimal()+ NoLegend()



# Figure 5e (Right; Pdcd1 vs Cd226) -------------------------------------------------------
ProcessingTable   <- FetchData(CD8T_FourClusters, vars = c("Pdcd1", "Cd226", "seurat_clusters"))
set.seed(42)
# Because NA value do not displayed in contour plot, replace NA value to a small random number ranged from 0 to 0.01
ProcessingTable[,1] <- ifelse(is.na(ProcessingTable[,1]),0,ProcessingTable[,1]) + runif(n = nrow(ProcessingTable), min = 0, max = 0.01)
ProcessingTable[,2] <- ifelse(is.na(ProcessingTable[,2]),0,ProcessingTable[,2]) + runif(n = nrow(ProcessingTable), min = 0, max = 0.01)
colnames(ProcessingTable) <- c("Pdcd1", "Cd226", "seurat_clusters")

Group_lists <- data.table(Group = paste(ifelse(ProcessingTable[,"Pdcd1"] > mean_cutoff(Gene = "Pdcd1", Seurat = CD8T_FourClusters),"Pdcd1+","Pdcd1-"),
                                        ifelse(ProcessingTable[,"Cd226"]  > mean_cutoff(Gene = "Cd226",  Seurat = CD8T_FourClusters),"Cd226+","Cd226-"), sep = ""),
                          Cluster = ProcessingTable[,"seurat_clusters"])

Group_numbers <- lapply(Clusters_to_Select, function(EachCluster){
  Cluster_Annotated = Group_lists[Cluster == EachCluster]$Group %>% sort %>% table %>% as.data.table 
  colnames(Cluster_Annotated)[2] = paste0("Cluster.",EachCluster)
  return(Cluster_Annotated)
}) %>% bind_cols()
colnames(Group_numbers)[1] = "Groups"
Group_numbers <- Group_numbers[,-grep("\\.\\.\\.", colnames(Group_numbers)), with=FALSE]

# Pdcd1 vs Cd226 Contour
ggplot(ProcessingTable, aes(x = Pdcd1, y = Cd226)) + geom_density_2d(aes(color = seurat_clusters), size = 0.5, contour_var = "ndensity") +
  geom_vline(xintercept= mean_cutoff(Gene = "Pdcd1", Seurat = CD8T_FourClusters), linetype="dashed", color = "red", size=0.5)+
  geom_hline(yintercept= mean_cutoff(Gene = "Cd226", Seurat = CD8T_FourClusters), linetype="dashed", color = "red", size=0.5)+
  facet_wrap(~seurat_clusters) + labs(x = "Pdcd1", y = "Cd226") + scale_color_brewer(palette = "Dark2") + theme_minimal()+ NoLegend()


# Suppl. Figure S14a ------------------------------------------------------
library(viridis)
Idents(CD8T) <- factor(Idents(CD8T), levels = c(4,5,0,2,1,3))
SF14.Markers <- c("Lef1","Ccr7","Tcf7","Sell","Cxcr6","Nr4a1","Cd69","Runx3","Nr4a3","Tnfsf10","Gzma","Gzmb","Ifng","Cst7","Prf1","Nkg7","Ctla4","Tigit","Pdcd1","Lag3","Havcr2","Tnfrsf9","Icos","Cd28","Cd226","Hif1a","Tox","Id2","Hopx","Eomes","Tbx21","Zeb2","Ikzf2","Foxp3","Il2ra")
DoHeatmap(CD8T, features = SF14.Markers, draw.line = T, label = F, assay = "SCT", disp.min = -2, disp.max = 2, slot = "scale.data", lines.width = 10) + 
  scale_fill_viridis(na.value = "white")


# Suppl. Figure S16a ------------------------------------------------------
# We did not use pre-defined gene lists due to the absence of a logical framework to account for differences between WT and KO genotypes. 
# Instead, we performed an exploratory analysis of differentially expressed genes that represent WT-associated clusters (e.g., clusters 1 and 2) and KO-associated clusters (e.g., clusters 0 and 3).
# To achieve this, we compared each WT cluster against the pooled KO clusters and vice versa. 
# Genes were selected based on the following criteria: a fold change greater than 1.5, expression in more than 25% of cells, and an adjusted p-value less than 0.05.
CD8T_FourClusters %<>% PrepSCTFindMarkers(assay = "SCT", verbose = TRUE)

# Cluster 3 or 0 vs WT specific clusters (1 and 2)  ------------------------------------------------
# Cluster 3 vs Cluster 1+2
cd8_0vs12 <- FindMarkers(CD8T_FourClusters, ident.1 = 0, ident.2 = c(1,2), recorrect_umi = FALSE)
cd8_0vs12_Table <- data.table(symbol = rownames(cd8_0vs12), cd8_0vs12)
cd8_0vs12_Genes <- cd8_0vs12_Table %>% .[pct.1 > 0.25]

# Cluster 0 vs Cluster 1+2
cd8_3vs12 <- FindMarkers(sample.cd8, ident.1 = 3, ident.2 = c(1,2), recorrect_umi = FALSE)
cd8_3vs12_Table <- data.table(symbol = rownames(cd8_3vs12), cd8_3vs12)
cd8_3vs12_Genes <- cd8_3vs12_Table %>% .[pct.1 > 0.25]

# Suppl. Figure S16a left (at glance)
Table_0or3_verse_1and2 <- merge.data.table(cd8_0vs12_Genes, cd8_3vs12_Genes, by = "symbol", all.x = TRUE, all.y = TRUE)
plot(Table_0or3_verse_1and2$avg_log2FC.x, Table_0or3_verse_1and2$avg_log2FC.y, xlim= c(-2,2), ylim= c(-2,2))
write.table(Table_0or3_verse_1and2, "Table_0or3_verse_1and2.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


# cluster 1 or 2 vs VKO specific clusters (0 and 3) ------------------------------------------------
# Cluster 1 vs Cluster 0+3
cd8_1vs03 <- FindMarkers(CD8T_FourClusters, ident.1 = 1, ident.2 = c(0,3), recorrect_umi = FALSE)
cd8_1vs03_Table <- data.table(symbol = rownames(cd8_1vs03), cd8_1vs03)
cd8_1vs03_Genes <- cd8_1vs03_Table %>% .[pct.1 > 0.25]

# Cluster 2 vs Cluster 0+3
cd8_2vs03 <- FindMarkers(CD8T_FourClusters, ident.1 = 2, ident.2 = c(0,3), recorrect_umi = FALSE)
cd8_2vs03_Table <- data.table(symbol = rownames(cd8_2vs03), cd8_2vs03)
cd8_2vs03_Genes <- cd8_2vs03_Table %>% .[pct.1 > 0.25]

# Suppl. Figure S16a right (at glance)
Table_1or2_verse_0and3 <- merge.data.table(cd8_1vs03_Genes, cd8_2vs03_Genes, by = "symbol", all.x = TRUE, all.y = TRUE)
plot(Table_1or2_verse_0and3$avg_log2FC.x, Table_1or2_verse_0and3$avg_log2FC.y, xlim= c(-2,4), ylim= c(-4,4))
write.table(Table_1or2_verse_0and3, "Table_1or2_verse_0and3.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


# Suppl. Figure S16b ----------------------------------------------------
# Because cluster 4 and 5 have same freq between WT and Vsir KO, we choose cluster 0,1 (WT specific) and 2,3 (Vsir-/- specific)
CD8T_FourClusters <- subset(CD8T, subset = seurat_clusters %in% c(0,1,2,3))
Idents(CD8T_FourClusters) <- factor(CD8T_FourClusters$seurat_clusters, levels = c(1,2,0,3))