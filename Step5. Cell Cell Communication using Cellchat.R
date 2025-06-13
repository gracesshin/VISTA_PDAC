############################################################################
# File 5. Cell Cell Communication using Cellchat  --------------------------
# Estimate Cell-Cell communication (CCC) using CellChat R-package ----------
############################################################################

# Environment configuration
library(dplyr)
library(data.table)
library(CellChat)
library(magrittr)
future::plan("multisession", workers = 64)

# Loading Dataset
MoMx <- readRDS("Pan02_MoMx_WT_and_KO.RDS")
CD8T <- readRDS("Pan02_CD8T_WT_and_KO.RDS")
CellChatDB <- CellChatDB.mouse 

# Attach Cell ID to Subset/Seurat_clusters
MoMx %<>% AddMetaData(metadata = paste0("MoMx_c",MoMx$seurat_clusters), col.name = "groups")
CD8T %<>% AddMetaData(metadata = paste0("CD8T_c",CD8T$seurat_clusters), col.name = "groups")

MoMx.CD8T <- merge(x = MoMx, y = CD8T)
MoMx.CD8T$groups <- as.factor(MoMx.CD8T$groups)

# Create & Adjust CellChat object
CellChat <- createCellChat(object = MoMx.CD8T, group.by = "groups", assay = "SCT")
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
CellChat@DB <- CellChatDB.use

# Processing CellChat
CellChat %<>% subsetData()  # Required process
CellChat %<>% identifyOverExpressedGenes()
CellChat %<>% identifyOverExpressedInteractions()
CellChat %<>% computeCommunProb(type = "triMean") 
CellChat %<>% filterCommunication(min.cells = 10) 
CellChat %<>% computeCommunProbPathway()

# > CellChat@idents %>% levels
# "CD8T_c0" "CD8T_c1" "CD8T_c2" "CD8T_c3" "CD8T_c4" "CD8T_c5" "MoMx_c0" "MoMx_c1" "MoMx_c2" "MoMx_c3" "MoMx_c4" "MoMx_c5" "MoMx_c6" "MoMx_c7"

# CD8+ T cell 
# CD8T_c1 CD8T_c2 [WT specific CD8+T cluster]       : Level 2,3
# CD8T_c0 CD8T_c3 [Vsir-/- specific CD8+T cluster]  : Level 1,4

# Macrophage 
# MoMx_c1 MoMx_c3 [WT specific macrophage cluster]  : Level 8,10
# MoMx_c2 [Vsir-/- specific macrophage cluster]     : Level 9

# Extract Genotype specific CCC with P-value < 0.05
Macrophage_to_CD8 = subsetCommunication(CellChat, slot.name = "net", sources.use = 8:10, targets.use = 1:4, thresh = 1) %>% as.data.table %>% .[pval < 0.05]
CD8_to_Macrophage = subsetCommunication(CellChat, slot.name = "net", sources.use = 1:4, targets.use = 8:10, thresh = 1) %>% as.data.table %>% .[pval < 0.05]

# Finalize interaction-probability dataset
# Low expressed genes were removed after this process

Macrophage_to_CD8_Export <- list(
  data.table(Genotype = "VKO-specific", Macrophage_to_CD8[source == "MoMx_c2" & target == "CD8T_c0"]),
  data.table(Genotype = "VKO-specific", Macrophage_to_CD8[source == "MoMx_c2" & target == "CD8T_c3"]),
  data.table(Genotype = "WT-specific", Macrophage_to_CD8[source == "MoMx_c1" & target == "CD8T_c1"]),
  data.table(Genotype = "WT-specific", Macrophage_to_CD8[source == "MoMx_c1" & target == "CD8T_c2"]),
  data.table(Genotype = "WT-specific", Macrophage_to_CD8[source == "MoMx_c3" & target == "CD8T_c1"]),
  data.table(Genotype = "WT-specific", Macrophage_to_CD8[source == "MoMx_c3" & target == "CD8T_c2"])
) %>% rbindlist
write.table(Macrophage_to_CD8_Export, "Macrophage_to_CD8_Export.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

CD8_to_Macrophage_Export <- list(
  data.table(Genotype = "VKO-specific", CD8_to_Macrophage[source == "CD8T_c0" & target == "MoMx_c2"]),
  data.table(Genotype = "VKO-specific", CD8_to_Macrophage[source == "CD8T_c3" & target == "MoMx_c2"]),
  data.table(Genotype = "WT-specific", CD8_to_Macrophage[source == "CD8T_c1" & target == "MoMx_c1"]),
  data.table(Genotype = "WT-specific", CD8_to_Macrophage[source == "CD8T_c2" & target == "MoMx_c1"]),
  data.table(Genotype = "WT-specific", CD8_to_Macrophage[source == "CD8T_c1" & target == "MoMx_c3"]),
  data.table(Genotype = "WT-specific", CD8_to_Macrophage[source == "CD8T_c2" & target == "MoMx_c3"])
) %>% rbindlist
write.table(CD8_to_Macrophage_Export, "CD8_to_Macrophage_Export.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


