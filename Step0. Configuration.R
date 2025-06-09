############################################################################
# File 0. Environment Configuration ----------------------------------------
############################################################################

##### Environment configuration #####
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(magrittr)
library(tidyverse)
library(Seurat.utils)

# Custom scripts
# 1. Convert Seurat to Monocle CDS object
Monocle3_seurat_convert3 = function(population_seurat, resolution){
  require(SeuratWrappers)
  require(monocle3)
  require(dplyr)
  require(Seurat)
  require(magrittr)
  
  input_subset_seurat = population_seurat
  cds_monocle3 = SeuratWrappers::as.cell_data_set(input_subset_seurat)
  set.seed(39)
  cds_monocle3 <- cluster_cells(cds_monocle3, resolution= resolution, reduction_method = "UMAP")
  cds_monocle3 <- learn_graph(cds_monocle3, use_partition = TRUE, close_loop = FALSE, verbose = TRUE)
  cds_monocle3 <- order_cells(cds_monocle3) # pop an window and set the root range of cells
  return(cds_monocle3)
  
}

# 2. Calculate gene expression level according to pseudotime points
Genexpression_timepoints =  function(seuratobjt, genesymbol, divider){
  timegene = data.table(time = floor(seuratobjt$pseudotime/divider), gene = seuratobjt[genesymbol, ]@assays$SCT$data)
  timelvls = sort(timegene$time) %>% unique
  timetabl =  lapply(timelvls, function(lvl){
    time_ealvl = timegene[time == lvl]
    time_summy = summary(time_ealvl$gene)
    time_mean  = mean(time_ealvl$gene)
    time_se    = sd(time_ealvl$gene)/sqrt(length(time_ealvl$gene))
    return(data.table(mean = time_mean, up_2sd = time_mean+time_se*2, dn_2sd = time_mean-time_se*2))
  }) %>% rbindlist
  return(timetabl)
}

# 3. Calculate frequency/number for given Seurat object according to genotype using cut-off value
Calculate_TwoGene_Expressing = function(Original, Seurat, Gene1, Gene2, method = NULL, Mode){
  Orig_Gene1 = Original[Gene1,]@assays$SCT$data
  Orig_Gene2 = Original[Gene2,]@assays$SCT$data
  
  Expr_Gene1  = Seurat[Gene1,]@assays$SCT$data
  Expr_Gene2  = Seurat[Gene2,]@assays$SCT$data
  
  Label.A = paste0(Gene1,"+")
  Label.B = paste0(Gene1,"-")
  Label.C = paste0(Gene2,"+")
  Label.D = paste0(Gene2,"-")
  
  if(is.null(method)){
    Expr_Table  = paste0(ifelse(Expr_Gene1 > 0, Label.A, Label.B), ifelse(Expr_Gene2  > 0, Label.C, Label.D)) %>% table
  }else{
    if(method == "Median"){
      Expr_Table  = paste0(ifelse(Expr_Gene1 > median(Orig_Gene1), Label.A, Label.B), ifelse(Expr_Gene2  > median(Orig_Gene2),Label.C, Label.D)) %>% table
    }else if(method == "Mean"){
      Expr_Table  = paste0(ifelse(Expr_Gene1 > mean(Orig_Gene1),Label.A, Label.B),   ifelse(Expr_Gene2  > mean(Orig_Gene2),Label.C, Label.D)) %>% table
    }else{
      stop("@@ Type proper in methods; Median, Mean, or NULL")
    }
  }
  if(Mode == "Freq"){
    return(Expr_Table/sum(Expr_Table))
  }else if(Mode == "Number"){
    return(Expr_Table)
  }else if(Mode == "GeneRatio"){
    count_Gene1 = sum(Expr_Table[names(Expr_Table) %in% c(paste0(Label.A,Label.C), paste0(Label.A,Label.D))])
    count_Gene2 = sum(Expr_Table[names(Expr_Table) %in% c(paste0(Label.A,Label.C), paste0(Label.B,Label.C))])
    return(count_Gene1/count_Gene2)
  }
}

# 3. Calculate frequency/number for given Seurat object according to genotype
Calculate_TwoGenotype = function(Seurat, Mode){
  genotype_cluster = Seurat$sample %>% table
  if(Mode == "Freq"){
    return(genotype_cluster/sum(genotype_cluster))
  }else if(Mode == "Number"){
    return(genotype_cluster)
  }
}


# 4. Calculate cut-off value for contour graph
mean_cutoff = function(Gene, Seurat){
  Genes_in_Seurat = Seurat[Gene,]@assays$SCT@data %>% as.numeric
  return(round(mean(Genes_in_Seurat),1))
}
