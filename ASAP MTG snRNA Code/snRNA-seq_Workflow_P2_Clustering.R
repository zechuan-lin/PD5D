

library(Seurat)
library(cowplot)
library(ggplot2)
library(ggsci)
library(dplyr)
library(psych)
library(pheatmap)
library(harmony)
library(DOSE)
library(GOSemSim)
library(enrichplot)
library(stringr)
library(reshape2)
library(sciplot)

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegration_Oct2023_FinalOSR_6HC_MTG_Part1.rds")

SeuratObject <- FindNeighbors(SeuratObject, reduction = "harmony", dims = 1:60)
SeuratObject <- FindClusters(SeuratObject, resolution = 1.5, algorithm = 4, method = "igraph")
SeuratObject <- RunUMAP(SeuratObject, reduction = "harmony", dims = 1:60)

SeuratObject <- SeuratObject[-grep("MT-",rownames(SeuratObject@assays$RNA)),]

metadataDF <- SeuratObject@meta.data

write.table(metadataDF, file = "Files/PostFilteringMetadata_Table.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

cells_per_cluster_table <- as.data.frame(table(SeuratObject$seurat_clusters))

colnames(cells_per_cluster_table) <- c("Cluster","Number of Cells")

write.table(cells_per_cluster_table, file = "Files/Cells_Per_Unassigned_Cluster_Table.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

samples_per_cluster_table <- group_by(SeuratObject@meta.data, seurat_clusters) %>% summarise(Sample_Count = length(unique(sample_id)))

colnames(samples_per_cluster_table) <- c("Cluster","Number of Samples")

write.table(samples_per_cluster_table, file = "Files/Samples_Per_Unassigned_Cluster_Table.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

saveRDS(SeuratObject,"/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegration_Oct2023_FinalOSR_6HC_MTG_Part2.rds")

SeuratObject_UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", label = TRUE, pt.size = 0.01, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(SeuratObject_UMAP_Clusters, filename = "Figures/SeuratObject_UMAP_Clusters.pdf", device = "pdf", width = 6, height = 4, units = "in")

SeuratObject_Case_Group_UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", group.by = "case", pt.size = 0.1, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(SeuratObject_Case_Group_UMAP_Clusters, filename = "Figures/SeuratObject_Case_Group_UMAP_Clusters.pdf", device = "pdf", width = 6, height = 4, units = "in")

SeuratObject_Region_Group_UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", group.by = "batch",pt.size = 0.1, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(SeuratObject_Region_Group_UMAP_Clusters, filename = "Figures/SeuratObject_Region_Group_UMAP_Clusters.pdf", device = "pdf", width = 6, height = 4, units = "in")

SeuratObject_Case_Split_UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", split.by = "case", label = TRUE, ncol = 1, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(SeuratObject_Case_Split_UMAP_Clusters, filename = "Figures/SeuratObject_Case_Split_UMAP_Clusters.pdf", device = "pdf", width = 6, height = 8, units = "in")

SeuratObject_Region_Split_UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", split.by = "batch", label = TRUE, ncol = 1, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(SeuratObject_Region_Split_UMAP_Clusters, filename = "Figures/SeuratObject_Region_Split_UMAP_Clusters.pdf", device = "pdf", width = 6, height = 32, units = "in")

MarkerGenes <- c("ENO2","RBFOX3","SLC17A6","SLC17A7","SLC32A1","GAD1","GAD2","AQP4","GFAP","PLP1","MBP","OLIG1","VCAN","BCAN","CX3CR1","P2RY12","FLT1","CLDN5","IL7R","CD96","CD8A","RELN","CALB2","CNR1","PTPRC","CD27","CD28","CCR7","PRF1")

data_barplot <- FetchData(SeuratObject, vars = c("ident",MarkerGenes), slot = "data")

data_barplot_melt <- melt(data_barplot)

data_barplot_melt$ident <- as.vector(data_barplot_melt$ident)
data_barplot_melt$variable <- as.vector(data_barplot_melt$variable)
data_barplot_melt$value <- as.numeric(as.vector(data_barplot_melt$value))

data_barplot_melt_sum <- group_by(data_barplot_melt,ident,variable) %>% summarise(mean = mean(value), SE = se(value))

data_barplot_melt_sum$ident <- factor(data_barplot_melt_sum$ident, levels = unique(data_barplot_melt_sum$ident))

data_barplot_melt_sum$variable <- factor(data_barplot_melt_sum$variable, levels = unique(MarkerGenes))

Batch1to22_barchart <- ggplot(data_barplot_melt_sum, aes(x = ident, y = mean, fill = ident)) + 
  geom_bar(aes(x = ident, y = mean), stat = "identity", alpha = 1) + 
  geom_errorbar(aes(x = ident, ymin = mean-SE, ymax = mean+SE, colour = ident), width = 0.4, alpha = 0.9, size = 0.5) + 
  ggplot2::facet_grid(rows = vars(variable), scales = "free_y", switch = "y") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 12, angle = 45, face = "bold", vjust = 0.5),
        axis.text.y = element_blank(), axis.ticks = element_blank(), panel.background = element_blank(),
        strip.background = element_blank(), strip.placement = "outside", 
        strip.text.y = element_text(size = 12, angle = 180, face = "bold"),
        strip.text.y.left = element_text(angle = 0)) + NoLegend()


ggsave(Batch1to22_barchart,filename = "Files/Marker_Barchart", device = "pdf", width = 12, height = 16, units = "in")


Markerggplots <- function(SeurObj,Genes){
  for (i in Genes) {
    TempViolin <- VlnPlot(SeurObj, features = i ,pt.size = 0)
    ggsave(TempViolin, filename = paste("Figures/",i,"_VlnPlot.pdf",sep = ""), device = "pdf", width = 12, height = 4, units = "in")
  }}

Markerggplots(SeuratObject,MarkerGenes)



Markerggplotspt1 <- function(SeurObj,Genes){
  for (i in Genes) {
    TempViolin <- VlnPlot(SeurObj, features = i ,pt.size = 1)
    ggsave(TempViolin, filename = paste("Figures/",i,"_Pt1_VlnPlot.pdf",sep = ""), device = "pdf", width = 12, height = 4, units = "in")
  }}

Markerggplotspt1(SeuratObject,MarkerGenes)


