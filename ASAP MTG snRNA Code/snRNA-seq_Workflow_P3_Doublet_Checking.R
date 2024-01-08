library(Seurat)
library(cowplot)
library(ggplot2)
library(ggsci)
library(dplyr)
library(psych)
library(pheatmap)
library(harmony)
#library(clusterProfiler)
#library(org.Hs.eg.db)
library(DOSE)
library(GOSemSim)
library(enrichplot)
library(stringr)
library(reshape2)
library(sciplot)
library(tidyr)
library(scDblFinder)

dir.create("Figures/Cluster_QC")

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegration_Oct2023_FinalOSR_6HC_MTG_Part2.rds")

CellNumberQC <- read.delim("Files/Cells_Per_Unassigned_Cluster_Table.tsv")

SampleNumberQC <- read.delim("Files/Samples_Per_Unassigned_Cluster_Table.tsv")

CellNumberQCFilter <- CellNumberQC[CellNumberQC$Number.of.Cells >= 200,]

SampleNumberQCFilter <- SampleNumberQC[SampleNumberQC$Number.of.Samples >= 20,]

SeuratObject <- subset(SeuratObject, subset = seurat_clusters %in% CellNumberQCFilter$Cluster & seurat_clusters %in% SampleNumberQCFilter$Cluster)

SeuratObject@meta.data$seurat_clusters <- factor(as.numeric(as.vector(SeuratObject@meta.data$seurat_clusters)), levels = sort(unique(as.numeric(as.vector(SeuratObject@meta.data$seurat_clusters)))))

SeuratObject.sce <- as.SingleCellExperiment(SeuratObject)
SeuratObject.sce <- scDblFinder(SeuratObject.sce, samples = "sample_id",clusters = "seurat_clusters")

barcode_doublet_Table <- as.data.frame(cbind(as.vector(rownames(SeuratObject.sce@colData)),as.vector(SeuratObject.sce@colData$scDblFinder.class),as.vector(SeuratObject.sce@colData$scDblFinder.score)))

SeuratObject@meta.data$class <- barcode_doublet_Table$V2

SeuratObject@meta.data$dblscore <- barcode_doublet_Table$V3

Doublet_UMAP <- DimPlot(SeuratObject, label = FALSE, repel = TRUE, pt.size = 0, label.size = 2.5, group.by = "class") + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(Doublet_UMAP, filename = "Figures/Cluster_QC/Unassigned_Doublet_UMAPclusters_scRNA_seq.pdf", device = "pdf", width = 6, height = 4, units = "in")

singlet_barcodes <- rownames(SeuratObject@meta.data)[as.vector(SeuratObject@meta.data$class) %in% "singlet"]

Doublets_Removed_UMAP <- DimPlot(SeuratObject, label = TRUE, repel = TRUE, pt.size = 0, label.size = 2.5, cells = singlet_barcodes) + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(Doublets_Removed_UMAP, filename = "Figures/Cluster_QC/Unassigned_DoubletsRemoved_UMAPclusters_scRNA_seq.pdf", device = "pdf", width = 6, height = 4, units = "in")

metadata <- SeuratObject@meta.data

doublet_counts <- metadata %>% count(seurat_clusters, class)

cor_doublet_counts <- data.frame(matrix(ncol = 3, nrow = 0))

colnames(cor_doublet_counts) <- colnames(doublet_counts)

for (i in unique(doublet_counts$seurat_clusters)) {
  temptable <- doublet_counts[doublet_counts$seurat_clusters %in% i,]
  if (length(temptable$seurat_clusters) == 2){
    cor_doublet_counts <- rbind(cor_doublet_counts, temptable)
  }
  else {
    cell_class <- temptable$class[1]
    if (cell_class == "singlet") {
      temptable2 <- temptable
      temptable2[2] <- "doublet"
      temptable2[3] <- 0
      cor_doublet_counts <- rbind(cor_doublet_counts, temptable, temptable2)
    }
    else if (cell_class == "doublet") {
      temptable2 <- temptable
      temptable2[2] <- "doublet"
      temptable2[3] <- 0
      cor_doublet_counts <- rbind(cor_doublet_counts, temptable, temptable2)
  }
 }
}

doublet_counts <- cor_doublet_counts

doublet_percentage <- doublet_counts %>% group_by(seurat_clusters) %>% mutate(doublet_percent = n[which(class == "doublet")]/sum(n)*100, total_cells = sum(n))

doublet_percentage <- doublet_percentage[doublet_percentage$class == "doublet",]

dblscoretable <- metadata %>% group_by(seurat_clusters) %>% summarise(avgdblscore = mean(as.numeric(dblscore)))

doublet_percentage$avgdblscore <- dblscoretable$avgdblscore

write.table(doublet_percentage, file = "Files/Doublet_Percentage_Per_Cluster_Table.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

singlet_clusters <- as.numeric(as.vector(doublet_percentage$seurat_clusters[doublet_percentage$doublet_percent < 30]))

SeuratObjectSinglets <- subset(SeuratObject, subset = class == "singlet" & seurat_clusters %in% singlet_clusters)

cells_per_cluster_table <- as.data.frame(table(SeuratObjectSinglets$seurat_clusters))

colnames(cells_per_cluster_table) <- c("Cluster","Number.of.Cells")

write.table(cells_per_cluster_table, file = "Files/Cells_Per_Unassigned_Cluster_Table_Singlets.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

samples_per_cluster_table <- group_by(SeuratObjectSinglets@meta.data, seurat_clusters) %>% summarise(Sample_Count = length(unique(sample_id)))

colnames(samples_per_cluster_table) <- c("Cluster","Number.of.Samples")

write.table(samples_per_cluster_table, file = "Files/Samples_Per_Unassigned_Cluster_Table_Singlets.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

CellNumberQCFilter <- cells_per_cluster_table[cells_per_cluster_table$Number.of.Cells >= 200,]

SampleNumberQCFilter <- samples_per_cluster_table[samples_per_cluster_table$Number.of.Samples >= 20,]

SeuratObjectSinglets <- subset(SeuratObjectSinglets, subset = seurat_clusters %in% CellNumberQCFilter$Cluster & seurat_clusters %in% SampleNumberQCFilter$Cluster)

SeuratObjectSinglets@meta.data$seurat_clusters <- factor(as.numeric(as.vector(SeuratObjectSinglets@meta.data$seurat_clusters)), levels = sort(unique(as.numeric(as.vector(SeuratObjectSinglets@meta.data$seurat_clusters)))))

Singlet_UMAP <- DimPlot(SeuratObjectSinglets, label = TRUE, repel = TRUE, pt.size = 0, label.size = 2.5) + 
    theme(axis.text = element_text(size=8),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 8),
          title = element_text(size = 12),
          legend.key.size = unit(0.4,"cm"))

ggsave(Singlet_UMAP, filename = "Figures/Cluster_QC/Unassigned_Singlets_UMAPclusters_scRNA_seq.pdf", device = "pdf", width = 6, height = 4, units = "in")

ncount_per_cluster <- VlnPlot(SeuratObject, features = c("nCount_RNA"), split.by = "seurat_clusters")

ggsave(ncount_per_cluster, filename = "Figures/Cluster_QC/nCount_per_Cluster_ViolinPlot.pdf", device = "pdf", width = 8, height = 6, units = "in")
mean_ncounts <- metadata %>% group_by(seurat_clusters) %>% summarise(mean = mean(nCount_RNA), SE = se(nCount_RNA))

median_ncounts <- metadata %>% group_by(seurat_clusters) %>% summarise(median = median(nCount_RNA), SE = se(nCount_RNA))
VlnPlot(SeuratObject, features = c("nFeature_RNA"), ncol = 3, split.by = "seurat_clusters")

nfeature_per_cluster <- VlnPlot(SeuratObject, features = c("nFeature_RNA"), split.by = "seurat_clusters")

ggsave(nfeature_per_cluster, filename = "Figures/Cluster_QC/nFeature_per_Cluster_ViolinPlot.pdf", device = "pdf", width = 8, height = 6, units = "in")

ClusterTable <- as.data.frame(unique(SeuratObjectSinglets$seurat_clusters))

ClusterTable <- as.data.frame(ClusterTable[order(ClusterTable),])

colnames(ClusterTable) <- "Seurat_Clusters"

write.table(ClusterTable, file = "Files/PostQCFiltering_SeuratClustersTable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

cluster_assignment_table <- as.data.frame(cbind(rownames(SeuratObjectSinglets@meta.data),as.numeric(as.vector(SeuratObjectSinglets@meta.data$seurat_clusters))))

colnames(cluster_assignment_table) <- c("Barcode","Cluster")

cluster_assignment_table$Cluster <- as.numeric(cluster_assignment_table$Cluster)

write.table(cluster_assignment_table, file = "Files/Cell_to_Cluster_Assignment_Table_Initial.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

saveRDS(SeuratObjectSinglets,"/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegration_Oct2023_FinalOSR_6HC_MTG_Part3.rds")


