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
library(tidyr)
library(scDblFinder)
library(optparse)

option_list = list(
  make_option(c("-f", "--seurat"), type="character", default=NULL,
              help="Path to seurat object, outfile from script 4", metavar="character"),
    make_option(c("-d", "--dims"), type="integer", default=60,
              help="dimensions of reduction to used for FindNeighbours in previous script", metavar="integer"),
    make_option(c("-c", "--mincells"), type="integer", default=200,
              help="minimum number of cells required for retention of cluster, default based on ~600,000 cells, scale accordingly", metavar="integer"),
    make_option(c("-s", "--minsamples"), type="integer", default=20,
              help="minimum number of samples required for retention of cluster, default based on 20 samples, scale accordingly", metavar="integer"),
    make_option(c("-p", "--pthreshold"), type="integer", default=30,
              help="doublet percentage threshold - clusters that have a doublet cell percentage that equals or exceeds this threshold will be excluded", metavar="integer"),
    make_option(c("-o", "--outfile"), type="character", default="NULL",
              help="name of output SeuratObject", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Assigning variables - read through script before assigning to ensure correct value assignment

seuratpath <- opt$seurat
dimsval <- opt$dims
mincells <- opt$mincells
minsamples <- opt$minsamples
pthreshold <- opt$pthreshold
outfile <- opt$outfile

dir.create("Figures/Cluster_QC")

SeuratObject <- readRDS(seuratpath)

CellNumberQC <- read.delim("Files/Cells_Per_Unassigned_Cluster_Table.tsv")

SampleNumberQC <- read.delim("Files/Samples_Per_Unassigned_Cluster_Table.tsv")

#Cluster filtering - clusters need to represent at least [mincells] cells and more than [minsamples] samples to be retained - default cut offs in Assigning Variables section based on approx. 600,000 cells and ~100 samples, scale for size ofown dataset accordingly

CellNumberQCFilter <- CellNumberQC[CellNumberQC$Number.of.Cells >= mincells,]

SampleNumberQCFilter <- SampleNumberQC[SampleNumberQC$Number.of.Samples >= minsamples,]

SeuratObject <- subset(SeuratObject, subset = seurat_clusters %in% CellNumberQCFilter$Cluster & seurat_clusters %in% SampleNumberQCFilter$Cluster)

SeuratObject@meta.data$seurat_clusters <- factor(as.numeric(as.vector(SeuratObject@meta.data$seurat_clusters)), levels = sort(unique(as.numeric(as.vector(SeuratObject@meta.data$seurat_clusters)))))

#Converting to single cell experiment object to run find doublets using the scDblFinder package using our current set of clusters

SeuratObject.sce <- as.SingleCellExperiment(SeuratObject)
SeuratObject.sce <- scDblFinder(SeuratObject.sce, samples = "sample_id",clusters = "seurat_clusters")

barcode_doublet_Table <- as.data.frame(cbind(as.vector(rownames(SeuratObject.sce@colData)),as.vector(SeuratObject.sce@colData$scDblFinder.class),as.vector(SeuratObject.sce@colData$scDblFinder.score)))

#Adding doublet information to Seurat Object metadata and making UMAPs

SeuratObject@meta.data$class <- barcode_doublet_Table$V2

SeuratObject@meta.data$dblscore <- barcode_doublet_Table$V3

Doublet_UMAP <- DimPlot(SeuratObject, label = FALSE, repel = TRUE, pt.size = 0, label.size = 2.5, group.by = "class") + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(Doublet_UMAP, filename = "Figures/Cluster_QC/Unassigned_Doublet_UMAPclusters_scRNA_seq.pdf", device = "pdf", width = 6, height = 4, units = "in")

metadata <- SeuratObject@meta.data

doublet_counts <- metadata %>% dplyr::count(seurat_clusters, class)

#Building count table of doublets and singlets and calculating doublet % and mean doublet score for each cluster

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
      temptable2[2] <- "singlet"
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

#Subsetting Seurat Object to retain only singlet cells AND those clusters where less than [pthreshold] of constituent cells are doublets (default - 30%) and writing tables of cells/samples per clusters after filtering

singlet_clusters <- as.numeric(as.vector(doublet_percentage$seurat_clusters[doublet_percentage$doublet_percent < pthreshold]))

SeuratObjectSinglets <- subset(SeuratObject, subset = class == "singlet" & seurat_clusters %in% singlet_clusters)

#Reapplying cluster QC criteria as above after removal of doublets/doublet clusters, saving QC related figures and tables and finally the cleaned up Seurat Object

cells_per_cluster_table <- as.data.frame(table(SeuratObjectSinglets$seurat_clusters))

colnames(cells_per_cluster_table) <- c("Cluster","Number.of.Cells")

write.table(cells_per_cluster_table, file = "Files/Cells_Per_Unassigned_Cluster_Table_PostClusterQC.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

samples_per_cluster_table <- group_by(SeuratObjectSinglets@meta.data, seurat_clusters) %>% summarise(Sample_Count = length(unique(sample_id)))

colnames(samples_per_cluster_table) <- c("Cluster","Number.of.Samples")

write.table(samples_per_cluster_table, file = "Files/Samples_Per_Unassigned_Cluster_Table_PostClusterQC.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

CellNumberQCFilter <- cells_per_cluster_table[cells_per_cluster_table$Number.of.Cells >= mincells,]

SampleNumberQCFilter <- samples_per_cluster_table[samples_per_cluster_table$Number.of.Samples >= minsamples,]

SeuratObjectSinglets <- subset(SeuratObjectSinglets, subset = seurat_clusters %in% CellNumberQCFilter$Cluster & seurat_clusters %in% SampleNumberQCFilter$Cluster)

SeuratObjectSinglets@meta.data$seurat_clusters <- factor(as.numeric(as.vector(SeuratObjectSinglets@meta.data$seurat_clusters)), levels = sort(unique(as.numeric(as.vector(SeuratObjectSinglets@meta.data$seurat_clusters)))))

ncount_per_cluster <- VlnPlot(SeuratObject, features = c("nCount_RNA"), split.by = "seurat_clusters")

ggsave(ncount_per_cluster, filename = "Figures/Cluster_QC/nCount_per_Cluster_ViolinPlot.pdf", device = "pdf", width = 8, height = 6, units = "in")

nfeature_per_cluster <- VlnPlot(SeuratObject, features = c("nFeature_RNA"), split.by = "seurat_clusters")

ggsave(nfeature_per_cluster, filename = "Figures/Cluster_QC/nFeature_per_Cluster_ViolinPlot.pdf", device = "pdf", width = 8, height = 6, units = "in")

#Re-running UMAP and saving Seurat Object to clean up after doublet/doublet cluster removal

SeuratObjectSinglets <- RunUMAP(SeuratObjectSinglets, reduction = "harmony", dims = 1:dimsval)

meta <- SeuratObjectSinglets@meta.data

write.table(meta, file = "Files/PostClusterQCMetadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(SeuratObjectSinglets,outfile)


