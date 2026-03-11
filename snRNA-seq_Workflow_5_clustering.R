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
library(optparse)

option_list = list(
  make_option(c("-f", "--seurat"), type="character", default=NULL,
              help="Path to seurat object, outfile from previous script", metavar="character"),
    make_option(c("-d", "--dims"), type="integer", default=60,
              help="dimensions of reduction to use for FindNeighbours, based on PC cut off from elbow plot (script 3)", metavar="integer"),
    make_option(c("-r", "--resolution"), type="character", default="1.6",
              help="Resolution to pass to FindClusters - run FindClusters across range of resolutions script to find suitable value", metavar="character"),
    make_option(c("-o", "--outfile"), type="character", default="NULL",
              help="name of output SeuratObject", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Assigning variables - read through script before assigning to ensure correct value assignment

seuratpath <- opt$seurat
dimsval <- opt$dims
resval <- as.numeric(opt$resolution)
outfile <- opt$outfile

SeuratObject <- readRDS(seuratpath)

#dims selected based on PC cut off from elbow plot

SeuratObject <- FindNeighbors(SeuratObject, reduction = "harmony", dims = 1:dimsval)

#algorithm = 4 == Leiden algorithm. It is recommended that a resolution is selected after running clustering tests across a range of resolutions

SeuratObject <- FindClusters(SeuratObject, resolution = resval, algorithm = 4, method = "igraph")
SeuratObject <- RunUMAP(SeuratObject, reduction = "harmony", dims = 1:dimsval)

#Mitochondrially encoded genes removed from downstream analyses due to association of increased level of mitochondrial reads with poor quality cells - if there are differences in pathways associated with the mitochondria, we want to be be certain this is due to disease related changes and not because there are more damaged/poor quality cells in the disease state

SeuratObject <- SeuratObject[-grep("MT-",rownames(SeuratObject@assays$RNA)),]

#Write cells/samples represented by each cluster for cluster QC observations and save the Seurat Object

cells_per_cluster_table <- as.data.frame(table(SeuratObject$seurat_clusters))

colnames(cells_per_cluster_table) <- c("Cluster","Number of Cells")

write.table(cells_per_cluster_table, file = "Files/Cells_Per_Unassigned_Cluster_Table.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

samples_per_cluster_table <- group_by(SeuratObject@meta.data, seurat_clusters) %>% summarise(Sample_Count = length(unique(sample_id)))

colnames(samples_per_cluster_table) <- c("Cluster","Number of Samples")

write.table(samples_per_cluster_table, file = "Files/Samples_Per_Unassigned_Cluster_Table.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

saveRDS(SeuratObject,outfile)


