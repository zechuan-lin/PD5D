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
              help="Path to seurat object, output from previous script", metavar="character"),
  make_option(c("-m", "--markers"), type="character", default=NULL,
              help="Path to plain text file of marker genes to visualise expression of, no header, one per line", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Assigning variables - read through script before assigning to ensure correct value assignment

seuratpath <- opt$seurat
markerfile <- opt$markers

#Characterising the clusters after QC

SeuratObject <- readRDS(seuratpath)

markerdf <- read.delim(markerfile, header = FALSE)

MarkerGenes <- markerdf[[1]]

#Making UMAPs after cluster QC and re-running RunUMAP

SeuratObject_UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", label = TRUE, pt.size = 0.01, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(SeuratObject_UMAP_Clusters, filename = "Figures/SeuratObject_UMAP_Clusters_PostClusterQC.pdf", device = "pdf", width = 6, height = 4, units = "in")

#Splitting UMAP apart by case, case can be altered to variable of choice.

SeuratObject_Case_Split_UMAP_Clusters <- DimPlot(SeuratObject, reduction = "umap", split.by = "case", label = TRUE, ncol = 1, label.size=2.5, repel = TRUE) + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(SeuratObject_Case_Split_UMAP_Clusters, filename = "Figures/SeuratObject_Case_Split_UMAP_Clusters_PostClusterQC.pdf", device = "pdf", width = 6, height = 8, units = "in")

data_barplot <- FetchData(NormalizeData(SeuratObject, normalization.method = "RC",scale.factor = 10000), vars = c("ident",MarkerGenes), slot = "data")

data_barplot_melt <- melt(data_barplot)

data_barplot_melt$ident <- as.vector(data_barplot_melt$ident)
data_barplot_melt$variable <- as.vector(data_barplot_melt$variable)
data_barplot_melt$value <- as.numeric(as.vector(data_barplot_melt$value))

data_barplot_melt_sum <- group_by(data_barplot_melt,ident,variable) %>% summarise(mean = mean(value), SE = se(value))

data_barplot_melt_sum$ident <- factor(data_barplot_melt_sum$ident, levels = unique(data_barplot_melt_sum$ident))

data_barplot_melt_sum$variable <- factor(data_barplot_melt_sum$variable, levels = unique(MarkerGenes))

barchart <- ggplot(data_barplot_melt_sum, aes(x = ident, y = mean, fill = ident)) + 
  geom_bar(aes(x = ident, y = mean), stat = "identity", alpha = 1) + 
  geom_errorbar(aes(x = ident, ymin = mean-SE, ymax = mean+SE, colour = ident), width = 0.4, alpha = 0.9, size = 0.5) + 
  ggplot2::facet_grid(rows = vars(variable), scales = "free_y", switch = "y") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 12, angle = 45, face = "bold", vjust = 0.5),
        axis.text.y = element_blank(), axis.ticks = element_blank(), panel.background = element_blank(),
        strip.background = element_blank(), strip.placement = "outside", 
        strip.text.y = element_text(size = 12, angle = 180, face = "bold"),
        strip.text.y.left = element_text(angle = 0)) + NoLegend()

ggsave(barchart,filename = paste("Files/Marker_Barchart.pdf",sep = ""), device = "pdf", width = 12, height = 18, units = "in")

Markerggplots <- function(SeurObj,Genes){
  for (i in Genes) {
    TempViolin <- VlnPlot(SeurObj, features = i ,pt.size = 0)
    ggsave(TempViolin, filename = paste("Figures/",i,"_VlnPlot.pdf",sep = ""), device = "pdf", width = 18, height = 4, units = "in")
  }}

Markerggplots(SeuratObject,MarkerGenes)

#Seurat clusters can now be assigned major cell type identities based on the expression of the marker genes using the RenameIdents() function, allowing for downstream analyses and visualisations on named clusters 


