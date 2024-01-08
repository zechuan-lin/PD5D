
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
library(glmpca)
library(SeuratWrappers)
library(Azimuth)
library(sciplot)

reference <- LoadReference(path = "~/databases/Azimuth_References/human_motorcortex")

SeuratObject <- readRDS("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegration_Oct2023_FinalOSR_6HC_MTG_PostClusterQCUMAP.rds")

query <- SCTransform(
  object = SeuratObject,
  assay = "RNA",
  new.assay.name = "refAssay",
  residual.features = rownames(x = reference$map),
  reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
  method = 'glmGamPoi',
  ncells = 2000,
  n_genes = 2000,
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)

MTGanchors <- FindTransferAnchors(
  reference = reference$map,
  query = query,
  k.filter = NA,
  reference.neighbors = "refdr.annoy.neighbors",
  reference.assay = "refAssay",
  query.assay = "refAssay",
  reference.reduction = "refDR",
  normalization.method = "SCT",
  features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
  dims = 1:50,
  n.trees = 20,
  mapping.score.k = 100
)

refdata <- lapply(X = "subclass", function(x) {
  reference$map[[x, drop = TRUE]]
})
names(x = refdata) <- "subclass"
if (FALSE) {
  refdata[["impADT"]] <- GetAssayData(
    object = reference$map[['ADT']],
    slot = 'data'
  )
}

query <- TransferData(
  reference = reference$map,
  dims = 1:50,
  query = query,
  anchorset = MTGanchors,
  refdata = refdata,
  n.trees = 20,
  store.weights = TRUE
)

query <- AddMetaData(
  object = query,
  metadata = MappingScore(anchors = MTGanchors),
  col.name = "mapping.score"
)

query <- RenameIdents(query, `1` = "Oligodendrocytes", `2` = "GLU_Neurons_1", `3` = "GLU_Neurons_2", `4` = "GLU_Neurons_3", `5` = "GLU_Neurons_4", `6` = "Astrocytes", `7` = "GLU_Neurons_5", `8` = "GLU_Neurons_6",`9` = "GLU_Neurons_7", `10` = "Microglia", `11` = "GABA_Neurons_1", `12` = "GLU_Neurons_8",`13` = "GLU_Neurons_9", `14` = "OPCs", `16` = "GABA_Neurons_2", `17`= "GABA_Neurons_3", `18`="GABA_Neurons_4", `19`="Endothelial_Cells", `20`="GABA_Neurons_5", `21`="GABA_Neurons_6", `22` = "GLU_Neurons_10",`23` = "GABA_Neurons_7", `24` = "GLU_Neurons_11", `25` = "GLU_Neurons_12", `26` = "GLU_Neurons_13", `27` = "GLU_Neurons_14", `29` = "GLU_Neurons_15", `32` = "GABA_Neurons_8", `33` = "GABA_Neurons_9", `34`= "GABA_Neurons_10", `35`="GLU_Neurons_16", `36`="GLU_Neurons_17", `38`="Unknown_Cluster_38", `39` = "GABA_Neurons_11", `40`="Pericytes", `41` = "GABA_Neurons_12",`42` = "GLU_Neurons_18", `45` = "GABA_Neurons_13", `46` = "GABA_Neurons_14", `47` = "GABA_Neurons_15", `48`= "GLU_Neurons_19", `49`="GABA_Neurons_16", `50`="GLU_Neurons_20", `52`="Immune_Cells", `54`="GABA_Neurons_17", `67`="Unknown_Cluster_67")


query@meta.data$CellSubtypes <- Idents(query)

# First predicted metadata field, change to visualize other predicted metadata
id <- "subclass"[1]
predicted.id <- paste0("predicted.", id)

# DimPlot of the reference
#DimPlot(object = reference$plot, reduction = "refUMAP", group.by = id, label = TRUE) + NoLegend()

# DimPlot of the query, colored by predicted cell type

DimPlot(query, reduction = "umap", group.by = predicted.id, label = TRUE, repel = TRUE, pt.size = 0, label.size = 2.5) +               
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"),
        legend.position = "none")

UMAP_Predicted_ID <- DimPlot(query, reduction = "umap", group.by = predicted.id, label = TRUE, repel = TRUE, pt.size = 0, label.size = 2.5) +               
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"),
        legend.position = "none")

ggsave(UMAP_Predicted_ID, filename = "Figures/Azimuth_UMAP_GroupedByPredictedID_scRNA_seq.pdf", device = "pdf", width = 6, height = 6, units = "in")

UMAP_Predicted_ID_legend <- DimPlot(query, reduction = "umap", group.by = predicted.id, label = TRUE, repel = TRUE, pt.size = 0, label.size = 2.5) +
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        title = element_text(size = 12),
        legend.key.size = unit(0.4,"cm"))

ggsave(UMAP_Predicted_ID_legend, filename = "Figures/Azimuth_UMAP_GroupedByPredictedID_scRNA_seq_legend.pdf", device = "pdf", width = 7, height = 6, units = "in")

# Plot the score for the predicted cell type of the query
FeaturePlot(object = query, features = paste0(predicted.id, ".score"), reduction = "umap")

PredictedIDScore_AssignedClusters_FeaturePlot <- FeaturePlot(object = query, features = paste0(predicted.id, ".score"), reduction = "umap")

ggsave(PredictedIDScore_AssignedClusters_FeaturePlot, filename = "Figures/PredictedIDScore_AssignedClusters_FeaturePlot.pdf", device = "pdf", width = 6, height = 6, units = "in")

VlnPlot(object = query, features = paste0(predicted.id, ".score"), group.by = "CellSubtypes", pt.size = 0) + NoLegend()

PredictedIDScore_AssignedClusters_VlnPlot <- VlnPlot(object = query, features = paste0(predicted.id, ".score"), group.by = "CellSubtypes", pt.size = 0) + NoLegend()

ggsave(PredictedIDScore_AssignedClusters_VlnPlot, filename = "Figures/PredictedIDScore_AssignedClusters_VlnPlot.pdf", device = "pdf", width = 12, height = 6, units = "in")

predicted_id_score_meantable <- group_by(query@meta.data,CellSubtypes) %>% summarise(mean = mean(predicted.subclass.score), SE = se(predicted.subclass.score))

ggplot(predicted_id_score_meantable, aes(x = CellSubtypes, y = mean, fill = CellSubtypes)) + 
  geom_bar(aes(x = CellSubtypes, y = mean), stat = "identity", alpha = 1) + 
  geom_errorbar(aes(x = CellSubtypes, ymin = mean-SE, ymax = mean+SE, colour = CellSubtypes), width = 0.4, alpha = 0.9, size = 0.5) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 12, angle = 90, face = "bold",hjust=0.95,vjust=0.2), axis.ticks = element_blank(),
        strip.background = element_blank(), strip.placement = "outside", 
        strip.text.y = element_text(size = 12, angle = 180, face = "bold"),
        strip.text.y.left = element_text(angle = 0)) + NoLegend()

MeanPredictedIDScore_AssignedClusters_Barchart <- ggplot(predicted_id_score_meantable, aes(x = CellSubtypes, y = mean, fill = CellSubtypes)) + 
  geom_bar(aes(x = CellSubtypes, y = mean), stat = "identity", alpha = 1) + 
  geom_errorbar(aes(x = CellSubtypes, ymin = mean-SE, ymax = mean+SE, colour = CellSubtypes), width = 0.4, alpha = 0.9, size = 0.5) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 12, angle = 90, face = "bold",hjust=0.95,vjust=0.2), axis.ticks = element_blank(),
        strip.background = element_blank(), strip.placement = "outside", 
        strip.text.y = element_text(size = 12, angle = 180, face = "bold"),
        strip.text.y.left = element_text(angle = 0)) + NoLegend()

ggsave(MeanPredictedIDScore_AssignedClusters_Barchart, filename = "Figures/MeanPredictedIDScore_AssignedClusters_Barchart.pdf", device = "pdf", width = 12, height = 6, units = "in")

write.table(predicted_id_score_meantable, file = "Files/predicted_id_score_meantable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

# Plot the score for the predicted cell type of the query
FeaturePlot(object = query, features = "mapping.score", reduction = "umap")

MappingScore_AssignedClusters_FeaturePlot <- FeaturePlot(object = query, features = "mapping.score", reduction = "umap")

ggsave(MappingScore_AssignedClusters_FeaturePlot, filename = "Figures/MappingScore_AssignedClusters_FeaturePlot.pdf", device = "pdf", width = 6, height = 6, units = "in")

VlnPlot(object = query, features = "mapping.score", group.by = "CellSubtypes", pt.size = 0) + NoLegend()

MappingScore_AssignedClusters_VlnPlot <- VlnPlot(object = query, features = "mapping.score", group.by = "CellSubtypes", pt.size = 0) + NoLegend()

ggsave(MappingScore_AssignedClusters_VlnPlot, filename = "Figures/MappingScore_AssignedClusters_VlnPlot.pdf", device = "pdf", width = 12, height = 6, units = "in")

mapping_score_meantable <- group_by(query@meta.data,CellSubtypes) %>% summarise(mean = mean(mapping.score), SE = se(mapping.score))

ggplot(mapping_score_meantable, aes(x = CellSubtypes, y = mean, fill = CellSubtypes)) + 
  geom_bar(aes(x = CellSubtypes, y = mean), stat = "identity", alpha = 1) + 
  geom_errorbar(aes(x = CellSubtypes, ymin = mean-SE, ymax = mean+SE, colour = CellSubtypes), width = 0.4, alpha = 0.9, size = 0.5) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 12, angle = 90, face = "bold",hjust=0.95,vjust=0.2), axis.ticks = element_blank(),
        strip.background = element_blank(), strip.placement = "outside", 
        strip.text.y = element_text(size = 12, angle = 180, face = "bold"),
        strip.text.y.left = element_text(angle = 0)) + NoLegend()

MeanMappingScore_AssignedClusters_Barchart <- ggplot(mapping_score_meantable, aes(x = CellSubtypes, y = mean, fill = CellSubtypes)) + 
  geom_bar(aes(x = CellSubtypes, y = mean), stat = "identity", alpha = 1) + 
  geom_errorbar(aes(x = CellSubtypes, ymin = mean-SE, ymax = mean+SE, colour = CellSubtypes), width = 0.4, alpha = 0.9, size = 0.5) + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 12, angle = 90, face = "bold",hjust=0.95,vjust=0.2), axis.ticks = element_blank(),
        strip.background = element_blank(), strip.placement = "outside", 
        strip.text.y = element_text(size = 12, angle = 180, face = "bold"),
        strip.text.y.left = element_text(angle = 0)) + NoLegend()

ggsave(MeanMappingScore_AssignedClusters_Barchart, filename = "Figures/MeanMappingScore_AssignedClusters_Barchart.pdf", device = "pdf", width = 12, height = 6, units = "in")

write.table(mapping_score_meantable, file = "Files/mapping_score_meantable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

#barchart of predicted cell type composition for each assigned cluster

predictedcelltype_compositiontable <- group_by(query@meta.data,CellSubtypes, predicted.subclass) %>% summarise(length(predicted.subclass))

predictedcelltype_compositiontable <- ungroup(predictedcelltype_compositiontable)

colnames(predictedcelltype_compositiontable)[3] <- "predictedsubclass_size"

predictedcelltype_compositiontable <- group_by(predictedcelltype_compositiontable, CellSubtypes) %>% mutate(predictedsubclasspercentage=(predictedsubclass_size/sum(predictedsubclass_size))*100)

ggplot(predictedcelltype_compositiontable, aes(fill=predicted.subclass, y=predictedsubclasspercentage, x=CellSubtypes)) + 
  geom_bar(position="stack", stat="identity") +
  ylab("Proportion") +
  xlab("Cell Type") +
  scale_fill_discrete(name = "predicted.subclass") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6),
        axis.text.y = element_text(size = 6),
        legend.text=element_text(size=6)) +
  scale_y_continuous(breaks=c(0,25,50,75,100))

PredictedCellTypeComposition_Barchart <- ggplot(predictedcelltype_compositiontable, aes(fill=predicted.subclass, y=predictedsubclasspercentage, x=CellSubtypes)) + 
  geom_bar(position="stack", stat="identity") +
  ylab("Proportion") +
  xlab("Cell Type") +
  scale_fill_discrete(name = "predicted.subclass") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.text=element_text(size=12)) +
  scale_y_continuous(breaks=c(0,25,50,75,100))

ggsave(PredictedCellTypeComposition_Barchart, filename = "Figures/PredictedCellTypeComposition_Barchart.pdf", device = "pdf", width = 14, height = 8, units = "in")

write.table(predictedcelltype_compositiontable, file = "Files/predictedcelltype_compositiontable.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

write.table(query@meta.data, file = "Files/Azimuth_Metadata.tsv", quote = FALSE, row.names = FALSE, sep = "\t")


