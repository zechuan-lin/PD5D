library(Seurat)
library(cowplot)
library(ggplot2)
library(ggsci)
library(dplyr)
library(psych)
library(pheatmap)
library(harmony)
#library(clusterProfiler)
library(DOSE)
library(GOSemSim)
library(enrichplot)
library(stringr)
library(glmpca)
library(SeuratWrappers)
library(Matrix)

sparsedata <- as(as.matrix(read.delim("/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/merge_94_samples.txt.gz",sep = "\t",row.names = 1)),"sparseMatrix")

metadata <- read.delim("Files/preliminary_metadata_for_94_samples_postQCfiltering.tsv", sep = "\t", stringsAsFactors = FALSE)

SeuratObject <- CreateSeuratObject(counts = sparsedata,
                            project = "FullIntegration",
                            min.cells = 3)
rm(sparsedata)

SeuratObject@meta.data$sample_id <- metadata$library_id
SeuratObject@meta.data$case <- metadata$case
SeuratObject@meta.data$batch <- metadata$batch
SeuratObject@meta.data$sex <- metadata$sex
SeuratObject@meta.data$RIN <- metadata$RIN
SeuratObject@meta.data$PMI <- metadata$PMI
SeuratObject@meta.data$age <- metadata$age
age_bracket <- as.vector(cut(SeuratObject@meta.data$age, c(50,60,70,80,90,100,110)))
age_bracket <- gsub("\\(|]","",age_bracket)
age_bracket <- gsub(",","-",age_bracket)
age_bracket <- gsub("90-100","90-110",gsub("100-110","90-110",age_bracket)) 
SeuratObject@meta.data$age_bracket <- age_bracket
SeuratObject@meta.data$PMI_bracket <- droplevels(cut(as.numeric(SeuratObject@meta.data$PMI), seq(0,6,by=1)))
SeuratObject@meta.data$RIN_bracket <- droplevels(cut(as.numeric(SeuratObject@meta.data$RIN), seq(0,10,by=1)))

SeuratObject[["percent.mt"]] <- PercentageFeatureSet(SeuratObject, pattern = "^MT-")

saveRDS(SeuratObject,"/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegration_Oct2023_FinalOSR_6HC_MTG_Combined_Matrix.rds")

cellbarcodesdf <- data.frame(cbind(rownames(SeuratObject@meta.data)[200000:200200],colnames(SeuratObject@assays$RNA@counts)[200000:200200]))

colnames(cellbarcodesdf) <- c("metadata","countmatrix")

write.table(cellbarcodesdf, file = "Files/cell_barcodes_df.txt", sep = "\t", quote = FALSE)

SeuratObject <- NormalizeData(SeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)

SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(SeuratObject)

SeuratObject <- ScaleData(SeuratObject, features = all.genes, verbose = FALSE)

#finding the top 30 principal components for cells
SeuratObject <- RunGLMPCA(SeuratObject, features=SeuratObject@assays$RNA@var.features, L = 70)

SeuratObject <- RunHarmony(SeuratObject, group.by.vars = c("sample_id","batch","sex","age_bracket","RIN_bracket","PMI_bracket"), plot_convergence = TRUE, reduction = "glmpca", theta = c(0.4,0.4,0.4,0.4,0.4,0.4))

harmony_embeddings <- Embeddings(SeuratObject, 'harmony')

SeuratObject_ElbowPlot <- ElbowPlot(SeuratObject, reduction = "harmony", ndims = 70)

ggsave2("Figures/SeuratObject_ElbowPlot.pdf", SeuratObject_ElbowPlot, device = "pdf", width = 4, height = 4, units = "in")

saveRDS(SeuratObject,"/n/scratch3/users/j/jap0606/FullIntegration_OutlierSamplesRemoved/FullIntegration_Oct2023_FinalOSR_6HC_MTG_Part1.rds")

meta <- SeuratObject@meta.data

write.table(meta, file = "Files/Part1Metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)

