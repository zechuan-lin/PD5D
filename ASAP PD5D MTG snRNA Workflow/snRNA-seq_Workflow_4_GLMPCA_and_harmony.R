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
library(glmpca)
library(SeuratWrappers)
library(Matrix)
library(optparse)

option_list = list(
  make_option(c("-f", "--matrix"), type="character", default=NULL, 
              help="Path to filtered matrix, output from script 2", metavar="character"),
    make_option(c("-m", "--metadata"), type="character", default="NULL", 
              help="metadata where each row corresponds to column of --matrix, in the same order", metavar="character"),
    make_option(c("-t", "--theta"), type="character", default=0.4,
              help="theta value per covariate for harmony, ideal total should not exceed 2-2.5", metavar="character"),
    make_option(c("-c", "--ncovariates"), type="integer", default=6,
              help="number of harmony covariates", metavar="integer"),
    make_option(c("-l", "--ldglm"), type="integer", default=NULL,
              help="latent dimensions for GLM-PCA, should be large enough to be able to see the ideal cut-off on the elbow plot output for the size of dataset", metavar="integer"),
    make_option(c("-o", "--outfile"), type="character", default="NULL",
              help="name of output SeuratObject", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Assigning variables - read through script before assigning to ensure correct value assignment

matrixpath <- opt$matrix
metapath <- opt$metadata
thetaval <- as.numeric(opt$theta)
covn <- opt$ncovariates
ldglm <- opt$ldglm
outfile <- opt$outfile

#Making required folders

dir.create("Files")
dir.create("Figures")

#Reading in pre-filtered counts matrix
sparsedata <- as(as.matrix(read.delim(matrixpath,sep = "\t",row.names = 1)),"sparseMatrix")

#Reading in pre-aggregated metadata (columns: sample id, case, batch, sex, RIN, PMI, age)
metadata <- read.delim(metapath, sep = "\t", stringsAsFactors = FALSE)

#Creating Seurat object
SeuratObject <- CreateSeuratObject(counts = sparsedata,
                            project = "FullIntegration")

rm(sparsedata)

#Adding metadata to Seurat Object, generating rational bins for continous variables as Harmony prefers discreet, this will vary depending upon the metadata.
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

#Counts per ten thousand transformation followed by log1p() transformation
SeuratObject <- NormalizeData(SeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)

#Finding top 2000 most varables features
SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(SeuratObject)

#Scaling, although GLMPCA doesn't use scaled data, the subsequent functions require scaled data to be present
SeuratObject <- ScaleData(SeuratObject, features = all.genes, verbose = FALSE)

#Finding the top 30 principal components for cells, GLM-PCA.
SeuratObject <- RunGLMPCA(SeuratObject, features=SeuratObject@assays$RNA@var.features, L = ldglm)

#Run harmony to correct PCA embeddings for our clincial covariates. Total theta should be between 2-2.5. Increase/decrease theta value for less/more covariates accordingly

SeuratObject <- RunHarmony(SeuratObject, group.by.vars = c("sample_id","batch","sex","age_bracket","RIN_bracket","PMI_bracket"), plot_convergence = TRUE, reduction.use = "glmpca", theta = rep(thetaval,covn))

harmony_embeddings <- Embeddings(SeuratObject, 'harmony')

#Making elbow plot to determine the number of PCs used to cluster, ndims == L from runGLMPCA
SeuratObject_ElbowPlot <- ElbowPlot(SeuratObject, reduction = "harmony", ndims = ldglm)

#Saving Seurat Object, Elbow plot figure and inital metadata. Interpret Elbow plot to determine PCs to use when clustering, see https://satijalab.org/seurat/articles/pbmc3k_tutorial. Note: Elbow plots for large datasets can have "weak" elbows, ending up as more of a curve. Some optimisation may be required
ggsave2("Figures/SeuratObject_ElbowPlot.pdf", SeuratObject_ElbowPlot, device = "pdf", width = 4, height = 4, units = "in")

saveRDS(SeuratObject,outfile)

meta <- SeuratObject@meta.data

write.table(meta, file = "Files/InitialMetadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

