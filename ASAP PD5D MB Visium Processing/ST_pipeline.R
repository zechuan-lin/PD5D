library(Seurat)
library(reticulate)
leidenalg<- import("leidenalg")
library(hdf5r)
library(SeuratWrappers)
library(harmony)
library(data.table)
library(glmpca)
library(ggplot2)
random_seed = 1234
job_name <- "20231102"
resolutions_list <- c(0.1, 0.2, 0.3, 0.4, 0.5)

# remove duplicate columns from a sparse matrix, source: https://stackoverflow.com/a/51457395
duplicated.dgCMatrix <- function (dgCMat, MARGIN, include.all.zero.vectors = TRUE) {
    MARGIN <- as.integer(MARGIN)
    J <- rep(1:ncol(dgCMat), diff(dgCMat@p))
    I <- dgCMat@i + 1
    x <- dgCMat@x
    if (MARGIN == 1L) {
        ## check duplicated rows
        names(x) <- J
        if (include.all.zero.vectors) {
            RowLst <- split(x, factor(I, levels = 1:nrow(dgCMat)))
        } else {
            RowLst <- split(x, I)  ## will do `factor(I)` internally in `split`
        }
        result <- duplicated.default(RowLst)
    } else if (MARGIN == 2L) {
        ## check duplicated columns
        names(x) <- I
        if (include.all.zero.vectors) {
            ColLst <- split(x, factor(J, levels = 1:ncol(dgCMat)))
        } else {
            ColLst <- split(x, J)  ## will do `factor(J)` internally in `split`
        }
        result <- duplicated.default(ColLst)
    } else {
        warning("invalid MARGIN; return NULL")
        result <- NULL
    }
    result
}

# remove variablefeatures if they appear too rarely in the spatial samples
filter_variablefeatures <- function(comb_seurat, gene_list, nzero_spot_thresh=20,frac_sample_thresh=0.3, bin_size=50) {
    num_samples <- length(unique(comb_seurat@meta.data$orig.ident))
    sample_thresh <- floor(num_samples * frac_sample_thresh)
    print(paste("Genes must be nonzero in >", nzero_spot_thresh,
                " spots in >", sample_thresh, " samples", sep=""))
    sample_ids <- comb_seurat@meta.data$orig.ident
    pass_genes <- c()
    ii <- 0
    
    for (gene_chunk in split(gene_list, ceiling(seq_along(gene_list)/bin_size))) {
        gene_matrix <- comb_seurat@assays$Spatial@counts[gene_chunk, ]
        gene_matrix_gt0 <- gene_matrix > 0
        for (row in 1:nrow(gene_matrix_gt0)) {
            gene <- rownames(gene_matrix_gt0)[row]
            idxs <- gene_matrix_gt0[row, ] # slow step
            nzero_spots <- sample_ids[idxs]
            nzero_spots <- table(nzero_spots)
            pass_samples <- names(nzero_spots)[nzero_spots > nzero_spot_thresh]
            ii <- ii+1
            if (length(pass_samples) > sample_thresh) {
                pass_genes <- c(pass_genes, gene)
                print(paste(length(pass_genes), "/", ii, "Passed:", gene, sep=" "))
            } else {
                print(paste(length(pass_genes), "/", ii, "Failed:", gene, sep=" "))
            }
        }
    }
    return(pass_genes)
}
cluster_pipeline <- function(comb_seurat, job_name,
                        harmony_vars = c("orig.ident", "batch", "sex", "expired_age"),
                        harmony_thetas = c(0.4, 0.4, 0.4, 0.4),
                        resolutions_list = c(0.5, 1.5),
                        random_seed = 1234,
                        remove_mito_ribo = FALSE) {

    dir.create(job_name)
    checkpoints_file <- file.path(job_name, "checkpoints.txt")
    write(paste(Sys.time(), job_name),file=checkpoints_file,append=FALSE)
    write(paste(Sys.time(), "Harmony covars:"),file=checkpoints_file,append=TRUE)
    write(harmony_vars,file=checkpoints_file,append=TRUE)
    write(harmony_thetas,file=checkpoints_file,append=TRUE)
    write(paste(Sys.time(), "random seed:", random_seed),file=checkpoints_file,append=TRUE)
    write(paste(Sys.time(), "cluster res:"),file=checkpoints_file,append=TRUE)
    write(resolutions_list,file=checkpoints_file,append=TRUE)
    set.seed(random_seed)   
    comb_seurat <- NormalizeData(comb_seurat,normalization.method = "LogNormalize",scale.factor = 10000)
    comb_seurat <- FindVariableFeatures(comb_seurat,selection.method = "vst",nfeatures = 4000)
    comb_seurat <- ScaleData(comb_seurat,features = rownames(comb_seurat))
    saveRDS(comb_seurat, file.path(job_name, "comb_seurat_scale.rds"))
    gene_list <- VariableFeatures(comb_seurat)
    pass_genes <- filter_variablefeatures(comb_seurat, gene_list,nzero_spot_thresh=20,frac_sample_thresh=0.3, bin_size=50)
    VariableFeatures(comb_seurat) <- pass_genes

    # filter out spots that are all-zero in variable genes, or duplicates
    varfeat_counts <- comb_seurat@assays$Spatial@counts[VariableFeatures(comb_seurat), ]
    dup_spots <- which(colSums(varfeat_counts) == 0 |           duplicated.dgCMatrix(varfeat_counts, 2))
    total_zero_count <- length(which(colSums(varfeat_counts) == 0))
    print(sprintf("Variable features passing filtering: %s", length(pass_genes)))
    print(sprintf("Total spots all zero in filtered variable features: %s", total_zero_count))
    print(sprintf("Total dup spots removed: %s", length(dup_spots)))
    comb_seurat <- comb_seurat[, -dup_spots]
    
    # To save memory, remove all genes which are not in VariableFeatures
    comb_seurat <- subset(comb_seurat, features = VariableFeatures(comb_seurat))
    gc(verbose = TRUE)

    set.seed(random_seed)
    print("Running GLMPCA")
    comb_seurat <- RunGLMPCA(comb_seurat,
                       features = pass_genes, L = 30,
                       minibatch = "stochastic", batch_size = 10000)

    pdf(file.path(job_name, "glmpca_out.pdf"))
    dm_pt <- DimPlot(object = comb_seurat, reduction = "glmpca",pt.size = .1, group.by = "orig.ident") + NoLegend()
    print(dm_pt)
    vn_pt <- VlnPlot(object = comb_seurat, features = "GLMPC_1",group.by = "orig.ident", pt.size = .1) + NoLegend()
    print(vn_pt)
    dev.off()
    saveRDS(comb_seurat, file.path(job_name, "comb_seurat_glmpca.rds")) 
    set.seed(random_seed)
    print("Running Harmony")
    pdf(file.path(job_name, "harmony_out.pdf"))
    comb_seurat <- readRDS(file.path(job_name, "comb_seurat_harmony_O2.rds"))
    dm_pt <- DimPlot(object = comb_seurat, reduction = "harmony",
                 pt.size = .1, group.by = "orig.ident") + NoLegend()
    print(dm_pt)
    vn_pt <- VlnPlot(object = comb_seurat, features = "harmony_1",group.by = "orig.ident", pt.size = .1) + NoLegend()
    print(vn_pt)
    dev.off()
    saveRDS(comb_seurat, file.path(job_name, "comb_seurat_harmony.rds"))
comb_seurat <- readRDS(file.path(job_name, "comb_seurat_harmony.rds"))

for (res in resolutions_list) {
	set.seed(random_seed)
        print(paste("Finding Clusters with res=", res, sep=''))
        column_name <- paste("Spatial_snn_res.", res, sep='')
        comb_seurat <- FindNeighbors(comb_seurat,reduction = "harmony", dims = 1:30)
	comb_seurat <- FindClusters(comb_seurat,resolution = res, algorithm=4, method="igraph")
	comb_seurat <- RunUMAP(comb_seurat,reduction = "harmony", dims = 1:30)
	saveRDS(comb_seurat, file.path(job_name, paste("comb_seurat_umap",res,".rds", sep='')))
        pdf(file.path(job_name,
                      paste("leiden_out_res", res, ".pdf", sep='')))
        dm_pt <- DimPlot(comb_seurat, reduction = "umap",
                     label = TRUE, pt.size = .1,
                     group.by = column_name) + ggtitle(column_name)
        print(dm_pt)
        dev.off()
}
}

merged_seurat <- "raw_Matrix.rds"
combined <- readRDS(merged_seurat)
job_name <- "20231102"
harmony_vars <- c("orig.ident", "sex", "expired_age","PMI","RIN","batch")
harmony_thetas <- c(0.4, 0.4, 0.4, 0.4,0.4,0.4)
cluster_pipeline(combined, job_name,harmony_vars = harmony_vars,harmony_thetas = harmony_thetas,resolutions_list = resolutions_list,random_seed = 1234)
