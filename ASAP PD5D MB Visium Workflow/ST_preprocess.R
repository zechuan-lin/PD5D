library(Seurat)
library(hdf5r)
library(SeuratWrappers)
library(data.table)
library(plyr)

merge_spatial_files <- function(input_file_list, controls_only=FALSE, age_batch_size=10) {
    st_input_files <- read.csv(input_file_list, comment.char = '#')
    data.list <- list();count <- list();j=1
    for (i in c(1:nrow(st_input_files))) {
        row_vals = st_input_files[i, ]
        if (controls_only && row_vals$diagnosis != "HC") {
            print(paste("Skipping non-control:", row_vals$sample_spatial))
            next
        }
	print(paste("Processing:", row_vals$sample_spatial))
        if (row_vals$fold_removed) {
            print(paste(i, row_vals$sample_name, "fold removed"), sep=', ')
            st_data <- readRDS(row_vals$dir)
	    names(st_data@images) <- make.names(row_vals$sample_spatial);
            st_data <- subset(st_data, subset=selected_spot == FALSE)
        } else {
            print(paste(i, row_vals$sample_name, "no fold removed"), sep=', ')
            st_data <- Load10X_Spatial(row_vals$dir,slice = make.names(row_vals$sample_spatial))
        }
        st_data@meta.data["orig.ident"] <- row_vals$sample_spatial
	st_data@meta.data["sample_id"] <- row_vals$sample_namew
        st_data@meta.data["batch"] <- row_vals$batch
        st_data@meta.data["sex"] <- row_vals$sex
        st_data@meta.data["expired_age"] <- row_vals$expired_age
        st_data@meta.data["diagnosis"] <- row_vals$diagnosis
        st_data@meta.data["PMI"] <- row_vals$PMI
        st_data@meta.data["RIN"] <- row_vals$RIN

	st_data[["percent.mt"]] <- PercentageFeatureSet(st_data, pattern = "^MT-")
	data.list[[j]] = st_data
	j=j+1
	df <- as.data.frame(AggregateExpression(st_data,slot="counts")$Spatial)
	colnames(df)<- row_vals$sample_spatial;df[,1]<- 100000000*df[,1]/sum(df[,1])
	df$gene <- rownames(df);count[[row_vals$sample_spatial]] <- df
    }
    combined <- merge(data.list[[1]], y = data.list[2:length(data.list)])
    combined$batch <- as.factor(combined$batch)
    combined$sex <- as.factor(combined$sex)
    combined@meta.data$expired_age <- droplevels(cut(combined@meta.data$expired_age, seq(0,130,by=age_batch_size)))
    combined$diagnosis <- as.factor(combined$diagnosis)
    combined@meta.data$PMI <- droplevels(cut(as.numeric(combined@meta.data$PMI), seq(0,10,by=1)))
    combined@meta.data$RIN <- droplevels(cut(as.numeric(combined@meta.data$RIN), seq(0,10,by=1)))
    write.table(combined@meta.data, file = "meta.txt", sep = "\t", quote = FALSE)
    cc <- join_all(count,by="gene",type="full");cc[is.na(cc)] <- 0; rownames(cc) <- cc[,2];cc <- cc[,-2]
    write.table(cc, file = "mean.ST_sample_counts", sep = "\t", quote = FALSE)
    return(combined)
}

SeuratObject <- merge_spatial_files("infiles.csv", controls_only=FALSE, age_batch_size=10)
saveRDS(SeuratObject, "raw_Matrix.rds")

SeuratObject <- readRDS("raw_Matrix.rds")
SeuratObject[["percent.mt"]] <- PercentageFeatureSet(SeuratObject, pattern = "^MT-")
SeuratObject <- SplitObject(SeuratObject, split.by = "ident")
data.list <- list()
for (sample in names(SeuratObject)) {
    st <- SeuratObject[[sample]]
    sub <- subset(st, subset = nFeature_Spatial < 8000 & nFeature_Spatial > 1000 & nCount_Spatial < 50000 & percent.mt < 10)
    print(paste(sample,":","Filter out", ncol(st) - ncol(sub), "samples because of the outlier QC metrics, with", ncol(sub),"samples left."))
    pdf(paste(sample, "afterQC.pdf",sep=""))
    p3 <- SpatialFeaturePlot(sub, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"))
    print(p3);dev.off()
    data.list[[sample]] = sub
}
combined <- merge(data.list[[1]], y = data.list[2:length(data.list)])
saveRDS(combined,"raw_Matrix.rds")

