
##### 2023-04-21 ######

##########################################################
#####            batch 1-7                      ##########
#### 20 samples in total; output from cellranger-arc #####
##########################################################


##set up

library(ArchR)
set.seed(1)

##########################
###Step 0: SET UP#########
##########################

##set the default number of threads for parallelized operations

addArchRThreads(threads = 10)


###input data


inDir <- "/data/neurogen/ASAP/Multiomics/cellranger_multiome/"

fragFiles <- list.files(inDir,pattern="atac_fragments.tsv.gz$",full.names=TRUE,recursive=TRUE)

## get rid of cellranger_arc_aggr output for each batch and _combine output for batch3

fragFiles <- fragFiles[grep("_aggr|_combine",fragFiles, invert=TRUE)]

## get sample ID

##write a function to extract the 3rd to last element

last3rd <- function(invect){
	last3 <- invect[length(invect)-2]
	return(last3)
}

names(fragFiles) <- sapply(strsplit(fragFiles,"/"),last3rd)

####create Arrow Files

#library(here)

##here() starts at /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/archr

setwd("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/01_archr")

dir.create(here("2023-04_batch1-7"))


## add a reference genome annotation
addArchRGenome("hg38")

## setwd so all files will be under 2023-04_batch1-7

setwd(here("2023-04_batch1-7"))

#####################################
####Step 1: create arrow files#######
#####################################

ArrowFiles <- createArrowFiles(
  inputFiles = fragFiles,
  sampleNames = names(fragFiles),
  QCDir = "QualityControl",
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

################################
#####Step 2:Inferring Doublets##
################################

## Skip this as 10/20 with UMAP projection R square < 0.9
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)


############################################
####Step 3: Creating an ArchRProject #####
############################################


###Step 3.1: Creating an ArchRProject
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "First20samples",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)


## add batch information

Batch <- rep(1,length(proj$Sample))

Batch[which(proj$Sample %in% c("BN1814SN","BN1939SN","BN1959SN"))] <- 2
Batch[which(proj$Sample %in% c("BN0415SN","BN0464SN","BN1719SN"))] <- 3
Batch[which(proj$Sample %in% c("BN1750SN","BN1805SN","BN1812SN"))] <- 4
Batch[which(proj$Sample %in% c("BN1822SN","BN1827SN","BN1848SN"))] <- 5
Batch[which(proj$Sample %in% c("BN1849SN","BN1862SN","BN1872SN"))] <- 6
Batch[which(proj$Sample %in% c("BN1902SN","BN1957SN","BN2003SN"))] <- 7

proj$Batch <- as.factor(Batch)


##save our original proj using saveArchRProject from ArchR.

saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj", load = FALSE)

#proj <- loadArchRProject(path = "Save-Proj")

##Skip this step due to low R^2 for some samples in inferring doublets
##filter putative doublets 
#proj <- filterDoublets(ArchRProj = proj)

###Step 3.2: Plotting QC metrics - log10(Unique Fragments) vs TSS enrichment score

df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))


p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10(Number of Unique Fragments)",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = proj, addDOC = FALSE)


###Step 3.3: Plotting Sample Statistics

##Make a ridge plot for each sample for the TSS enrichment scores

p1 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )

###Make a violin plot for each sample for the TSS enrichment scores
p2 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

###Make a ridge plot for each sample for the log10(unique nuclear fragments).

p3 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )

###Make a violin plot for each sample for the log10(unique nuclear fragments).
p4 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj, addDOC = FALSE, width = 4, height = 4)

####Step 3.4: Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.

p1 <- plotFragmentSizes(ArchRProj = proj)

p2 <- plotTSSEnrichment(ArchRProj = proj)

plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


###Step 4: Dimensionality Reduction and Clustering
##implements an iterative LSI dimensionality reduction
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

#### Batch Effect Correction

proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

##Step 5.1: clustering
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")


table(proj$Clusters)

#  C1   C10   C11   C12   C13   C14   C15   C16    C2    C3    C4    C5    C6
#   45  2198  9011  6224  6309 10400  6570  6301   770  1016  4282   532   811
#   C7    C8    C9
#  909  3859  2087

###see which sample resides in which clusters

cM <- confusionMatrix(paste0(proj$Clusters), paste0(proj$Sample))


###plot this confusion matrix as a heatmap


library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)

plotPDF(p, name = "Cluster_confusion_matrix_across_sample.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

## see which batch resides in which clusters

cM2 <- confusionMatrix(paste0(proj$Clusters), paste0(proj$Batch))


###plot this confusion matrix as a heatmap


cM2 <- cM2 / Matrix::rowSums(cM2)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM2),
    color = paletteContinuous("whiteBlue"),
    border_color = "black"
)

plotPDF(p, name = "Cluster_confusion_matrix_across_batch.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

####Step 6.1: Visualizing in a 2D UMAP Embedding
##add a UMAP embedding 
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")



##visualize various attributes
##color by clusters

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

#plotPDF(p, name = "Plot-UMAP-Cell-Clusters.pdf",
#	ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

###by sample
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")


# by batch

p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Batch", embedding = "UMAP")

#plotPDF(p, name = "Plot-UMAP-Cell-Sample.pdf",
#        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

## to visualize two plots side by side
#ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2,p3, name = "Plot-UMAP-Sample-Clusters-Batch.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


###Step 6.2: t-Stocastic Neighbor Embedding (t-SNE)

proj <- addTSNE(
    ArchRProj = proj, 
    reducedDims = "IterativeLSI", 
    name = "TSNE", 
    perplexity = 30
)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")

p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "TSNE")

p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Batch", embedding = "TSNE")

#ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2,p3, name = "Plot-TSNE-Sample-Clusters-Batch.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

###Step 6.3: Dimensionality Reduction After Harmony

###UMAP
proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)


p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")

p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")

p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Batch", embedding = "UMAPHarmony")

#ggAlignPlots(p3, p4, type = "h")

plotPDF(p3,p4,p5, name = "Plot-UMAP2Harmony-Sample-Clusters-Batch.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

####t-SNE

proj <- addTSNE(
    ArchRProj = proj, 
    reducedDims = "Harmony", 
    name = "TSNEHarmony", 
    perplexity = 30
)

p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "TSNEHarmony")

p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "TSNEHarmony")

p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Batch", embedding = "TSNEHarmony")

#ggAlignPlots(p3, p4, type = "h")

plotPDF(p3,p4,p5, name = "Plot-TSNE2Harmony-Sample-Clusters-Batch.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)



####Step 7: Identifying Marker Genes
##identify marker genes based on gene scores
markersGS <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

##get a list of DataFrame objects, one for each of our clusters, containing the relevant marker features using the getMarkers() function:
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C1

###visualize all markers

## from XT

#marker <- list(c("AQP4","GFAP"),c("TH","VMAT2","DAT"),c("MOBP","PLP1","MBP"),c("VCAN","OLIG1","OLIG2"),c("CLDN5","FLT1","ENG"),c("PDGFRB"),c("ABCA9"),c("CX3CR1","P2RY12","ARHGAP15"),c("VGLUT2","VGLUT1"),c("GAD1","GAD2","VGAT"),c("MRC1","CD163"),c("CD8A","CD96"))
#names(marker) <- c("astrocytes","DA neurons","oligodendrocytes","OPCs","Endothelial cell","Pericytes","Fibroblasts","Microglia","Glutamate neurons","GABA neurons","Monocytes","Memory_CD8_T_Cells")

#features <- c("ENO2","RBFOX3","NEFM","NEFL","TH", "DAT", "VMAT2","VGLUT1","VGLUT2","VGAT","GAD1","GAD2","AQP4","GFAP","PLP1", "OLIG1", "VCAN","CX3CR1","P2RY12","FLT1","PDGFRB","ABCA9")

#The first 4 genes are pan-neuronal markers 

markerGenes  <- c(
	"ENO2","RBFOX3","NEFM","NEFL", #pan-neuronal markers
	"AQP4","GFAP", #astrocytes
	"TH","VMAT2","DAT", #DA neurons
	"MOBP","PLP1","MBP", #oligodendrocytes
	"VCAN","OLIG1","OLIG2",#OPCs
	"CLDN5","FLT1","ENG", # Endothelial cell
	"PDGFRB", #Pericytes
	"ABCA9", #Fibroblasts
	"CX3CR1","P2RY12","ARHGAP15", #Microglia
	"VGLUT2","VGLUT1", #Glutamate neurons
	"GAD1","GAD2","VGAT", #GABA neurons
	"MRC1","CD163", #Monocytes
	"CD8A","CD96" #Memory_CD8_T_Cells
  )

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

## from JP

markerGenesJP <- c(
	"ENO2","RBFOX3",# – General Neurons Markers
	"TH", "SLC6A3", "SLC18A2", # – Dopaminergic Neurons Markers
	"SLC17A6", "SLC17A7", # – GLU Neuron Markers
	"SLC32A1", "GAD1", "GAD2", # – GABA Neuron Markers
	"AQP4", "GFAP", # – Astrocyte Markers
	"PLP1", "MBP", "OLIG1", # – Oligodendrocyte Markers
	"VCAN", "BCAN", # – OPC Markers
	"CX3CR1", "P2RY12", # – Microglia Markers
	"FLT1", "CLDN5", # – Endothelial Cell Markers
	"CPED1" # – Pericytes
)



##plot heatmap


plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

##print out marker gene for each cluster
dir.create("First20samples/MarkerGene")

outList <- list()
for (j in 1:length(markerList)){
	cluster <- paste0("C",j)
	outList[[cluster]] <- markerList[[j]]
}

library(openxlsx)

write.xlsx(outList, file = "First20samples/MarkerGene/2023-04_Batch1-7_MarkerGenePerCluster.xlsx")


###Step 7.4: Visualizing Marker Genes on an Embedding

## subset to prevent plotEmbedding failture due to missing of some genes
## Error: FeatureName (VMAT2,DAT,VGLUT2,VGLUT1,VGAT) does not exist! See getFeatures

markerGenes1  <- c(
        "AQP4","GFAP", #astrocytes
        "TH",#"VMAT2","DAT", #DA neurons
        "MOBP","PLP1","MBP", #oligodendrocytes
        "VCAN","OLIG1","OLIG2",#OPCs
        "CLDN5","FLT1","ENG", # Endothelial cell
        "PDGFRB", #Pericytes
        "ABCA9", #Fibroblasts
        "CX3CR1","P2RY12","ARHGAP15", #Microglia
        #"VGLUT2","VGLUT1", #Glutamate neurons
        "GAD1","GAD2",#"VGAT", #GABA neurons
        "MRC1","CD163", #Monocytes
        "CD8A","CD96" #Memory_CD8_T_Cells
  )

##UMAP
p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes1, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

###TSNE

p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = markerGenes1,
    embedding = "TSNE",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)

plotPDF(plotList = p,
    name = "Plot-TSNE-Marker-Genes-WO-Imputation.pdf",
    ArchRProj = proj,
    addDOC = FALSE, width = 5, height = 5)

#####Step 7.5 Marker Genes Imputation with MAGIC

proj <- addImputeWeights(proj)

##UMAP
p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes1, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

###TSNE
p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = markerGenes1,
    embedding = "TSNE",
    imputeWeights = getImputeWeights(proj)
)

plotPDF(plotList = p,
    name = "Plot-TSNE-Marker-Genes-W-Imputation.pdf",
    ArchRProj = proj,     
    addDOC = FALSE, width = 5, height = 5)

### Step 7.6 Track Plotting with ArchRBrowser

p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000
)


plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

##save our original proj using saveArchRProject from ArchR.

saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj-01-noAnnotation", load = FALSE)

######Step 8: Defining Cluster Identity (SKIP for now)

####Step 9: Pseudo-bulk Replicates in ArchR


proj2 <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters")


saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-Proj2", load = FALSE)


# Save all your objects in My_Object.RData
save.image(file = "AllSkin_08032021.RData")

####Step 10: 10 Calling Peaks with ArchR

##need run: module load MACS2 first

###!!!!!module load python/2.7.3; module load macs2/2.1.1.20160309
##doesnot work

###conda install -c bioconda macs2
#pathToMacs2 <- findMacs2()


###############################
####All Skin: cellranger-atac
### only 1299/76924 (~1.7%) cells left
###############################


##re-load
setwd("/data/rama/labMembers/mao20/Saladi/scatac/results/ArchR/Skin")
#proj <- loadArchRProject(path = "Analysis")

# Save all your objects in My_Object.RData
#save.image(file = "AllSkin_08032021.RData")



## Loading the workspace

load("AllSkin_08032021.RData")

###Step 10.2: peak calling
pathToMacs2 <- findMacs2()

proj2 <- addReproduciblePeakSet(
    ArchRProj = proj2, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2
    #pathToMacs2="/data/rama/labMembers/mao20/miniconda3/envs/Signac/bin/macs2"
)

###error
## might due to module load old version macs2?
## try without module load, use macs2 in conda environment
### worked!!!

getPeakSet(proj2)


##Step 10.4: Add Peak Matrix

##save our original proj2 using the saveArchRProject() function. This ArchRProject contains the MACS2-derived merged peak set.
saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-Proj2", load = FALSE)

##prepare for downstream analyses, we can create a new ArchRProject

proj3 <- addPeakMatrix(proj2)

getAvailableMatrices(proj3)


#####Step 11:  Identifying Marker Peaks with ArchR

##Marker features are features that are unique to a specific cell grouping. These can be very useful in understanding cluster- or cell type-specific biology. 


markersPeaks <- getMarkerFeatures(
    ArchRProj = proj3, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markersPeaks

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList

###save markerPeaks
dir.create(paste0("./Analysis/MarkerPeak"))

outList <- list()
for (j in 1:length(markerList)){
        cluster <- paste0("C",j)
        outList[[cluster]] <- markerList[[j]]
}

sample <- "allSkin"

write.xlsx(outList, file = paste0("./Analysis/MarkerPeak/",sample,"_MarkerPeakPerCluster.xlsx"))

###Step 11.2: Plotting Marker Peaks in ArchR

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = TRUE
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj3, addDOC = FALSE)


##skip 11.2.2: Marker Peak MA and Volcano Plots
## 11.2.3: Marker Peaks in Browser Tracks

p <- plotBrowserTrack(
    ArchRProj = proj3, 
    groupBy = "Clusters", 
    geneSymbol = c("TEAD3"),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)["C2"],
    upstream = 50000,
    downstream = 50000
)

plotPDF(p, name = "Plot-Tracks-With-Features-TEAD3", width = 5, height = 5, ArchRProj = proj3, addDOC = FALSE)


##skip 11.2.4: Pairwise Testing Between Groups

#####Step 12: Motif and Feature Enrichment with ArchR

devtools::install_github("GreenleafLab/chromVARmotifs")
library(chromVARmotifs)

proj3 <- addMotifAnnotations(ArchRProj = proj3, motifSet = "cisbp", name = "Motif")


###Step 12.2: Motif Enrichment in Marker Peaks

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj3,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.01 & Log2FC >= 1"
  )


enrichMotifs


heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj3, addDOC = FALSE)

## Step 12.3: ArchR Enrichment
##Step 12.3.1 Encode TF Binding Sites

proj3 <- addArchRAnnotations(ArchRProj = proj3, collection = "EncodeTFBS")


enrichEncode <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj3,
    peakAnnotation = "EncodeTFBS",
    cutOff = "FDR <= 0.01 & Log2FC >= 1"
  )


enrichEncode

heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)

plotPDF(heatmapEncode, name = "EncodeTFBS-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj3, addDOC = FALSE)


###Step 13: ChromVAR Deviatons Enrichment with ArchR

##Step 13.1: 13.1 Motif Deviations

##sample background peaks

proj3 <- addBgdPeaks(proj3)


proj3 <- addDeviationsMatrix(
  ArchRProj = proj3, 
  peakAnnotation = "Motif",
  force = TRUE
)


#saveArchRProject(ArchRProj = proj3, outputDirectory = "Save-Proj3", load = FALSE)

plotVarDev <- getVarDeviations(proj3, name = "MotifMatrix", plot = TRUE)

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj3, addDOC = FALSE)


###extract a subset of motifs for downstream analysis



###Step 14: Footprinting with ArchR

##Motif Footprinting

motifPositions <- getPositions(proj3)

motifPositions

motifs <- c("JUNB", "JUN", "JUND", "FOS", "FOSB", "TEAD3")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs

##proj3 <- addGroupCoverages(ArchRProj = proj3, groupBy = "Clusters")


seFoot <- getFootprints(
  ArchRProj = proj3, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters"
)

###Step 14.2 Normalization of Footprints for Tn5 Bias

##14.2.1 Subtracting the Tn5 Bias

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj3, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

##14.2.2 Dividing by the Tn5 Bias

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj3, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)


##14.2.3 Footprinting Without Normalization for Tn5 Bias

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj3, 
  normMethod = "None",
  plotName = "Footprints-No-Normalization",
  addDOC = FALSE,
  smoothWindow = 5
)


##### 2023-07-09 ######

##########################################################
#####            all batches                      ##########
#### 101 samples in total; output from cellranger-arc #####
##########################################################


##set up
## load R session in /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/01_archr
library(ArchR)
set.seed(1)

##########################
###Step 0: SET UP#########
##########################

##set the default number of threads for parallelized operations

addArchRThreads(threads = 10)


###input data


inDir <- "/data/neurogen/ASAP/Multiomics/cellranger_multiome/2023_summer_cellranger-arc_count"

fragFiles <- list.files(inDir,pattern="atac_fragments.tsv.gz$",full.names=TRUE,recursive=TRUE)
## exclude old version for three resequenced samples: BN1610SNR_old, BN1862SN_old, and BN1959SN_old
## exclude cellranger_arc_aggr output for each batch and _combine output for batch3

fragFiles <- fragFiles[grep("_aggr|_combine|_old",fragFiles, invert=TRUE)]

## get sample ID

##write a function to extract the 3rd to last element

last3rd <- function(invect){
	last3 <- invect[length(invect)-2]
	return(last3)
}

names(fragFiles) <- sapply(strsplit(fragFiles,"/"),last3rd)

####create Arrow Files

library(here)

#dir.create(here("2023-07_allBatch"))
dir.create(here("2023-08-v3"))

## add a reference genome annotation
addArchRGenome("hg38")

## setwd so all files will be under 2023-07_allBatch

#setwd(here("2023-07_allBatch"))

## setwd so all files will be under 2023-08-v3

setwd(here("2023-08-v3"))

#####################################
####Step 1: create arrow files#######
#####################################

ArrowFiles <- createArrowFiles(
  inputFiles = fragFiles,
  sampleNames = names(fragFiles),
  QCDir = "QualityControl",
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

### remove BN0655SN.arrow with extremely low cells (157 in cellranger-arc; 14 passed filters in ArchR)

ArrowFiles <- ArrowFiles[!ArrowFiles=="BN0655SN.arrow"]

################################
#####Step 2:Inferring Doublets##
################################

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1 
)


############################################
####Step 3: Creating an ArchRProject #####
############################################


###Step 3.1: Creating an ArchRProject
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "AllSamples",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

#numberOfCells(1): 355307
#medianTSS(1): 6.646
#medianFrags(1): 10612

## Step 3.1.1 add metadata
library(openxlsx)
metadata <- read.xlsx("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/data/metadata/2023-07-26-midbrain-scATAC-metadata.xlsx")

for(meta in names(metadata)){
        if(meta != "Sample"){
                proj@cellColData[[meta]] <- metadata[match(as.character(proj@cellColData$Sample),metadata$Sample),meta]
        }
}




##save our original proj using saveArchRProject from ArchR.

saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj1", load = FALSE)



#proj <- loadArchRProject(path = "Save-Proj1")


###Step 3.2: Plotting QC metrics - log10(Unique Fragments) vs TSS enrichment score

df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))


p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10(Number of Unique Fragments)",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = proj, addDOC = FALSE)


###Step 3.3: Plotting Sample Statistics

##Make a ridge plot for each sample for the TSS enrichment scores

p1 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )

###Make a violin plot for each sample for the TSS enrichment scores
p2 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

###Make a ridge plot for each sample for the log10(unique nuclear fragments).

p3 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )

###Make a violin plot for each sample for the log10(unique nuclear fragments).
p4 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj, addDOC = FALSE, width = 4, height = 4)

####Step 3.4: Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.

p1 <- plotFragmentSizes(ArchRProj = proj)

p2 <- plotTSSEnrichment(ArchRProj = proj)

plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

## 3.5 filtering doublets

proj <- filterDoublets(proj)

### Filtering 16754 cells from ArchRProject!

#numberOfCells(1): 338553
#medianTSS(1): 6.652
#medianFrags(1): 10382

## save proj with doublets removed using saveArchRProject from ArchR.

saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj1-Doublets-filtered", load = FALSE)



### Step 3.6 filtering low quality samples

#samples2rm <- c("BN0662SN","BN0415SN","BN1730SN","BN1865SN","BN1412SN","BN0855SN","BN1737SN","BN0952SN","BN0655SN","BN9950SN","BN0644SN","BN0651SN","BN0615SN")
samples2rm <- c("BN0662SN","BN0415SN","BN1865SN","BN1412SN","BN0855SN","BN1737SN","BN0952SN","BN0655SN","BN9950SN","BN0644SN","BN0651SN","BN0615SN")
idxSample <- BiocGenerics::which(!proj$Sample %in% samples2rm)
cellsSample <- proj$cellNames[idxSample]
proj1 <- proj[cellsSample, ]

## 8 samples that are at the borderline: "BN1730SN"  "BN1128SN"  "BN0329SN"  "BN1910SN"  "BN1317SN"  "BN1351SN"  "BN1864SNR" "BN2037SN" 

## 89 samples retained

#numberOfCells(1): 331548
#medianTSS(1): 6.684
#medianFrags(1): 10528


## save proj with low quality samples removed using saveArchRProject from ArchR.

saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-Proj1-low-quality-samples-filtered", load = FALSE)


###Step 4: Dimensionality Reduction and Clustering
##implements an iterative LSI dimensionality reduction
proj1 <- addIterativeLSI(ArchRProj = proj1, useMatrix = "TileMatrix", name = "IterativeLSI")

#### Batch Effect Correction

proj1 <- addHarmony(
    ArchRProj = proj1,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "BatchLib"
)

##Step 5.1: clustering
### default resolution: 0.8
#proj1 <- addClusters(input = proj1, reducedDims = "IterativeLSI",nOutlier=20)

## resolution 0.8

proj1 <- addClusters(input = proj1, reducedDims = "IterativeLSI", names="Clusters",nOutlier=50, force=TRUE)

table(proj1$Clusters)


###see which sample resides in which clusters

cM <- confusionMatrix(paste0(proj1$Clusters), paste0(proj1$Sample))


###plot this confusion matrix as a heatmap


library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)

plotPDF(p, name = "Cluster_confusion_matrix_across_sample.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

## see which batch resides in which clusters

cM2 <- confusionMatrix(paste0(proj1$Clusters), paste0(proj1$BatchLib))


###plot this confusion matrix as a heatmap


cM2 <- cM2 / Matrix::rowSums(cM2)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM2),
    color = paletteContinuous("whiteBlue"),
    border_color = "black"
)

plotPDF(p, name = "Cluster_confusion_matrix_across_batch.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

####Step 6.1: Visualizing in a 2D UMAP Embedding
##add a UMAP embedding 
proj1 <- addUMAP(ArchRProj = proj1, reducedDims = "IterativeLSI")

##visualize various attributes
##color by clusters

p1 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

###by sample
p2 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Sample", embedding = "UMAP")


# by batch

p3 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "BatchLib", embedding = "UMAP")

plotPDF(p1,p2,p3, name = "Plot-UMAP-Sample-Clusters-Batch.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)


####Step 6.2: t-Stocastic Neighbor Embedding (t-SNE)
#
#proj1 <- addTSNE(
#    ArchRProj = proj1, 
#    reducedDims = "IterativeLSI", 
#    name = "TSNE", 
#    perplexity = 30
#)
#
#p1 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
#
#p2 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
#
#p3 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Batch", embedding = "TSNE")
#
##ggAlignPlots(p1, p2, type = "h")
#
#plotPDF(p1,p2,p3, name = "Plot-TSNE-Sample-Clusters-Batch.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

###Step 6.3: Dimensionality Reduction After Harmony

###UMAP
proj1 <- addUMAP(
    ArchRProj = proj1, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)


p3 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")

p4 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")

p5 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "BatchLib", embedding = "UMAPHarmony")

#ggAlignPlots(p3, p4, type = "h")

plotPDF(p3,p4,p5, name = "Plot-UMAP2Harmony-Sample-Clusters-Batch.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

####t-SNE

#proj <- addTSNE(
#    ArchRProj = proj, 
#    reducedDims = "Harmony", 
#    name = "TSNEHarmony", 
#    perplexity = 30
#)
#
#p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "TSNEHarmony")
#
#p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "TSNEHarmony")
#
#p5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Batch", embedding = "TSNEHarmony")
#
##ggAlignPlots(p3, p4, type = "h")
#
#plotPDF(p3,p4,p5, name = "Plot-TSNE2Harmony-Sample-Clusters-Batch.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)



####Step 7: Identifying Marker Genes
##identify marker genes based on gene scores
markersGS <- getMarkerFeatures(
    ArchRProj = proj1,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

##get a list of DataFrame objects, one for each of our clusters, containing the relevant marker features using the getMarkers() function:
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C1

###visualize all markers

## from XT

#marker <- list(c("AQP4","GFAP"),c("TH","VMAT2","DAT"),c("MOBP","PLP1","MBP"),c("VCAN","OLIG1","OLIG2"),c("CLDN5","FLT1","ENG"),c("PDGFRB"),c("ABCA9"),c("CX3CR1","P2RY12","ARHGAP15"),c("VGLUT2","VGLUT1"),c("GAD1","GAD2","VGAT"),c("MRC1","CD163"),c("CD8A","CD96"))
#names(marker) <- c("astrocytes","DA neurons","oligodendrocytes","OPCs","Endothelial cell","Pericytes","Fibroblasts","Microglia","Glutamate neurons","GABA neurons","Monocytes","Memory_CD8_T_Cells")

#features <- c("ENO2","RBFOX3","NEFM","NEFL","TH", "DAT", "VMAT2","VGLUT1","VGLUT2","VGAT","GAD1","GAD2","AQP4","GFAP","PLP1", "OLIG1", "VCAN","CX3CR1","P2RY12","FLT1","PDGFRB","ABCA9")

#The first 4 genes are pan-neuronal markers 

### replace (VMAT2,DAT,VGLUT2,VGLUT1,VGAT)  with SLCxxxx from JP's list

markerGenes  <- c(
	"ENO2","RBFOX3","NEFM","NEFL", #pan-neuronal markers
	"AQP4","GFAP", #astrocytes
	#"TH","VMAT2","DAT", #DA neurons
	"TH", "SLC6A3", "SLC18A2",#DA neurons
	"MOBP","PLP1","MBP", #oligodendrocytes
	"VCAN","OLIG1","OLIG2",#OPCs
	"CLDN5","FLT1","ENG", # Endothelial cell
	"PDGFRB", #Pericytes
	"ABCA9", #Fibroblasts
	"CX3CR1","P2RY12","ARHGAP15", #Microglia
	#"VGLUT2","VGLUT1", #Glutamate neurons
	"SLC17A6", "SLC17A7",#Glutamate neurons
	#"GAD1","GAD2","VGAT", #GABA neurons
	"SLC32A1", "GAD1", "GAD2", #GABA neurons
	"MRC1","CD163", #Monocytes
	"CD8A","CD96" #Memory_CD8_T_Cells
  )

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

### from JP
#
#markerGenesJP <- c(
#	"ENO2","RBFOX3",# – General Neurons Markers
#	"TH", "SLC6A3", "SLC18A2", # – Dopaminergic Neurons Markers
#	"SLC17A6", "SLC17A7", # – GLU Neuron Markers
#	"SLC32A1", "GAD1", "GAD2", # – GABA Neuron Markers
#	"AQP4", "GFAP", # – Astrocyte Markers
#	"PLP1", "MBP", "OLIG1", # – Oligodendrocyte Markers
#	"VCAN", "BCAN", # – OPC Markers
#	"CX3CR1", "P2RY12", # – Microglia Markers
#	"FLT1", "CLDN5", # – Endothelial Cell Markers
#	"CPED1" # – Pericytes
#)
#


##plot heatmap


plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj1, addDOC = FALSE)

##print out marker gene for each cluster
dir.create("MarkerGene")

outList <- list()
for (j in 1:length(markerList)){
	cluster <- paste0("C",j)
	outList[[cluster]] <- markerList[[j]]
}

library(openxlsx)

write.xlsx(outList, file = "MarkerGene/2023-08_MarkerGenePerCluster.xlsx")


###Step 7.4: Visualizing Marker Genes on an Embedding

##UMAP
p <- plotEmbedding(
    ArchRProj = proj1, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
    ArchRProj = proj1, 
    addDOC = FALSE, width = 5, height = 5)

###TSNE

#p <- plotEmbedding(
#    ArchRProj = proj1,
#    colorBy = "GeneScoreMatrix",
#    name = markerGenes1,
#    embedding = "TSNE",
#    quantCut = c(0.01, 0.95),
#    imputeWeights = NULL
#)
#
#plotPDF(plotList = p,
#    name = "Plot-TSNE-Marker-Genes-WO-Imputation.pdf",
#    ArchRProj = proj1,
#    addDOC = FALSE, width = 5, height = 5)
#
#####Step 7.5 Marker Genes Imputation with MAGIC

proj1 <- addImputeWeights(proj1)

##UMAP
p <- plotEmbedding(
    ArchRProj = proj1, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj1)
)

plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
    ArchRProj = proj1, 
    addDOC = FALSE, width = 5, height = 5)


##UMAPHarmony
p <- plotEmbedding(
    ArchRProj = proj1,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAPHarmony",
    imputeWeights = getImputeWeights(proj1)
)

plotPDF(plotList = p,
    name = "Plot-UMAPHarmony-Marker-Genes-W-Imputation.pdf",
    ArchRProj = proj1,
    addDOC = FALSE, width = 5, height = 5)



####TSNE
#p <- plotEmbedding(
#    ArchRProj = proj1,
#    colorBy = "GeneScoreMatrix",
#    name = markerGenes1,
#    embedding = "TSNE",
#    imputeWeights = getImputeWeights(proj1)
#)
#
#plotPDF(plotList = p,
#    name = "Plot-TSNE-Marker-Genes-W-Imputation.pdf",
#    ArchRProj = proj1,     
#    addDOC = FALSE, width = 5, height = 5)
#
### Step 7.6 Track Plotting with ArchRBrowser

p <- plotBrowserTrack(
    ArchRProj = proj1, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000
)


plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes.pdf", 
    ArchRProj = proj1, 
    addDOC = FALSE, width = 5, height = 5)

##save our original proj1 using saveArchRProject from ArchR.

saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-Proj-01-woAnnotation", load = FALSE)

######Step 8: Defining Cluster Identity (label transfer from scRNA annotation)


### Step 8.1: annotation I - clustersManual: marker gene based manual annotation


labelNew <- c("Oligodendrocytes","Astrocytes","Astrocytes","Astrocytes","DA_neurons","GABA_neurons","Oligodendrocytes","Oligodendrocytes","Oligodendrocytes","Oligodendrocytes","Oligodendrocytes","Oligodendrocytes","Oligodendrocytes","Oligodendrocytes","Oligodendrocytes","Oligodendrocytes","Oligodendrocytes","Oligodendrocytes","OPCs","Endothelial_cells","Pericytes","Fibroblasts","Microglia")

names(labelNew) <- c("C1","C2","C17","C18","C22","C23",paste0("C",3:14),"C19","C16","C21","C20","C15")

## assign each cell identity based on their majority of votes in the cluster

proj1$ClustersManual <- mapLabels(proj1$Clusters, newLabels = labelNew, oldLabels = names(labelNew))


p <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "ClustersManual", embedding = "UMAPHarmony")
plotPDF(p, name = "Plot-UMAP2Harmony-ClustersManual-MarkerGene.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

#### Step 8.2: annotation II - copy annotation for overlapping cells and assign others based on majority of votes in the cluster
# from XT

scRNAanno <- read.table("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/data/annotation/2023_08_Midbrain_snRNA_annotation_xt.txt",header=TRUE, sep="\t")

##reformat "BN1730SN_BN1730SN_AAACCGCGTCCACAAA-1" in column cell to "BN1730SN#AAACCGCGTCCACAAA-1" in ArchR project

subset23 <- function(x){
	sub("_","#",substring(x,unlist(gregexpr('_',x))[1]+1,nchar(x)))
}

scRNAanno$Sample <- sapply(scRNAanno$cell,subset23)


### check overlapping cells between scRNA and scATAC

# Load library
library(VennDiagram)

# Generate 2 sets of cells
cell_RNA <- scRNAanno$Sample
cell_ATAC <- rownames(proj1@cellColData)

# Chart
venn.diagram(
  x = list(cell_RNA, cell_ATAC),
  category.names = c("Ncell(RNA)" , "Ncell(ATAC)"),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("blue", "red"),

  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",

  ## output
  print.mode=c("raw","percent"),
  filename = here("2023-08-v3/AllSamples/Plots","Overlapping_cells_between_RNA_and_ATAC.png"),
  output = TRUE
)



## add scRNA Cluster

proj1$Clusters_scRNA <- scRNAanno[match(rownames(proj1@cellColData),scRNAanno$Sample),"anno"]

## create a confusion matrix between scATAC clusters and scRNA-seq clusters

cM <- confusionMatrix(proj1$Clusters, proj1$Clusters_scRNA)


## for each of our scATAC-seq clusters, we identify the cell type from scRNAanno which best defines that cluster.

labelNew <- colnames(cM)[2:ncol(cM)][apply(cM[,2:ncol(cM)], 1, which.max)]

names(labelNew) <- rownames(cM)

## assign each cell identity based on their majority of votes in the cluster

proj1$Clusters2 <- mapLabels(proj1$Clusters, newLabels = labelNew, oldLabels = names(labelNew))

## plot by new clusters

p3 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Diagnosis", embedding = "UMAPHarmony")

p4 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters_scRNA", embedding = "UMAPHarmony")

p5 <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters2", embedding = "UMAPHarmony")

#ggAlignPlots(p3, p4, type = "h")

plotPDF(p3,p4,p5, name = "Plot-UMAP2Harmony-Diagnosis-NewClusters-scRNA-transfer-NA.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

##### Step 8.3: Annotation III
### refine cell/cluster annotation: for overlapping cells, copy scRNA annotation; for non-overlapping cells, perform query-reference cross-platform linkage

## create a confusion matrix between scATAC clusters and scRNA-seq clusters

cM <- confusionMatrix(proj2$Clusters, proj2$Clusters_scRNA)

###plot this confusion matrix as a heatmap

cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM),
    color = paletteContinuous("whiteBlue"),
    border_color = "black"
)

plotPDF(p, name = "Cluster_confusion_matrix_across_scRNACluster_res0.8.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

###
##> table(proj1$Clusters_scRNA,useNA="always")
#
#         Astrocytes          DA neurons DA-Glu-GABA neurons   Endothelial cells
#              20129                   5                  63                1227
#        Fibroblasts        GABA neurons    Glu-GABA neurons   Glutamate neurons
#                150                1262                 329                   7
# Memory_CD8_T_Cells           Microglia           Monocytes    Oligodendrocytes
#                272               17742                 787              163964
#               OPCs           Pericytes                <NA>
#               9509                 951              115151
#


## Step 8.3: Cross-platform linkage of scATAC-seq cells with scRNA-seq cells for non-overlapping cells

projt <- subsetCells(ArchRProj = proj1, cellNames = proj1$cellNames[which(is.na(proj1$Clusters_scRNA))])


seRNA <- readRDS("/data/neurogen/Xufei/scRNA.Midbrain/backupO2/Midbrain_sc_diet.rds")

## Unconstrained Integration

projt <- addGeneIntegrationMatrix(
    ArchRProj = projt, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "anno",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

p <- plotEmbedding(ArchRProj = projt, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAPHarmony")


plotPDF(p, name = "Plot-UMAP2Harmony-scRNA-prediction-nonoverlapping.pdf", ArchRProj = projt, addDOC = FALSE, width = 5, height = 5)

cM <- as.matrix(confusionMatrix(projt$Clusters, projt$predictedGroup_Un))

### add annotation back to proj1 for non-overlapping cells


proj1$Clusters_anno3 <- proj1$Clusters_scRNA

dft <- projt@cellColData[,c("Sample","predictedGroup_Un")]

proj1$Clusters_anno3[which(is.na(proj1$Clusters_anno3))] <- projt@cellColData[match(rownames(proj1@cellColData)[which(is.na(proj1$Clusters_anno3))],rownames(projt@cellColData)),"predictedGroup_Un"]

## plot before refinement
p <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters_anno3", embedding = "UMAPHarmony")


plotPDF(p, name = "Plot-UMAP2Harmony-Clusters_anno3-before-refinement.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)


## create a confusion matrix between scATAC clusters and scRNA-seq clusters

cM <- confusionMatrix(proj1$Clusters, proj1$Clusters_anno3)

###plot this confusion matrix as a heatmap

cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM),
    color = paletteContinuous("whiteBlue"),
    border_color = "black"
)

plotPDF(p, name = "Cluster_confusion_matrix_across_scRNACluster-anno3.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


## merge DA neurons DA-Glu-GABA neurons and Glutamate neurons due to low abundance

proj1$Clusters_annoF <- proj1$Clusters_anno3
proj1$Clusters_annoF[which(proj1$Clusters_annoF %in% c("DA neurons", "DA-Glu-GABA neurons", "Glutamate neurons"))] <- "other_neurons"

table(proj1$Clusters_annoF)

## save proj before generating pseudo replicates using saveArchRProject from ArchR.

saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-Proj1-withAnnotation", load = FALSE)


## Step 8.4: Cross-platform linkage of scATAC-seq cells with scRNA-seq cells for non-overlapping cells use downsampled balanced reference dataset

projt <- subsetCells(ArchRProj = proj1, cellNames = proj1$cellNames[which(is.na(proj1$Clusters_scRNA))])


seRNA <- readRDS("/data/neurogen/Xufei/scRNA.Midbrain/backupO2/Midbrain_sc_diet.rds")

## check compostion and subset all cell types to the largest common number of cells

table(seRNA$anno)

#         Astrocytes          DA neurons DA-Glu-GABA neurons   Endothelial cells 
#              29894                 487                 478                3333 
#        Fibroblasts        GABA neurons    Glu-GABA neurons   Glutamate neurons 
#                369                2604                3425                 332 
# Memory_CD8_T_Cells           Microglia           Monocytes    Oligodendrocytes 
#                690               25973                1308              254751 
#               OPCs           Pericytes 
#              17334                2144 

### choose 332 cells from each clusters

### downsample

Idents(object = seRNA) <- "anno"


# Downsample the number of cells per identity class
seRNA332 <- subset(x = seRNA, downsample = 332)


## Unconstrained Integration

projt <- addGeneIntegrationMatrix(
    ArchRProj = projt,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA332,
    addToArrow = FALSE,
    groupRNA = "anno",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

p <- plotEmbedding(ArchRProj = projt, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAPHarmony")


plotPDF(p, name = "Plot-UMAP2Harmony-scRNA-prediction-nonoverlapping-downsample332.pdf", ArchRProj = projt, addDOC = FALSE, width = 5, height = 5)

cM <- as.matrix(confusionMatrix(projt$Clusters, projt$predictedGroup_Un))

### add annotation back to proj1 for non-overlapping cells


proj1$Clusters_anno4 <- proj1$Clusters_scRNA

dft <- projt@cellColData[,c("Sample","predictedGroup_Un")]

proj1$Clusters_anno4[which(is.na(proj1$Clusters_anno4))] <- projt@cellColData[match(rownames(proj1@cellColData)[which(is.na(proj1$Clusters_anno4))],rownames(projt@cellColData)),"predictedGroup_Un"]

## plot before refinement
p <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters_anno4", embedding = "UMAPHarmony")


plotPDF(p, name = "Plot-UMAP2Harmony-Clusters_anno4-downsample332-before-refinement.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

## save proj before generating pseudo replicates using saveArchRProject from ArchR.

saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-Proj1-withAnnotation-downsample", load = FALSE)

### CCA-based scRNA label transfer relies on reference data cell composition, is not reliable

### need to do real multiome analysis on the data




#### Step 8.5: check annotation concordance for overlapping cells only

##proj1 <- loadArchRProject(path = "Save-Proj1-withAnnotation")
## so changes will be saved to "Save-Proj1-withAnnotation"



projo <- subsetCells(ArchRProj = proj1, cellNames = proj1$cellNames[which(!is.na(proj1$Clusters_scRNA))])


###UMAP: for subset otherwise plotEmbedding will return error and plot all cells in original project
projo <- addUMAP(
    ArchRProj = projo, 
    reducedDims = "Harmony",
    name = "UMAPHarmony",
    nNeighbors = 30,
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
)
colorATAC <- c("Astrocytes"="#D51F26","DA_neurons"="#272E6A","Endothelial_cells"="#208A42","Fibroblasts"="#89288F","GABA_neurons"="#F47D2B","Microglia"="#FEE500","Oligodendrocytes"="#8A9FD1","OPCs"="#C06CAB","Pericytes"="#E6C2DC")
p3 <- plotEmbedding(ArchRProj = projo, colorBy = "cellColData", name = "ClustersManual", embedding = "UMAPHarmony",labelAsFactors=TRUE,pal=colorATAC)

plotPDF(p3, name = "Plot-UMAP2Harmony-Cluster-annotation-for-Overlapping-scRNAvsATAC-labelTest.pdf", ArchRProj = projo, addDOC = FALSE, width = 5, height = 5)

## plot by new clusters

p3 <- plotEmbedding(ArchRProj = projo, colorBy = "cellColData", name = "ClustersManual", embedding = "UMAPHarmony")

p4 <- plotEmbedding(ArchRProj = projo, colorBy = "cellColData", name = "Clusters_scRNA", embedding = "UMAPHarmony")

plotPDF(p3,p4, name = "Plot-UMAP2Harmony-Cluster-annotation-for-Overlapping-scRNAvsATAC.pdf", ArchRProj = projo, addDOC = FALSE, width = 5, height = 5)

### specify same color for same cell types

colorATAC <- c("Astrocytes"="#D51F26","DA_neurons"="#272E6A","Endothelial_cells"="#208A42","Fibroblasts"="#89288F","GABA_neurons"="#F47D2B","Microglia"="#FEE500","Oligodendrocytes"="#8A9FD1","OPCs"="#C06CAB","Pericytes"="#E6C2DC")
p3 <- plotEmbedding(ArchRProj = projo, colorBy = "cellColData", name = "ClustersManual", embedding = "UMAPHarmony",pal=colorATAC, colorTitle=names(colorATAC))

colorRNA <- c("Astrocytes"="#D51F26","DA neurons"="#272E6A","Endothelial cells"="#208A42","Fibroblasts"="#89288F","GABA neurons"="#F47D2B","Microglia"="#FEE500","Oligodendrocytes"="#8A9FD1","OPCs"="#C06CAB","Pericytes"="#E6C2DC","DA-Glu-GABA neurons"="#90D5E4","Glu-GABA neurons"="#89C75F","Glutamate neurons"="#F37B7D","Memory_CD8_T_Cells"="#9983BD","Monocytes"="#D24B27")
p4 <- plotEmbedding(ArchRProj = projo, colorBy = "cellColData", name = "Clusters_scRNA", embedding = "UMAPHarmony",pal=colorRNA,colorTitle=names(colorRNA))

##convert Clusters_scRNA to level and the first nine corresponding to those in ClustersManual
#projo$Clusters_scRNA_level <- as.factor(projo$Clusters_scRNA)
#projo$Clusters_scRNA_level <- factor(projo$Clusters_scRNA_level,levels=names(colorRNA))

#p4 <- plotEmbedding(ArchRProj = projo, colorBy = "cellColData", name = "Clusters_scRNA_level", embedding = "UMAPHarmony",pal=colorRNA)
## does not work

plotPDF(p3,p4, name = "Plot-UMAP2Harmony-Cluster-annotation-for-Overlapping-scRNAvsATAC-consistentColor.pdf", ArchRProj = projo, addDOC = FALSE, width = 5, height = 5)


## create a confusion matrix between scATAC manual clusters and scRNA-seq clusters for overlapping cells

cM <- confusionMatrix(projo$ClustersManual, projo$Clusters_scRNA)

###plot this confusion matrix as a heatmap

cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM),
    color = paletteContinuous("whiteBlue"),
    border_color = "black"
)

plotPDF(p, name = "Cluster_confusion_matrix_Overlapping-scRNAvsATAC.pdf", ArchRProj = projo, addDOC = FALSE, width = 5, height = 5)

### compare clustering agreement with Adjusted Rand index

library(pdfCluster)

adj.rand.index(projo$ClustersManual, projo$Clusters_scRNA)

## [1] 0.9739967

### It is high but largely driven by dominant cell types such as Astrocytes,Oligodendrocytes,OPCs,Microglia
### exclude cells that with identical annotation of these major cell types

C2 <- data.frame(ClustersManual=projo$ClustersManual, Clusters_scRNA=projo$Clusters_scRNA)

C3 <- C2[!((C2$ClustersManual==C2$Clusters_scRNA) & (C2$ClustersManual %in% c("Astrocytes","Oligodendrocytes","OPCs","Microglia"))),]

adj.rand.index(C3$ClustersManual, C3$Clusters_scRNA)

#[1] 0.4948691



#### Step 8.6: subset ATAC-unique cells and redo all steps to see if new clustering helps

proj1 <- loadArchRProject(path="Save-Proj1-withAnnotation")

projo <- subsetCells(ArchRProj = proj1, cellNames = proj1$cellNames[which(is.na(proj1$Clusters_scRNA))])

projo <- addIterativeLSI(ArchRProj = projo, useMatrix = "TileMatrix", name = "IterativeLSI2")

#### Batch Effect Correction

projo <- addHarmony(
    ArchRProj = projo,
    reducedDims = "IterativeLSI2",
    name = "Harmony2",
    groupBy = "BatchLib"
)

## clustering
projo <- addClusters(input = projo, reducedDims = "IterativeLSI2", name="ClustersA")


table(projo$ClustersA)

#   C1   C10   C11   C12   C13   C14   C15   C16   C17   C18   C19    C2    C3 
# 5816   168    57  1254  1960  8583   921  2102  2448  1650  1054  3113 18359 
#   C4    C5    C6    C7    C8    C9 
#11664  8510 18210  1003 25802  2477 

###see which sample resides in which clusters


cM <- confusionMatrix(paste0(projo$ClustersA), paste0(projo$Sample))


###plot this confusion matrix as a heatmap


library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM),
    color = paletteContinuous("whiteBlue"),
    border_color = "black"
)

plotPDF(p, name = "ATAC-unique-cell-Cluster_confusion_matrix_across_sample.pdf", ArchRProj = projo, addDOC = FALSE, width = 5, height = 5)

## see which batch resides in which clusters

cM2 <- confusionMatrix(paste0(projo$ClustersA), paste0(projo$BatchLib))


###plot this confusion matrix as a heatmap


cM2 <- cM2 / Matrix::rowSums(cM2)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM2),
    color = paletteContinuous("whiteBlue"),
    border_color = "black"
)

plotPDF(p, name = "ATAC-unique-cell-Cluster_confusion_matrix_across_batch.pdf", ArchRProj = projo, addDOC = FALSE, width = 5, height = 5)

####Visualizing in a 2D UMAP Embedding
##add a UMAP embedding 
projo <- addUMAP(ArchRProj = projo, reducedDims = "IterativeLSI2", name="UMAPA")



p1 <- plotEmbedding(ArchRProj = projo, colorBy = "cellColData", name = "ClustersA", embedding = "UMAPA")

###by sample
p2 <- plotEmbedding(ArchRProj = projo, colorBy = "cellColData", name = "Sample", embedding = "UMAPA")


# by batch

p3 <- plotEmbedding(ArchRProj = projo, colorBy = "cellColData", name = "BatchLib", embedding = "UMAPA")

plotPDF(p1,p2,p3, name = "ATAC-unique-cell-Plot-UMAP-Sample-Clusters-Batch.pdf", ArchRProj = projo, addDOC = FALSE, width = 5, height = 5)



### Dimensionality Reduction After Harmony

###UMAP
projo <- addUMAP(
    ArchRProj = projo,
    reducedDims = "Harmony2",
    name = "UMAPHarmony2",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
)


p3 <- plotEmbedding(ArchRProj = projo, colorBy = "cellColData", name = "ClustersA", embedding = "UMAPHarmony2")

p4 <- plotEmbedding(ArchRProj = projo, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony2")

p5 <- plotEmbedding(ArchRProj = projo, colorBy = "cellColData", name = "BatchLib", embedding = "UMAPHarmony2")

#ggAlignPlots(p3, p4, type = "h")

plotPDF(p3,p4,p5, name = "ATAC-unique-cell-Plot-UMAP2Harmony-Sample-Clusters-Batch.pdf", ArchRProj = projo, addDOC = FALSE, width = 5, height = 5)


### marker genes



markerGenes  <- c(
        "ENO2","RBFOX3","NEFM","NEFL", #pan-neuronal markers
        "AQP4","GFAP", #astrocytes
        #"TH","VMAT2","DAT", #DA neurons
        "TH", "SLC6A3", "SLC18A2",#DA neurons
        "MOBP","PLP1","MBP", #oligodendrocytes
        "VCAN","OLIG1","OLIG2",#OPCs
        "CLDN5","FLT1","ENG", # Endothelial cell
        "PDGFRB", #Pericytes
        "ABCA9", #Fibroblasts
        "CX3CR1","P2RY12","ARHGAP15", #Microglia
        #"VGLUT2","VGLUT1", #Glutamate neurons
        "SLC17A6", "SLC17A7",#Glutamate neurons
        #"GAD1","GAD2","VGAT", #GABA neurons
        "SLC32A1", "GAD1", "GAD2", #GABA neurons
        "MRC1","CD163", #Monocytes
        "CD8A","CD96" #Memory_CD8_T_Cells
  )


### Visualizing Marker Genes on an Embedding

##UMAP
p <- plotEmbedding(
    ArchRProj = projo,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAPA",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)

plotPDF(plotList = p,
    name = "ATAC-unique-cell-Plot-UMAP-Marker-Genes-WO-Imputation.pdf",
    ArchRProj = projo,
    addDOC = FALSE, width = 5, height = 5)

##### Marker Genes Imputation with MAGIC

projo <- addImputeWeights(projo,reducedDims = "IterativeLSI2")

##UMAP
p <- plotEmbedding(
    ArchRProj = projo,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAPA",
    imputeWeights = getImputeWeights(projo)
)

plotPDF(plotList = p,
    name = "ATAC-unique-cell-Plot-UMAP-Marker-Genes-W-Imputation.pdf",
    ArchRProj = projo,
    addDOC = FALSE, width = 5, height = 5)


##UMAPHarmony
p <- plotEmbedding(
    ArchRProj = projo,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAPHarmony2",
    imputeWeights = getImputeWeights(projo)
)

plotPDF(plotList = p,
    name = "ATAC-unique-cell-Plot-UMAPHarmony-Marker-Genes-W-Imputation.pdf",
    ArchRProj = projo,
    addDOC = FALSE, width = 5, height = 5)



### annotation

labelNew <- c("Microglia","Microglia",rep("Oligodendrocytes",7),rep("Doublets",2),"DA neurons","GABA neurons",rep("Astrocytes",3),"OPCs","Doublets","Endothelial cells")

names(labelNew) <- c("C1","C2",paste0("C",3:9),"C10","C11","C12","C13",paste0("C",14:16),"C17","C18","C19")


## assign each cell identity

projo$ClustersManualA <- mapLabels(projo$ClustersA, newLabels = labelNew, oldLabels = names(labelNew))


p <- plotEmbedding(ArchRProj = projo, colorBy = "cellColData", name = "ClustersManualA", embedding = "UMAPHarmony2")
plotPDF(p, name = "ATAC-unique-cell-Plot-UMAP2Harmony-ClustersManualA-MarkerGene.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)

#### subset neurons and get finer resolution

projN <- subsetCells(ArchRProj = projo, cellNames = projo$cellNames[which(projo$ClustersManualA %in% c("DA neurons","GABA neurons"))])

projN <- addIterativeLSI(ArchRProj = projN, useMatrix = "TileMatrix", name = "IterativeLSI3")

#### Batch Effect Correction

projN <- addHarmony(
    ArchRProj = projN,
    reducedDims = "IterativeLSI3",
    name = "Harmony3",
    groupBy = "BatchLib"
)

## clustering
projN <- addClusters(input = projN, reducedDims = "IterativeLSI3", name="ClustersAN",maxClusters=10)


table(projN$ClustersAN)



####Visualizing in a 2D UMAP Embedding
##add a UMAP embedding 
projN <- addUMAP(ArchRProj = projN, reducedDims = "IterativeLSI3", name="UMAPAN")



p1 <- plotEmbedding(ArchRProj = projN, colorBy = "cellColData", name = "ClustersAN", embedding = "UMAPAN")

###by sample
p2 <- plotEmbedding(ArchRProj = projN, colorBy = "cellColData", name = "Sample", embedding = "UMAPAN")


# by batch

p3 <- plotEmbedding(ArchRProj = projN, colorBy = "cellColData", name = "BatchLib", embedding = "UMAPAN")

plotPDF(p1,p2,p3, name = "ATAC-unique-neuron-Plot-UMAP-Sample-Clusters-Batch.pdf", ArchRProj = projN, addDOC = FALSE, width = 5, height = 5)



### Dimensionality Reduction After Harmony

###UMAP
projN <- addUMAP(
    ArchRProj = projN,
    reducedDims = "Harmony3",
    name = "UMAPHarmony3",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine"
)


p3 <- plotEmbedding(ArchRProj = projN, colorBy = "cellColData", name = "ClustersAN", embedding = "UMAPHarmony3")

p4 <- plotEmbedding(ArchRProj = projN, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony3")

p5 <- plotEmbedding(ArchRProj = projN, colorBy = "cellColData", name = "BatchLib", embedding = "UMAPHarmony3")

#ggAlignPlots(p3, p4, type = "h")

plotPDF(p3,p4,p5, name = "ATAC-unique-neuron-Plot-UMAP2Harmony-Sample-Clusters-Batch.pdf", ArchRProj = projN, addDOC = FALSE, width = 5, height = 5)

### Visualizing Marker Genes on an Embedding

##UMAP
p <- plotEmbedding(
    ArchRProj = projN,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAPAN",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)

plotPDF(plotList = p,
    name = "ATAC-unique-neuron-Plot-UMAP-Marker-Genes-WO-Imputation.pdf",
    ArchRProj = projN,
    addDOC = FALSE, width = 5, height = 5)

##### Marker Genes Imputation with MAGIC

projN <- addImputeWeights(projN,reducedDims = "IterativeLSI3")

##UMAP
p <- plotEmbedding(
    ArchRProj = projN,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAPAN",
    imputeWeights = getImputeWeights(projN)
)

plotPDF(plotList = p,
    name = "ATAC-unique-neuron-Plot-UMAP-Marker-Genes-W-Imputation.pdf",
    ArchRProj = projN,
    addDOC = FALSE, width = 5, height = 5)


##UMAPHarmony
p <- plotEmbedding(
    ArchRProj = projN,
    colorBy = "GeneScoreMatrix",
    name = markerGenes,
    embedding = "UMAPHarmony3",
    imputeWeights = getImputeWeights(projN)
)

plotPDF(plotList = p,
    name = "ATAC-unique-neuron-Plot-UMAPHarmony-Marker-Genes-W-Imputation.pdf",
    ArchRProj = projN,
    addDOC = FALSE, width = 5, height = 5)


### Track Plotting with ArchRBrowser

p <- plotBrowserTrack(
    ArchRProj = projN,
    groupBy = "ClustersAN",
    geneSymbol = c("TH", "SLC6A3", "SLC18A2","SLC17A6", "SLC17A7","SLC32A1", "GAD1", "GAD2"),
    upstream = 50000,
    downstream = 50000
)


plotPDF(plotList = p,
    name = "Plot-Tracks-Marker-Neuron-Genes-ATACuniqueNeuron.pdf",
    ArchRProj = projN,
    addDOC = FALSE, width = 5, height = 5)


### annotation

labelNew <- c("Doublets?","DA/GLU neurons","DA/GLU neurons","Glu-GABA neurons?","GABA neurons")

names(labelNew) <- paste0("C",1:5)


## assign each cell identity

projN$ClustersManualAN <- mapLabels(projN$ClustersAN, newLabels = labelNew, oldLabels = names(labelNew))


p <- plotEmbedding(ArchRProj = projN, colorBy = "cellColData", name = "ClustersManualAN", embedding = "UMAPHarmony3")
plotPDF(p, name = "ATAC-unique-neuron-Plot-UMAP2Harmony-ClustersManualAN-MarkerGene.pdf", ArchRProj = projN, addDOC = FALSE, width = 5, height = 5)



## save ArchR project

saveArchRProject(ArchRProj = projo, outputDirectory = "Save-Proj1-withAnnotation-ATACUniqueAll-manual", load = FALSE)

saveArchRProject(ArchRProj = projN, outputDirectory = "Save-Proj1-withAnnotation-ATACUniqueNeuron-manual", load = FALSE)

### add annotation back to proj1 for non-overlapping cells


proj1$Clusters_anno5 <- proj1$Clusters_scRNA

dfo <- projo@cellColData[,c("Sample","ClustersManualA")]

proj1$Clusters_anno5[which(is.na(proj1$Clusters_anno5))] <- projo@cellColData[match(rownames(proj1@cellColData)[which(is.na(proj1$Clusters_anno5))],rownames(projo@cellColData)),"ClustersManualA"]

## plot before refinement
p <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters_anno5", embedding = "UMAPHarmony")


plotPDF(p, name = "Plot-UMAP2Harmony-Clusters_anno5-manual-anno-nonoverlapping.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)


### majory vote for each cluster

cM <- confusionMatrix(proj1$Clusters, proj1$Clusters_anno5)


## for each of our scATAC-seq clusters, we identify the cell type from scRNAanno which best defines that cluster.

labelNew <- colnames(cM)[apply(cM, 1, which.max)]

names(labelNew) <- rownames(cM)

## assign each cell identity based on their majority of votes in the cluster

proj1$ClustersManualAMajor <- mapLabels(proj1$Clusters, newLabels = labelNew, oldLabels = names(labelNew))


p <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "ClustersManualAMajor", embedding = "UMAPHarmony")
plotPDF(p, name = "ATAC-unique-cell-Plot-UMAP2Harmony-ClustersManualA-MarkerGene-MajorVote.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)


### add annotation back to proj1 for non-overlapping neurons

proj1$Clusters_anno6 <- proj1$Clusters_scRNA

dfo <- projo@cellColData[,c("Sample","ClustersManualA")]

proj1$Clusters_anno6[which(is.na(proj1$Clusters_anno6))] <- projo@cellColData[match(rownames(proj1@cellColData)[which(is.na(proj1$Clusters_anno6))],rownames(projo@cellColData)),"ClustersManualA"]

proj1$Clusters_anno6[which(rownames(proj1@cellColData) %in% rownames(projN@cellColData))] <- projN@cellColData[match(rownames(proj1@cellColData)[which(rownames(proj1@cellColData) %in% rownames(projN@cellColData))],rownames(projN@cellColData)),"ClustersManualAN"]

## plot before refinement
p <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters_anno6", embedding = "UMAPHarmony")


plotPDF(p, name = "Plot-UMAP2Harmony-Clusters_anno6-double-manual-anno-nonoverlapping.pdf", ArchRProj = proj1, addDOC = FALSE, width = 5, height = 5)


## save ArchR project

saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-Proj1-withAnnotation-ATACUnique-doublemanualAnno", load = FALSE)


### Track Plotting with ArchRBrowser

p <- plotBrowserTrack(
    ArchRProj = proj1,
    groupBy = "ClustersManualAMajor",
    geneSymbol = markerGenes,
    upstream = 50000,
    downstream = 50000
)


plotPDF(plotList = p,
    name = "Plot-Tracks-Marker-Genes-ClustersManualAMajor.pdf",
    ArchRProj = proj1,
    addDOC = FALSE, width = 5, height = 5)



## change ranges

#proj1 <- loadArchRProject(path="Save-Proj1-withAnnotation")

### Track Plotting with ArchRBrowser

p <- plotBrowserTrack(
    ArchRProj = proj1,
    groupBy = "ClustersManual",
    geneSymbol = markerGenes,
    tileSize = 100,
    upstream = 2000,
    downstream = 2000
)


plotPDF(plotList = p,
    name = "Plot-Tracks-Marker-Genes-ClustersManual.pdf",
    ArchRProj = proj1,
    addDOC = FALSE, width = 5, height = 5)



### Step 8.7: check how many cells in scRNA retained in ATAC-seq

## load scRNA seurat project
seRNA <- readRDS("/data/neurogen/Xufei/scRNA.Midbrain/backupO2/Midbrain_sc_diet.rds")

## plot by annotation


library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)

p <- DimPlot(seRNA,raster=FALSE)

ggsave("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/01_archr/2023-08-v3/AllSamples/Plots/midbrain-scRNA-UMAP.pdf")


###add a column for the existence of cells in ATAC-seq

### get overlapping cells


### change the format of cell in ATAC ("BN1752SN#AGTCCTGAGCGATACT-1") to match those in seRNA ("BN1730SN_BN1730SN_AAACCGCGTCCACAAA-1")
## for samples in c("BN1814SN","BN1864SN","BN1430SN","BN1585SN","BN1610SN","BN1305SN"): replace SNR with SN in ATAC
## BN1308SNR and BN1960SNR are labeled same in both

formatCell <- function(x){
       sample <-  substring(x,1,unlist(gregexpr('#',x))[1]-1)
       cell <- substring(x,unlist(gregexpr('#',x))[1]+1,nchar(x))
       sample <- ifelse(sample %in%  c("BN1814SNR","BN1864SNR","BN1430SNR","BN1585SNR","BN1610SNR","BN1305SNR"),substring(sample,1,nchar(sample)-1),sample)	
       paste0(sample,"_",sample,"_",cell)      
}


ATACcell <- sapply(rownames(proj1@cellColData),formatCell)


###check how many cells are included

length(which(colnames(seRNA) %in% ATACcell))
#[1] 233496
## increased from 215553


##check if the shared cells are less in those six sample with unknown resequence time

#table(seRNA@meta.data$sample_id[which(rownames(seRNA@meta.data) %in% ATACcell)])

#BN0329SN  BN0339SN  BN0341SN  BN0347SN  BN0348SN  BN0452SN  BN0464SN  BN0536SN
#     1525      4840      2037      6167      1308       592      1435      2647
# BN0737SN  BN0746SN  BN0934SN  BN1076SN  BN1128SN  BN1144SN  BN1160SN  BN1204SN
#     3601      1063      1666      1681      2871      1881       533       764
# BN1206SN  BN1221SN  BN1266SN  BN1305SN BN1308SNR  BN1317SN  BN1339SN  BN1340SN
#     4391      2017       471      4012      2157      1761      3628      4572
# BN1351SN  BN1361SN  BN1424SN  BN1430SN  BN1504SN  BN1506SN  BN1518SN  BN1535SN
#      239      1606      1012      3784      4778      2913      4934      3678
# BN1546SN  BN1554SN  BN1567SN  BN1578SN  BN1585SN  BN1610SN  BN1614SN  BN1644SN
#     3776      5112       981      3698      4943      3625      4690       751
# BN1719SN  BN1722SN  BN1726SN  BN1730SN  BN1741SN  BN1747SN  BN1749SN  BN1750SN
#     2759      6454      4754       395      5215      2159      3694      2045
# BN1752SN  BN1756SN  BN1762SN  BN1805SN  BN1809SN  BN1812SN  BN1814SN  BN1817SN
#     2958      2569      1630      2779      3442      2151         5      6372
# BN1822SN  BN1827SN  BN1836SN  BN1839SN  BN1842SN  BN1848SN  BN1849SN  BN1855SN
#     1969      1818      2342      3576      5210      3024      2037      3534
# BN1862SN  BN1864SN  BN1867SN  BN1872SN  BN1878SN  BN1902SN  BN1910SN  BN1934SN
#     1017       730       293      1227      1518      3232       867      1770
# BN1939SN  BN1957SN  BN1959SN BN1960SNR  BN1974SN  BN1975SN  BN2003SN  BN2015SN
#     2455       683      6349      1582      1282      3424      1124      3184
# BN2029SN  BN2037SN  BN2050SN  BN2057SN  BN2149SN  BN2152SN  BN9944SN  BN9947SN
#     3031       666      3003      1921      1516      1509      5353      2353
# BN9966SN
#     2376

### only 1814 is suspicious

### let's check the agreement of these cells in these six samples

sharedCellinATACName <- names(ATACcell[(which(ATACcell %in% colnames(seRNA)))])

sharedCellinATACin6SNR <- rownames(proj1@cellColData)[which(rownames(proj1@cellColData) %in% sharedCellinATACName & as.character(proj1@cellColData$Sample) %in% c("BN1814SNR","BN1864SNR","BN1430SNR","BN1585SNR","BN1610SNR","BN1305SNR"))]

## subset only shared cells in these 6 samples

proj6 <- subsetCells(ArchRProj = proj1, cellNames = sharedCellinATACin6SNR)


### need to copy annotations for these cells only

scRNAanno <- read.table("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/data/annotation/2023_08_Midbrain_snRNA_annotation_xt.txt",header=TRUE, sep="\t")


subset23r <- function(x){
        gsub("SN#","SNR#",sub("_","#",substring(x,unlist(gregexpr('_',x))[1]+1,nchar(x))))
}

scRNAanno$Sample <- sapply(scRNAanno$cell,subset23r)

## add scRNA Cluster

proj6$Clusters_scRNA <- scRNAanno[match(rownames(proj6@cellColData),scRNAanno$Sample),"anno"]




###UMAP: for subset otherwise plotEmbedding will return error and plot all cells in original project
proj6 <- addUMAP(
    ArchRProj = proj6,
    reducedDims = "Harmony",
    name = "UMAPHarmony",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    force = TRUE
)

###plot by Manunal
p <- plotEmbedding(ArchRProj = proj6, colorBy = "cellColData", name = "ClustersManual", embedding = "UMAPHarmony")
plotPDF(p, name = "six-SNR-sample-shared-Plot-UMAP2Harmony-ClustersManual.pdf", ArchRProj = proj6, addDOC = FALSE, width = 5, height = 5)

## plot by scRNAanno
p <- plotEmbedding(ArchRProj = proj6, colorBy = "cellColData", name = "Clusters_scRNA", embedding = "UMAPHarmony")
plotPDF(p, name = "six-SNR-sample-shared-Plot-UMAP2Harmony-scRNA.pdf", ArchRProj = proj6, addDOC = FALSE, width = 5, height = 5)

## seem to agree with each other well
tdf <- proj6@cellColData[,c("Sample","Clusters_scRNA","ClustersManual")]

#table(tdf$Clusters_scRNA,tdf$ClustersManual)
#
#                      Astrocytes DA_neurons Endothelial_cells Fibroblasts
#  Astrocytes                1190          0                 0           0
#  DA neurons                   3          0                 0           0
#  DA-Glu-GABA neurons          0          1                 0           0
#  Endothelial cells            0          0               130           0
#  Fibroblasts                  0          0                 8           0
#  GABA neurons                 0          2                 0           0
#  Glu-GABA neurons             0          3                 0           0
#  Memory_CD8_T_Cells           0          0                12           0
#  Microglia                    3          0                 0           0
#  Monocytes                    0          0                 1           0
#  Oligodendrocytes            23          0                 3           0
#  OPCs                         3          2                 1           2
#  Pericytes                    0          0                72           0
#
#                      GABA_neurons Microglia Oligodendrocytes  OPCs
#  Astrocytes                     0         1                2     3
#  DA neurons                     0         0                1     0
#  DA-Glu-GABA neurons            0         0                0     0
#  Endothelial cells              1         8               14     0
#  Fibroblasts                    0         0                0     0
#  GABA neurons                  71         0                0     0
#  Glu-GABA neurons               3         0                0     0
#  Memory_CD8_T_Cells             0         2                2     0
#  Microglia                      0       884                7     0
#  Monocytes                      0        34                0     0
#  Oligodendrocytes               0         2            13871     1
#  OPCs                           1         1               10   704
#  Pericytes                      1         9                6     1
#
#####check existence in scRNA seurat project

seRNA$grouping <- ifelse(rownames(seRNA@meta.data) %in% ATACcell,"shared","RNA-unique")

p <- DimPlot(seRNA,group.by="grouping", raster=FALSE)

ggsave("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/01_archr/2023-08-v3/AllSamples/Plots/midbrain-scRNA-UMAP-omicExistence.pdf")

p <- DimPlot(seRNA,split.by="grouping", raster=FALSE)

ggsave("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/01_archr/2023-08-v3/AllSamples/Plots/midbrain-scRNA-UMAP-splitbyomicExistence.pdf")


## plot stacked bar plot for each cell type


p1 <- ggplot(seRNA@meta.data, aes(x=factor(anno,level=c("DA neurons","Glutamate neurons","GABA neurons","Glu-GABA neurons","DA-Glu-GABA neurons","Astrocytes","Oligodendrocytes","OPCs","Microglia","Monocytes","Memory_CD8_T_Cells","Endothelial cells","Pericytes","Fibroblasts")), fill=grouping)) + geom_bar() + theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90))  +
  theme(panel.grid = element_blank(),panel.border = element_blank())
ggsave("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/01_archr/2023-08-v3/AllSamples/Plots/midbrain-scRNA-stackedBarplot-byomicExistence.pdf")

p1 <- ggplot(seRNA@meta.data, aes(x=factor(anno,level=c("DA neurons","Glutamate neurons","GABA neurons","Glu-GABA neurons","DA-Glu-GABA neurons","Astrocytes","Oligodendrocytes","OPCs","Microglia","Monocytes","Memory_CD8_T_Cells","Endothelial cells","Pericytes","Fibroblasts")), fill=grouping)) + geom_bar(position = "fill") + theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90))  +
  theme(panel.grid = element_blank(),panel.border = element_blank())
ggsave("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/01_archr/2023-08-v3/AllSamples/Plots/midbrain-scRNA-stackedBarplot-byomicExistence-percent.pdf")


### Step 8.8: final cell type annotation

## update the number of overlapping cells

### previous version should be correct; as for those six samples BN1814SNR","BN1864SNR","BN1430SNR","BN1585SNR","BN1610SNR","BN1305SNR, they are either new tissue block or fill in samples

## check the performance of inferred-gene score-based annotation (Step 8.5)
## annotate the ATAC-unique cells (Step 8.6)
## annotate the ATAC-unique neurons (Step 8.6)

### Do neurons have lower TSS enrichment scores compared to other cell stypes? 

dfo <- as.data.frame(proj1@cellColData[which(!is.na(proj1$Clusters_scRNA)),]) 

library(ggplot2)

# grouped boxplot TSS
ggplot(dfo, aes(x=Clusters_scRNA, y=TSSEnrichment, fill=Diagnosis)) + 
    geom_boxplot()+ theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90))

ggsave("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/01_archr/2023-08-v3/AllSamples/Plots/midbrain-TSS-comparison-by-CellType-Diagnosis.pdf")

# grouped boxplot nFrags
ggplot(dfo, aes(x=Clusters_scRNA, y=nFrags, fill=Diagnosis)) +
    geom_boxplot()+ theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90))

ggsave("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/01_archr/2023-08-v3/AllSamples/Plots/midbrain-nFrags-comparison-by-CellType-Diagnosis.pdf")

###########################################################
###### finalize cell type annotation on 2023-09-26#########
###########################################################

library(ArchR)
set.seed(1)
addArchRThreads(threads = 10) 
library(here) 
addArchRGenome("hg38") 
setwd(here("2023-08-v3")) 
proj1 <- loadArchRProject(path="Save-Proj1-withAnnotation")


####plot cell numbers in ascending order per cell type using manual annotation

library(ggpubr)

dfPC <- as.data.frame(table(proj1$ClustersManual))

colnames(dfPC) <- c("cellType","count")

dfPC2 <- dfPC

dfPC2$cellType <- factor(dfPC2$cellType, levels=dfPC2$cellType[order(dfPC2$count,decreasing=TRUE)])


#ggbarplot(dfPC2, "cellType", y="count", fill = "cellType", label=TRUE,lab.vjust=0.5,lab.hjust=0.3,lab.size=3,orientation = "horiz") + theme(legend.position = "none")

ggbarplot(dfPC2, "cellType", y="count", fill = "cellType", label=TRUE,orientation = "horiz") + theme(legend.position = "none")

ggsave("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/01_archr/2023-08-v3/AllSamples/Plots/midbrain-cell-count-by-CellType-sorted.pdf")



####plot cell numbers in ascending order per cell type & diagnosis using manual annotation

proj1$ClustersManualDiagnosis <- paste0(proj1$Diagnosis,"_",proj1$ClustersManual)

dfPC <- as.data.frame(table(proj1$ClustersManualDiagnosis))

colnames(dfPC) <- c("cellType","count")


dfPC2 <- dfPC

dfPC2$cellType <- factor(dfPC2$cellType, levels=dfPC2$cellType[order(dfPC2$count,decreasing=TRUE)])


#ggbarplot(dfPC2, "cellType", y="count", fill = "cellType", label=TRUE,lab.vjust=0.5,lab.hjust=0.3,lab.size=3,orientation = "horiz") + theme(legend.position = "none")

ggbarplot(dfPC2, "cellType", y="count", fill = "cellType", label=TRUE,orientation = "horiz") + theme(legend.position = "none")

ggsave("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/01_archr/2023-08-v3/AllSamples/Plots/midbrain-cell-count-by-CellType-Diagnosis-sorted.pdf")


#### 



##reformat "BN1730SN_BN1730SN_AAACCGCGTCCACAAA-1" in column cell to "BN1730SN#AAACCGCGTCCACAAA-1" in ArchR project
## and change SN to SNR for 6 samples

scRNAanno <- read.table("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/data/annotation/2023_08_Midbrain_snRNA_annotation_xt.txt",header=TRUE, sep="\t")



subset23adj <- function(x){
        sample <-  substring(x,1,unlist(gregexpr('_',x))[1]-1)
        cell <- substring(x,unlist(gregexpr('_',x))[2]+1,nchar(x))
        sample <- ifelse(sample %in%  c("BN1814SN","BN1864SN","BN1430SN","BN1585SN","BN1610SN","BN1305SN"),paste0(sample,"R"),sample)
        paste0(sample,"#",cell)
}


scRNAanno$Sample <- sapply(scRNAanno$cell,subset23adj)

### check overlapping cells between scRNA and scATAC
# Load library
library(VennDiagram)
# Generate 2 sets of cells
cell_RNA <- scRNAanno$Sample
cell_ATAC <- rownames(proj1@cellColData)

# Chart
venn.diagram(
  x = list(cell_RNA, cell_ATAC),
  category.names = c("Ncell(RNA)" , "Ncell(ATAC)"),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("blue", "red"),

  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",

  ## output
  print.mode=c("raw","percent"),
  filename = here("2023-08-v3/AllSamples/Plots","Overlapping_cells_between_RNA_and_ATAC_IDadjusted.png"),
  output = TRUE
)



## add scRNA Cluster

proj1$Clusters_scRNA <- gsub(" ","_",scRNAanno[match(rownames(proj1@cellColData),scRNAanno$Sample),"anno"])

### re-run UMAP as the embeddings seem to change over the time 
###UMAP
proj1 <- addUMAP(
    ArchRProj = proj1,
    reducedDims = "Harmony",
    name = "UMAPHarmony",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    force = TRUE
)


#### get cells passed QC for each sample for each cell type and generate a text file per sample per cell type

df <- data.frame(Sample=proj1$Sample,Cell=paste0("CB:Z:",gsub("^BN\\S+#","",rownames(proj1@cellColData))),CellType=gsub(" ","_",proj1$ClustersManual))

for (celltype in unique(df$CellType)){
        ctDir <- paste0("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/02_peakcalling/cell/",celltype)
        dir.create(ctDir)
        ## input for subsetting cells from filtered bams for each sample
        for (sample in unique(df$Sample[which(df$CellType==celltype)])){
                dft <- data.frame(df[which(df$Sample==sample & df$CellType==celltype),"Cell"])
                write.table(dft, file=paste0(ctDir,"/",celltype,"-",sample,"-cell.txt"), quote=F, sep="\t", row.names=F, col.names=F)
        }
}

###subset bam file per sample per cell type and merge for each celltype



### ?exclude "Fibroblasts" (C20) and "Pericytes"(C21) at peak calling step

###finalize cluster: remove fibroblast and pericytes
### generate celltype-sample
### run peakcalling


### get cells passed QC for each sample for each cell type for each group and generate a text file per sample per cell type per group
## No need to re-genenrate cell.txt and can use the bams file from precious celltype level analysis
## but do need to generate file for each celltype for each group


df <- data.frame(Sample=proj1$Sample,Cell=paste0("CB:Z:",gsub("^BN\\S+#","",rownames(proj1@cellColData))),CellType=gsub(" ","_",proj1$ClustersManual), Diagnosis=proj1$Diagnosis)
ctDir <- "/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/02_peakcalling/diagnosis"

for (case in unique(df$Diagnosis)){
	
	df0 <- data.frame(df[which(df$Diagnosis==case),c("Diagnosis","CellType","Sample")])
	for (celltype in unique(df0$CellType)){
		## get all samples for each CellType per group
		dft <- unique(data.frame(df0[which(df0$CellType==celltype),]))
		write.table(dft, file=paste0(ctDir,"/",case,"-",celltype,"-sample.txt"), quote=F, sep="\t", row.names=F, col.names=F)
	}
 
}



### plot scatter plot for TH and Vglut2
## subset only DA neurons
projDA <-  subsetCells(ArchRProj = proj1, cellNames = proj1$cellNames[which(proj1$ClustersManual=="DA_neurons")])
 
gsMatrix <- getMatrixFromProject(ArchRProj = projDA,useMatrix = "GeneScoreMatrix")

THscore <- gsMatrix@assays@data@listData$GeneScoreMatrix[which(gsMatrix@elementMetadata$name == "TH"),]
VGLUT2score <- gsMatrix@assays@data@listData$GeneScoreMatrix[which(gsMatrix@elementMetadata$name == "SLC17A6"),]
VGLUT1score <- gsMatrix@assays@data@listData$GeneScoreMatrix[which(gsMatrix@elementMetadata$name == "SLC17A7"),]

DAdf <- as.data.frame(cbind(TH=THscore,VGLUT2=VGLUT2score,VGLUT1=VGLUT1score))

### use pairs

pdf(here("2023-08-v3/AllSamples/Plots/TH-VGLUT1-VGLUT2-genescore-pairs.pdf"))

pairs(DAdf)

dev.off()




### ArchR peak calling

proj2 <- addGroupCoverages(ArchRProj = proj1,minReplicates=3,maxReplicates=5, groupBy = "ClustersManual",force=TRUE)


pathToMacs2 <- findMacs2()

proj2 <- addReproduciblePeakSet(
    ArchRProj = proj2,
    groupBy = "ClustersManual",
    force=TRUE, 
    #threads =20,### temporary to speed up peak calling
    pathToMacs2 = pathToMacs2
)

####save  ArchRProject containing MACS2-derived merged peak set.

saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-Proj1-withPeak-ClustersManual", load = FALSE)



###

## Add peak matrix
#proj1 <- loadArchRProject(path = "Save-Proj1-withPeak")
proj3 <- addPeakMatrix(proj2)

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj3, 
    useMatrix = "PeakMatrix", 
    groupBy = "ClustersManual",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

#### get marker peaks

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

####Plotting Marker Peaks 

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap-ClustersManual", width = 8, height = 6, ArchRProj = proj3, addDOC = FALSE)


##### Step 11: Motif and Feature Enrichment 

#### motif annotation
proj3 <- addMotifAnnotations(ArchRProj = proj3, motifSet = "cisbp", name = "Motif")


####Motif Enrichment in Marker Peaks

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj3,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

## plot 

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap-ClustersManual", width = 8, height = 6, ArchRProj = proj3, addDOC = FALSE)

#### ENCODE enrichment

## ENCODE TFBS annotation

proj3 <- addArchRAnnotations(ArchRProj = proj3, collection = "EncodeTFBS")


## motif enrichment
enrichEncode <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj3,
    peakAnnotation = "EncodeTFBS",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

##plot

heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)

plotPDF(heatmapEncode, name = "EncodeTFBS-Enriched-Marker-Heatmap-ClustersManual", width = 8, height = 6, ArchRProj = proj3, addDOC = FALSE)



### ChromVAR Deviatons Enrichment

## add background peak sets

proj3 <- addBgdPeaks(proj3)

### compute per-cell deviation for cisbp

proj3 <- addDeviationsMatrix(
  ArchRProj = proj3, 
  peakAnnotation = "Motif",
  force = TRUE
)

## plot 

plotVarDev <- getVarDeviations(proj3, name = "MotifMatrix", plot = TRUE)

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores-ClustersManual", width = 5, height = 5, ArchRProj = proj3, addDOC = FALSE)

saveArchRProject(ArchRProj = proj3, outputDirectory = "Save-Proj2-motif-enrichment-ClustersManual", load = FALSE)







#### new cluster by integrating ClustersManual and Diagnosis


proj1$ClustersManualDiagnosis <- paste0(proj1$ClustersManual,"_",proj1$Diagnosis)

### ArchR peak calling

proj4 <- addGroupCoverages(ArchRProj = proj1,minReplicates=3,maxReplicates=5, groupBy = "ClustersManualDiagnosis",force=TRUE)


pathToMacs2 <- findMacs2()

proj4 <- addReproduciblePeakSet(
    ArchRProj = proj4,
    groupBy = "ClustersManualDiagnosis",
    force=TRUE, 
    #threads =20,### temporary to speed up peak calling
    pathToMacs2 = pathToMacs2
)

####save  ArchRProject containing MACS2-derived merged peak set.

#saveArchRProject(ArchRProj = proj4, outputDirectory = "Save-Proj1-withPeak-ClustersManualDiagnosis", load = FALSE)



###

## Add peak matrix
#proj1 <- loadArchRProject(path = "Save-Proj1-withPeak")
proj5 <- addPeakMatrix(proj4)

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj5, 
    useMatrix = "PeakMatrix", 
    groupBy = "ClustersManualDiagnosis",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

#### get marker peaks

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

####Plotting Marker Peaks 

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap-ClustersManualDiagnosis", width = 8, height = 6, ArchRProj = proj5, addDOC = FALSE)


##### Step 11: Motif and Feature Enrichment 

#### motif annotation
proj5 <- addMotifAnnotations(ArchRProj = proj5, motifSet = "cisbp", name = "Motif")


####Motif Enrichment in Marker Peaks

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

## plot 

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap-ClustersManualDiagnosis", width = 8, height = 6, ArchRProj = proj5, addDOC = FALSE)

### ChromVAR Deviatons Enrichment

## add background peak sets

proj5 <- addBgdPeaks(proj5)

### compute per-cell deviation for cisbp

proj5 <- addDeviationsMatrix(
  ArchRProj = proj5, 
  peakAnnotation = "Motif",
  force = TRUE
)

## plot 

plotVarDev <- getVarDeviations(proj5, name = "MotifMatrix", plot = TRUE)

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores-ClustersManualDiagnosis", width = 5, height = 5, ArchRProj = proj5, addDOC = FALSE)


####save  ArchRProject containing motif enrichment analysis.

saveArchRProject(ArchRProj = proj5, outputDirectory = "Save-Proj2-motif-enrichment-ClustersManualDiagnosis", load = FALSE)








####Step 9: Calling Peaks w/ Macs2


## Step 9.1: collect cells per sample per celltype for peakcalling

metaCellType <- as.data.frame(proj1@cellColData[,c("Sample","Clusters2")])
subsetCell <- function(x){
        substring(x,unlist(gregexpr('#',x))[1]+1,nchar(x))
}

metaCellType$cell <- sapply(rownames(metaCellType),subsetCell)

write.table(metaCellType,"/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/02_peakcalling/2023-08-28-cells-passed-QC.txt",quote=FALSE,row.names=FALSE,sep="\t") 


### group by cell type
#for(cType in unique(metaCellType$Clusters2)){
#       cTypeNew <- sub(" ","_",cType)
#       cTDir <- paste0("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/02_peakcalling/",cTypeNew)
#       dir.create(cTDir)
#       tdf <- metaCellType[grep(cType,metaCellType$Clusters2),]
#       write.table(tdf,paste0(cTDir,"/",cTypeNew,"_cellList.txt"),quote=FALSE,row.names=FALSE,sep="\t")        
#}
#




## Step 9.2: ArchR peak calling

proj1$Clusters4 <- paste0(proj1$Clusters_annoF,"_",proj1$Diagnosis)

#> table(proj1$Clusters4)
#
#          Astrocytes_HC         Astrocytes_ILBD           Astrocytes_PD
#                   9237                    9900                   11715
#   Endothelial cells_HC  Endothelial cells_ILBD    Endothelial cells_PD
#                   1053                     956                     829
#         Fibroblasts_HC        Fibroblasts_ILBD          Fibroblasts_PD
#                     64                      55                      31
#        GABA neurons_HC       GABA neurons_ILBD         GABA neurons_PD
#                   1054                    1317                     873
#    Glu-GABA neurons_HC   Glu-GABA neurons_ILBD     Glu-GABA neurons_PD
#                    514                     655                     394
#  Memory_CD8_T_Cells_HC Memory_CD8_T_Cells_ILBD   Memory_CD8_T_Cells_PD
#                     86                     103                      83
#           Microglia_HC          Microglia_ILBD            Microglia_PD
#                   9625                    7975                    7858
#           Monocytes_HC          Monocytes_ILBD            Monocytes_PD
#                    230                     296                     261
#    Oligodendrocytes_HC   Oligodendrocytes_ILBD     Oligodendrocytes_PD
#                  96700                   76131                   75065
#                OPCs_HC               OPCs_ILBD                 OPCs_PD
#                   6992                    4690                    4636
#       other_neurons_HC      other_neurons_ILBD        other_neurons_PD
#                     38                      32                       5
#           Pericytes_HC          Pericytes_ILBD            Pericytes_PD
#                    817                     666                     612

proj1 <- addGroupCoverages(ArchRProj = proj1,minReplicates=3,maxReplicates=5, groupBy = "Clusters4")
## for Save-Proj1-withPeakLarge
# proj1 <- addGroupCoverages(ArchRProj = proj1,minReplicates=8,maxReplicates=15, groupBy = "Clusters4")
pathToMacs2 <- findMacs2()

proj1 <- addReproduciblePeakSet(
    ArchRProj = proj1, 
    groupBy = "Clusters4", 
    #threads =20,### temporary to speed up peak calling
    pathToMacs2 = pathToMacs2
)

####save  ArchRProject containing MACS2-derived merged peak set.

saveArchRProject(ArchRProj = proj1, outputDirectory = "Save-Proj1-withPeak", load = FALSE)


### plot by number of cells for each subtype

#proj1 <- loadArchRProject(path = "Save-Proj1-withPeak")


## Step 9.3: MACS2-pseudobulk-fragments(all samples merged together)

## extract fragments

allFrags <- getFragmentsFromProject(ArchRProj=proj1)

##back to  0-BASED
## the original fragments.tsv from cellranger-arc count is 0-based, changed to 1-based in CreateArrow(): {data.table(V2=.$V2 + 1, V3=.$V3, V4=.$V4)}

###put into loop due to memory issue


for(i in 1:length(allFrags)){
  dftm <- data.frame(seqnames=seqnames((allFrags[[i]])),
  starts=start((allFrags[[i]]))-1,
  ends=end((allFrags[[i]])),
  names=c(rep(".", length((allFrags[[i]])))),
  scores=c(rep(".", length((allFrags[[i]])))),
  strands=strand((allFrags[[i]])))

  write.table(dftm, file=paste0("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/02_peakcalling/bed/",names(allFrags[i]),"-fragments.bed"), quote=F, sep="\t", row.names=F, col.names=F)
}
 ### failed as no strand information included in the original fragments output from cellranger-arc which is required by MACS2



## Step 9.4: MACS2-pseudobulk-bams(all samples merged together)

#### get cells passed QC for each sample and generate a text file per sample

proj1 <- loadArchRProject(path="Save-Proj1-withAnnotation")

###for each sample get cell identities

df <- data.frame(Sample=proj1$Sample,Cell=paste0("CB:Z:",gsub("^BN\\S+#","",rownames(proj1@cellColData))))

for (sample in unique(df$Sample)){
	## input for subsetting cells from filtered bams for each sample
	dft <- data.frame(df[which(df$Sample==sample),"Cell"])
	write.table(dft, file=paste0("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/02_peakcalling/cell/",sample,"-cell-passed-QC.txt"), quote=F, sep="\t", row.names=F, col.names=F)
}

#### get cells passed QC for each sample for each cell type and generate a text file per sample per cell type

proj1 <- loadArchRProject(path="Save-Proj1-withAnnotation")

###for each cell type each sample get cell identities
projt <- subsetCells(ArchRProj = proj1, cellNames = proj1$cellNames[which(!is.na(proj1$Clusters_scRNA))])


df <- data.frame(Sample=projt$Sample,Cell=paste0("CB:Z:",gsub("^BN\\S+#","",rownames(projt@cellColData))),CellType=gsub(" ","_",projt$Clusters_scRNA))

for (celltype in unique(df$CellType)){
	ctDir <- paste0("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/02_peakcalling/cell/",celltype)
	dir.create(ctDir)
        ## input for subsetting cells from filtered bams for each sample
	for (sample in unique(df$Sample[which(df$CellType==celltype)])){
        	dft <- data.frame(df[which(df$Sample==sample & df$CellType==celltype),"Cell"])
        	write.table(dft, file=paste0(ctDir,"/",celltype,"-",sample,"-cell.txt"), quote=F, sep="\t", row.names=F, col.names=F)
	}
}

###subset bam file per sample per cell type and merge for each celltype









#### Step 10: Identifying Marker Peaks

## Add peak matrix
#proj1 <- loadArchRProject(path = "Save-Proj1-withPeak")
proj2 <- addPeakMatrix(proj1)

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj2, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters4",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

#### get marker peaks

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

####Plotting Marker Peaks 

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj2, addDOC = FALSE)


##### Step 11: Motif and Feature Enrichment 

#### motif annotation
proj2 <- addMotifAnnotations(ArchRProj = proj2, motifSet = "cisbp", name = "Motif")


####Motif Enrichment in Marker Peaks

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj2,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

## plot 

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj2, addDOC = FALSE)

#### ENCODE enrichment

## ENCODE TFBS annotation

proj2 <- addArchRAnnotations(ArchRProj = proj2, collection = "EncodeTFBS")


## motif enrichment
enrichEncode <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj2,
    peakAnnotation = "EncodeTFBS",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

##plot

heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)

plotPDF(heatmapEncode, name = "EncodeTFBS-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj2, addDOC = FALSE)


### Step 12: ChromVAR Deviatons Enrichment

## add background peak sets

proj2 <- addBgdPeaks(proj2)

### compute per-cell deviation for cisbp

proj2 <- addDeviationsMatrix(
  ArchRProj = proj2, 
  peakAnnotation = "Motif",
  force = TRUE
)

## plot 

plotVarDev <- getVarDeviations(proj2, name = "MotifMatrix", plot = TRUE)

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj2, addDOC = FALSE)

### compute per-cell deviation for ENCODE TFBS

proj2 <- addDeviationsMatrix(
  ArchRProj = proj2, 
  peakAnnotation = "EncodeTFBS",
  force = TRUE
)


### plot

plotVarDev <- getVarDeviations(proj2, plot = TRUE, name = "EncodeTFBSMatrix")

plotPDF(plotVarDev, name = "Variable-EncodeTFBS-Deviation-Scores", width = 5, height = 5, ArchRProj = proj2, addDOC = FALSE)

####save  ArchRProject containing motif enrichment analysis.

saveArchRProject(ArchRProj = proj2, outputDirectory = "Save-Proj2-motif-enrichment", load = FALSE)




#### Step 13: Footprinting ########



#### SUPPLEMENT #####

###############################################################################################
######2023-12-14 get peak count matrix from MACS2 called peak set for each main cell type #####
###############################################################################################

library(ArchR)
set.seed(1)
addArchRThreads(threads = 10)
library(here)
addArchRGenome("hg38")
setwd("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/01_archr/2023-08-v3"))
proj1 <- loadArchRProject(path="Save-Proj1-withAnnotation")


### get bed file for peak set for main cell type #####

peakDF <- read.table("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/02_peakcalling/peaks/2023-09-26-celltype/intersection/2023-11-30-peak-intersection-using-comprehensive-peakset-halfOfEither-7-major-celltypes.txt",header=FALSE, sep="\t")


peakDF <- peakDF[,c(1,2,3,11,12)]
colnames(peakDF) <- c("chr","start","end","celltype","overlap")
peakDF$peak <- paste0(peakDF$chr,"_",peakDF$start,"_",peakDF$end)

### remove cell types that don't overlap with the peak
peakDF3 <- peakDF[which(peakDF$overlap>0),]

library(dplyr)


countDF3 <- peakDF3 %>% dplyr::count(peak)

### get celltype-specific peaks

ctPeak <- peakDF3[which(peakDF3$peak %in% countDF3$peak[which(countDF3$n==1)]),c("chr","start","end","celltype")]

## split by data frame

ctList <- split(ctPeak,f=ctPeak$celltype)

### save one bed file for each cell type's specific peaks

inDir <- "/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/02_peakcalling/peaks/2023-09-26-celltype/intersection"
for( i in 1:length(ctList))
  write.table(ctList[[i]],paste0(inDir,"/../uniquePeaks/cellTypeSpecific/2023-12-14-",names(ctList[i]),'-cell-type-specific-peak.bed'),quote=FALSE, row.names=FALSE,sep="\t")


#### 2023-12-15 ####
#### update bed file with other peak information in narrowPeak and share with Rosan

#### may not need for now #####


#> names(ctList)
#[1] "Astrocytes"        "DA_neurons"        "Endothelial_cells"
#[4] "GABA_neurons"      "Microglia"         "Oligodendrocytes" 
#[7] "OPCs"

#### add peak sets

#### add all celltype-specific peak at once

allCTPeak <- makeGRangesFromDataFrame(ctPeak)

proj2 <- addPeakSet(proj1, peakSet = allCTPeak, force=TRUE)

### add peak matrix

proj2 <- addPeakMatrix(proj2)

## get peak matrix
peakmtx = getMatrixFromProject(
  ArchRProj = proj2,
  useMatrix = "PeakMatrix",
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

saveRDS(peakmtx,paste0(inDir,"/../uniquePeaks/cellTypeSpecific/2023-12-14-all-celltype-specific.peak.mtx.rds"))

##save ArchR project

saveArchRProject(ArchRProj = proj2, outputDirectory = "2023-12-14-Save-Proj2-withAllCellTypePeakMatrix", load = FALSE)


### get peak pseudobulk count matrix for each cell type
### https://jef.works/blog/2020/04/06/quickly-creating-pseudobulks/
## pseudobulk accessibility at peak per sample 
## celltype in the function are samples
getPseudobulk <- function(mat, samples) {
	mat.summary <- do.call(cbind, lapply(levels(samples), function(sp) {
	#sps <- names(samples)[samples==sp] 
	sps <- which(samples==sp)
	pseudobulk <- Matrix::rowSums(as.matrix(mat[, sps]))
	return(pseudobulk)
	}))
	colnames(mat.summary) <- levels(samples)
	return(mat.summary)
}

### get sparse matrix

ctPeakMtx <- peakmtx@assays@data@listData$PeakMatrix

rownames(ctPeakMtx) <- paste(as.data.frame(peakmtx@rowRanges)$seqnames,as.data.frame(peakmtx@rowRanges)$start,as.data.frame(peakmtx@rowRanges)$end,sep="_")

###get per cell type
for (cellType in names(ctList)){
	## get peak subset for this cell type
	tmpDF <- as.data.frame(ctList[cellType])
	tmpPeak <- paste(tmpDF[,1],tmpDF[,2],tmpDF[,3],sep="_")
	## get cells for this cell type
	tmpCell <- proj2$cellNames[which(proj2$ClustersManual == cellType)]
	## subset the matrix for this cell type
	tmpMtx <- ctPeakMtx[tmpPeak,tmpCell]
	## replace the colnames with samples
	colnames(tmpMtx) <- gsub("\\#.*","",colnames(tmpMtx))
	## get pseudobulk peak matrix for this cell type
	samples <- as.factor(colnames(tmpMtx))
	##
	names(samples) <- colnames(tmpMtx)
	mat.summary <- getPseudobulk(tmpMtx, samples)		
	write.xlsx(data.frame(mat.summary),paste0(inDir,"/../uniquePeaks/cellTypeSpecific/2023-12-15-",cellType,"-specific-peak-pseudubulk-count-matrix.xlsx"),rowNames=TRUE)
	saveRDS(mat.summary,paste0(inDir,"/../uniquePeaks/cellTypeSpecific/2023-12-15-",cellType,"-specific-peak-pseudubulk-count-matrix.rds"))
}





#### save cell counts per sample per cell type

ct.summary <- data.frame(unclass(table(proj2$Sample,proj2$ClustersManual)))

## get rid of Fibroblasts and Pericytes

ct.summary <- ct.summary[,which(!colnames(ct.summary) %in% c("Fibroblasts","Pericytes"))]

ct.summary$Sample <- rownames(ct.summary)
ct.summary <- ct.summary[,c(8,1:7)]
##save output

library(openxlsx)
write.xlsx(ct.summary,paste0(inDir,"/../uniquePeaks/cellTypeSpecific/2023-12-14-sample-by-celltype-summary-table.xlsx"))





#####exclude Fibroblasts and Pericytes


cellsSample <- proj2$cellNames[which( proj2$ClustersManual %ni% c("Fibroblasts","Pericytes"))]

proj3 <- subsetArchRProject(
  ArchRProj = proj2,
  cells = cellsSample,
  outputDirectory = "2023-12-17-Save-Proj3-7major-cell-type",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)

###plot clustering with annotation

p1 <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "ClustersManual", embedding = "UMAPHarmony")

plotPDF(p1, name = "Plot-UMAPHarmony-ClustersManual.pdf", ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)


### plot their scRNA_label

## remove cells not existing in scRNA-seq
cellsSample <- proj3$cellNames[which(!proj3$Clusters_scRNA=="NA")]


###subset proj3 with only common cells in scRNA

projo <- subsetCells(ArchRProj = proj3, cellNames = proj3$cellNames[which(!is.na(proj3$Clusters_scRNA))])


###UMAP: for subset otherwise plotEmbedding will return error and plot all cells in original project
projo <- addUMAP(
    ArchRProj = projo,
    reducedDims = "Harmony",
    name = "UMAPHarmony",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    force = TRUE
)

###failed due to the 81th sample's PeakMatrix
#proj3_NAscRNA <- subsetArchRProject(
#  ArchRProj = proj3,
#  cells = cellsSample,
#  outputDirectory = "2023-12-18-Save-Proj3-7major-cell-type-scRNA-NA-removed",
#  dropCells = TRUE,
#  logFile = NULL,
#  threads = getArchRThreads(),
#  force = FALSE
#)



### plot manual annotation and scRNA annotation

#colorATAC <- c("Astrocytes"="#D51F26","DA_neurons"="#272E6A","Endothelial_cells"="#208A42","Fibroblasts"="#89288F","GABA_neurons"="#F47D2B","Microglia"="#FEE500","Oligodendrocytes"="#8A9FD1","OPCs"="#C06CAB","Pericytes"="#E6C2DC")
colorATAC <- c("Astrocytes"="#D51F26","DA_neurons"="#272E6A","Endothelial_cells"="#208A42","GABA_neurons"="#F47D2B","Microglia"="#FEE500","Oligodendrocytes"="#8A9FD1","OPCs"="#C06CAB")
p3 <- plotEmbedding(ArchRProj = projo, colorBy = "cellColData", name = "ClustersManual", embedding = "UMAPHarmony",pal=colorATAC, colorTitle=names(colorATAC))

#colorRNA <- c("Astrocytes"="#D51F26","DA neurons"="#272E6A","Endothelial cells"="#208A42","Fibroblasts"="#89288F","GABA neurons"="#F47D2B","Microglia"="#FEE500","Oligodendrocytes"="#8A9FD1","OPCs"="#C06CAB","Pericytes"="#E6C2DC","DA-Glu-GABA neurons"="#90D5E4","Glu-GABA neurons"="#89C75F","Glutamate neurons"="#F37B7D","Memory_CD8_T_Cells"="#9983BD","Monocytes"="#D24B27")
colorRNA <- c("Astrocytes"="#D51F26","DA neurons"="#272E6A","Endothelial cells"="#208A42","Fibroblasts"="#89288F","GABA neurons"="#F47D2B","Microglia"="#FEE500","Oligodendrocytes"="#8A9FD1","OPCs"="#C06CAB","Pericytes"="#E6C2DC","DA-Glu-GABA neurons"="#90D5E4","Glu-GABA neurons"="#89C75F","Glutamate neurons"="#F37B7D","Memory_CD8_T_Cells"="#9983BD","Monocytes"="#D24B27")
p4 <- plotEmbedding(ArchRProj = projo, colorBy = "cellColData", name = "Clusters_scRNA", embedding = "UMAPHarmony",pal=colorRNA,colorTitle=names(colorRNA))

plotPDF(p3,p4, name = "2023-12-18-Plot-UMAP2Harmony-Cluster-annotation-for-Overlapping-scRNAvsATAC-consistentColor.pdf", ArchRProj = projo, addDOC = FALSE, width = 5, height = 5)

### plot two in one

p <- ggAlignPlots(p3, p4, type = "h", draw=FALSE)

plotPDF(p, name = "2023-12-18-Plot-UMAP2Harmony-Cluster-annotation-for-Overlapping-scRNAvsATAC-consistentColor-2in1.pdf", ArchRProj = projo, addDOC = FALSE, width = 5, height = 5)





##################################
### 2023-12-18 aQTL analysis #####
## https://github.com/KellisLab/AD_regulome_analysis/tree/main/snATAC.processing/7.aQTL.calling

### Step a1: get the counts of fragments for each sample

inDir <- "/data/neurogen/ASAP/Multiomics/cellranger_multiome/2023_summer_cellranger-arc_count"

fragFiles <- list.files(inDir,pattern="atac_fragments.tsv.gz$",full.names=TRUE,recursive=TRUE)
## exclude old version for three resequenced samples: BN1610SNR_old, BN1862SN_old, and BN1959SN_old
## exclude cellranger_arc_aggr output for each batch and _combine output for batch3

fragFiles <- fragFiles[grep("_aggr|_combine|_old",fragFiles, invert=TRUE)]

## get sample ID

##write a function to extract the 3rd to last element

last3rd <- function(invect){
        last3 <- invect[length(invect)-2]
        return(last3)
}

names(fragFiles) <- sapply(strsplit(fragFiles,"/"),last3rd)

###count the number of fragments in each sample

depthTable <- data.frame()

for (i in 1:length(fragFiles)){
	#tmpDF <- read.table(gzfile(frag),)
	frag <- fragFiles[i]
	cmd <- paste0('zless -S ', frag, ' | grep -v "#" | wc -l')
	#print(names(fragFiles)[i])
	count <- system(cmd, intern=TRUE)	
	tmpDF <- data.frame(SAMPLE <- names(fragFiles)[i],COUNT=count)
	depthTable <- rbind(depthTable,tmpDF)
}
colnames(depthTable) <- c("SAMPLE","COUNT")

outDir <- "/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/03_aQTL"

write.table(depthTable,paste0(outDir,"/2023-12-20-ASAP-midbrain-fragments-count-per-sample.txt"),quote=FALSE,sep="\t",row.names=FALSE)

### Step a2: prepare covariate file

library(openxlsx)
metadata <- read.xlsx("/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/data/metadata/2023-07-26-midbrain-scATAC-metadata.xlsx")

metadata <- metadata[grep("SNR|repeat",metadata$Sample, invert=TRUE),c(2,4,5,6,8)]


rownames(metadata) <- metadata$Sample

metadata$Sex <- ifelse(metadata$Sex=="F",0,1)

metadata <- metadata[,2:5]

tmeta <- data.frame(t(metadata))

### remove BN1610

tmeta <- tmeta[,which(colnames(tmeta) != "BN1610SN")]

## write.table

write.table(data.frame("id"=rownames(tmeta),tmeta),paste0(outDir,"/2023-12-20-ASAP-midbrain-covariate.txt"), row.names=FALSE,quote=FALSE)


### Step a3: prepare genotype file
#gtFile <- read.table(paste0(outDir,"/2023-12-20-ASAP-midbrain-genotype.012"),header=FALSE,sep="\t")
#https://www.biostars.org/p/270984/

#!/usr/bin/Rscript

snps<-read.delim(paste0(outDir,"/2023-12-20-ASAP-midbrain-genotype.012"),header=F,row.names=1)
pos<-read.delim(paste0(outDir,"/2023-12-20-ASAP-midbrain-genotype.012.pos"),header=F)
indv<-read.delim(paste0(outDir,"/2023-12-20-ASAP-midbrain-genotype-sample.txt"),header=F)
colnames(snps)<-paste(pos[,1],pos[,2],sep=':')
rownames(snps)<-indv[,1]
snps<-as.matrix(snps)
snps.convert<-t(snps)
    
write.table(snps.convert, file= "snps.convert_all_pruned", col.names = FALSE, row.names = TRUE)

### file too large; need to split by chromosome


### Step a4: get peak file

### normalize by cell number and sequencing depth
### done at local due to cannot access to compute node on 2023-12-21
## /Users/me73/Dropbox (Partners HealthCare)/projects/2023-ASAP/local_Rscript/2023-12-21-peak-matrix-normalization.R


### Step a5: change sample names and order



