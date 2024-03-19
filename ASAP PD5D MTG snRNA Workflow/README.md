# **PD5D snRNA-seq Workflow**

Code applied to snRNA-seq data from the middle temporal gyrus, part of the ASAP PD5D Parkinson's Cell Atlas project. Analysis steps are split into scripts that should be run in numerical order. The workflow is designed to be used with outputs from the 10X Cell Ranger pipeline cellranger count. An aggregated matrix of all samples is required for the sample quality control script (snRNA-seq_Workflow_1_sample_QC_1.R). snRNA-seq_Workflow_2_cell_QC_for_count_matrix.R performs cell QC for each sample individually, which can be submitted for all samples by the script Batch_run_snRNA-seq_Workflow_2_cell_QC_for_count_matrix.pl, providing the cellranger count output matrix folder for every sample is in the run directory. After QC the sample matrices are then aggregated by snRNA-seq_Workflow_3_merge_count_matrices.pl. snRNA-seq_Workflow_4_GLMPCA_and_harmony.R then takes the aggregate matrix as an input and outputs a Seurat Object. snRNA-seq_Workflow_4_GLMPCA_and_harmony.R also reads in a set of bespoke metadata so this section will require editing by other users to account for their own set of metadata. Again, it should be noted that later code may depend upon specific metadata columns in the Seurat Object, and so will need to be manually adjusted by other users. For small datasets it may be easier to skip scripts 2 and 3 and do the cell quality control and data aggregation in R itself as part of a Seurat workflow by customizing script 4. The proceeding scripts (5-7) then all take the Seurat Object from the preceding script as their input.

## **snRNA-seq_Workflow_1_sample_QC_1.R**

•	PC plot
 
•	Sex concordance check

•	Sample clustering

•	Sample gene expression RLE plot

•	Outliers manually selected based on outputs

**Command line**

	Rscript single_cell_sample_QC.R --expression [Bulk_sample_expression_matrix].txt --sex [sex_information_table].txt --outfile [output_prefix]

	•	--expression/-e: The bulk sample gene count matrix where rows are genes and columns are samples
	•	--sex/-x: Sex information of samples, with rows are sample and a column for sex information, file header is required
	•	--outfile/-o: Prefix of the output file

**Outputs**

Files

•	Sample PC values (PC.[output_prefix].txt)

Figures

•	Sample PC plot, sample Hierarchical Clustering plot, sample sex concordant check plot ([output_prefix].pdf)


## **Batch_run_snRNA-seq_Workflow_2_cell_QC_for_count_matrix.pl**

•	This script runs single cell QC for count matrix for each sample

•	Put this script and the R script snRNA-seq_Workflow_2_cell_QC_for_count_matrix.R (see below) within a folder containing the Cell Ranger count output matrices folders for all your samples in order to run the R QC script across all your samples. You should rename each sample output matrices folder to be the sample id so the outputs from the R QC script - filtered count matrices - can be distinguished from one another.

**Command line**

	perl Batch_run_QC_with_Single_cell_QC_for_count_matrix.pl

## **snRNA-seq_Workflow_2_cell_QC_for_count_matrix.R**

•	Genes expressed in less than 3 cells removed

•	Cells with >= 5% mitochondrial counts removed

•	Cells with <= 200 or >= median+3*median average deviation unique genes removed (MAD calculated prior)

**Command line**

	Rscript Single_cell_QC_for_count_matrix.R [CellRanger_output_folder_of_a_sample]

**Outputs**

Files
	
•	A count matrix for the relevant sample after QC (exMatrix_[CellRanger_output_folder_of_a_sample].txt)
 
•	A QC summary table with the number of cells and genes that were removed during QC, and the final number of cells and genes retained in the data (Summary_qc_with_seurat.txt)

## **snRNA-seq_Workflow_3_merge_count_matrices.pl**

•	Merges all of the filtered count matrices for each sample (output from the previous script)

•	The following line “my @file = @ARGV = <exMatrix_*>;” should automatically detect the default outputs from the previous script in run directory of this script, but if the output file name pattern for the previous script is changed, section between the <> should be altered to pattern match the new output file names

**Command line**

	perl snRNA-seq_Workflow_3_merge_count_matrices.pl [gene_list_with_gene_each_row].txt [output_file]

**Outputs**

Files

•	combined matrix of all samples ([output_file])

## **snRNA-seq_Workflow_4_GLMPCA_and_harmony**

•	Seurat Object assembly

•	Library normalization and log1p transformation

•	Variable feature selection - 2000 features

•	Dimensionality reduction - GLM-PCA (2019)

•	Batch correction using Harmony (2019)

**Command line**

	Rscript snRNA-seq_Workflow_4_GLMPCA_and_harmony.R --matrix [matrix_path] --metadata [metadata_path] --theta [theta_for_harmony] --ncovariates [number_of_harmony_covariates] --ldglm [glm_pca_dimensions] --outfile [outfile]

	•	--matrix/-f: Path to filtered matrix, output from previous script
	•	--metadata/-m: Metadata where each row corresponds to a column of --matrix, in the same order
	•	--theta/-t: Theta value per covariate for harmony, ideal total should not exceed 2-2.5 (default = 0.4)
	•	--ncovariates/-c: Number of harmony covariates (default = 6)
	•	--ldglm/-l: Latent dimensions for GLM-PCA, should be large enough to be able to see the ideal cut-off on the elbow plot output
	•	--outfile/-o: name of output Seurat Object

**Outputs**

Files

•	Seurat Object ([outfile])

•	Seurat Object metadata (InitialMetadata.tsv)

## **snRNA-seq_Workflow_5_clustering.R**

•	Clustering using Leiden algorithm

•	Remove mitochondrial-encoded genes to prevent confounding effects with cell quality in downstream analyses

**Command line**

	Rscript snRNA-seq_Workflow_5_clustering.R --seurat [seuratpath] --dims [dimensions] --resolution [resolution] --outfile [outfile]

	•	--seurat/-f: Path to Seurat Object, output from previous script
	•	--dims/-d: Dimensions of reduction to use for FindNeighbours, based on PC cut off from elbow plot output from previous script
	•	--resolution/-r: Resolution to pass to FindClusters - run FindClusters across range of resolutions script to find suitable value
	•	--outfile/-o: name of output Seurat Object

**Outputs**

Files

•	Seurat Object ([outfile])

•	Cells per cluster table (Cells_Per_Unassigned_Cluster_Table.tsv)

•	Samples per cluster table (Samples_Per_Unassigned_Cluster_Table.tsv)

## **snRNA-seq_Workflow_6_cluster_and_doublet_qc.R**

•	Doublet identification and removal using the scDblFinder package (2021)

•	Removal of clusters that are constituted by >= user defined % (default = 30%) doublet cells

•	Removal of clusters representing low numbers of cell/samples (user defined parameters)

•	Re-run UMAP after cluster QC to clean up visualization

**Command line** 

	Rscript snRNA-seq_Workflow_6_cluster_and_doublet_qc.R --seurat [seuratpath] --dims [dimensions] --mincells [minimumcells] --minsamples [minimumsamples] --pthreshold [percentagethreshold] --outfile [outfile]

	•	--seurat/-f: Path to Seurat Object, output from previous script
	•	--dims/-d: Dimensions of reduction used for FindNeighbours in previous script
	•	--mincells/-c: minimum number of cells required for retention of cluster, default (200) based on ~600,000 cells, scale accordingly
	•	--minsamples/-s: minimum number of samples required for retention of cluster, default (20) based on ~100 samples, scale accordingly
	•	--pthreshold/-p: doublet percentage threshold - clusters that have a doublet cell percentage that equals or exceeds this threshold will be excluded (default = 30)
	•	--outfile/-o: name of output Seurat Object

**Outputs**

Files

•	Seurat Object ([outfile])
 
•	Doublet percentage per cluster table (Doublet_Percentage_Per_Cluster_Table.tsv)

•	Cells per cluster table (Cells_Per_Unassigned_Cluster_Table_PostClusterQC.tsv)

•	Samples per cluster table (Samples_Per_Unassigned_Cluster_Table_PostClusterQC.tsv)

•	Seurat Object metadata (PostClusterQCMetadata.tsv)

Figures

•	UMAP colored by doublets and singlets (pre-doublet removal) (Unassigned_Doublet_UMAPclusters_scRNA_seq.pdf)
	
•	Violin plot of total counts per cell, split by cluster (post-doublet removal) (nCount_per_Cluster_ViolinPlot.pdf)
	
•	Violin plot of total unique features per cell, split by cluster (post-doublet removal)(nFeature_per_Cluster_ViolinPlot.pdf)

## **snRNA-seq_Workflow_7_UMAP_and_marker_genes.R**

•	Generate UMAPs and marker gene expression figures to characterize clusters

**Command line**

	Rscript snRNA-seq_Workflow_7_UMAP_and_marker_genes.R --seurat [seuratpath]

	•	--seurat/-f: Path to Seurat Object, output from previous script
	•	--markers/-m: Path to plain text file of marker genes to visualise expression of, no header, one per line

**Outputs**

Figures

•	UMAP (SeuratObject_UMAP_Clusters_PostClusterQC.pdf)
	
•	Split UMAP - default: by case (disease status) (SeuratObject_Case_Split_UMAP_Clusters_PostClusterQC.pdf)

•	Barchart of marker gene expression across all clusters (Marker_Barchart.pdf)

•	Violin plot of each marker gene across all cluster ([genename]_VlnPlot.pdf)


## **Key Packages/Requirements**

	•	R (4.2.3)
	•	Seurat (4.3.0)
	•	SeuratWrappers (0.3.0)
	•	SeuratObject (4.1.3)
	•	glmpca (0.2.0)
	•	harmony (1.2.0)
	•	SingleCellExperiment (1.20.1)
	•	scDblFinder (1.12.0)
	•	tidyverse packages (ggplot2 (3.4.4), dplyr (1.0.9), tidyr (1.2.0), stringr (1.4.0))
	•	reshape2 (1.4.4)
	•	leidenalg (0.9.1, via conda-forge)
	•	optparse (1.7.4)


## **Contributions and Maintenance**

This code is developed and maintained by Dr. Jacob Parker and Dr. Zechuan Lin

Queries regarding scripts 1 and 2: zechuan.lin@yale.edu
Queries regarding scripts 2 - 4: Jacob.parker@yale.edu

## **License**

Copyright 2024 Jacob Parker and Zechuan Lin
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
