# **PD5D Midbrain Visium Workflow**

This processing pipeline consists of a pair of scripts that perform aggregation, QC (ST_preprocess.R), integration and clustering (ST_pipeline.R) on midbrain visium spatial transcriptomic data as part of the ASAP PD5D Parkinson's Cell Atlas project.

The input to ST_preprocess.R is a csv file (infiles.csv) with the following column names and information:

| Column name | Contents |
| :---------- | :----------: |
| sample_name | sample identifier, duplicates tolerated to allow for replicates and multiple runs covering different regions of the same tissue slice |
| dir | If manual removal of unwanted regions was performed on sample (e.g. regions containing tissue folds), field contains a path to a Seurat object (.rds file) of the amended data. If no amendments were made, field contains the path to the spaceranger count output directory for the sample |
| diagnosis | clinical diagnosis (HC – healthy control, ILB – incidental Lewy body disease, PD – Parkinson’s disease) |
| fold_removed | Boolean value denoting whether the manual removed of unwanted regions was performed (e.g. tissue folds). Tells script whether to read in an amended Seurat object (.rds file) or use Load10X_Spatial() to read in unprocessed spaceranger count output files |
| batch | batch number, numeric |
| sex | subject sex, either “M” or “F” |
| expired_age | expired age, numeric |
| PMI | post-mortem interval, numeric |
| RIN | RNA integrity number, numeric |
| sample_spatial | unique identifier consisting of sample identifier in sample_name column + spatial region information (medial, lateral, ventral, dorsal) |

The input to ST_pipeline.R is the output from ST_preprocess.R, so the scripts can be run consecutively without manual interference.

ST_pipeline.R performs sample integration using the Harmony package and runs the Seurat function FindClusters() across a range of resolutions that can be specified in-script. For each resolution the script outputs a Seurat object (.rds) containing the cluster annotations and a UMAP plot to enable optimal resolution selection.

Code written by Dr. Xufei Teng, modified and uploaded by Dr. Jacob Parker
