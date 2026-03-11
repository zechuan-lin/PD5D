#!usr/bin/R

##### This script perform quality control for each single cell sequencing sample ##############

library(sctransform)
library(Seurat)
library(harmony)

args = commandArgs()
input = args[6]
mydata = Read10X(data.dir = paste(input))

##### Next to create Seurat object, and control minimum number of cells per gene, minimum expressed gene per cell and proportion of MT reads

min_cell = 3
UF_min = 200#min of unique features/cell
UF_max = 9372.98 #pre compute max unqiue features, which is equals to three times of Median Absolute Deviation + median across all samples
mit_pro = 5#mt%

seuobj1 = CreateSeuratObject(counts = mydata, project = paste(input), min.cells = min_cell)# set min cells and min features (genes)
seuobj1[["percent.mt"]] = PercentageFeatureSet(seuobj1, pattern = "^MT-")#mt proportion


##### Conducting QC: just for cells since we want to combine gene count across cells, conducting gene QC may set counts below the cutoff (eg gene found > 3 cells) into zero during combining across cells ###########

qcseu1 = subset(seuobj1, subset = percent.mt < mit_pro)#filtering mt% < 5%
qcseu2 = subset(qcseu1, subset = nFeature_RNA > UF_min & nFeature_RNA < UF_max)#filtering with UMI min


mydata = as.data.frame(mydata)
fit1 = mydata[,colnames(qcseu2)]#get filtered cells based on barcode
colnames(fit1) = paste(input,"_",colnames(fit1),sep = "")

write.table(fit1,file = paste("exMatrix_",input,".txt",sep = ""),quote = F, row.names = T, col.names = T, sep = "\t")#print out the result after QC

qsum = matrix(nrow = 1, ncol = 9)
qsum[1,1] = paste(input)
qsum[1,2] = nrow(mydata) - nrow(seuobj1)
qsum[1,3] = ncol(mydata) - ncol(seuobj1)
qsum[1,4] = abs(nrow(qcseu1) - nrow(seuobj1))
qsum[1,5] = abs(ncol(qcseu1) - ncol(seuobj1))
qsum[1,6] = nrow(qcseu1) - nrow(qcseu2)
qsum[1,7] = ncol(qcseu1) - ncol(qcseu2)
qsum[1,8] = nrow(qcseu2)
qsum[1,9] = ncol(qcseu2)


write.table(qsum,file = "Summary_qc_with_seurat.txt",quote = F, append = T, sep = "\t",col.names = F, row.names = F)#output the QC summary of the file
