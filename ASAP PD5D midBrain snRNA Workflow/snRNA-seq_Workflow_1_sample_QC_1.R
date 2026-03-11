#!usr/bin/R

#### This script using bulk sample gene expression to generate PC plot, sex concordant check, sample clustering and sample gene expression RLE plot for sample quality control purpose ############
#### Require package: optparse, factoextra, ggplot2 and reshape2 ############################################################################

library("optparse")

opt_list = list(make_option(c("-e", "--expression"), action="store", default=NA,help="The bulk sample gene count matrix with rows are genes and columns are samples",type="character"),
        make_option(c("-x", "--sex"), action="store", default=NA,help = "Sex information of samples, with rows are sample and a column for sex information, file header is required",type="character"),
        make_option(c("-o","--outfile"),action = "store",default = NA,help = "Prefix of the output file", type = "character")
        #make_option(c("-p","--refpos"),action = "store", default = NA, help = "Gene position on reference genome, with first column is gene name, second is chromosome, thrid and fourth is start and stop position.", type = "character")
        )


library(factoextra)
library(ggplot2)


mydata = read.table(paste(opt$input),head = T, row.names = 1, sep = "\t")#cols are individuals, raw are features

pc_ana = function(x){
    pc = prcomp(as.matrix(x),scale = T)
    return(pc)
}

##### Conduct PCA analysis ###############
mypc = pc_ana(mydata)

pdf(paste(input,".pdf",sep = ""),onefile = T)
par(mar=c(7,5,5,2))
fviz_eig(mypc)
#fviz_pca_ind(mypc,repel = TRUE,col.ind = c("#00AFBB"))#visulization contrubuting factors
pc_val = mypc$rotation
pc_val = as.data.frame(pc_val)
ggplot(data = pc_val,aes(x = pc_val[,1], y = pc_val[,2])) + geom_point(col = "darkred") + geom_text(aes(label = row.names(pc_val)),size = 1.8,col = c("#00AFBB")) + xlab("PC1") + ylab("PC2") + theme(plot.margin = unit(c(2,2,2,2), "cm"))
#plot(pc_val[,1],pc_val[,2],xlab = "PC1",ylab = "PC2",cex.lab = 1.5,cex.axis = 1.5,pch = 19,labels=row.names(pc_val))#here plot contribution of genes to PCs, too slow!!! disable it unless necessary

##### Next to clustering samples #########

no_norm = mydata
mydata = scale(mydata)

mydis = dist(t(mydata),method = "euclidean")
par(mar=c(6,5,5,2))
myclu = hclust(mydis,method = "average")
plot(myclu,cex = 0.4)

##### Next to plot sex genes: XIST for female and RPS4Y1 for male #############

sex_info = read.table(paste(opt$sex), sep = "", head = T)
par(mar=c(7,7,7,4))
plot(mydata[grep("XIST",row.names(mydata)),],mydata[grep("RPS4Y1",row.names(mydata)),],xlab = "XIST", ylab = "RPS4Y1",cex.lab = 1.5,col = as.factor(sex_info[,2]),pch = 19)
legend("topright",c("Female","Male"),col= c("black","red"),pch = 19)


##### Next to calculate RLE for genes expression ########################

library(reshape2)
no_norm[no_norm < 0] = 0
no_norm = no_norm[rowSums(no_norm) > 0,]
log_exp = log10(no_norm + 1e-4)
rle=log_exp-apply(log_exp, 1, median)
rle=melt(cbind(ID=rownames(rle), rle), variable.name = "Sample",value.name ="FPKM", id="ID");
bymedian <- with(rle, reorder(Sample, FPKM, IQR))
op=par(mar=c(7,3,3,1))
boxplot(FPKM ~ bymedian, data=rle, outline=F, las=2, boxwex=1, col='gray', cex.axis=0.3, main="Relative Log Expression", xlab="", ylab="RLE", frame=F,ylim=c(-4,2))

dev.off()

write.table(mypc$rotation,file = paste("PC.",args[6],".txt",sep = ""),row.names = TRUE, col.names = TRUE, quote = F, sep ="\t")
