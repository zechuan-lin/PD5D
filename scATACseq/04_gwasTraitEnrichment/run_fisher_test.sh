
GWAS=/data/bioinformatics/external_data/externalData/GWAS_catalog_20220810
wd=/PHShome/rw552


for input in maoxuan/*.counts;
do
Rscript $wd/fisher_test.R $input $GWAS/GWAS_20220810.v1.02.counts.v2 $input.pdf
done