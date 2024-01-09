#!/bin/bash


## create libraries.csv for cellranger-arc for each sample 

## bash this.sh fastqDir outputDir sample.txt lsfDir

## Usage: bash writeLibraryCSV4CellrangerArc.sh /data/neurogen/ASAP/Multiomics/cellranger_multiome /data/neurogen/ASAP/Multiomics/cellranger_multiome/2023_summer_cellranger-arc_count  /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/data/metadata/midbrain_multiome_last_80_sample_label.txt /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/src/lsf/01_cellranger-arc

fastqDir=$1
targetDir=$2 
input=$3
lsfDir=$4

cd $targetDir

while read -r label sample;
do
	# find fastq files for GEX
	echo "$sample"
	echo "fastqs,sample,library_type" > ${sample}_libraries.csv
	paths=($(find $fastqDir -type f -name "${sample}_*.fastq.gz" | xargs dirname| sort | uniq))
	for path0 in ${paths[*]}
	do
		echo "$path0,$sample,Gene Expression" >> ${sample}_libraries.csv
	done 
	# find fastq files for ATAC
	echo "$label"
	paths=($(find $fastqDir -type f -name "${label}_*.fastq.gz" | xargs dirname| sort | uniq))
        for path0 in ${paths[*]}
        do
                echo "$path0,$label,Chromatin Accessibility" >> ${sample}_libraries.csv
        done	
	##submit jobs
	echo -e "
	#!/bin/bash
	#BSUB -J ${sample}
	#BSUB -o $lsfDir/${sample}.out
	#BSUB -e $lsfDir/${sample}.err
	
	source /PHShome/me73/miniconda3/etc/profile.d/conda.sh

	conda activate ArchR;
	
	cellranger-arc count --id=${sample} --reference=/data/bioinformatics/referenceGenome/Homo_sapiens/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 --libraries=${targetDir}/${sample}_libraries.csv --localcores=8 --localmem=64
	
	conda deactivate	
	" > $lsfDir/${sample}.lsf
	 bsub -q bigmem < $lsfDir/${sample}.lsf
done < "$input"
