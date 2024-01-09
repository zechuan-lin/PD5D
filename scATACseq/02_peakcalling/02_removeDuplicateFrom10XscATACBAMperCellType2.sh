#!/bin/bash

## prepare BAM file for each cell type annotated in ArchR
## first submit lsf jobs to extract read for all cells from corresponding sample's filtered-deduplicated-cellsPassedQC BAM
## after this: merge per cell type


## Usage: bash /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/src/bash/removeDuplicateFrom10XscATACBAMperCellType2.sh  /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/02_peakcalling /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/src/lsf/02_peakcalling/bamCellType

targetDir=$1/bam 
cellIn=$1/cell
lsfDir=$2

cd $targetDir

for d in $cellIn/* 
do	
	echo $d
        if [ -d $d ]; then
		cellType=$(basename $d)
		echo $cellType
		mkdir -p $cellType
		for file in $d/*-cell.txt
		do	
			(sfile=$(basename $file) &&
			sample=$(echo $sfile |cut -d'-' -f2 ) &&
			echo $sfile &&
			echo $sample &&
			bam="${targetDir}/${sample}_06_filtered.bam"
			if [ ! -f {lsfDir}/${cellType}_${sample}_bam.err ];then 	
				echo "
				#!/bin/bash
				#BSUB -J ${lsfDir}/${cellType}_${sample}_bam
				#BSUB -o ${lsfDir}/${cellType}_${sample}_bam.out
				#BSUB -e ${lsfDir}/${cellType}_${sample}_bam.err
				module load samtools
				cd $targetDir
				#bam="\""${targetDir}/${sample}_06_filtered.bam"\""
				#echo $bam
				##extract header
				samtools view -H $bam > ${cellType}/${sample}_SAM_header
				# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
				samtools view $bam | LC_ALL=C grep -F -f $file > ${cellType}/${sample}_SAM_body
				cat ${cellType}/${sample}_SAM_header ${cellType}/${sample}_SAM_body > ${cellType}/${cellType}_${sample}.sam
				samtools view -b ${cellType}/${cellType}_${sample}.sam > ${cellType}/${cellType}_${sample}.bam 
				rm ${cellType}/${cellType}_${sample}.sam ${cellType}/${sample}_SAM_body ${cellType}/${sample}_SAM_header
				samtools index ${cellType}/${cellType}_${sample}.bam
				samtools sort -o ${cellType}/${cellType}_${sample}_sorted.bam ${cellType}/${cellType}_${sample}.bam 
				rm ${cellType}/${cellType}_${sample}.bam ${cellType}/${cellType}_${sample}.bam.bai
				samtools index ${cellType}/${cellType}_${sample}_sorted.bam " > ${lsfDir}/${cellType}_${sample}_bam.lsf &&
				bsub -q normal < ${lsfDir}/${cellType}_${sample}_bam.lsf
			fi)
		done
	fi

done
