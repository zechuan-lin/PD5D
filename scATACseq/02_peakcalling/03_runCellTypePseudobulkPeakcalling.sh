#!/bin/bash

## celltype specific peak calling

## input: CELLTYPE_SAMPLE.bed from /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/src/bash/removeDuplicateFrom10XscATACBAMperCellType2.sh

## first: merge per cell type; sort bams

## second: convert to bed

## third: peak calling

## fourth: blacklist filtering

## Usage: bash /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/src/bash/runCellTypePseudobulkPeakcalling.sh  /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/02_peakcalling  /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/src/lsf/02_peakcalling/peakcalling


inDir=$1/bam 
peakDir=$1/peaks
lsfDir=$2

cd $inDir

for d in $inDir/* 
do	
	echo $d
        if [ -d $d ]; then
		cellType=$(basename $d)
		echo $cellType
		## merge bams for each cell type
		if [ ! -f {lsfDir}/${cellType}_peakcalling.err ];then

			(echo "
			#!/bin/bash
                        #BSUB -J ${lsfDir}/${cellType}_peakcalling
                        #BSUB -o ${lsfDir}/${cellType}_peakcalling.out
                        #BSUB -e ${lsfDir}/${cellType}_peakcalling.err
                        
			module load samtools
			module load bedtools
			
			## merge all cells per celltype
			cd $inDir/${cellType} 
			samtools merge ${cellType}_merged.bam ${cellType}_*_sorted.bam
			samtools index ${cellType}_merged.bam
			samtools sort -o ${cellType}_final.bam ${cellType}_merged.bam
			samtools index ${cellType}_final.bam
			rm ${cellType}_merged.bam ${cellType}_merged.bam.bai
			
			## shift reads
			samtools view -bf 0x2 ${cellType}_final.bam | bedtools bamtobed -i stdin | awk -v OFS="\""\t"\"" '{if (\$6=="\""+"\""){print \$1,\$2+4,\$3+4,\$4,\$5,\$6} else if (\$6=="\""-"\""){print \$1,\$2-5,\$3-5,\$4,\$5,\$6}}' > ${cellType}_fragments.bed 
			
			## peak calling
			cd $peakDir
			macs2 callpeak -f BED -t $inDir/${cellType}/${cellType}_fragments.bed -n ${cellType}-peakcalling  --outdir 2023-09-26-celltype --nomodel --shift -100 --extsize 200 -g hs --keep-dup all -B
			## filter blacklist regions
			bedtools intersect -v -a $peakDir/2023-09-26-celltype/${cellType}-peakcalling_peaks.narrowPeak -b /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/src/lsf/02_peakcalling/ENCFF356LFX.bed.gz | awk 'BEGIN{OFS="\""\t"\""} {if (\$5>1000) \$5=1000; print \$0}'| grep -P 'chr[0-9XY]+(?!_)' | gzip -nc > $peakDir/2023-09-26-celltype/${cellType}-peakcalling_peaks.narrowPeak.blacklist.filtered.bed.gz

			zless -S $peakDir/2023-09-26-celltype/${cellType}-peakcalling_peaks.narrowPeak.blacklist.filtered.bed.gz > $peakDir/2023-09-26-celltype/${cellType}-peakcalling_peaks.narrowPeak.blacklist.filtered.bed 
			" > ${lsfDir}/${cellType}_peakcalling.lsf &&
                        bsub -q normal < ${lsfDir}/${cellType}_peakcalling.lsf)
		fi
	fi

done
