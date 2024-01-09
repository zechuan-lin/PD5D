#!/bin/bash


### detect cell type unique peaks using bedtools


## Usage: bash /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/src/bash/runBedtools4CellTypeUniquePeaks.sh  /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/02_peakcalling/peaks/2023-09-26-celltype 

inDir=$1

cd $inDir


beds=($(find . -type f -name '*-peakcalling_peaks.narrowPeak.blacklist.filtered.bed'))

for bed in "${beds[@]}"
do	
	echo $bed
	echo ${beds[@]/$bed}
	#array=("${bedss[@]/$bed}")
	cellType=$(basename $bed -peakcalling_peaks.narrowPeak.blacklist.filtered.bed)
	echo $cellType
	
	bedtools intersect -wa -wb -a ${cellType}-peakcalling_peaks.narrowPeak.blacklist.filtered.bed -b ${beds[@]/$bed}  tissue-peakcalling_peaks.narrowPeak.blacklist.filtered.bed -v > uniquePeaks/${cellType}-unique-peaks.bed	
done
