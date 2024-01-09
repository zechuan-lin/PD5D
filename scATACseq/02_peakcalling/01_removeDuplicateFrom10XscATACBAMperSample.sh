#!/bin/bash

## prepare BAM file for each sample
## first filter low quality reads
## second remove duplicates within the sample
## third filter by cell tags
## fourth convert to bed

## Usage: bash /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/src/bash/removeDuplicateFrom10XscATACBAMperSample.sh  /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/results/02_peakcalling/bam /data/neurogen/ASAP/Multiomics/cellranger_multiome/2023_summer_cellranger-arc_count /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/src/lsf/02_peakcalling/bambed

targetDir=$1 
bamIn=$2
lsfDir=$3

cd $targetDir

for d in $bamIn/* 
do	
	echo $d
	if [ -d $d ]; then
		sample=$(basename $d)
		#echo $sample
		bam="${bamIn}/${sample}/outs/atac_possorted_bam.bam"
		#echo $bam
	 	#if [ ! -f ${lsfDir}/${sample}_bambed.err ] && [ $sample != "BN0347SN" ] && [ $sample != "BN0339SN" ];then
		if [ ! -f ${lsfDir}/${sample}_bambed.err ];then	
			(echo $sample &&
			echo "
			#!/bin/bash
                        #BSUB -J ${lsfDir}/${sample}_bambed
                        #BSUB -o ${lsfDir}/${sample}_bambed.out
                        #BSUB -e ${lsfDir}/${sample}_bambed.err
                        module load samtools
			module load bedtools
			### set environment TMPDIR as a trial for solving 'cannot create temp file for here-document: No space left on device'
			##export TMPDIR=/data/bioinformatics/projects/donglab/ASAP_scATAC_2023/tmp
			##systemd-path temporary
			##doesnot seem to work
			#echo $TMPDIR
			##Step 1: filter low quality reads
			samtools view -b -F 0x4 -f 0x2 -q 30 -o ${sample}_01_filterLow.bam $bam
			
			##Step 2: remove duplicates
			## I. SAMTOOLS time consuming
			## sort by name
			samtools collate -o ${sample}_02_namecollate.bam ${sample}_01_filterLow.bam /data/bioinformatics/projects/donglab/ASAP_scATAC_2023/tmp/${sample}
			## Add ms and MC tags for markdup to use later
			samtools fixmate -m ${sample}_02_namecollate.bam ${sample}_03_fixmate.bam
			## Markdup needs position order
			samtools sort -o ${sample}_04_positionsort.bam ${sample}_03_fixmate.bam
			##Finally remove duplicates
			samtools markdup -r ${sample}_04_positionsort.bam  ${sample}_05_dupRM.bam
			samtools index ${sample}_05_dupRM.bam
			### clean house
			rm ${sample}_01_filterLow.bam ${sample}_02_namecollate.bam ${sample}_03_fixmate.bam ${sample}_04_positionsort.bam 	
			## II. PICARD MarkDuplicates
			#picard MarkDuplicates -I ${sample}_01_filterLow.bam  -O ${sample}_marked_duplicates.bam -M ${sample}_marked_dup_metrics.txt 
			#picard MarkDuplicates -I ${sample}_01_filterLow.bam -O ${sample}_dedup.bam -M ${sample}_dedup_metrics.txt -REMOVE_DUPLICATES true
			#samtools index ${sample}_dedup.bam
			#rm ${sample}_01_filterLow.bam
			###SWITCH BACK to SAMTOOLS due to complains from PICARD about the format of 10x bam files
			
			##Step 3: extract cells passed QC in ArchR
			# Save the header lines
			samtools view -H ${sample}_05_dupRM.bam > ${sample}_SAM_header
			# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
			samtools view ${sample}_05_dupRM.bam | LC_ALL=C grep -F -f $targetDir/../cell/${sample}-cell-passed-QC.txt > ${sample}_SAM_body
			# Combine header and body
			cat ${sample}_SAM_header ${sample}_SAM_body > ${sample}_06_filtered.sam
			# Convert filtered.sam to BAM format
			samtools view -b ${sample}_06_filtered.sam > ${sample}_06_filtered.bam
			samtools index ${sample}_06_filtered.bam
			rm ${sample}_SAM_body ${sample}_SAM_header ${sample}_06_filtered.sam
			# does not seem to be necessary
			# the BAM file should be sorted by read name beforehand
			#samtools sort -n -T ${sample}_aln.sorted -o ${sample}_07_sorted.bam ${sample}_06_filtered.bam
			#samtools index ${sample}_07_sorted.bam

			##Step 4:shift reads
			# The bedtools command should extract the paired-end alignments as bedpe format, then the awk command should shift the fragments as needed
			##One can easily use samtools and bamToBed together as part of a UNIX pipe. In this example, we will only convert properly-paired (FLAG == 0x2) reads to BED format.
			### need to escape $VARIABLE within awk in shell script, otherwise it will interpret by shell not awk
			samtools view -bf 0x2 ${sample}_06_filtered.bam | bedtools bamtobed -i stdin | awk -v OFS="\""\t"\"" '{if (\$6=="\""+"\""){print \$1,\$2+4,\$3+4,\$4,\$5,\$6} else if (\$6=="\""-"\""){print \$1,\$2-5,\$3-5,\$4,\$5,\$6}}' > ${sample}_fragments.bed " > ${lsfDir}/${sample}_bambed.lsf &&
                        bsub -q normal < ${lsfDir}/${sample}_bambed.lsf)
		fi
	fi 
done

