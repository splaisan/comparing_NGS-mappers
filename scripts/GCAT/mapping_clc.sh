#! /usr/bin/env bash
## script: 'mapping_clc.sh'
## ©SP-BITS, 2015 v1.0
# last edit: 2015-01-07

# required:
# CLCGWB
# Picard Version (broad build): 1.127

## define global variables
reffasta=$STAR_INDEXES/hg19_24.fa
refgtf=$STAR_INDEXES//hg19_24.gtf
readlen=100
refidx=$STAR_INDEXES/index_${readlen}

# create new folder
outfolder=hg19_CLC-mapping
mkdir -p ${outfolder}

## reads
infolder=reads
f=gcat_set_041

fq1=${infolder}/${f}_1.fq.gz
fq2=${infolder}/${f}_2.fq.gz

clcpe=${f}_clc-pe

# using 'nthr' processors in parallel (again limited by our RAM!)
# mem is more demanding and needs more than 3Gb per thread
nthr=8
maxmem="48G"

# log all to file
logfile=${outfolder}/clc_mapping-log.txt
exec &> >(tee -a "${logfile}")

# report run duration
starttime=$(date +%s)

## MAPPING DONE USING CLCGWN CLI

########################### post process results ##############
## Add RG tag in header
rgstring="@RG\tID:NA12878\tSM:gcat_set_041-CLC\tPL:ILLUMINA"

cmd="bam polishBam --RG \"${rgstring}\" \
	-i ${outfolder}/gcat_set_041_clc.bam \
	-o ${outfolder}/gcat_set_041_clc-rg.bam \
	-l ${outfolder}/polishBam_${clcpe}.log && \
	rm ${outfolder}/${clcpe}_Aligned.out.bam"

# execute the command
echo
echo "# ${cmd}"
eval ${cmd}

##### convert to sorted bam, index & cleanup
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar SortSam \
	I=${outfolder}/gcat_set_041_clc-rg.bam \
	O=${outfolder}/${clcpe}.bam \
	SO=coordinate \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=LENIENT"

echo
echo "# ${cmd}"
eval ${cmd}

## ValidateSamFile ##
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar ValidateSamFile \
	I=${outfolder}/${clcpe}.bam \
	O=${outfolder}/${clcpe}-mappings_ValidateSam.txt \
	M=SUMMARY"

echo
echo "# ${cmd}"
eval ${cmd}

##### get basic stats from the resulting bam file
cmd="samtools flagstat \
	${outfolder}/${clcpe}.bam \
	>${outfolder}/${clcpe}-mappings_flagstat.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## samtools stats ##
cmd="samtools stats \
	${outfolder}/${clcpe}.bam \
	>${outfolder}/${clcpe}-mappings_stats.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## samtools plot-bamstats ##
cmd="plot-bamstats \
	${outfolder}/${clcpe}-mappings_stats.txt \
	-p ${outfolder}/${clcpe}-bamstats_plots/clc"

echo
echo "# ${cmd}"
eval ${cmd}

## plot mapping quality score distribution
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar \
	QualityScoreDistribution \
	CHART=${outfolder}/QS_distribution.pdf \
	I=${outfolder}/${clcpe}.bam \
	O=${outfolder}/${clcpe}-mappings_QS.txt"

echo
echo "# ${cmd}"
eval ${cmd}

#=====================
# report run duration
#=====================

endtime=$(date +%s)
dur=$(echo ${endtime}-${starttime} | bc)
echo
echo "# ended at $date: $endtime"
echo "Done in ${dur} sec"
