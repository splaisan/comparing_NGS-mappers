#! /usr/bin/env bash
## script: 'mapping_hisat.sh'
## ©SP-BITS, 2015 v1.0
# last edit: 2015-04-28

# required:
# HISAT version b0.6
# Picard Version (broad build): 1.120

## define global variables
reffasta=$HISAT_INDEXES/hg19_24.fa
refidx=$HISAT_INDEXES/hg19_24

# create new folder
outfolder=hg19_hisat-mapping
mkdir -p ${outfolder}

## reads
infolder=reads
f=gcat_set_041

fq1=${infolder}/${f}_1.fq.gz
fq2=${infolder}/${f}_2.fq.gz

hisatpe=${f}_hisat-pe

# using 'nthr' processors in parallel (again limited by our RAM!)
# mem is more demanding and needs more than 3Gb per thread
nthr=8
maxmem="48G"

# log all to file
logfile=${outfolder}/hisat_mapping-log.txt
exec &> >(tee -a "${logfile}")

# report run duration
starttime=$(date +%s)

# build HISAT command for zipped paired fastq files
## compute length-1
idxlen=$( echo ${readlen}-1 | bc)

cmd="hisat \
	-p ${nthr} \
	-x ${refidx} \
	-q \
	--phred33 \
	-1 ${fq1} \
	-2 ${fq2} \
	-S ${outfolder}/${hisatpe}.sam"

echo
echo "# ${cmd}"
eval ${cmd}

########################### post process results ##############

## Add RG tag in header
rgstring="@RG\tID:NA12878\tSM:gcat_set_041-HISAT\tPL:ILLUMINA"

cmd2="bam polishBam \
	--RG \"${rgstring}\" \
	--CO \"@CO\t${cmd}\" \
	-i ${outfolder}/${hisatpe}.sam \
	-o ${outfolder}/${hisatpe}_q.bam \
	-l ${outfolder}/polishBam_${hisatpe}.log"

# execute the command
echo
echo "# ${cmd2}"
eval ${cmd2}

##### convert to sorted bam, index & cleanup
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar SortSam \
	I=${outfolder}/${hisatpe}_q.bam \
	O=${outfolder}/${hisatpe}.bam \
	SO=coordinate \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=LENIENT"

echo
echo "# ${cmd}"
eval ${cmd}

## ValidateSamFile ##
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar ValidateSamFile \
	I=${outfolder}/${hisatpe}.bam \
	O=${outfolder}/${hisatpe}-mappings_ValidateSam.txt \
	M=SUMMARY"

echo
echo "# ${cmd}"
eval ${cmd}

##### get basic stats from the resulting bam file
cmd="samtools flagstat \
	${outfolder}/${hisatpe}.bam \
	>${outfolder}/${hisatpe}-mappings_flagstat.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## samtools stats ##
cmd="samtools stats \
	${outfolder}/${hisatpe}.bam \
	>${outfolder}/${hisatpe}-mappings_stats.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## samtools plot-bamstats ##
cmd="plot-bamstats \
	${outfolder}/${hisatpe}-mappings_stats.txt \
	-p ${outfolder}/${hisatpe}-bamstats_plots/hisat"

echo
echo "# ${cmd}"
eval ${cmd}

## plot mapping quality score distribution
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar \
	QualityScoreDistribution \
	CHART=${outfolder}/QS_distribution.pdf \
	I=${outfolder}/${hisatpe}.bam \
	O=${outfolder}/${hisatpe}-mappings_QS.txt"

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