#! /usr/bin/env bash
## script: 'mapping_hisat.sh'
## ©SP-BITS, 2015 v1.0
# last edit: 2015-04-29

# required:
# HISAT version HISAT b0.6
# Picard Version (broad build): 1.120

## define global variables
reffasta=$STAR_INDEXES/hg19_24.fa
refidx=$HISAT_INDEXES/hg19_24

# create new folder
outfolder=hg19_hisat-mapping
mkdir -p ${outfolder}

## reads
infolder=reads
f=shuffled_PE_NA18507_GAIIx_100_chr21
fq1=${infolder}/${f}_2_1.fq.gz
fq2=${infolder}/${f}_2_2.fq.gz

hisatpe=${f}_hisat-pe

# using 'nthr' processors in parallel
nthr=8
maxmem="48G"

# log all to file
logfile=${outfolder}/hisat_mapping-log.txt
exec &> >(tee -a "${logfile}")

# report run duration
starttime=$(date +%s)

# build HISAT command for zipped paired fastq files
mapcmd="hisat \
	-p ${nthr} \
	-x ${refidx} \
	-q \
	--phred33 \
	-1 ${fq1} \
	-2 ${fq2} \
	-S ${outfolder}/${hisatpe}.sam"

echo
echo "# ${mapcmd}"
eval ${mapcmd}

########################### post process results ##############
## Add RG tag in header
rgstring='@RG\tID:NA18507\tLB:lib-NA18507\tPU:unkn-0.0\tPL:ILLUMINA\tSM:GAIIx-chr21-HISAT'

cmd="bam polishBam \
	--RG \"${rgstring}\" \
	--CO \"@CO\t${mapcmd}\" \
	-i ${outfolder}/${hisatpe}.sam \
	-o ${outfolder}/${hisatpe}_q.bam \
	-l ${outfolder}/polishBam_${hisatpe}.log"

# execute the command
echo
echo "# ${cmd}"
eval ${cmd}

##### convert to sorted bam, index & cleanup
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar SortSam \
	I=${outfolder}/${hisatpe}_q.bam \
	O=${outfolder}/${hisatpe}-mappings.bam \
	SO=coordinate \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=LENIENT"

echo
echo "# ${cmd}"
eval ${cmd}

## extract chr21 mappings for comparison
chr="chr21"
thr=8
cmd="samtools view -bh \
	${outfolder}/${hisatpe}-mappings.bam ${chr} | \
	samtools sort -@ ${thr} - ${outfolder}/${chr}_hisat-mappings \
	&& samtools index ${outfolder}/${chr}_hisat-mappings.bam"

echo
echo "# ${cmd}"
eval ${cmd}

## ValidateSamFile ##
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar \
	ValidateSamFile \
	I=${outfolder}/${chr}_hisat-mappings.bam \
	O=${outfolder}/${chr}_hisat-mappings_ValidateSam.txt \
	M=SUMMARY"

echo
echo "# ${cmd}"
eval ${cmd}

##### get basic stats from the resulting bam file
cmd="samtools flagstat \
	${outfolder}/${chr}_hisat-mappings.bam \
	>${outfolder}/${chr}_hisat-mappings_flagstat.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## samtools stats ##
cmd="samtools stats \
	${outfolder}/${chr}_hisat-mappings.bam \
	>${outfolder}/${chr}_hisat-mappings_stats.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## samtools plot-bamstats ##
cmd="plot-bamstats \
	${outfolder}/${chr}_hisat-mappings_stats.txt \
	-p ${outfolder}/${chr}_hisat-bamstats_plots/hisat"

echo
echo "# ${cmd}"
eval ${cmd}

## plot mapping quality score distribution
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar \
	QualityScoreDistribution \
	CHART=${outfolder}/QS_distribution.pdf \
	I=${outfolder}/${chr}_hisat-mappings.bam \
	O=${outfolder}/${chr}_hisat-mappings_QS.txt"

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
