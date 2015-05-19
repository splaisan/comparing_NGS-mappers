#! /usr/bin/env bash
## script: 'mapping_star.sh'
## ©SP-BITS, 2015 v1.0
# last edit: 2015-01-07

# required:
# rnaSTAR version STAR_2.4.0h1
# Picard Version (broad build): 1.127

## define global variables
reffasta=$STAR_INDEXES/hg19_24.fa
refgtf=$STAR_INDEXES//hg19_24.gtf
readlen=100
refidx=$STAR_INDEXES/index_${readlen}

# create new folder
outfolder=hg19_star-mapping
mkdir -p ${outfolder}

## reads
infolder=reads
f=gcat_set_041

fq1=${infolder}/${f}_1.fq.gz
fq2=${infolder}/${f}_2.fq.gz

starpe=${f}_star-pe

# using 'nthr' processors in parallel (again limited by our RAM!)
# mem is more demanding and needs more than 3Gb per thread
nthr=8
maxmem="48G"

# log all to file
logfile=${outfolder}/star_mapping-log.txt
exec &> >(tee -a "${logfile}")

# report run duration
starttime=$(date +%s)

# build STAR command for zipped paired fastq files, keep unmapped reads
## compute length-1
idxlen=$( echo ${readlen}-1 | bc)

cmd="STAR \
	--runMode alignReads \
	--runThreadN ${nthr} \
	--genomeDir ${refidx} \
	--genomeFastaFiles ${reffasta} \
	--sjdbGTFfile ${refgtf} \
	--sjdbOverhang ${idxlen} \
	--readFilesCommand zcat \
	--readFilesIn ${fq1} ${fq2} \
	--outSAMunmapped Within \
	--outFileNamePrefix ${outfolder}/${starpe}_\
	--outSAMtype BAM Unsorted"

echo
echo "# ${cmd}"
eval ${cmd}

########################### post process results ##############
## Add RG tag in header
rgstring="@RG\tID:NA12878\tSM:gcat_set_041-STAR\tPL:ILLUMINA"

cmd="bam polishBam --RG \"${rgstring}\" \
	-i ${outfolder}/${starpe}_Aligned.out.bam \
	-o ${outfolder}/${starpe}_Aligned.out-rg.bam \
	-l ${outfolder}/polishBam_${starpe}.log && \
	rm ${outfolder}/${starpe}_Aligned.out.bam"

# execute the command
echo
echo "# ${cmd}"
eval ${cmd}

##### convert to sorted bam, index & cleanup
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar SortSam \
	I=${outfolder}/${starpe}_Aligned.out-rg.bam \
	O=${outfolder}/${starpe}.bam \
	SO=coordinate \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=LENIENT"

echo
echo "# ${cmd}"
eval ${cmd}

## ValidateSamFile ##
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar ValidateSamFile \
	I=${outfolder}/${starpe}.bam \
	O=${outfolder}/${starpe}-mappings_ValidateSam.txt \
	M=SUMMARY"

echo
echo "# ${cmd}"
eval ${cmd}

##### get basic stats from the resulting bam file
cmd="samtools flagstat \
	${outfolder}/${starpe}.bam \
	>${outfolder}/${starpe}-mappings_flagstat.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## samtools stats ##
cmd="samtools stats \
	${outfolder}/${starpe}.bam \
	>${outfolder}/${starpe}-mappings_stats.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## samtools plot-bamstats ##
cmd="plot-bamstats \
	${outfolder}/${starpe}-mappings_stats.txt \
	-p ${outfolder}/${starpe}-bamstats_plots/star"

echo
echo "# ${cmd}"
eval ${cmd}

## plot mapping quality score distribution
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar \
	QualityScoreDistribution \
	CHART=${outfolder}/QS_distribution.pdf \
	I=${outfolder}/${starpe}.bam \
	O=${outfolder}/${starpe}-mappings_QS.txt"

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
