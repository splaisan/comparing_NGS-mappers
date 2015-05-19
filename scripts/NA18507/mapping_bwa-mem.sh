#! /usr/bin/env bash
## script: 'mapping_bwa-mem.sh'
## ©SP-BITS, 2013 v1.0
# last edit: 2015-01-05

# required:
# bwa version: 0.7.10-r789
# Picard Version (broad build): 1.127
# bamutils Version: 1.0.12b | https://github.com/statgen/bamUtil

# report run duration
starttime=$(date +%s)

## define global variables
refgen=$BWA_INDEXES/hg19.fa

## your picard sample
infolder=reads
f=shuffled_PE_NA18507_GAIIx_100_chr21

# create new folder
outfolder=hg19_bwa-mapping
mkdir -p ${outfolder}

# log all to file
logfile=${outfolder}/bwa_mapping-log.txt
exec &> >(tee -a "${logfile}")

# common
fq1=${infolder}/${f}_2_1.fq.gz
fq2=${infolder}/${f}_2_2.fq.gz

# marking secondary hits with -M to comply with Picard
# using 'nthr' processors in parallel (again limited by our RAM!)
# mem is more demanding and needs more than 3Gb per thread
nthr=8
maxmem="48G"

rgstring='@RG\tID:NA18507\tLB:lib-NA18507\tPU:unkn-0.0\tPL:ILLUMINA\tSM:GAIIx-chr21-BWA.mem'

echo
echo "# mapping paired reads with **bwa mem**"

## default numeric settings are left unchanged as:
#  -k INT     minimum seed length [19]
#  -w INT     band width for banded alignment [100]
#  -d INT     off-diagonal X-dropoff [100]
#  -r FLOAT   look for internal seeds inside a seed longer than {-k}
#             * FLOAT [1.5]
#  -c INT     skip seeds with more than INT occurrences [10000]
#  -A INT     score for a sequence match [1]
#  -B INT     penalty for a mismatch [4]
#  -O INT     gap open penalty [6]
#  -E INT     gap extension penalty; a gap of size k cost {-O} + {-E}*k [1]
#  -L INT     penalty for clipping [5]
#  -U INT     penalty for an unpaired read pair [17]
#  -T INT     minimum score to output [30]
# '-p': first query file consists of interleaved paired-end sequences
#
# '-M': mark shorter split hits as secondary (for Picard/GATK compatibility)

bwape=${f}_bwa-mem-pe

# store the full command line to include it in the next part
cmd="bwa mem -R \"${rgstring}\" \
	-M \
	-t ${nthr} \
	${refgen} \
	${fq1} ${fq2} > ${outfolder}/${bwape}.sam"

# execute the command
echo
echo "# ${cmd}"
eval ${cmd}

########################### post process results ##############
#### add @PG group to results with **bamUtil 'bam polishBam' **
# important formatting issues here
# $'\t' add a 'tab' character
# $'\'' add single quotes around ${cmd} to avoid it being interpreted
# finally, $pgstring is called quoted to be replaced by its building variables

# we first detect the current version of BWA to store it
bwaver=$(expr "$(bwa 2>&1)" : '.*Version:\ \(.*\)Contact.*')

# we then create a '@PG' line to store our action and include the 'BWA version number' AND the 'full command'
# clean extra spaces with sed and the POSIX character class [:space:]
# some explanation: [:space:] is functionally identical to [ \t\r\n\v\f]
#cli=$(echo ${cmd} | sed -e "s/[[:space:]]\+/ /g")

#pgstring=@PG$'\t'ID:BWA$'\t'PN:bwa$'\t'VN:${bwaver}$'\t'CL:$'\''${cli}$'\''
pgstring="@PG\tID:BWA\tPN:bwa\tVN:${bwaver}"

# clean up initial file after success to save space
cmd="bam polishBam --PG \"${pgstring}\" \
	-i ${outfolder}/${bwape}.sam \
	-o ${outfolder}/${bwape}_pg.sam \
	-l ${outfolder}/polishBam_${bwape}.log && rm ${outfolder}/${bwape}.sam"

# execute the command
echo
echo "# ${cmd}"
eval ${cmd}

##### convert to sorted bam, index & cleanup
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar SortSam\
	I=${outfolder}/${bwape}_pg.sam \
	O=${outfolder}/${bwape}.bam \
	SO=coordinate \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=LENIENT"

# execute the command
echo
echo "# ${cmd}"
eval ${cmd}

## extract chr21 mappings for comparison
chr="chr21"
thr=8
cmd="samtools view -bh \
	${outfolder}/${bwape}.bam ${chr} | \
	samtools sort -@ ${thr} - ${outfolder}/${chr}_bwa-mappings && \
	samtools index ${outfolder}/${chr}_bwa-mappings.bam"

echo
echo "# ${cmd}"
eval ${cmd}

## ValidateSamFile ##
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar ValidateSamFile \
	I=${outfolder}/${chr}_bwa-mappings.bam \
	O=${outfolder}/${chr}_bwa-mappings_ValidateSam.txt \
	M=SUMMARY"

echo
echo "# ${cmd}"
eval ${cmd}

##### get basic stats from the resulting bam file
cmd="samtools flagstat \
	${outfolder}/${chr}_bwa-mappings.bam \
	>${outfolder}/${chr}_bwa-mappings_flagstat.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## samtools stats ##
cmd="samtools stats \
	${outfolder}/${chr}_bwa-mappings.bam \
	>${outfolder}/${chr}_bwa-mappings_stats.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## plot mapping quality score distribution
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar QualityScoreDistribution \
	CHART=${outfolder}/QS_distribution.pdf \
	I=${outfolder}/${chr}_bwa-mappings.bam \
	O=${outfolder}/${chr}_bwa-mappings_QS.txt"

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
