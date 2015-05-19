#! /usr/bin/env bash
## script: 'mapping_bowtie2.sh'
## ©SP-BITS, 2013 v1.0
# last edit: 2015-01-05

# required:
# bowtie2 version: 2.2.4
# Picard Version (broad build): 1.127
# bamutils Version: 1.0.12b | https://github.com/statgen/bamUtil

# report run duration
starttime=$(date +%s)

## define global variables
refgen=$BOWTIE2_INDEXES/hg19

## your picard sample
infolder=reads
f=shuffled_PE_NA18507_GAIIx_100_chr21

# create new folder
outfolder=hg19_bowtie2-mapping
mkdir -p ${outfolder}

# log all to file
logfile=${outfolder}/bowtie2_mapping-log.txt
exec &> >(tee -a "${logfile}")

# common
fq1=${infolder}/${f}_2_1.fq.gz
fq2=${infolder}/${f}_2_2.fq.gz

# marking secondary hits with -M to comply with Picard
# using 'nthr' processors in parallel (again limited by our RAM!)
# mem is more demanding and needs more than 3Gb per thread
nthr=8
maxmem="48G"

echo
echo "# mapping paired reads with **bowtie2**"

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

bowpe=${f}_bowtie2-pe

# store the full command line to include it in the next part
cmd="bowtie2 \
	-x  $BOWTIE2_INDEXES/hg19 \
	-q \
	-1 ${fq1} \
	-2 ${fq2} \
	--phred33 \
	--sensitive \
	--sensitive-local \
	-p ${nthr} \
	-S ${outfolder}/${bowpe}.sam \
	--rg-id NA18507 \
	--rg SM:GAIIx-chr21-bowtie2 \
	--rg PL:ILLUMINA"

# execute the command
echo "# ${cmd}"
eval ${cmd}

########################### post process results ##############
##### convert to sorted bam, index & cleanup
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar SortSam\
	I=${outfolder}/${bowpe}.sam \
	O=${outfolder}/${bowpe}_querysrt.bam \
	SO=queryname \
	VALIDATION_STRINGENCY=LENIENT"

# execute the command
echo "# ${cmd}"
eval ${cmd}

## FixMateInformation
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar \
	FixMateInformation \
	I=${outfolder}/${bowpe}_querysrt.bam \
	O=${outfolder}/${bowpe}.bam \
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
	${outfolder}/${bowpe}.bam ${chr} | \
	samtools sort -@ ${thr} - ${outfolder}/${chr}_bowtie2-mappings && \
	samtools index ${outfolder}/${chr}_bowtie2-mappings.bam"

echo
echo "# ${cmd}"
eval ${cmd}

## ValidateSamFile ##
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar ValidateSamFile \
	I=${outfolder}/${chr}_bowtie2-mappings.bam \
	O=${outfolder}/${chr}_bowtie2-mappings_ValidateSam.txt \
	M=SUMMARY"

echo
echo "# ${cmd}"
eval ${cmd}

##### get basic stats from the resulting bam file
cmd="samtools flagstat \
	${outfolder}/${chr}_bowtie2-mappings.bam \
	>${outfolder}/${chr}_bowtie2-mappings_flagstat.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## samtools stats ##
cmd="samtools stats \
	${outfolder}/${chr}_bowtie2-mappings.bam \
	>${outfolder}/${chr}_bowtie2-mappings_stats.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## plot mapping quality score distribution
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar \
	QualityScoreDistribution \
	CHART=${outfolder}/QS_distribution.pdf \
	I=${outfolder}/${chr}_bowtie2-mappings.bam \
	O=${outfolder}/${chr}_bowtie2-mappings_QS.txt"

echo
echo "# ${cmd}"
eval ${cmd}

#=====================
# report run duration
#=====================

endtime=$(date +%s)
dur=$(echo ${endtime}-${starttime} | bc)
echo "# ended at $date: $endtime"
echo "Done in ${dur} sec"
