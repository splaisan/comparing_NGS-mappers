#! /usr/bin/env bash
## script: 'mapping_rtg.sh'
## ©SP-BITS, 2013 v1.0
# last edit: 2015-01-05

# required:
# rtg version: v3.4 build 75c3e66 (2014-12-23)
# Picard Version (broad build): 1.127

# report run duration
starttime=$(date +%s)

## define global variables
refgen=$RTG_INDEXES/hg19_sdf

## your picard sample
infolder=reads
f=shuffled_PE_NA18507_GAIIx_100_chr21

# create new folder
outfolder=hg19_rtg-mapping
mkdir -p ${outfolder}

# log all to file
logfile=${outfolder}/rtg_mapping-log.txt
exec &> >(tee -a "${logfile}")

# common
fq1=${infolder}/${f}_2_1.fq.gz
fq2=${infolder}/${f}_2_2.fq.gz

# marking secondary hits with -M to comply with Picard
# using 'nthr' processors in parallel (again limited by our RAM!)
# mem is more demanding and needs more than 3Gb per thread
nthr=8
maxmem="48G"

rgstring="@RG\tID:NA18507\tSM:GAIIx-chr21-RTG\tPL:ILLUMINA"

echo
echo "# mapping paired reads with **RTG**"

## default numeric settings are left unchanged as:
# --aligner-band-width=FLOAT   aligner indel band width scaling factor,
#                                    fraction of read length allowed as an indel
#                                    (Default is 0.5)
# --aligner-mode=STRING        pick the aligner to use (Must be one of
#                                    [auto, table, general]) (Default is auto)
# --bed-regions=FILE           restrict calibration to mappings falling
#                                    within the supplied BED regions
# --gap-extend-penalty=INT     penalty for a gap extension during alignment
#                                    (Default is 1)
# --gap-open-penalty=INT       penalty for a gap open during alignment
#                                    (Default is 19)
# -c, --indel-length=INT           guaranteed number of positions that will be
#                                    detected in a single indel (Default is 1)
# -b, --indels=INT                 guaranteed minimum number of indels which
#                                    will be detected (Default is 1)
# -M, --max-fragment-size=INT      maximum permitted fragment size when mating
#                                    paired reads (Default is 1000)
# -m, --min-fragment-size=INT      minimum permitted fragment size when mating
#                                    paired reads (Default is 0)
# --mismatch-penalty=INT       penalty for a mismatch during alignment
#                                    (Default is 9)
# -d, --orientation=STRING         orientation for proper pairs (Must be one of
#                                    [fr, rf, tandem, any]) (Default is any)
# --pedigree=FILE              genome relationships pedigree containing sex
#                                    of sample
# --repeat-freq=INT            maximum repeat frequency (Default is 90%)
# --sex=SEX                    sex of sample (Must be one of [male, female,
#                                    either])
# --soft-clip-distance=INT     soft clip alignments if indels occur INT bp
#                                    from either end (Default is 5)
# -s, --step=INT                   step size (Default is word size)
# -a, --substitutions=INT          guaranteed minimum number of substitutions
#                                    which will be detected (Default is 1)
# --unknowns-penalty=INT       penalty for unknown nucleotides during
#                                    alignment (Default is 5)
# -w, --word=INT                   word size (Default is 22, or read length /
#                                    2, whichever is smaller)

rtgpe=${f}_rtg-pe

# store the full command line to include it in the next part
cmd="rtg RTG_MEM=${maxmem} map \
	-T ${nthr} \
	-o ${outfolder}/results \
	--sam-rg="\"${rgstring}\"" \
	-t ${refgen} \
	-F fastq \
	-q sanger \
	-l ${fq1} \
	-r ${fq2} \
	--read-names \
	--no-unmapped \
	--no-unmated \
	-a 1 \
	-b 1 \
	-c 1 \
	-d fr"

# execute the command
echo
echo "# ${cmd}"
eval ${cmd}

########################### post process results ##############
##### convert to sorted bam, index & cleanup
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar SortSam \
	I=${outfolder}/results/alignments.bam \
	O=${outfolder}/results/${rtgpe}_querysrt.bam \
	SO=queryname \
	VALIDATION_STRINGENCY=LENIENT"

# execute the command
echo
echo "# ${cmd}"
eval ${cmd}

## FixMateInformation
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar FixMateInformation \
	I=${outfolder}/results/${rtgpe}_querysrt.bam \
	O=${outfolder}/${rtgpe}.bam \
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
	${outfolder}/${rtgpe}.bam ${chr} | \
	samtools sort -@ ${thr} - ${outfolder}/${chr}_rtg-mappings && \
	samtools index ${outfolder}/${chr}_rtg-mappings.bam"

echo
echo "# ${cmd}"
eval ${cmd}

## ValidateSamFile ##
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar ValidateSamFile \
	I=${outfolder}/${chr}_rtg-mappings.bam \
	O=${outfolder}/${chr}_rtg-mapping_ValidateSam.txt \
	M=SUMMARY"

echo
echo "# ${cmd}"
eval ${cmd}

##### get basic stats from the resulting bam file
cmd="samtools flagstat \
	${outfolder}/${chr}_rtg-mappings.bam \
	>${outfolder}/${chr}_rtg-mapping_flagstat.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## samtools stats ##
cmd="samtools stats \
	${outfolder}/${chr}_rtg-mappings.bam \
	>${outfolder}/${chr}_rtg-mapping_stats.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## plot mapping quality score distribution
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar QualityScoreDistribution \
	CHART=${outfolder}/QS_distribution.pdf \
	I=${outfolder}/${chr}_rtg-mappings.bam \
	O=${outfolder}/${chr}_rtg-mapping_QS.txt"

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
