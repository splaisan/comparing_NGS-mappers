#!/bin/sh

# gene-level count from BAM data
bam=$1
pref=$(basename $1 ".bam")

gff=$BOWTIE2_INDEXES/hg19_ensGene.gff

# only take correctly paired mappings from flag
# replace new cigar flags by 'M' to comply with HTseq
samtools view -f 0x0002 ${bam} | \
        awk -F '\t' 'BEGIN{OFS="\t"} NF >= 6 {gsub(/[=X]/, "M", $6);} {print}' | \
        htseq-count -f sam -r pos -s no -a 10 \
        -t exon -i gene_id -m intersection-strict - ${gff} \
         > ${bam%%.bam}-counts.txt

# summarize results
# compute descriptive stats using qstats
echo "# processing ${bam}" > ${pref}-counts_summary.txt
grep "^_" ${bam%%.bam}-counts.txt >> ${pref}-counts_summary.txt
echo "# coverage stats" >> ${pref}-counts_summary.txt
awk 'BEGIN{FS="\t"; OFS="\t"}{if($1~/ENSG/ && $2>0) print $2}' \
        ${bam%%.bam}-counts.txt | qstats -s >> ${pref}-counts_summary.txt

