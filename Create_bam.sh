#!/bin/bash

REF="control"
PATH_TO_PICARD=$1
PATH_TO_GATK=$2
FILE=$3
OUT_DIR=$4
TMP_DIR=$5
pp=$6

cd $TMP_DIR

if [ $pp == "p" ]
then
	bwa mem -t 10 ${OUT_DIR}/${REF}.fa ${OUT_DIR}/${FILE}.1.fq ${OUT_DIR}/${FILE}.2.fq > ${FILE}.sam
    samtools view -Su ${FILE}.sam | samtools sort - | samtools rmdup /dev/stdin ${FILE}.sr.rm.bam
else
	bwa mem -t 10 ${OUT_DIR}/${REF}.fa ${OUT_DIR}/${FILE}.fq > ${FILE}.sam
    samtools view -Su ${FILE}.sam | samtools sort - | samtools rmdup -s /dev/stdin ${FILE}.sr.rm.bam
fi

java -jar ${PATH_TO_PICARD}/picard.jar AddOrReplaceReadGroups \
    I= ${FILE}.sr.rm.bam \
    O= ${FILE}.sr.rm.hd.bam \
    SORT_ORDER=coordinate \
    RGID=foo \
    RGLB=bar \
    RGPL=illumina \
    RGSM=Sample1 \
    RGPU=L001 \
    CREATE_INDEX=True

java -jar ${PATH_TO_GATK}/GenomeAnalysisTK.jar -T RealignerTargetCreator \
-R ${OUT_DIR}/${REF}.fa -I ${FILE}.sr.rm.hd.bam -o ${FILE}_IndelRealigner.intervals

java -jar ${PATH_TO_GATK}/GenomeAnalysisTK.jar -T IndelRealigner \
-R ${OUT_DIR}/${REF}.fa -I ${FILE}.sr.rm.hd.bam -targetIntervals ${FILE}_IndelRealigner.intervals \
-o ${OUT_DIR}/${FILE}.final.bam