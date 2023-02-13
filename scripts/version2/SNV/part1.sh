
#!/bin/bash
source settings_snv.sh

echo -e "Please enter sample directory name:"
SAMPLE_DIR_NAME=$1
FASTQ="${ALL_SAMPLES_FASTQ_DIR}/${SAMPLE_DIR_NAME}"
SAMPLE_ID=$(ls ${FASTQ}/*_R1_*.fastq.gz | awk -F"_R1" '{print $1}'| awk -F"/" '{print $(NF)}')
echo "sample ID : $SAMPLE_ID"
SAMPLE_NAME=$(echo $SAMPLE_ID | awk -F "_" '{print $1}')
OUTPUTS="${ALL_SAMPLES_OUTPUT_DIR}/${SAMPLE_DIR_NAME}"
FASTQ1=${FASTQ}/${SAMPLE_ID}_R1_001
FASTQ2=${FASTQ}/${SAMPLE_ID}_R2_001
SAMPLE_ID_PATH=${OUTPUTS}/${SAMPLE_ID}


echo "sample name: $SAMPLE_NAME"
echo "sample IF PATH: $SAMPLE_ID_PATH"

mkdir -p $OUTPUTS

 

echo -e "\n\n\n************************** Align to the reference genome and sort the BAM file **************************\n"
bwa mem -M -t $NCORE $REF <(gunzip -c ${FASTQ1}.fastq.gz) <(gunzip -c ${FASTQ2}.fastq.gz) | samtools view -bS | samtools sort -@ 32 -m 6G > ${SAMPLE_ID_PATH}.sorted.bam


echo -e "\n\n\n************************** Make index of BAM file **************************\n"
samtools index ${SAMPLE_ID_PATH}.sorted.bam ${SAMPLE_ID_PATH}.sorted.bai

echo -e "\n\n\n ************************** Replace all read groups with a single new read group and assign all reads to this read group ************************** \n"
picard-tools AddOrReplaceReadGroups I=${SAMPLE_ID_PATH}.sorted.bam O=${SAMPLE_ID_PATH}.grp.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$SAMPLE_NAME

echo -e "\n\n\n************************** Make index of fixed group file **************************\n"
samtools index ${SAMPLE_ID_PATH}.grp.bam

echo -e "\n\n\n ************************** Mark duplicates **************************\n"
picard-tools MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 METRICS_FILE=${SAMPLE_ID_PATH}.metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT INPUT=${SAMPLE_ID_PATH}.grp.bam OUTPUT=${SAMPLE_ID_PATH}.dedup.bam


echo -e "\n\n\n ************************** Make index of deduplicated bam file **************************\n"
samtools index ${SAMPLE_ID_PATH}.dedup.bam

echo -e "\n\n\n ************************** Base Quality Score Recalibration 1 **************************\n"
$GATK --java-options "-Xmx7g" BaseRecalibrator -R $REF -I ${SAMPLE_ID_PATH}.dedup.bam --known-sites $DBSNP -O ${SAMPLE_ID_PATH}.bqrecal.table


echo -e "\n\n\n ************************** Base Quality Score Recalibration 2 **************************\n"
$GATK --java-options "-Xmx7g" ApplyBQSR -R $REF -I ${SAMPLE_ID_PATH}.dedup.bam --bqsr-recal-file ${SAMPLE_ID_PATH}.bqrecal.table -O  ${SAMPLE_ID_PATH}.bqrecal.bam


echo -e "\n\n\n ************************** Germline variant Calling by HaplotypeCaller **************************\n"
$GATK --java-options "-Xmx4g" HaplotypeCaller -R $REF -I ${SAMPLE_ID_PATH}.bqrecal.bam -ERC BP_RESOLUTION  -L $INTERVAL_LIST -O ${SAMPLE_ID_PATH}.raw.g.vcf














