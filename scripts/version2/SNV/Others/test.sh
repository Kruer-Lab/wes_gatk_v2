#!/bin/bash
source settings_snv.sh

#FASTQ1=/media/sheetalshetty/Data2/Pipeline_Data2/Fastqs/Set2/M_F504-001-U/M_F504-001-U_BH3LNLDSXX_L003_R1_001
#FASTQ2=/media/sheetalshetty/Data2/Pipeline_Data2/Fastqs/Set2/M_F504-001-U/M_F504-001-U_BH3LNLDSXX_L003_R2_001
#OUTPUTS=/media/sheetalshetty/Data/Pipeline_V2/Data/Outputs/M_F504-001-U
#SAMPLE_ID="M_F504-001-U"
#SAMPLE_NAME="M_F504-001-U"
#SAMPLE_ID_PATH=${OUTPUTS}/${SAMPLE_ID}
#SAMPLE_ID_PATH=/media/sheetalshetty/Data/Pipeline_V2/Data/Outputs/M_F504-001-U/M_F504-001-U

echo -e "Please enter sample directory name:"
SAMPLE_DIR_NAME=$1
FASTQ="${ALL_SAMPLES_FASTQ_DIR}/${SAMPLE_DIR_NAME}"
SAMPLE_ID=$(ls ${FASTQ}/*_R1_*.fastq.gz | awk -F"-DNA" '{print $1}'| awk -F"/" '{print $(NF)}')
echo "sample ID : $SAMPLE_ID"
SAMPLE_NAME=$(echo $SAMPLE_ID | awk -F "_" '{print $1}')
OUTPUTS="${ALL_SAMPLES_OUTPUT_DIR}/${SAMPLE_DIR_NAME}"
FASTQ1=${FASTQ}/${SAMPLE_ID}_R1_001
FASTQ2=${FASTQ}/${SAMPLE_ID}_R2_001
SAMPLE_ID_PATH=${OUTPUTS}/${SAMPLE_ID}


echo "sample name: $SAMPLE_NAME"
echo "sample name: $SAMPLE_ID_PATH"

mkdir -p $OUTPUTS
